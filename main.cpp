#include <iostream>
#include <cstdlib>
#include <chrono>
#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>

#include <ultimaille/all.h>
#define STB_IMAGE_IMPLEMENTATION
#include "stb/stb_image.h"
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb/stb_image_write.h"
#include "bin2mesh.h"
#include "voronoi.h"

using namespace UM;
//#include "bin2mesh.h"

#define FOR(i, n) for(int i = 0; i < n; i++)



void multiplyBy2(std::vector<std::vector<bool>>& image) {
    int height = image.size();
    if (height == 0) return;
    int width = image[0].size();
    std::vector<std::vector<bool>> newImage(2 * height, std::vector<bool>(2 * width, false));
    for (int y = 0; y < height; y++) {
        for (int x = 0; x < width; x++) {
            FOR(Y, 2) FOR(X, 2)
                newImage[2 * y + Y][2 * x + X] = image[y][x];
        }
    }
    image = newImage;
}

void writefile(const std::vector<std::vector<bool>>& image, const std::string& filename) {
    int height = image.size();
    if (height == 0) return;
    int width = image[0].size();
    unsigned char* img = (unsigned char*)malloc(height * width);
    FOR(y, height) {
        FOR(x, width) {
            img[y * width + x] = 255 * (int)!image[y][x];
        }
    }
    int channels = 1;
    stbi_write_png(filename.c_str(), width, height, channels, img, width * channels);
    free(img);
}

void readfile(const std::string& filename, std::vector<std::vector<bool>>& image) {
    image.clear();
    int width, height, channels;
    unsigned char *img = stbi_load(filename.c_str(), &width, &height, &channels, 1);
    if(img == NULL) {
        printf("Error in loading the image\n");
        exit(1);
    }
    printf("Loaded image with a width of %dpx, a height of %dpx and %d channels\n", width, height, channels);

    const int threshold = 128;
    
    FOR(y, height) {
        image.push_back(std::vector<bool>(width, false));
        FOR(x, width) {
            image[y][x] = !((int)img[y * width + x] >= threshold);
        }
    }
    suppress_mono_nei(image);
    multiplyBy2(image);
    stbi_image_free(img);
}

void suppressRedondance(std::vector<std::vector<bool>>& contour) {
    int height = contour.size();
    if (height == 0) return;
    int width = contour[0].size();
    bool square[3][3];
    FOR(y, height) {
        FOR(x, width) {
            int count_deltay = 0;
            int count_deltax = 0;
            if (!contour[y][x]) continue;
            for (int i = -1; i <= 1; i++) {
                for (int j = -1; j <= 1; j++) {
                    square[i + 1][j + 1] = contour[y + i][x + j];
                    if (abs(i) + abs(j) == 1) {
                        if (abs(i)) count_deltay++;
                        if (abs(j)) count_deltax++;
                    }
                }
            }
            square[1][1] = false;
            bool is_disconnected = false;
            for (int i = -1; i <= 1; i++) {
                for (int j = -1; j <= 1; j++) {
                    if(!square[i + 1][j + 1]) continue;
                    bool flag = true;
                    for (int I = -1; I <= 1; I++) {
                        for (int J = -1; J <= 1; J++) {
                            if(abs(I) + abs(J) == 0) continue;
                            if ((0 <= i + I + 1 && i + I + 1 < 3) && (0 <= j + J + 1 && j + J + 1 < 3) && square[i + I + 1][j + J + 1]) {
                                flag = false;
                            }
                        }
                    }
                    if (flag) is_disconnected = true;
                }
            }
            if (count_deltax && count_deltay && !is_disconnected) contour[y][x] = false;
        }
    }
    //suppress_mono_nei(contour);
}

void detectContour(const std::vector<std::vector<bool>>& image, std::vector<std::vector<bool>>& contour) {
    contour.clear();
    int height = image.size();
    if (height == 0) return;
    int width = image[0].size();
    FOR(y, height) {
        contour.push_back(std::vector<bool>(width, false));
        FOR(x, width) {
            if (!image[y][x]) continue;
            for (int i = -1; i <= 1; i++) {
                for (int j = -1; j <= 1; j++) {
                    if (abs(i) + abs(j) != 1) continue;

                    if ((0 <= y + i && y + i < height) && (0 <= x + j && x + j < width) && !image[y + i][x + j]) {
                        contour[y][x] = true;
                    }
                }
            }
        }
    }
}

void poly_centroid_area(const std::vector<vec2> &polygon, vec2 &centroid, double &area) {
    centroid = vec2(0,0);
    area = 0;
    int n = polygon.size();
    for (int i=n-1, j=0; j<n; i=j, j++) {
        double ai = polygon[i].x*polygon[j].y - polygon[j].x*polygon[i].y;
        centroid = centroid + (polygon[i]+polygon[j])*ai;
        area += ai;
    }
    centroid = centroid/(3*area);
    area /= 2;
}

int main(int argc, char** argv) {
    std::string filename = "filename";
    if (argc >= 2) filename = argv[1];
    std::vector<std::vector<bool>> image;
    readfile(filename, image);
    std::vector<std::vector<bool>> contour;
    detectContour(image, contour);
    suppressRedondance(contour);
    writefile(contour, filename + "_rest_contour.png");
    PolyLine polyline;
    int echantillonnage = 5;
    mesh_from_contour(contour, polyline, echantillonnage, filename + "_output.geogram");

    { // colocate vertices TODO: David, demerde-toi pour ne pas avoir de doublons :)
        std::vector<bool> to_kill(polyline.nverts(), true);
        { // TODO: ugly hack because PolyLine lacks delete_vertices
            std::vector<int> old2new;
            colocate(*polyline.points.data, old2new, 1e-3);
            for (int v : range(polyline.nverts()))
                to_kill[old2new[v]] = false;
        }
        std::vector<int> old2new;
        polyline.points.delete_points(to_kill, old2new);

        for (int s : range(polyline.nsegments()))
            for (int lv : range(2))
                polyline.vert(s, lv) = old2new[polyline.vert(s, lv)];
    }

    /*
    { // scale
        const double boxsize = 10.;
        const double shrink  = 1.3;
        vec3 bbmin, bbmax;
        polyline.points.util.bbox(bbmin, bbmax);
        float maxside = std::max(bbmax.x-bbmin.x, bbmax.y-bbmin.y);
        for (vec3 &p : polyline.points)
            for (int d : range(2))
                p[d] = (p[d] - (bbmax[d]+bbmin[d])/2.)*boxsize/(shrink*maxside) + boxsize/2.;
        write_geogram("scaled.geogram", polyline);
    }
    */

    PointAttribute<bool> lock(polyline.points);
    for (int s : range(polyline.nsegments()))
        for (int lv : range(2))
            lock[polyline.vert(s, lv)] = true;
    write_geogram("cleaned.geogram", polyline, {{{"selection", lock.ptr}}, {}});

    Polygons mvoro;
    { // compute Voronoi diagram
        vec3 bbmin, bbmax;
        polyline.points.util.bbox(bbmin, bbmax);
        std::vector<vec2> seeds;
        for (vec3 &p : polyline.points) {
            seeds.emplace_back(p.x, p.y);
        }

        for (int lloyditer : range(10)) {
            std::vector<std::vector<vec2> > vorverts(seeds.size());
            KNN<2> knn(seeds);

#pragma omp parallel for
            for (size_t i=0; i<seeds.size(); i++) {
                ConvexCell cc({bbmin.x, bbmin.y}, {bbmax.x, bbmax.y});
                bool res = voronoi_cell(knn, cc, vorverts[i], seeds, i, 16);
                assert(res);
            }

            for (size_t i=0; i<seeds.size(); i++) {
                if (lock[i]) continue;
                double site_area;
                vec2 site_centroid;
                poly_centroid_area(vorverts[i], site_centroid, site_area);
                seeds[i] = site_centroid;
                polyline.points[i].x = site_centroid.x;
                polyline.points[i].y = site_centroid.y;
            }

            mvoro = Polygons();
            int off = 0;
            int cnt = 0;
            for (auto vcell : vorverts) {
                for (vec2 &p : vcell)
                    mvoro.points.push_back({p.x, p.y, 0});
                off = mvoro.create_facets(1, vcell.size());
                for (int i : range(vcell.size())) {
                    mvoro.vert(off, i) = cnt++;
                }
            }
        }
    }

    { // glue vertices of the Voronoi diagram
        { // colocate vertices of the diagram
            std::vector<int> old2new;
            colocate(*mvoro.points.data, old2new, 1e-9);

            for (int t : facet_iter(mvoro)) {
                for (int lv : range(mvoro.facet_size(t))) {
                    mvoro.vert(t, lv) = old2new[mvoro.vert(t, lv)];
                }
            }
        }

        { // remove duplicates and degenerate facets
            Polygons mvoro2;
            mvoro2.points = mvoro.points;
            for (int f : facet_iter(mvoro)) {
                std::vector<int> new_facet;
                for (int lv : range(mvoro.facet_size(f))) {
                    if (mvoro.vert(f, lv)!=mvoro.vert(f, (lv+1)%mvoro.facet_size(f)))
                        new_facet.push_back(mvoro.vert(f, lv));
                }
                if (new_facet.size()<3) continue;
                int off_f = mvoro2.create_facets(1, new_facet.size());
                for (int lv : range(new_facet.size()))
                    mvoro2.vert(off_f, lv) = new_facet[lv];
            }
            mvoro.facets = mvoro2.facets;
            mvoro.offset = mvoro2.offset;
        }

        { // remove isolated vertices
            std::vector<bool> to_kill(mvoro.nverts(), true);
            for (int t : facet_iter(mvoro)) {
                for (int v : facet_vert_iter(mvoro, t)) {
                    to_kill[v] = false;
                }
            }
            mvoro.delete_vertices(to_kill);
        }
    }
    write_geogram("voronoi.geogram", mvoro);

    Polygons mdel;
    {  // compute a Delaunay-ish triangulation of the point set
        SurfaceConnectivity fec(mvoro);
        mdel.points = polyline.points;

        for (int v : vert_iter(mvoro)) {
            if (fec.is_boundary_vert(v)) continue;
            int cnt = 0;
            { // count incident Voronoi facets
                int cir = fec.v2c[v];
                do {
                    cir = fec.next_around_vertex(cir);
                    cnt++;
                } while (cir != fec.v2c[v]);
            }
            int off_f = mdel.create_facets(1, cnt);
            { // create dual facet
                int cnt = 0;
                int cir = fec.v2c[v];
                do {
                    mdel.vert(off_f, cnt++) = fec.c2f[cir];
                    cir = fec.next_around_vertex(cir);
                } while (cir != fec.v2c[v]);
            }
        }

    }

    { // no guarantees on the quality, simple fan-triangulate
        int nb_triangles = 0;
        for (int f : facet_iter(mdel))
            nb_triangles += (mdel.facet_size(f) - 2);

        std::vector<int> new_facets;
        new_facets.reserve(nb_triangles*3);

        for (int f : facet_iter(mdel)) {
            int v0 = mdel.vert(f, 0);
            for (int lv=1; lv+1<mdel.facet_size(f); lv++) {
                new_facets.push_back(v0);
                new_facets.push_back(mdel.vert(f, lv));
                new_facets.push_back(mdel.vert(f, lv+1));
            }
        }
        mdel.facets = new_facets;
        mdel.offset.resize(nb_triangles+1);
        for (int t=0; t<nb_triangles+1; t++)
            mdel.offset[t] = t*3;
    }

    write_geogram("del.geogram", mdel);
    return 0;
}
