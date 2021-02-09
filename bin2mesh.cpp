#include "bin2mesh.h"
#include <map>
#include <chrono>
#include <ultimaille/all.h>
#include <algorithm>

#define FOR(i, n) for(int i = 0; i < n; i++)

using namespace UM;

void suppress_mono_nei(std::vector<std::vector<bool>> & image) {
    int height = image.size();
    if (height == 0) return;
    int threshold = 128;
    int width = image[0].size();
    bool flag;
    do {
        flag = false;
        for (int y = 0; y < height; y++) {
            for (int x = 0; x < width; x++) {
                if (image[y][x]) {
                    int count = 0;
                    for (int i = -1; i <= 1; i++) {
                        for (int j = -1; j <= 1; j++) {
                            if ((0 <= y + i && y + i < height) && (0 <= x + j && x + j < width) && image[y + i][x + j]) {
                                count++;
                            }
                        }
                    }
                    if (count <= 2) {
                        image[y][x] = false;
                        flag = true;
                    }
                }
            }
        }
    } while (flag == true);
}

void mesh_from_contour(std::vector<std::vector<bool>> & contour, PolyLine& polyline, int nb_pixels_per_edge, std::string meshfilename) {
    polyline.clear();
    int height = contour.size();
    if (height == 0) return;
    int width = contour[0].size();
    int start_x = 0, start_y = 0;

    suppress_mono_nei(contour);

    std::vector<std::vector<bool>> done(height, std::vector<bool>(width, false));
    bool flag = false;
    do{
        flag = false;
        {//detect point to start with
            for(int y = start_y; y < height; y ++) {
                if (flag) break;
                for (int x = 0; x < width; x++) {
                    if (contour[y][x] && !done[y][x]) {
                        start_x = x;
                        start_y = y;
                        flag = true;
                        break;
                    }
                }
            }
        }
        if (!flag) break;

        {//create point set with boundary starting by (start_x, start_y)
            int old_x = start_x, old_y = start_y;
            int cur_x = start_x, cur_y = start_y;
            int count_act = 0;
            int tot_count = 0;

            int id_startpoint = polyline.points.size();
            bool first_time = true;
            do {
                done[cur_y][cur_x] = true;
                bool flag = false;
                for (int i = -1; i <= 1; i++) {
                    if (flag) break;
                    for (int j = -1; j <= 1; j++) {
                        int new_x = cur_x + j, new_y = cur_y + i;
                        if (!(0 <= new_y && new_y < height) || !(0 <= new_x && new_x < width)) {
                            continue;
                        }
                        if (contour[new_y][new_x] && (new_x != cur_x || new_y != cur_y) && (new_x != old_x || new_y != old_y)) {
                            flag = true;
                            old_x = cur_x; old_y = cur_y;
                            cur_x = new_x; cur_y = new_y;
                            count_act++;
                            if (count_act >= nb_pixels_per_edge) {
                                polyline.points.push_back(vec3(new_x, height - 1 - new_y, 0));
                                if (!first_time) {
                                    int offset = polyline.create_segments(1);
                                    polyline.vert(offset, 0) = polyline.points.size() - 2;
                                    polyline.vert(offset, 1) = polyline.points.size() - 1;
                                }
                                else {
                                    first_time = false;
                                }
                                count_act = 0;
                            }
                            break;
                        }
                    }
                }
            } while (cur_x != start_x || cur_y != start_y);
            int offset = polyline.create_segments(1);
            polyline.vert(offset, 0) = polyline.points.size() - 1;
            polyline.vert(offset, 1) = id_startpoint;
        }
    } while (flag);

    std::cout << polyline.segments.size() << "\n";

    testNuagePoint(polyline, 2./3. * (double)nb_pixels_per_edge);
    write_geogram(meshfilename.c_str(), polyline);
}


int permx[32] = { 11,3,13,28,14,8,17,7,5,11,0,1,6,19,18,15,2,16,27,31,26,29,20,4,12,24,9,25,21,23,30,22 };
int permy[32] = { 20,31,25,9,15,23,18,27,30,29,21,17,28,14,24,11,26,16,4,20,5,7,12,10,3,0,1,8,22,13,2,6 };
vec3 permuted_grid(int i, int j, int n, double size) {
    return vec3(i / (double)n - (n - 1 - permy[j]) / (double)(n * n), j / (double)n - permx[i] / (double)(n * n), 0) * size;
}

bool intersectionSegment(vec2 A, vec2 B, vec2 C, vec2 D) {
    vec2 I = B - A; //vecteur AB
    vec2 J = D - C; //vecteur CD
    double k, m; // a determiner

    double denominateur = (I.x * J.y - I.y * J.x);
    if (denominateur < 0.00001) return false;
    //m = -(-I.x * A.y + I.x * C.y + I.y * A.x - I.y * C.x) / denominateur;
    m = (I.x * (A.y - C.y) - I.y * (A.x - C.x)) / denominateur;
    //k = -(A.x * J.y - C.x * J.y - J.x * A.y + J.x * C.y) / denominateur;
    k = (J.x * (A.y - C.y) - J.y * (A.x - C.x)) / denominateur;
    if (0 < m && m < 1 && 0 < k && k < 1) {
        return true;
    }return false;
}

void testNuagePoint(PolyLine& polyline, double distmin) {

    const int grid_size = 32;
    std::vector<vec3> seeds;
    
    const int sizeNuage = 200;
    int minX = polyline.points[0].x, minY = polyline.points[0].y, maxX = polyline.points[0].x, maxY = polyline.points[0].y;
    FOR(i, polyline.points.size()) {
        if (minX > polyline.points[i].x) minX = polyline.points[i].x;
        if (maxX < polyline.points[i].x) maxX = polyline.points[i].x;
        if (minY > polyline.points[i].y) minY = polyline.points[i].y;
        if (maxY < polyline.points[i].y) maxY = polyline.points[i].y;
        
    }
    vec2 pInfini(10000, 10000);
    //vec3 startPoint = polyline.points[0] - vec3(sizeNuage / 2, sizeNuage / 2, 0);
    for (int y = minY; y < maxY; y += sizeNuage) {
        for (int x = minX; x < maxX; x += sizeNuage) {
            vec3 startPoint(x, y, 0);
            std::cout << startPoint << "\n";
            for (int j = 0; j < grid_size; j++)
                for (int i = 0; i < grid_size; i++)
                    seeds.push_back(startPoint + permuted_grid(i, j, grid_size, sizeNuage));

            PointSet& ps = polyline.points;
            FOR(i, seeds.size()) {
                bool trop_proche = false;
                int count = 0;
                FOR(s, polyline.nsegments()){
                    vec3 P03D = polyline.points[polyline.vert(s, 0)];
                    vec3 P13D = polyline.points[polyline.vert(s, 1)];
                    vec2 P0(P03D.x, P03D.y);
                    vec2 P1(P13D.x, P13D.y);
                    if (intersectionSegment(P0, P1, vec2(seeds[i].x, seeds[i].y), pInfini) || intersectionSegment(P1, P0, vec2(seeds[i].x, seeds[i].y), pInfini)) {
                        count++;
                    }
                    double dist = std::min((vec2(seeds[i].x, seeds[i].y) - P1).norm(), (vec2(seeds[i].x, seeds[i].y) - P0).norm());
                    if (dist < distmin) trop_proche = true;
                }
                if (count % 2 != 0 && !trop_proche) {
                    ps.push_back(seeds[i]);
                }
            }
        }
    }

}
