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

    return 0;
}
