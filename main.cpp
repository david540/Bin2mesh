#include <iostream>
#include <cstdlib>
#include <chrono>
#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>

#include <ultimaille/io/medit.h>
#include <ultimaille/io/geogram.h>
#define STB_IMAGE_IMPLEMENTATION
#include "stb/stb_image.h"


//#include "bin2mesh.h"

#define FOR(i, n) for(int i = 0; i < n; i++)



void readfile(const std::string& filename, std::vector<std::vector<bool>>& image) {
    int width, height, channels;
    unsigned char *img = stbi_load(filename.c_str(), &width, &height, &channels, 1);
    if(img == NULL) {
        printf("Error in loading the image\n");
        exit(1);
    }
    printf("Loaded image with a width of %dpx, a height of %dpx and %d channels\n", width, height, channels);

    const int threshold = 128;
    ;
    FOR(y, height) {
        image.push_back(std::vector<bool>());
        FOR(x, width) {
            image[y].push_back((int)img[y * width + x] >= threshold);
        }
    }
}

int main(int argc, char** argv) {
    std::string filename = "filename";
    if (argc >= 2) filename = argv[1];
    std::vector<std::vector<bool>> image;
    readfile(filename, image);

  

    return 0;
}
