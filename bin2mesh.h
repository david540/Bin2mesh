#pragma once
#include <array>
#include <vector>
#include <ultimaille/pointset.h>
#include <ultimaille/surface.h>
#include <ultimaille/attributes.h>
#include <ultimaille/range.h>
#include <ultimaille/hboxes.h>

using namespace UM;
void mesh_from_contour(std::vector<std::vector<bool>>& contour, PolyLine& polyline, int nb_pixels_per_edge, std::string meshfilename);
void suppress_mono_nei(std::vector<std::vector<bool>>& contour);
void testNuagePoint(PolyLine& polyline, double distmin);