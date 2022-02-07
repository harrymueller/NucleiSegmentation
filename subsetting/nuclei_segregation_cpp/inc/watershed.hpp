#ifndef WATERSHED_
#define WATERSHED_

// cv
#include <opencv2/opencv.hpp>
#include <opencv2/core.hpp>
#include <opencv2/imgproc.hpp>
using namespace cv;

// image class
#include "image.hpp"

namespace Watershed {
    Image get_markers(Image, Image);
    Image apply_watershed(Image, Image);
}
#endif