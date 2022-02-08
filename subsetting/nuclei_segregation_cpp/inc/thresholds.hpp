#ifndef THRESHOLDS_
#define THRESHOLDS_

#include <stdio.h>

// cv
#include <opencv2/opencv.hpp>
#include <opencv2/imgproc.hpp>
using namespace cv;

// image class
#include "image.hpp"

/*
    Wrappers for threshold functions in OpenCV
*/
namespace Thresholds {
    Image global_threshold(Image, double, double maxValue = 255, int type = THRESH_BINARY);
    Image gaussian_threshold(Image, int, double, double maxValue = 255, 
        int adaptiveMethod = ADAPTIVE_THRESH_GAUSSIAN_C, int thresholdType = THRESH_BINARY);
}
#endif