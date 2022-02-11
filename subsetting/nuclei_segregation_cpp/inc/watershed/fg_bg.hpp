#ifndef FG_BG_
#define FG_BG_

// cv
#include <opencv2/opencv.hpp>
#include <opencv2/core.hpp>
#include <opencv2/imgproc.hpp>
using namespace cv;

// image class
#include "image.hpp"

// thresholding namespace
#include "thresholds.hpp"

namespace FG_BG {
    Image get_sure_fg(Image, int distanceType = DIST_L2, int maskSize = DIST_MASK_PRECISE );
    Image get_sure_bg(Image, Mat kernel = Mat());
}
#endif