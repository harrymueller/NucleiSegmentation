#include "fg_bg.hpp"

Image FG_BG::get_sure_fg(Image thresholdedI, int distanceType, int maskSize)
{
    Mat dest;
    distanceTransform(thresholdedI.get_im(), dest, distanceType, maskSize, CV_8UC1);
   
    // normalise
    Mat norm;
    normalize(dest, norm, 255.0, 0.0, NORM_MINMAX);
    dest.~Mat();

    // convert back to CV_8UC1 image (shouldn't be needed?)
    Mat normDest;
    norm.convertTo(normDest, CV_8UC1);
    norm.~Mat();
    
    Image distTransform = Image(normDest);

    // gaussian threshold to get peaks
    Image thresholdedDist = Thresholds::gaussian_threshold(distTransform, 21, -30);
    thresholdedDist = thresholdedDist.applyMask(distTransform.get_im());

    return thresholdedDist;
}

Image FG_BG::get_sure_bg(Image thresholdedI, Mat kernel)
{
    Mat dest;
    dilate(thresholdedI.get_im(), dest, kernel, Point(-1, -1), 5);
    
    Image bg = Image(Mat(dest.size(), CV_8UC1, Scalar(255)));
    bg = bg.applyMask(dest);
    return bg;
}