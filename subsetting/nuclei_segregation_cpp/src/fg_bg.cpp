#include "fg_bg.hpp"

/*
    Get sure foreground by performing a euclidean dist transform, 
        then a gaussian threshold to get the peaks
*/
Image FG_BG::get_sure_fg(Image thresholdedI, int distanceType, int maskSize)
{
    Mat distM;
    // perform distance transformation
    distanceTransform(thresholdedI.get_im(), distM, distanceType, maskSize, CV_8UC1);
   
    // normalise
    normalize(distM, distM, 255.0, 0.0, NORM_MINMAX);

    // convert back to CV_8UC1 image (shouldn't be needed?)
    Image distTransform = Image(distM);
    distTransform.convertType(CV_8UC1);
    //thresholdedI.display();
    // gaussian threshold to get peaks
    //Image thresholdedDist = distTransform.duplicate();
    // apply gaussian threshold to distTransform
    Image gaussThres = Thresholds::gaussian_threshold(distTransform, 21, -30); 
    // use distTransform (bin image) to get dist peaks
    //distTransform.applyMask(gaussThres);

    return gaussThres;
}

/*
    Get the sure background by dilating the thresholded image
*/
Image FG_BG::get_sure_bg(Image thresholdedI, Mat kernel)
{
    Mat dest;
    dilate(thresholdedI.get_im(), dest, kernel, Point(-1, -1), 5);
    
    Image bg = Image(Mat(dest.size(), CV_8UC1, Scalar(255)));
    bg.applyMask(dest);
    return bg;
}