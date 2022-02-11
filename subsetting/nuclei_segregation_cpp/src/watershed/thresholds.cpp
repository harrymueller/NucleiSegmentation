#include "thresholds.hpp"

/*##############################
    FUNCTIONS
##############################*/
Image Thresholds::global_threshold(Image i, double thresh, double maxValue, int type)
{
    Mat dest;
    int ret = threshold(i.get_im(), dest, thresh, maxValue, type);

    return Image(dest);
}

Image Thresholds::gaussian_threshold(Image i, int blockSize, double C, double maxValue, int adaptiveMethod, int thresholdType)
{
    Mat dest;
    adaptiveThreshold(i.get_im(), dest, maxValue, adaptiveMethod, thresholdType, blockSize, C);

    return Image(dest);
}
