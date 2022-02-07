#include "watershed.hpp"

Image Watershed::get_markers(Image unknown, Image sure_fg)
{
    Mat markers;
    int ret = connectedComponents(sure_fg.get_im(), markers);
    //ret, markers = cv.connectedComponents(sure_fg)

    //# Makes known background label = 1, not 0
    markers += Scalar(1);

    Mat markers0;
    markers.convertTo(markers0, CV_32SC1);

    //# set label of unknown pixels to 0
    Image markersIm = Image(markers0);
    markersIm = markersIm.applyMask(unknown.get_im());

    return markersIm;
}

Image Watershed::apply_watershed(Image i, Image markers)
{
    //base_image = mask_image(im, gaussian_thresh == 255)
    //base_image = cv.cvtColor(base_image, cv.COLOR_GRAY2BGR) # CV_8UC3
    //markers = watershed(base_image, markers, im, get_filepath("watershed.png"))
    return i;
}