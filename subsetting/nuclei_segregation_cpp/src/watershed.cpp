#include "watershed.hpp"
#include <stdio.h>

/*
    Get the markers, or seed positions, for the watershed algorithm
*/
Image Watershed::get_markers(Image sure_fg)
{
    Mat markers;
    int ret = connectedComponents(sure_fg.get_im(), markers);
    
    return Image(markers);
}

/*
    Apply the watershed algo to segment i used seeds given by markers
*/
Image Watershed::apply_watershed(Image i, Image markers)
{
    // convert i from grey to RGB    
    i.convertColour(COLOR_GRAY2BGR);
    
    // ensure markers of correct format
    Image watershed_markers = markers.duplicate();
    watershed_markers.convertType(CV_32SC1);

    // apply watershed
    watershed(i.get_im(), watershed_markers.get_im());

    // reset colour of i
    i.convertColour(COLOR_BGR2GRAY);

    return watershed_markers;
}

/*
    Function to pass to Image::loop for Watershed::outline
*/
void watershed_loop(Mat im1, Mat im2, int i, int j) {
    int val = im2.at<int>(i,j);
    if (val == -1) {
        im1.at<Vec3b>(i,j) = Vec3b(0,0,255);
    }
}

/*
    Save the image original with outlines in red of the segments, as determined by watershed
*/
Image Watershed::outline(Image original, Image watershed) 
{
    // duplicate first
    Image outlined = original.duplicate();

    // ensure type is CV_8UC3
    outlined.convertColour(COLOR_GRAY2BGR);

    outlined.loop(watershed, &watershed_loop);
    return outlined;
}