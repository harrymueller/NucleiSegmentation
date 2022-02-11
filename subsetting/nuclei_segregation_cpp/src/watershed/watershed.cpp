#include "watershed.hpp"
#include <stdio.h>

void add_one(Mat im, int i, int j)
{
    im.at<int>(i,j) += 1;
}

void set_unknown(Mat markers, Mat unknown, int i, int j)
{
    Scalar val = unknown.at<uchar>(i,j);
    if (val != Scalar(0)) {
        markers.at<int>(i,j) = 0;
    }
}

/*
    Get the markers, or seed positions, for the watershed algorithm
*/
Image Watershed::get_markers(Image sure_fg, Image unknown)
{
    Mat markers;
    int ret = connectedComponents(sure_fg.get_im(), markers);
    
    // add 1 to all -> Makes known background label = 1, not 0
    Image markersI = Image(markers);
    markersI.apply(add_one);

    // unknown == 0
    markersI.apply(unknown, set_unknown);

    return markersI;
}

/*
    Function to correct outside border - 
        checks all contour pixels for surrounding segment, sets to background if not found
*/
void correct_border(Mat segments, int i, int j)
{
    int val = segments.at<int>(i,j);
    if (val == -1) {
        if (i > 0 && segments.at<int>(i-1,j) > 1) return;
        else if (i < segments.rows - 1 && segments.at<int>(i+1,j) > 1) return;
        else if (j > 0 && segments.at<int>(i,j-1) > 1) return;
        else if (j < segments.cols - 1 && segments.at<int>(i, j+1) > 1) return;
        else segments.at<int>(i,j) = 1;
    }
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

    watershed_markers.apply(correct_border);

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

    outlined.apply(watershed, &watershed_loop);
    return outlined;
}