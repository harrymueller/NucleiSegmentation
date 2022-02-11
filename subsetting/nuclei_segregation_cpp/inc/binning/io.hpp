#ifndef IO_
#define IO_

#include "image.hpp"

#include <stdio.h>
#include <iostream>
#include <fstream>

#include <opencv2/opencv.hpp>
using namespace cv;



namespace IO 
{
    Image read_markers(std::string);
    void sparse_segments(Mat, std::string);
    void bin_coords(Mat, std::string);
    int find_max_bin_id(Mat);
}

#endif