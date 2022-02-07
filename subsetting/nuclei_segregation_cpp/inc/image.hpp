#ifndef IMAGE_
#define IMAGE_

#include <stdio.h>

// cv
#include <opencv2/opencv.hpp>
#include <opencv2/core.hpp>
using namespace cv;

class Image {
    public:
        // attributes

        // constructor
        Image(std::string filename);
        Image(Mat im);

        // destructor
        ~Image();

        // accessor methods
        Mat get_im();
        void set_im(Mat);
        
        // methods
        int display();
        int save(std::string);
        Image duplicate();

        Image applyMask(Mat);
        Image subtract(Image);

    private:
        // attributes
        Mat im; // image matrix
};

#endif