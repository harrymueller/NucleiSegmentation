#ifndef IMAGE_
#define IMAGE_

#include <stdio.h>

// cv
#include <opencv2/opencv.hpp>
#include <opencv2/core.hpp>
using namespace cv;

/*
    Image class:
        stores a cv::Mat obj, and provides many methods for the given matrix
*/
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

        void applyMask(Image);
        void applyMask(Mat);

        void subtract(Image);

        void convertColour(int);
        void convertType(int);

        void loop(Image, std::function<void(Mat, Mat, int, int)>);
    private:
        // attributes
        Mat im; // image matrix
};

#endif