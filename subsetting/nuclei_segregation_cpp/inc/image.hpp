#include <stdio.h>

// cv
#include <opencv2/opencv.hpp>
using namespace cv;

class Image {
    public:
        // attributes

        // constructor
        Image(std::string filename);
        Image(Mat im);

        // accessor methods
        Mat get_im();
        void set_im(Mat);
        
        // methods
        int display();
        int save(std::string);
        Image duplicate();


    private:
        // attributes
        Mat im;
};