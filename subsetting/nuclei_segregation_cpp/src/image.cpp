#include "image.hpp"


/*##############################
    CONSTRUCTORS
##############################*/
void isValid(Mat);

Image::Image(Mat im)
{
    this->im = im;
    isValid(this->im);
}

Image::Image(std::string filename)
{
    this->im = imread( filename, IMREAD_GRAYSCALE );
    isValid(this->im);
}

void isValid(Mat im) {
    if (!im.data) {
        printf("Invalid input image.\n");
        exit(EXIT_FAILURE);
    }
}

/*##############################
    DESTRUCTOR
##############################*/
Image::~Image()
{
    ~this->im;
}

/*##############################
    ACCESSOR METHODS
##############################*/
Mat Image::get_im()
{
    return this->im;
}

void Image::set_im(Mat im)
{
    this->im = im;
}

/*##############################
    METHODS
##############################*/
int Image::display() 
{
    namedWindow("Display Image", WINDOW_AUTOSIZE );
    imshow("Display Image", this->im);
    waitKey(0);
    return 0;
}

int Image::save(std::string filename) 
{
    bool x = imwrite(filename, this->im);
    std::cout << "Saved image to " << filename << "\n";
    return x ? 0 : 1;
}   

Image Image::duplicate() 
{
    return Image(this->im);
}

// Pixels where mask != 0 are copied to a new Image obj
Image Image::applyMask(Mat mask)
{
    Mat dest;
    copyTo(this->get_im(), dest, mask);
    return Image(dest);
}

// subtracts the provided image from this image
Image Image::subtract(Image other) 
{
    Mat dest;
    cv::subtract(this->get_im(), other.get_im(), dest);
    return Image(dest);
}