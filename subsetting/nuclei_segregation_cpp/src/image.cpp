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
// displays this image in a window
int Image::display() 
{
    namedWindow("Display Image", WINDOW_AUTOSIZE );
    imshow("Display Image", this->im);
    waitKey(0);
    return 0;
}

// saves this image to the filename given
int Image::save(std::string filename) 
{
    bool x = imwrite(filename, this->im);
    std::cout << "Saved image to " << filename << "\n";
    return x ? 0 : 1;
}   

// returns a copy of this obj 
Image Image::duplicate() 
{
    return Image(this->im.clone());
}

// Pixels where mask != 0 are copied to a new Image obj
void Image::applyMask(Image mask)
{
    Mat dest;
    this->im.copyTo(dest, mask.get_im());
    this->set_im(dest);
}
void Image::applyMask(Mat mask)
{
    Mat dest;
    this->im.copyTo(dest, mask);
    this->set_im(dest);
}


// subtracts the provided image from this image
void Image::subtract(Image other) 
{
    cv::subtract(this->im, other.get_im(), this->im);
}

// converts this to a new colour (e.g. COLOR_GRAY2BGR)
void Image::convertColour(int new_colour)
{
    cvtColor(this->im, this->im, new_colour);
}

// converts this image to a new type (e.g. CV_8UC1)
void Image::convertType(int new_type)
{
    this->im.convertTo(this->im, new_type);
}   

// for each pixel in this, call func at each point w/ <this->im, image.get_im(), i, j>
void Image::loop(Image image, std::function<void(Mat, Mat, int, int)> func)
{
    int i,j;

    for ( i = 0; i < this->im.rows; i++ ) {
        for ( j = 0; j < this->im.cols; j++ ) {
            func(this->im, image.get_im(), i, j);
        }
    }
}