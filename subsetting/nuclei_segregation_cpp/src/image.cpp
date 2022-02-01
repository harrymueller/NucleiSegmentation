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
    this->im = imread( filename, 1 );
    isValid(this->im);
}

void isValid(Mat im) {
    if (!im.data) {
        printf("Invalid input image.\n");
        exit(EXIT_FAILURE);
    }
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
