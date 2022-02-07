#include "nuclei_segregation.hpp"
//#include "global_vars.h"

int main(int argc, char *argv[])
{
    // params
    OPTS* opts = (OPTS*) calloc(1, sizeof(OPTS));
    if (get_params(argc, argv, opts) != 0) exit(EXIT_FAILURE);

    // convert to C++ strings
    std::string output_dir = std::string(opts->output_dir);

    // read image
    Image orig(opts->input_file);
    
    // thresholding
    std::cout << "Applying thresholds...\n";

    Image globalI = Thresholds::global_threshold(orig, 100);
    Image gaussianI = Thresholds::gaussian_threshold(orig, 41, 0.03);
    
    Image thresholded = orig.applyMask(globalI.get_im());
    thresholded = thresholded.applyMask(gaussianI.get_im());
    
    thresholded.save(output_dir + "/thresholded.png");

    globalI.~Image(); gaussianI.~Image(); // cleaning up

    // isolating known fg and bg
    Image sure_fg = FG_BG::get_sure_fg(thresholded);
    sure_fg.save(output_dir + "/sure_fg.png");

    Image sure_bg = FG_BG::get_sure_bg(thresholded);
    sure_bg.save(output_dir + "/sure_bg.png");

    Image unknown = sure_bg.subtract(sure_fg);

    // watershed
    Image markers = Watershed::get_markers(unknown, sure_fg); // TODO check that this works
    
    Image base_image = orig.applyMask(thresholded.get_im());
    Mat base_im;
    base_image.get_im().convertTo(base_im, CV_8UC3);
    base_image = Image(base_im);
    
    watershed(base_image.get_im(), markers.get_im());
    markers.display();
    // TODO erode seeds?
    // TODO commenting in functions
    // TODO need Image class? or just use as namespace or smth?
        // extend Mat class?

    exit(EXIT_SUCCESS);
}



