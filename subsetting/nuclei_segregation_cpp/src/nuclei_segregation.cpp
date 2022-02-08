#include "nuclei_segregation.hpp"

/* TODO erode seeds?
    Q. erode seeds?
    Improve algo
    Improve UI
    Q. WTF did unknown do previously???
    Rename seeds / markers to be consistant
*/
int main(int argc, char *argv[])
{
    // params
    OPTS* opts = (OPTS*) calloc(1, sizeof(OPTS));
    if (get_params(argc, argv, opts) != 0) exit(EXIT_FAILURE);

    // convert to C++ strings
    std::string output_dir = std::string(opts->output_dir);

    // read image
    Image orig(opts->input_file);
    orig.save(output_dir + "/1 original.png");

    // thresholding
    std::cout << "Applying thresholds...\n";
    Image thresholded = orig.duplicate();

    Image globalI = Thresholds::global_threshold(thresholded, 100);
    thresholded.applyMask(globalI);

    Image gaussianI = Thresholds::gaussian_threshold(thresholded, 41, 0.03);
    thresholded.applyMask(gaussianI);
    
    thresholded.save(output_dir + "/2 thresholded.png");

    globalI.~Image(); gaussianI.~Image(); // cleaning up

    // isolate known fg
    Image sure_fg = FG_BG::get_sure_fg(thresholded);
    sure_fg.save(output_dir + "/3 sure_fg.png");

    // get markers from fg
    Image markers = Watershed::get_markers(sure_fg);
    sure_fg.~Image();

    // apply watershed algo
    markers = Watershed::apply_watershed(thresholded, markers);

    // outline segments, display and save
    Image watershed_outlined = Watershed::outline(orig, markers);
    watershed_outlined.display();
    watershed_outlined.save(output_dir + "/4 segmented.png");


    exit(EXIT_SUCCESS);
}



