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
    Image i(opts->input_file);
    
    // thresholding
    // global_thresh = thresholding.global_thresholding(im, val = 100, filename = get_filepath("global_threshold.png"))
    // gaussian_thresh = thresholding.gaussian_thresholding(mask_image(im, global_thresh == 255), extra_thresh = global_thresh, blockSize = 41, C = 0.03, filename = get_filepath("gaussian_threshold.png"))
    std::cout << "Applying thresholds...\n";

    


    exit(EXIT_SUCCESS);
}



