#include "nuclei_segregation.hpp"

/* TODO erode seeds?
    Q. erode seeds?
    Improve algo
    Improve UI
    Q. WTF did unknown do previously???
    Rename seeds / markers to be consistant
    marker type is int -> may be problem if >2**31-1 nuclei
    improve apply image to take lots of variables and pass to function
    bin ids starting from 0?
    use either tsv or csv consistantly
    improve namespaces
*/

void all_watershed(OPTS* opts, std::string output_dir)
{
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

    // isolate known fg
    Image sure_fg = FG_BG::get_sure_fg(thresholded);
    sure_fg.save(output_dir + "/3 sure_fg.png");

    Image sure_bg = FG_BG::get_sure_bg(thresholded);
    sure_bg.subtract(sure_fg);
    sure_bg.save(output_dir + "/unknown.png");

    // get markers from fg
    Image markers = Watershed::get_markers(sure_fg, sure_bg);

    // apply watershed algo
    markers = Watershed::apply_watershed(thresholded, markers);
    markers.save_to_csv(output_dir + "/segments.csv");

    // outline segments, display and save
    Image watershed_outlined = Watershed::outline(orig, markers);
    //watershed_outlined.display();
    watershed_outlined.save(output_dir + "/4 segmented.png");
}


void bin(OPTS* opts, std::string output_dir)
{
    Image markers = IO::read_markers(opts->input_file);
    //IO::sparse_segments(markers.get_im(), output_dir + "/sparse_segments.txt");
    //IO::bin_coords(markers.get_im(), output_dir + "/coords.csv");
}

int main(int argc, char *argv[])
{
    // params
    OPTS* opts = (OPTS*) calloc(1, sizeof(OPTS));
    if (get_params(argc, argv, opts) != 0) exit(EXIT_FAILURE);

    // convert to C++ strings
    std::string output_dir = std::string(opts->output_dir);

    //all_watershed(opts, output_dir);
    
    bin(opts, output_dir);

    exit(EXIT_SUCCESS);
}



