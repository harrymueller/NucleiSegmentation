#include "nuclei_segregation.hpp"
//#include "global_vars.h"

int main(int argc, char *argv[])
{
    // params
    OPTS* opts = (OPTS*) calloc(1, sizeof(OPTS));
    if (get_params(argc, argv, opts) != 0) exit(EXIT_FAILURE);

    // convert to C++ strings
    std::string output_dir = std::string(opts->output_dir);

    // make image
    Image i(opts->input_file);
    
    i.display();
    Image i2 = i.duplicate();
    i2.save(output_dir + std::string("/temp.png"));

    exit(EXIT_SUCCESS);
}



