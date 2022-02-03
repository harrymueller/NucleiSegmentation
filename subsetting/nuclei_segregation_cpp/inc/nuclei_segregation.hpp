/*##############################
    LIBRARIES
##############################*/
#include <stdlib.h>

/*##############################
    PREPROCESSOR VARIABLES
##############################*/

/*##############################
    STRUCTS
##############################*/
typedef struct CliArgs {
    char *input_file;
    char *output_dir; 
} OPTS;

/*##############################
    FUNCTIONS
##############################*/
// params.c
extern "C" int get_params(int, char **, OPTS *);

// image
#include "image.hpp"