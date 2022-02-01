/*##############################
    LIBRARIES
##############################*/
#include <stdio.h>
#include <string.h>
#include <getopt.h>



/*##############################
    PREPROCESSOR VARIABLES
##############################*/
#define CLI_OPTIONS "hi:o:"

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
// output functions
// params.c
int get_params(int, char **, OPTS *);
void usage();

// others
extern char *strdup(const char *);