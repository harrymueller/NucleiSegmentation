#include "params.h"

/*
    Sorts through parameters, checking that all req params are supplied, and adds to opts param
*/
int get_params(int argc, char *argv[], OPTS* opts)
{
    int option;
    opterr = 0;

    if (argc == 1) {
        usage();
        return 1;
    }

    while((option = getopt(argc, argv, CLI_OPTIONS)) != -1) {
        // -h -i: -o:
        switch (option) {
            case 'h': 
                usage();
                return 0;
            case 'i':
                opts->input_path = strdup(optarg);
                break;
            case 'o':
                opts->output_path = strdup(optarg);
                break;
            default: 
                printf("Unknown argument supplied.\n\n");
                usage();
                return 1;
        }
    }

    if (opts->input_path == NULL) {
        printf("Please supply an input file.\n");
        return 1;
    }
    if (opts->output_path == NULL) {
        printf("Please supply an output path.\n");
        return 1;
    }
    
    return 0;
}

void usage()
{
    printf("Usage: nuclei_segregation [OPTIONS] -i input\n");
    printf("Uses the watershed algorithm to segment nuclei of a ssDNA stained image.\n\n");

    // options
    printf("  -h\t Display this help menu\n");
    printf("  -i\t Input image\n");
    printf("  -o\t Output directory\n");
}