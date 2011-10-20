/*
 * @author Ruaridh Thomson
 * Simple teapot rasteriser (rasterizer)
 */

// teapot.cpp*
#include "headers.h"

// main program
int main(int argc, char *argv[]) {
    Options options;
    vector<string> filenames;
    
    if (argc==1) {
      printf("usage: teapot [--verbose] [--help] <filename.obj> ...\n");
      return 0;
    }
    
    // Parse command-line arguments
    for (int i = 1; i < argc; ++i) {
        
        if (!strcmp(argv[i], "--quiet")) options.quiet = true;
        else if (!strcmp(argv[i], "--verbose")) options.verbose = true;
        else if (!strcmp(argv[i], "--help") || !strcmp(argv[i], "-h")) {
            printf("usage: teapot [--verbose] [--help] <filename.obj> ...\n");
            return 0;
        }
        else filenames.push_back(argv[i]);
    }

    // Print welcome
    if (!options.quiet) {
        printf("teapot version %s of %s at %s\n",
               TEAPOT_VERSION, __DATE__, __TIME__);
        printf("Copyright ©2011 Ruaridh Thomson.\n");
        printf("The source code to teapot is covered by the MIT license.\n");
        printf("See the file LICENSE.txt for the full license.\n");
        fflush(stdout);
    }
    
    return 0;
}


