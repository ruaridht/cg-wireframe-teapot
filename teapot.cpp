/*
 * @author Ruaridh Thomson
 * Simple teapot rasteriser (rasterizer)
 */

// teapot.cpp*
#include "headers.h"

Wireframe wire;

// main program
int main(int argc, char *argv[]) {
    Options options;
    char * filename;
    
    if (argc==1) {
      printf("usage: teapot [--quiet] [--verbose] [--help] <filename.obj> ...\n");
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
        else filename = argv[i];
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
    
    glutInit(&argc, argv);
		glutInitWindowSize(WIDTH, HEIGHT);
		glutCreateWindow("Teatpot Renderer");
		
		wire.loadMesh(filename);
		wire.draw();
		/*
		gluOrtho2D(-WIDTH/2, WIDTH/2, -(float)HEIGHT/2,  (float)HEIGHT/2);
		glutDisplayFunc(myDisplay);// Callback function
		glutMainLoop();// Display everything and wait
		*/
    return 0;
}


