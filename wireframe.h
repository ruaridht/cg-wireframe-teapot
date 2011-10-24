#ifndef TEAPOT_WIREFRAME_H
#define TEAPOT_WIREFRAME_H

#include "headers.h"

using namespace std;

class Wireframe;

class Wireframe {
  
  void loadMesh(char * filename);
  
  void draw(char * strategy);
  void myDisplay();
  
  void bresenham();
};

#endif // TEAPOT_WIREFRAME_H