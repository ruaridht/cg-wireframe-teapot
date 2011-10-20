#ifndef TEAPOT_HEADERS_H
#define TEAPOT_HEADERS_H

#ifdef __APPLE__
#  include <OpenGL/gl.h>
#  include <OpenGL/glu.h>
#  include <GLUT/glut.h>
#else
#  include <GL/gl.h>
#  include <GL/glu.h>
#  include <GL/glut.h>
#endif

#include <cmath>
#include <vector>
#include <iostream>
#include <fstream>
#include "math.h"

#include "teapot.h"
#include "parser.h"
#include "wireframe.h"
#include "triangle.h"

#endif // TEAPOT_HEADERS_H