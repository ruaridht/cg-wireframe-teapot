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

using namespace std;

class Vector3f;

class Vector3f {

	float _item[3];

	public:

	float & operator [] (int i) {
		return _item[i];
    	}

	Vector3f(float _x, float _y, float _z) 
	{  _item[0] = _x ; _item[1] = _y ; _item[2] = _z; };

	Vector3f() {};


	Vector3f & operator = (Vector3f & obj) 
	{
		_item[0] = obj[0];
		_item[1] = obj[1];
		_item[2] = obj[2];
		
		return *this;
	};

	Vector3f & operator += (Vector3f & obj) 
	{
		_item[0] += obj[0];
		_item[1] += obj[1];
		_item[2] += obj[2];

		return *this;
	};
};

ostream & operator << (ostream & stream, Vector3f & obj) 
{
	stream << obj[0] << ' ' << obj[1] << ' ' << obj[2] << ' ';
};

#endif // TEAPOT_HEADERS_H