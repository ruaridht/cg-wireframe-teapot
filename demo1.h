#ifndef _rt_H
#define _rt_H

#include <cmath>
#include <vector>

#ifdef __APPLE__
#  include <OpenGL/gl.h>
#  include <OpenGL/glu.h>
#  include <GLUT/glut.h>
#else
#  include <GL/gl.h>
#  include <GL/glu.h>
#  include <GL/glut.h>
#endif

#include "Color.h"

using namespace std;

#ifdef M_PI
#undef M_PI
#endif
#define M_PI 3.14159265358979323846f

// Global inline functions
inline float Radians(float deg) {
    return ((float)M_PI/180.f) * deg;
}

inline float Degrees(float rad) {
    return (180.f/(float)M_PI) * rad;
}

inline int Round(float n) {
  return (int)floor(n + 0.5);
}

inline float FracPart(float n) {
  return (float)(n - floor(n));
}

inline float RFracPart(float n) {
  return (float)(1.0 - FracPart(n));
}

inline int IntPart(float n) {
  return (int)floor(n);
}

inline int iround(double x)
{
  return (int)floor(x + 0.5);
}

// Classes

class Vector3f;
class Triangle;
class TriangleMesh;
class Edge;
class Span;

class Edge
{
	public:
		Color Color1, Color2;
		int X1, Y1, X2, Y2;

		Edge(const Color &color1, int x1, int y1, const Color &color2, int x2, int y2)
		{
		  if(y1 < y2) {
    		Color1 = color1;
    		X1 = x1;
    		Y1 = y1;
    		Color2 = color2;
    		X2 = x2;
    		Y2 = y2;
    	} else {
    		Color1 = color2;
    		X1 = x2;
    		Y1 = y2;
    		Color2 = color1;
    		X2 = x1;
    		Y2 = y1;
    	}
		}
};

class Span
{
	public:
		Color Color1, Color2;
		int X1, X2;

		Span(const Color &color1, int x1, const Color &color2, int x2)
		{
		  if(x1 < x2) {
    		Color1 = color1;
    		X1 = x1;
    		Color2 = color2;
    		X2 = x2;
    	} else {
    		Color1 = color2;
    		X1 = x2;
    		Color2 = color1;
    		X2 = x1;
    	}
		}
};

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
	
	Vector3f & operator *= (float t)
	{
    _item[0] *= t;
    _item[1] *= t;
    _item[2] *= t;
    
    return *this;
	};
};

ostream & operator << (ostream & stream, Vector3f & obj) 
{
	stream << obj[0] << ' ' << obj[1] << ' ' << obj[2] << ' ';
};

class Triangle {
friend class TriangleMesh;

	//int _vertex[3];
public:
  int _vertex[3];
  
	Triangle(int v1, int v2, int v3) 
	{
		_vertex[0] = v1;  _vertex[1] = v2;  _vertex[2] = v3;  
		
	};
};

class TriangleMesh 
{
	vector <Vector3f> _v;
	vector <Triangle> _trig;
	float _xmax, _xmin, _ymax, _ymin, _zmin, _zmax;

public: 
	TriangleMesh(char * filename) { loadFile(filename) ;};
	TriangleMesh() {};
	void loadFile(char * filename);

	int trigNum() { return _trig.size() ;};
	int vNum() { return _v.size();};
	Vector3f v(int i) { return _v[i];};

	void getTriangleVertices(int i, Vector3f & v1, Vector3f & v2, Vector3f & v3)
	{
		v1 = _v[_trig[i]._vertex[0]]; 
		v2 = _v[_trig[i]._vertex[1]]; 
		v3 = _v[_trig[i]._vertex[2]]; 
	}
	
	Triangle getTriangle(int i)
	{
    return _trig[i];
	}
  
  int getTriangleVertexIndex1(int i)
  {
    int i1 = _trig[i]._vertex[0];
    //int i2 = _trig[i]._vertex[1];
    //int i3 = _trig[i]._vertex[2];
    //int inds[] = { i1, i2, i3 };
    return i1;
  }
  
  int getTriangleVertexIndex2(int i)
  {
    int i2 = _trig[i]._vertex[1];
    return i2;
  }
  
  int getTriangleVertexIndex3(int i)
  {
    int i3 = _trig[i]._vertex[2];
    return i3;
  }
};

// blah
inline void printVector3f(Vector3f v)
{
	printf("%f,", v[0] );
	printf("%f,", v[1] );
	printf("%f  ", v[2] );
}

#endif //_rt_H