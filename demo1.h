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

// Classes

class Vector3f;
class Triangle;
class TriangleMesh;

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

class Triangle {
friend class TriangleMesh;

	int _vertex[3];
public:

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
			
};

#endif //_rt_H
