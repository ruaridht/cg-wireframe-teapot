#ifndef TEAPOT_TRIANGLE_H
#define TEAPOT_TRIANGLE_H

#include "headers.h"

using namespace std;

class Triangle;
class TriangleMesh;

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

#endif // TEAPOT_TRIANGLE_H