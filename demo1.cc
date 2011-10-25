#ifdef __APPLE__
#  include <OpenGL/gl.h>
#  include <OpenGL/glu.h>
#  include <GLUT/glut.h>
#else
#  include <GL/gl.h>
#  include <GL/glu.h>
#  include <GL/glut.h>
#endif

#include <iostream>
#include <fstream>
#include "demo1.h"
#include "math.h"

using namespace std;

int nRows = 640;
int nCols = 480; 

TriangleMesh trig;

void AALine(int x0, int y0, int x1, int y1)
{
  int addr = (y0*640+x0)*4;
  int dx = x1-x0;
  int dy = y1-y0;
  int u, du, dv, uincr, vincr, v;
  /* By switching to (u,v), we combine all eight octants */
  if (abs(dx) > abs(dy))
  {
  	/* Note: If this were actual C, these integers would be lost
  	 * at the closing brace.  That's not what I mean to do.  Do what
  	 * I mean. */
  	du = abs(dx);
    dv = abs(dy);
  	u = x1;
  	v = y1;
    uincr = 4;
  	vincr = 640*4;
  	if (dx < 0) uincr = -uincr;
  	if (dy < 0) vincr = -vincr;
  } else {
  	du = abs(dy);
  	dv = abs(dx);
  	u = y1;
  	v = x1;
  	uincr = 640*4;
  	vincr = 4;
  	if (dy < 0) uincr = -uincr;
  	if (dx < 0) vincr = -vincr;
  }
  int uend = u + 2 * du;
  int d = (2 * dv) - du;                          /* Initial value as in Bresenham's */
  int incrS = 2 * dv;	                            /* Δd for straight increments */
  int incrD = 2 * (dv - du);	                    /* Δd for diagonal increments */
  int twovdu = 0;	                                /* Numerator of distance; starts at 0 */
  double invD = 1.0 / (2.0*sqrt(du*du + dv*dv));  /* Precomputed inverse denominator */
  double invD2du = 2.0 * (du*invD);               /* Precomputed constant */
  
  while (u < uend) {
  	/* Note: this pseudocode doesn't ensure that the address is
  	 * valid, or that it even represents a pixel on the same side of
  	 * the screen as the adjacent pixel */
  	//DrawPixelD(addr, twovdu*invD);
  	//DrawPixelD(addr + vincr, invD2du - twovdu*invD);
  	//DrawPixelD(addr - vincr, invD2du + twovdu*invD);
    
    glVertex2i(addr, twovdu*invD);
  	glVertex2i(addr + vincr, invD2du - twovdu*invD);
  	glVertex2i(addr - vincr, invD2du + twovdu*invD);
    
  	if (d < 0) {
	    /* choose straight (u direction) */
	    twovdu = d + du;
	    d = d + incrS;
  	} else {
	    /* choose diagonal (u+v direction) */
	    twovdu = d - du;
	    d = d + incrD;
	    v = v+1;
	    addr = addr + vincr;
  	}
  	u = u+1;
  	addr = addr+uincr;
  } //while (u < uend);
}

void XwPlot(int x, int y, float c)
{
  glColor4f(1,1,1,c);
  glVertex2i(x,y);
  glColor4f(1,1,1,1);
}

int XwRound(int n)
{
  return (float)floor(n + 0.5);
}

float XwFpart(float n)
{
  return (float)(n - floor(n));
}

float XwRfpart(float n)
{
  return (float)(1.0 - XwFpart(n));
}

float XwIpart(float n)
{
  return (float)floor(n);
}

void XiaolinWu(int x1, int y1, int x2, int y2)
{
  int dx = x2-x1;
  int dy = y2-y1;
  int tempx1, tempx2, tempdx;
  float gradient;
  float xend, yend, xgap, xpx11, ypx11, xpx12, ypx12, intery;
  
  if (fabs(dx) < fabs(dy)) {
    tempx1 = x1;
    tempx2 = x2;
    x1 = y1;
    x2 = y2;
    y1 = tempx1;
    y2 = tempx2;
    tempdx = dx;
    dx = dy;
    dy = tempdx;
  }
  if (x2 < x1) {
    tempx1 = x1;
    tempx2 = y1;
    x1 = x2;
    x2 = tempx1;
    y1 = y2;
    y2 = tempx2;
  }
  gradient = (float)dy/dx;
  
  // handle first endpoint
  xend = (float)XwRound(x2);
  yend = (float)(y1 + gradient*(xend - x1));
  xgap = XwRfpart((float)(x1 + 0.5));
  xpx11 = xend;
  ypx11 = XwIpart((float)yend);
  XwPlot((int)xpx11,(int)ypx11, XwRfpart(yend)*xgap);
  XwPlot((int)xpx11,(int)ypx11+1,XwFpart(yend)*xgap);
  intery = yend + gradient;
  
  // handle second endpoint
  xend = XwRound(x2);
  yend = y2 + gradient * (xend-x2);
  xgap = XwFpart(x2+0.5);
  xpx12 = xend;
  ypx12 = XwIpart(yend);
  XwPlot(xpx12, ypx12, XwRfpart(yend)*xgap);
  XwPlot(xpx12, ypx12+1, XwFpart(yend)*xgap);
  
  // main loop
  for (float x=(xpx11+1); x<=(xpx12-1); x++) {
    XwPlot((int)x,(int)XwIpart(intery),XwRfpart(intery));
    XwPlot((int)x,(int)XwIpart(intery)+1,XwFpart(intery));
    intery += gradient;
  }
}

void Bresenham(int x0, int y0, int x1, int y1)
{
  int dx, dy, sx=-1, sy=-1, err, e2;
  
  dx = abs(x1-x0);
  dy = abs(y1-y0);
  
  if (x0 < x1) { sx = 1; }
  if (y0 < y1) { sy = 1; }
  
  err = dx-dy;
  
  while (true) {
    glVertex2i(x0,y0);
    
    if (x0==x1 && y0==y1) { break; }
    e2 = 2*err;
    
    if (e2 > -dy) {
      err = err - dy;
      x0 = x0 + sx;
    }
    if (e2 < dx) {
      err = err + dx;
      y0 = y0 + sy;
    }
  }
}

void MidpointLine(int x1, int y1, int x2, int y2)
{
  int dx = x2-x1;
  int dy = y2-y1;
  int d = 2*dy-dx;
  int increE = 2*dy;
  int incrNE = 2*(dy-dx);
  int x = x1;
  int y = y1;
  
  glVertex2i(x,y);
  
  while (x < x2) {
    if (d <= 0) {
      d += increE;
      x++;
    } else {
      d += incrNE;
      x++;
      y++;
    }
    
    glVertex2i(x,y);
  }
}

void TriangleMesh::loadFile(char * filename)
{
	ifstream f(filename);


	if (f == NULL) {
		cerr << "failed reading polygon data file " << filename << endl;
		exit(1);
	}

	char buf[1024];
	char header[100];
	float x,y,z;
	float xmax,ymax,zmax,xmin,ymin,zmin;
	int v1, v2, v3, n1, n2, n3;

	xmax =-10000; ymax =-10000; zmax =-10000;
	xmin =10000; ymin =10000; zmin =10000;
	Vector3f av;
	av[0] = av[1] = av[2] = 0.f;

	while (!f.eof()) {
		    f.getline(buf, sizeof(buf));
		    sscanf(buf, "%s", header);  

		    if (strcmp(header, "v") == 0) {
			sscanf(buf, "%s %f %f %f", header, &x, &y, &z);  

		//	x *= 1000; y *= 1000; z *= 1000;

			_v.push_back(Vector3f(x,y,z));


			av[0] += x; av[1] += y; av[2] += z;

			if (x > xmax) xmax = x;
			if (y > ymax) ymax = y;
			if (z > zmax) zmax = z;

			if (x < xmin) xmin = x;
			if (y < ymin) ymin = y;
			if (z < zmin) zmin = z;
		    }
		    else if (strcmp(header, "f") == 0) {
			sscanf(buf, "%s %d %d %d", header, &v1, &v2, &v3);
			
			Triangle trig(v1-1, v2-1, v3-1);
			_trig.push_back(trig);

		    }
 	}

	_xmin = xmin; _ymin = ymin; _zmin = zmin;
	_xmax = xmax; _ymax = ymax; _zmax = zmax;
  
	float range; 
	if (xmax-xmin > ymax-ymin) range = xmax-xmin;
	else range = ymax-ymin;

	for (int j = 0; j < 3; j++) av[j] /= _v.size();

	for (int i = 0; i < _v.size(); i++) 
	{
		for (int j = 0; j < 3; j++) _v[i][j] = (_v[i][j]-av[j])/range*400;  
	}
	cout << "trig " << _trig.size() << " vertices " << _v.size() << endl;
	f.close();
};


void RotateX(Vector3f &v, float angle) {
  angle = Radians(angle);
  float s = sinf(angle);
  float c = cosf(angle);
  float y, z;
  y = v[1]*c - v[2]*s;
  z = v[1]*s + v[2]*c;
}

void RotateY(Vector3f &v, float angle) {
  angle = Radians(angle);
  float s = sinf(angle);
  float c = cosf(angle);
  float x, z;
  x = v[0]*c + v[2]*s;
  z = -v[0]*s + v[2]*c;
  v[0] = x;
  v[2] = z;
}

void RotateZ(Vector3f &v, float angle) {
  angle = Radians(angle);
  float s = sinf(angle);
  float c = cosf(angle);
  float x, y;
  x = v[0]*c - v[1]*s;
  y = v[0]*s + v[1]*c;
  v[0] = x;
  v[1] = y;
}

void Rotate(Vector3f &v1, Vector3f &v2, Vector3f &v3, float ax, float ay, float az) {
  RotateX(v1,ax);
  RotateX(v2,ax);
  RotateX(v3,ax);
  
  RotateY(v1,ay);
  RotateY(v2,ay);
  RotateY(v3,ay);
  
  RotateZ(v1,az);
  RotateZ(v2,az);
  RotateZ(v3,az);
}

void DoBresenham(Vector3f v1, Vector3f v2, Vector3f v3) {
  Bresenham(v1[0], v1[1], v2[0], v2[1]);
  Bresenham(v1[0], v1[1], v3[0], v3[1]);
  Bresenham(v2[0], v2[1], v3[0], v3[1]);
}

void DoMidpoint(Vector3f v1, Vector3f v2, Vector3f v3) {
  MidpointLine(v1[0], v1[1], v2[0], v2[1]);
  MidpointLine(v1[0], v1[1], v3[0], v3[1]);
  MidpointLine(v2[0], v2[1], v3[0], v3[1]);
}

void DoXiaolinWu(Vector3f v1, Vector3f v2, Vector3f v3) {
  XiaolinWu(v1[0], v1[1], v2[0], v2[1]);
  XiaolinWu(v1[0], v1[1], v3[0], v3[1]);
  XiaolinWu(v2[0], v2[1], v3[0], v3[1]);
}

void DoAALine(Vector3f v1, Vector3f v2, Vector3f v3) {
  AALine(v1[0], v1[1], v2[0], v2[1]);
  AALine(v1[0], v1[1], v3[0], v3[1]);
  AALine(v2[0], v2[1], v3[0], v3[1]);
}

void myDisplay()
{
	glClear(GL_COLOR_BUFFER_BIT); // Clear OpenGL Window

	int trignum = trig.trigNum();
	Vector3f v1, v2, v3;
	
	float rotation = 0.0f;
	float ax = 0.4f;
	float ay = 0.6f;
	float az = 0.5f;
  
	glColor4f(1,1,1,1.0);  // change the colour of the pixel and set alpha

	// for all the triangles, get the location of the vertices,
	// project them on the xy plane, and color the corresponding pixel by white
	//glBegin(GL_POINTS);
	while (true) {
	for (int i = 0; i < trignum-1; i++) {
		trig.getTriangleVertices(i, v1,v2,v3);
    
    Rotate(v1, v2, v3, rotation*ax, rotation*ay, rotation*az);
    
		glBegin(GL_POINTS);
			glVertex2i((int)v1[0],(int)v1[1]);
			glVertex2i((int)v2[0],(int)v2[1]);
			glVertex2i((int)v3[0],(int)v3[1]);
			
			//DoMidpoint(v1,v2,v3);
			DoBresenham(v1,v2,v3);
			//DoXiaolinWu(v1,v2,v3);
			//DoAALine(v1,v2,v3);
		glEnd();
	}
	rotation += 1.0;
	}
	//glEnd();
	glFlush();// Output everything
}


int main(int argc, char **argv)
{
	if (argc >  1)  {
		trig.loadFile(argv[1]);
	}
	else {
		cerr << argv[0] << " <filename> " << endl;
		exit(1);
	}

	int width, height;
	glutInit(&argc, argv);
	glutInitWindowSize(nRows, nCols);
	glutCreateWindow("SimpleExample");
	gluOrtho2D(-nRows/2, nRows/2, -(float)nCols/2,  (float)nCols/2);
	glutDisplayFunc(myDisplay);// Callback function
	glutMainLoop();// Display everything and wait
}
