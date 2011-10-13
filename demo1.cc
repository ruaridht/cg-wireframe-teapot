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

void Bresenham(int x1, int y1, int x2, int y2)
{
  int slope;
  int dx, dy, incE, incNE, d, x, y;
  
  // Reverse lines where x1 > x2
  if (x1 > x2) {
    Bresenham(x2,y2,x1,y1);
    return;
  }
  
  dx = x2-x1;
  dy = y2-y1;
  
  // Adjust y-increment for negatively sloped lines
  if (dy < 0) {
    slope = -1;
    dy = -dy;
  } else {
    slope = 1;
  }
  
  // Bresenham constants
  incE = 2*dy;
  incNE = (2*dy)-(2*dx);
  d = (2*dy)-dx;
  y = y1;
  
  // Blit
  for (x=x1; x<=x2; x++) {
    glVertex2i(x,y);
    
    if (d <= 0) {
      d += incE;
    } else {
      d += incNE;
      y += slope;
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
  //WritePixel(x,y);
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
    //WritePixel(x,y);
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



void myDisplay()
{
	glClear(GL_COLOR_BUFFER_BIT); // Clear OpenGL Window

	/** drawing a line for test **/

	/*** clear the Zbuffer here ****/

	int trignum = trig.trigNum();
	Vector3f v1, v2, v3;

	glColor3f(1,1,1);  // change the colour of the pixel


	//
	// for all the triangles, get the location of the vertices,
	// project them on the xy plane, and color the corresponding pixel by white
	//

	for (int i = 0 ; i < trignum-1; i++)  
	{
		/*** do the rasterization of the triangles here using glRecti ***/
		trig.getTriangleVertices(i, v1,v2,v3);

		//
		// colouring the pixels at the vertex location 
		// (just doing parallel projectiion to the xy plane. 
		// only use glBegin(GL_POINTS) for rendering the scene  
		//
		glBegin(GL_POINTS);
		//glBegin(GL_LINES);	
			glVertex2i((int)v1[0],(int)v1[1]);
			glVertex2i((int)v2[0],(int)v2[1]);
			glVertex2i((int)v3[0],(int)v3[1]);
			
			//MidpointLine((int)v1[0],(int)v1[1],(int)v2[0],(int)v2[1]);
			//MidpointLine((int)v1[0],(int)v1[1],(int)v3[0],(int)v3[1]);
			//MidpointLine((int)v2[0],(int)v2[1],(int)v3[0],(int)v3[1]);
			//MidpointLine((int)v2[0],(int)v2[1],(int)v1[0],(int)v1[1]);
			//MidpointLine((int)v3[0],(int)v3[1],(int)v1[0],(int)v1[1]);
			//MidpointLine((int)v3[0],(int)v3[1],(int)v2[0],(int)v2[1]);
			
			//Bresenham((int)v1[0],(int)v1[1],(int)v2[0],(int)v2[1]);
			//Bresenham((int)v1[0],(int)v1[1],(int)v3[0],(int)v3[1]);
			//Bresenham((int)v2[0],(int)v2[1],(int)v3[0],(int)v3[1]);
			//Bresenham((int)v2[0],(int)v2[1],(int)v1[0],(int)v1[1]);
			//Bresenham((int)v3[0],(int)v3[1],(int)v1[0],(int)v1[1]);
			//Bresenham((int)v3[0],(int)v3[1],(int)v2[0],(int)v2[1]);
			
			//XiaolinWu((int)v1[0],(int)v1[1],(int)v2[0],(int)v2[1]);
			//XiaolinWu((int)v1[0],(int)v1[1],(int)v3[0],(int)v3[1]);
			//XiaolinWu((int)v2[0],(int)v2[1],(int)v3[0],(int)v3[1]);
			//XiaolinWu((int)v2[0],(int)v2[1],(int)v1[0],(int)v1[1]);
			//XiaolinWu((int)v3[0],(int)v3[1],(int)v1[0],(int)v1[1]);
			//XiaolinWu((int)v3[0],(int)v3[1],(int)v2[0],(int)v2[1]);
		glEnd();
		
	}

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
