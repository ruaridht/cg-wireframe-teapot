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
#include <cstdlib>
#include <fstream>
#include "demo1.h"
#include "math.h"

using namespace std;

int nRows = 640;
int nCols = 480; 
float rotation = 100.0f;
float rgb[] = { 1.0, 1.0, 1.0 };
bool rgbup = false;

TriangleMesh trig;

#define swap(a,b)           {a^=b; b^=a; a^=b;}
#define absolute(i,j,k)     ( (i-j)*(k = ( (i-j)<0 ? -1 : 1)))

/* non-zero flag indicates the pixels needing swap back. */
void plot(int x, int y, int flag) {
	if (flag)
		glVertex2i(y, x);
	else
		glVertex2i(x, y);
}

void symwuline(int a1, int b1, int a2, int b2) {
	int           dx, dy, incr1, incr2, D, x, y, xend, c, pixels_left;
	int           x1, y1;
	int           sign_x, sign_y, step, reverse, i;

	dx = absolute(a2, a1, sign_x);
	dy = absolute(b2, b1, sign_y);
	/* decide increment sign by the slope sign */
	if (sign_x == sign_y)
		step = 1;
	else
		step = -1;

	if (dy > dx) {		/* chooses axis of greatest movement (make
				 		 * dx) */
		swap(a1, b1);
		swap(a2, b2);
		swap(dx, dy);
		reverse = 1;
	} else
		reverse = 0;
	/* note error check for dx==0 should be included here */
	if (a1 > a2) {		/* start from the smaller coordinate */
		x = a2;
		y = b2;
		x1 = a1;
		y1 = b1;
	} else {
		x = a1;
		y = b1;
		x1 = a2;
		y1 = b2;
	}


	/* Note dx=n implies 0 - n or (dx+1) pixels to be set */
	/* Go round loop dx/4 times then plot last 0,1,2 or 3 pixels */
	/* In fact (dx-1)/4 as 2 pixels are already plotted */
	xend = (dx - 1) / 4;
	pixels_left = (dx - 1) % 4;	/* number of pixels left over at the
					 			 * end */
	plot(x, y, reverse);
	if ( pixels_left < 0 ) return ;	/* plot only one pixel for zero
							* length vectors */
	plot(x1, y1, reverse);	/* plot first two points */
	incr2 = 4 * dy - 2 * dx;
	if (incr2 < 0) {	/* slope less than 1/2 */
		c = 2 * dy;
		incr1 = 2 * c;
		D = incr1 - dx;

		for (i = 0; i < xend; i++) {	/* plotting loop */
			++x;
			--x1;
			if (D < 0) {
                  			/* pattern 1 forwards */
				plot(x, y, reverse);
				plot(++x, y, reverse);
                                /* pattern 1 backwards */
				plot(x1, y1, reverse);
				plot(--x1, y1, reverse);
				D += incr1;
			} else {
				if (D < c) {
					/* pattern 2 forwards */
					plot(x, y, reverse);
					plot(++x, y += step, reverse);
					/* pattern 2 backwards */
					plot(x1, y1, reverse);
					plot(--x1, y1 -= step, reverse);	
				} else {
				        /* pattern 3 forwards */
					plot(x, y += step, reverse);
					plot(++x, y, reverse);
					/* pattern 3 backwards */
					plot(x1, y1 -= step, reverse);
					plot(--x1, y1, reverse);
				}
				D += incr2;
			}
		}		/* end for */

		/* plot last pattern */
		if (pixels_left) {
			if (D < 0) {
				plot(++x, y, reverse);	/* pattern 1 */
				if (pixels_left > 1)
					plot(++x, y, reverse);
				if (pixels_left > 2)
					plot(--x1, y1, reverse);
			} else {
				if (D < c) {
					plot(++x, y, reverse);	/* pattern 2  */
					if (pixels_left > 1)
						plot(++x, y += step, reverse);
					if (pixels_left > 2)
						plot(--x1, y1, reverse);
				} else {
				  /* pattern 3 */
					plot(++x, y += step, reverse);
					if (pixels_left > 1)
						plot(++x, y, reverse);
					if (pixels_left > 2)
						plot(--x1, y1 -= step, reverse);
				}
			}
		}		/* end if pixels_left */
	}
	/* end slope < 1/2 */
	else {			/* slope greater than 1/2 */
		c = 2 * (dy - dx);
		incr1 = 2 * c;
		D = incr1 + dx;
		for (i = 0; i < xend; i++) {
			++x;
			--x1;
			if (D > 0) {
			  /* pattern 4 forwards */
				plot(x, y += step, reverse);
				plot(++x, y += step, reverse);
			  /* pattern 4 backwards */
				plot(x1, y1 -= step, reverse);
				plot(--x1, y1 -= step, reverse);
				D += incr1;
			} else {
				if (D < c) {
				  /* pattern 2 forwards */
					plot(x, y, reverse);
					plot(++x, y += step, reverse);

 				  /* pattern 2 backwards */
					plot(x1, y1, reverse);
					plot(--x1, y1 -= step, reverse);
				} else {
				  /* pattern 3 forwards */
					plot(x, y += step, reverse);
					plot(++x, y, reverse);
				  /* pattern 3 backwards */
					plot(x1, y1 -= step, reverse);
					plot(--x1, y1, reverse);
				}
				D += incr2;
			}
		}		/* end for */
		/* plot last pattern */
		if (pixels_left) {
			if (D > 0) {
				plot(++x, y += step, reverse);	/* pattern 4 */
				if (pixels_left > 1)
					plot(++x, y += step, reverse);
				if (pixels_left > 2)
					plot(--x1, y1 -= step, reverse);
			} else {
				if (D < c) {
					plot(++x, y, reverse);	/* pattern 2  */
					if (pixels_left > 1)
						plot(++x, y += step, reverse);
					if (pixels_left > 2)
						plot(--x1, y1, reverse);
				} else {
				  /* pattern 3 */
					plot(++x, y += step, reverse);
					if (pixels_left > 1)
						plot(++x, y, reverse);
					if (pixels_left > 2) {
						if (D > c) /* step 3 */
						   plot(--x1, y1 -= step, reverse);
						else /* step 2 */
							plot(--x1, y1, reverse);
                         		}
				}
			}
		}
	}
}

void XwPlot(int x, int y, float c)
{
  glColor4f(1.0,1.0,1.0,c);
  glVertex2i(x,y);
  //glColor4f(1,1,1,1);
}

int XwRound(float n)
{
  return (int)floor(n + 0.5);
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
  float dx = (float)x2-(float)x1;
  float dy = (float)y2-(float)y1;
  float gradient;
  float xend, yend, xgap, xpx11, ypx11, xpx12, ypx12, intery;
  /*
  if (fabs(dx) < fabs(dy)) {
    swap(x1,y1);
    swap(x2,y2);
    swap(dx,dy);
  }
  if (x2 < x1) {
    swap(x1,x2);
    swap(y1,y2);
  }
  */
  gradient = dy/dx;
  
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
  for (float i=(xpx11+1); i<=(xpx12-1); i++) {
    XwPlot((int)i,(int)XwIpart(intery),(float)XwRfpart(intery));
    XwPlot((int)i,(int)XwIpart(intery)+1,(float)XwFpart(intery));
    intery += gradient;
  }
}

void GuptaSproul(int x1, int x2, int y1, int y2) {
  int dx = x2-x1;
  int dy = y2-y1;
  int d = 2*dy-dx;
  int increE = 2*dy;
  int incrNE = 2*(dy-dx);
  int x = x1;
  int y = y1;
  
  float numerator, denominator, capD, dUpper, dLower;
  
  glVertex2i(x,y);
  
  while (x < x2) {
    if (d <= 0) {
      numerator = d + dx;
      d += increE;
      
      x++;
    } else {
      numerator = d - dx;
      d += incrNE;
      
      x++;
      y++;
    }
    
    denominator = 2*(sqrt(dx*dx + dy*dy));
    capD = (float)numerator / (float)denominator;
    
    dUpper = (2.0*dx - 2.0*numerator)/denominator;
    dLower = (2.0*dx + 2.0*numerator)/denominator;
    
    glColor4f(1.0,1.0,1.0,capD);
    glVertex2i(x,y);
    
    glColor4f(1.0,1.0,1.0,dUpper);
    glVertex2i(x,y+1);
    
    glColor4f(1.0,1.0,1.0,dLower);
    glVertex2i(x,y-1);
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

void Scale(Vector3f &v1, Vector3f &v2, Vector3f &v3, float sc) {
  float x,y,z;
  x = v1[0]*sc;
  y = v1[1]*sc;
  z = v1[2]*sc;
  v1[0] = x;
  v1[1] = y;
  v1[2] = z;
  
  x = v2[0]*sc;
  y = v2[1]*sc;
  z = v2[2]*sc;
  v2[0] = x;
  v2[1] = y;
  v2[2] = z;
  
  x = v3[0]*sc;
  y = v3[1]*sc;
  z = v3[2]*sc;
  v3[0] = x;
  v3[1] = y;
  v3[2] = z;
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

void DoSymwuline(Vector3f v1, Vector3f v2, Vector3f v3) {
  symwuline(v1[0], v1[1], v2[0], v2[1]);
  symwuline(v1[0], v1[1], v3[0], v3[1]);
  symwuline(v2[0], v2[1], v3[0], v3[1]);
}

void DoGS(Vector3f v1, Vector3f v2, Vector3f v3) {
  GuptaSproul(v1[0], v1[1], v2[0], v2[1]);
  GuptaSproul(v1[0], v1[1], v3[0], v3[1]);
  GuptaSproul(v2[0], v2[1], v3[0], v3[1]);
}

void myDisplay()
{
  glClearColor(0.0,0.0,0.0,0.0); // similarly
	glClear(GL_COLOR_BUFFER_BIT); // Clear OpenGL Window

	int trignum = trig.trigNum();
	Vector3f v1, v2, v3;
	
	float ax = 0.6f;
	float ay = 0.6f;
	float az = 0.6f;
	double alpha = 1.0;
	
	float scaleCoeff = 3.0f;
  
  if (rgb[0] >= 1.0) {
    rgbup = false;
  } else if (rgb[0] <= 0.0) {
    rgbup = true;
  }
  
	glColor4f(rgb[0],rgb[1],rgb[2],alpha);  // change the colour of the pixel and set alpha
  //glColor4f(1.0,1.0,1.0,alpha);
  
	// for all the triangles, get the location of the vertices,
	// project them on the xy plane, and color the corresponding pixel by white
	for (int i = 0; i < trignum-1; i++) {
		trig.getTriangleVertices(i, v1,v2,v3);
    
    //Rotate(v1, v2, v3, rotation*ax, rotation*ay, rotation*az);
    //Scale(v1, v2, v3, scaleCoeff);
    //Translate();
    
		glBegin(GL_POINTS);
		  //glColor4f(1.0,1.0,1.0,alpha);
			//glColor4f(rgb[0],rgb[1],rgb[2],alpha);
			
			glVertex2i((int)v1[0],(int)v1[1]);
			glVertex2i((int)v2[0],(int)v2[1]);
			glVertex2i((int)v3[0],(int)v3[1]);
			
			//DoMidpoint(v1,v2,v3);
			//DoBresenham(v1,v2,v3);
			//DoXiaolinWu(v1,v2,v3);
			//DoSymwuline(v1,v2,v3);
			DoGS(v1,v2,v3);
		glEnd();
	}
	rotation += 0.2;
	/*
	if (rgbup) {
  	rgb[0] += 0.004;
  	rgb[1] += 0.016;
  	rgb[2] += 0.010;
  } else {
    rgb[0] -= 0.010;
  	rgb[1] -= 0.016;
  	rgb[2] -= 0.004;
  }
	*/
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
	
	//glEnable(GL_ALPHA_TEST);
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA);
	//glAlphaFunc(GL_GREATER,0.1f);
	
	glutDisplayFunc(myDisplay);// Callback function
	glutIdleFunc(myDisplay); // Idling
	glutMainLoop();// Display everything and wait
}
