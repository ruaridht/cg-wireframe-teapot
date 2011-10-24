#include "headers.h"

using namespace std;

TriangleMesh trig;
char * strat;

void Wireframe::bresenLine(int x0, int y0, int x1, int y1) {
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

void Wireframe::bresenham() {
  this.bresenLine(v1[0], v1[1], v2[0], v2[1]);
  this.bresenLine(v1[0], v1[1], v3[0], v3[1]);
  this.bresenLine(v2[0], v2[1], v3[0], v3[1]);
}

void Wireframe::myDisplay()
{
	glClear(GL_COLOR_BUFFER_BIT); // Clear OpenGL Window

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
		
		Rotate(v1, rotateX, rotateY, rotateZ);
		Rotate(v2, rotateX, rotateY, rotateZ);
		Rotate(v3, rotateX, rotateY, rotateZ);

		glBegin(GL_POINTS);
			glVertex2i((int)v1[0],(int)v1[1]);
			glVertex2i((int)v2[0],(int)v2[1]);
			glVertex2i((int)v3[0],(int)v3[1]);
			
			if (strat=="bresenham") {
		    this.bresenham(v1,v2,v3);
		  } else if (strat=="midpoint") {
		    //this.midpoint(v1,v2,v3);
		  }
		glEnd();
	}
	glFlush();// Output everything
}

void Wireframe::loadMesh(char * filename) {
  ParseFile(filename,trig);
}

void Wireframe::draw(char * strategy) {
  strat = strategy;
  gluOrtho2D(-WIDTH/2, WIDTH/2, -(float)HEIGHT/2,  (float)HEIGHT/2);
  glutDisplayFunc(this.myDisplay);// Callback function
  glutMainLoop();// Display everything and wait
}

