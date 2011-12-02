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
#include <cstring>
#include <fstream>
#include "demo1.h"
#include "math.h"

using namespace std;

int nRows = 640;
int nCols = 480;
float rotation = 2.0f;
float rx = 0.0f;
float ry = 0.0f;
float rz = 0.0f;
float deltaAngle = 0.0f;
float angle = 0.0f;
int mButton = -1;
int prevX = 0;
int prevY = 0;

TriangleMesh trig;

// For shading
vector<Vector3f> face_normals;
vector<Vector3f> vertex_normals;

Vector3f lx; // Result colour
Vector3f Lx; // Light colour
Vector3f Ax; // Ambient colour
Vector3f Dx; // Diffuse colour
Vector3f Sx; // Specular colour
float    Ka; // Ambient coefficient (intensity)
float    Kd; // Diffuse coefficient
float    Ks; // Specular coefficient
float    Att; // Attenuation coefficient
int      small_n; // Shine/roughness
Vector3f big_N; // Surface normal
Vector3f light; // Light vector
Vector3f reflection; // Reflection vector
Vector3f view; // View vector

#define swap(a,b)           {a^=b; b^=a; a^=b;}
#define absolute(i,j,k)     ( (i-j)*(k = ( (i-j)<0 ? -1 : 1)))
#define max(a,b)            ( (a<b) ? b : a)

/**********************************/
/*           Vector3f             */
/**********************************/

Vector3f Vector3f::Cross(Vector3f& v)
{
  float _x = _item[1] * v[2] - _item[2] * v[1];
  float _y = _item[2] * v[0] - _item[0] * v[2];
  float _z = _item[0] * v[1] - _item[1] * v[0];
  /*
  const float _x = y * v.z - z * v.y;
  const float _y = z * v.x - x * v.z;
  const float _z = x * v.y - y * v.x;
  */
  return Vector3f(_x, _y, _z);
}

Vector3f& Vector3f::Normalize()
{
  const float Length = sqrtf(_item[0]*_item[0] + _item[1]*_item[1] + _item[2]*_item[2]);
  float x = _item[0]/Length;
  float y = _item[1]/Length;
  float z = _item[2]/Length;
  _item[0] = x;
  _item[1] = y;
  _item[2] = z;
  /*
  const float Length = sqrtf(x * x + y * y + z * z);
  x /= Length;
  y /= Length;
  z /= Length;
  */
  return *this;
}

/**********************************/
/*           Load OBJ             */
/**********************************/

void TriangleMesh::loadFile(char * filename) {
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
      
      //  x *= 1000; y *= 1000; z *= 1000;
      _v.push_back(Vector3f(x,y,z));
      av[0] += x; av[1] += y; av[2] += z;

      if (x > xmax) xmax = x;
      if (y > ymax) ymax = y;
      if (z > zmax) zmax = z;

      if (x < xmin) xmin = x;
      if (y < ymin) ymin = y;
      if (z < zmin) zmin = z;
    } else if (strcmp(header, "f") == 0) {
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

  for (int i = 0; i < _v.size(); i++) {
    for (int j = 0; j < 3; j++) _v[i][j] = (_v[i][j]-av[j])/range*400;  
  }
  cout << "trig " << _trig.size() << " vertices " << _v.size() << endl;
  f.close();
};

/**********************************/
/*        Draw the pixels         */
/**********************************/
void plot(int x, int y, int flag) {
  // non-zero flag indicates the pixels needing to be swapped
  if (flag)
    glVertex2i(y, x);
  else
    glVertex2i(x, y);
}

void putpixel(int x, int y, float alpha) {
  glColor4f(1.0,1.0,1.0,alpha);
  glVertex2i(x,y);
}

/**********************************/
/*         Rasterization          */
/**********************************/



/**********************************/
/*        Line algorithms         */
/**********************************/
void EFLA(int x, int y, int x2, int y2) {
  bool yLonger=false;
  int incrementVal, endVal;

  int shortLen=y2-y;
  int longLen=x2-x;
  if (abs(shortLen)>abs(longLen)) {
    int swap=shortLen;
    shortLen=longLen;
    longLen=swap;
    yLonger=true;
  }
  
  endVal=longLen;
  if (longLen<0) {
    incrementVal=-1;
    longLen=-longLen;
  } else incrementVal=1;

  double decInc;
  if (longLen==0) decInc=(double)shortLen;
  else decInc=((double)shortLen/(double)longLen);
  double j=0.0;
  if (yLonger) {
    for (int i=0;i!=endVal;i+=incrementVal) {
      glVertex2i(x+(int)j,y+i);
      j+=decInc;
    }
  } else {
    for (int i=0;i!=endVal;i+=incrementVal) {
      glVertex2i(x+i,y+(int)j);
      j+=decInc;
    }
  }
}

void WuLine(int a1, int b1, int a2, int b2) {
  int dx, dy, incr1, incr2, D, x, y, xend, c, pixels_left;
  int x1, y1;
  int sign_x, sign_y, step, reverse, i;

  dx = absolute(a2, a1, sign_x);
  dy = absolute(b2, b1, sign_y);
  /* decide increment sign by the slope sign */
  if (sign_x == sign_y)
    step = 1;
  else
    step = -1;

  if (dy > dx) {  /* chooses axis of greatest movement (make * dx) */
    swap(a1, b1);
    swap(a2, b2);
    swap(dx, dy);
    reverse = 1;
  } else
    reverse = 0;
  /* note error check for dx==0 should be included here */
  if (a1 > a2) {  /* start from the smaller coordinate */
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
  
  xend = (dx - 1) / 4;
  pixels_left = (dx - 1) % 4;  /* number of pixels left over at the end */
  plot(x, y, reverse);
  if ( pixels_left < 0 ) return ;  /* plot only one pixel for zero length vectors */
  plot(x1, y1, reverse);  /* plot first two points */
  incr2 = 4 * dy - 2 * dx;
  if (incr2 < 0) {  /* slope less than 1/2 */
    c = 2 * dy;
    incr1 = 2 * c;
    D = incr1 - dx;

    for (i = 0; i < xend; i++) {  /* plotting loop */
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
    }    /* end for */

    /* plot last pattern */
    if (pixels_left) {
      if (D < 0) {
        plot(++x, y, reverse);  /* pattern 1 */
        if (pixels_left > 1)
          plot(++x, y, reverse);
        if (pixels_left > 2)
          plot(--x1, y1, reverse);
      } else {
        if (D < c) {
          plot(++x, y, reverse);  /* pattern 2  */
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
    }  /* end if pixels_left */
  }
  /* end slope < 1/2 */
  else { /* slope greater than 1/2 */
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
    }  /* end for */
    /* plot last pattern */
    if (pixels_left) {
      if (D > 0) {
        plot(++x, y += step, reverse);  /* pattern 4 */
        if (pixels_left > 1)
          plot(++x, y += step, reverse);
        if (pixels_left > 2)
          plot(--x1, y1 -= step, reverse);
      } else {
        if (D < c) {
          plot(++x, y, reverse);  /* pattern 2  */
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

void AALine(int x1, int y1, int x2, int y2, int linewidth) {
  int dx = fabs(x1 - x2); 
  int dy = fabs(y1 - y2);
  int error, sign, tmp;
  float ipix;
  int step = linewidth;
  
  // Check whether the slope is more vertical or horizontal.
  // Vertical means we anti-alias the sides, horizontal is top n' bottom
  if (dx >= dy) {
    if (x1 > x2) {
      tmp = x1;
      x1 = x2;
      x2 = tmp;
      tmp = y1;
      y1 = y2;
      y2 = tmp;
    }
    error = dx / 2;
    if (y2 > y1)
      sign = step;
    else
      sign = -step;
    putpixel(x1,y1,1.0);
    
    while (x1 < x2) {
      if ((error -= dy) < 0) {
        y1 += sign;
        error += dx;
      }
      x1 += step;
      ipix = (float)error / dx;
      
      if (sign == step)
        ipix = 1.0 - ipix;
      
      putpixel(x1, y1, 1.0);
      putpixel(x1, y1 - step, (1.0 - ipix));
      putpixel(x1, y1 + step, ipix);
    }
    putpixel(x2, y2, 1.0);
  } else {
    if (y1 > y2) {
      tmp = x1;
      x1 = x2;
      x2 = tmp;
      tmp = y1;
      y1 = y2;
      y2 = tmp;
    }
    error = dy / 2;
  
    if (x2 > x1)
      sign = step;
    else
      sign = -step;
      
    putpixel(x1, y1, 1.0);
    while (y1 < y2) {
      if ((error -= dx) < 0) {
        x1 += sign;
        error += dy;
      }
      y1 += step;
      ipix = (float)error / dy;
  
      if (sign == step)
        ipix = 1 - ipix;
  
      putpixel(x1, y1, 1.0);
      putpixel(x1 - step, y1, (1 - ipix));
      putpixel(x1 + step, y1, ipix);
    }
    putpixel(x2, y2, 1.0);
  }
}

void GuptaSproull(int x0, int y0, int x1, int y1) {
  int dx, dy, sx=-1, sy=-1, err, e2;
  float numerator, denominator, capD, dUpper, dLower;
  
  dx = fabs(x1-x0);
  dy = fabs(y1-y0);
  
  if (x0 < x1) { sx = 1; }
  if (y0 < y1) { sy = 1; }
  
  err = dx-dy;
  denominator = 2.0*(sqrt(dx*dx + dy*dy));
  dUpper = fabs((2.0*dx - 2.0*err)/(2.0*denominator));
  dLower = fabs((2.0*dx + 2.0*err)/(2.0*denominator));
  
  do {
    putpixel(x0,y0,1.0);
    putpixel(x0,y0+1,dUpper);
    putpixel(x0,y0-1,dLower);
    
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
    
    denominator = 2.0*(sqrt(dx*dx + dy*dy));
    dUpper = (2.0*dx - 2.0*err)/(2.0*denominator);
    dLower = (2.0*dx + 2.0*err)/(2.0*denominator);
  } while (true);
}

void Bresenham(int x0, int y0, int x1, int y1) {
  int dx, dy, sx=-1, sy=-1, err, e2;
  
  dx = abs(x1-x0);
  dy = abs(y1-y0);
  
  if (x0 < x1) { sx = 1; }
  if (y0 < y1) { sy = 1; }
  
  err = dx-dy;
  
  while (true) {
    putpixel(x0,y0,1.0);
    
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

void DDALine(int x1, int y1, int x2, int y2) {
  float x, y;
  int dx = x2-x1, dy = y2-y1;
  int n = max(fabs(dx),fabs(dy));
  
  float dt = n, dxdt = dx/dt, dydt = dy/dt;
  x = x1;
  y = y1;
  while( n-- ) {
    putpixel(Round(x),Round(y),1.0);
    x += dxdt;
    y += dydt;
  }
}

/**********************************/
/*            Shear               */
/**********************************/

void Shear(Vector3f &v1, Vector3f &v2, Vector3f &v3, float sx, float sy) {
  v1[0] += (sx*v1[2]);
  v1[1] += (sy*v1[2]);
  
  v2[0] += (sx*v2[2]);
  v2[1] += (sy*v2[2]);
  
  v3[0] += (sx*v3[2]);
  v3[1] += (sy*v3[2]);
}

/**********************************/
/*          Translation           */
/**********************************/

void Translate(Vector3f &v1, Vector3f &v2, Vector3f &v3, int tx, int ty, int tz) {
  v1[0] += tx;
  v1[1] += ty;
  v1[2] += tz;
  
  v2[0] += tx;
  v2[1] += ty;
  v2[2] += tz;
  
  v3[0] += tx;
  v3[1] += ty;
  v3[2] += tz;
}

/**********************************/
/*           Rotation             */
/**********************************/

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

/**********************************/
/*             Scale              */
/**********************************/

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

/**********************************/
/*          Line Drawing          */
/**********************************/
void DoDDA(Vector3f v1, Vector3f v2, Vector3f v3) {
  DDALine(v1[0], v1[1], v2[0], v2[1]);
  DDALine(v1[0], v1[1], v3[0], v3[1]);
  DDALine(v2[0], v2[1], v3[0], v3[1]);
}

void DoBresenham(Vector3f v1, Vector3f v2, Vector3f v3) {
  Bresenham(v1[0], v1[1], v2[0], v2[1]);
  Bresenham(v1[0], v1[1], v3[0], v3[1]);
  Bresenham(v2[0], v2[1], v3[0], v3[1]);
}

void DoWuLine(Vector3f v1, Vector3f v2, Vector3f v3) {
  WuLine(v1[0], v1[1], v2[0], v2[1]);
  WuLine(v1[0], v1[1], v3[0], v3[1]);
  WuLine(v2[0], v2[1], v3[0], v3[1]);
}

void DoGS(Vector3f v1, Vector3f v2, Vector3f v3) {
  GuptaSproull(v1[0], v1[1], v2[0], v2[1]);
  GuptaSproull(v1[0], v1[1], v3[0], v3[1]);
  GuptaSproull(v2[0], v2[1], v3[0], v3[1]);
}

void DoAAL(Vector3f v1, Vector3f v2, Vector3f v3, int width) {
  AALine(v1[0], v1[1], v2[0], v2[1], width);
  AALine(v1[0], v1[1], v3[0], v3[1], width);
  AALine(v2[0], v2[1], v3[0], v3[1], width);
}

void DoEFLA(Vector3f v1, Vector3f v2, Vector3f v3) {
  EFLA(v1[0], v1[1], v2[0], v2[1]);
  EFLA(v1[0], v1[1], v3[0], v3[1]);
  EFLA(v2[0], v2[1], v3[0], v3[1]);
}

/**********************************/
/*            Display             */
/**********************************/

void myRefresh()
{
  glClear(GL_COLOR_BUFFER_BIT); // Clear OpenGL Window
  
  int trignum = trig.trigNum();
  Vector3f v1, v2, v3;
  
  // Translate by: 
  int tx = 0;
  int ty = 0;
  int tz = 0;
  
  // Shear by:
  float sx = 1.0;
  float sy = 1.0;
  
  float scaleCoeff = 1.0f;
  
  glColor4f(1.0,1.0,1.0,1.0); // Needed for every algorithm except anti-aliasing
  
  // for all the triangles, get the location of the vertices,
  // project them on the xy plane, and color the corresponding pixel by white
  for (int i = 0; i < trignum-1; i++) {
    trig.getTriangleVertices(i, v1,v2,v3);
    
    // Let's manipulate the vertices!
    Rotate(v1, v2, v3, rx, ry, rz);
    //Scale(v1, v2, v3, scaleCoeff);
    //Translate(v1, v2, v3, tx, ty, tz);
    //Shear(v1,v2,v3,sx,sy); // To keep things simple, only shear x,y
    
    glBegin(GL_POINTS);
      glVertex2i((int)v1[0],(int)v1[1]);
      glVertex2i((int)v2[0],(int)v2[1]);
      glVertex2i((int)v3[0],(int)v3[1]);
      
      // Let's draw some lines!
      //DoDDA(v1,v2,v3);           // Simple DDA lines
      DoBresenham(v1,v2,v3);     // Standard bresenham/midpoint lines
      //DoWuLine(v1,v2,v3);     // Wu patterns
      //DoGS(v1,v2,v3);            // Gupta-Sproul + anti-aliasing
      //DoAAL(v1,v2,v3,1); // Anti-aliased lines, hard to see
      //DoEFLA(v1,v2,v3);          // An Extremely Fast Line Algorithm (found on the web)
    glEnd();
  }
  
  glFlush(); // Output everything
}

void myKB(unsigned char key, int x, int y)
{
  glClear(GL_COLOR_BUFFER_BIT); // Clear OpenGL Window
  
  int trignum = trig.trigNum();
  Vector3f v1, v2, v3;

  if (key == 'w') {
    rz += rotation;
  } else if (key == 's') {
    rz -= rotation;
  } else if (key == 'a') {
    ry -= rotation;
  } else if (key == 'd') {
    ry += rotation;
  }
  
  glColor4f(1.0,1.0,1.0,1.0); // Needed for every algorithm except anti-aliasing
  
  for (int i = 0; i < trignum-1; i++) {
    trig.getTriangleVertices(i, v1,v2,v3);
    
    // Let's manipulate the vertices!
    Rotate(v1, v2, v3, rx, ry, rz);
    
    glBegin(GL_POINTS);
      glVertex2i((int)v1[0],(int)v1[1]);
      glVertex2i((int)v2[0],(int)v2[1]);
      glVertex2i((int)v3[0],(int)v3[1]);
      
      DoBresenham(v1,v2,v3);
    glEnd();
  }
  
  glFlush(); // Output everything
}

void myMouse(int button, int state, int x, int y) {
  // only start motion if the left button is pressed
	if (button == GLUT_LEFT_BUTTON) {
		// when the button is released
		if (state == GLUT_UP) {
			angle += deltaAngle;
			mButton = -1;
		}
		else  {// state = GLUT_DOWN
			mButton = GLUT_LEFT_BUTTON;
			prevX = x;
		}
	}
	if (button == GLUT_RIGHT_BUTTON) {
    // when the button is released
		if (state == GLUT_UP) {
			angle += deltaAngle;
			mButton = -1;
		}
		else  {// state = GLUT_DOWN
			mButton = GLUT_RIGHT_BUTTON;
			prevX = x;
		}
	}
}

void myMotion(int x, int y) {
  // this will only be true when the left button is down
	if (mButton == GLUT_LEFT_BUTTON) {
    if (prevX > x) {
      ry += (prevX-x);
    } else {
      ry -= (x-prevX);
    }
    if (prevY > y) {
      rz += (prevY-y);
    } else {
      rz -= (y-prevY);
    }
	}
	
	prevX = x;
	prevY = y;
	
	myRefresh();
}

/**********************************/
/*          For Marker            */
/**********************************/

void myDisplay()
{
  glClear(GL_COLOR_BUFFER_BIT); // Clear OpenGL Window
  
  int trignum = trig.trigNum();
  Vector3f v1, v2, v3;
  
  // Translate by:
  // Change these to change the translation on each axis.
  int tx = 100;
  int ty = 100;
  int tz = 0;
  
  // Shear by:
  // Change these values to affect the shear amount (note: keep between 0.0 and 2.0)
  float sx = 1.0;  // Shear amount on x-axis
  float sy = 1.0;  // shear amount on y-axis
  
  float  scaleCoeff = 3.0f; // Scale the object by this
  
  glColor4f(1.0,1.0,1.0,1.0); // Needed for every algorithm except anti-aliasing
  
  for (int i = 0; i < trignum-1; i++) {
    trig.getTriangleVertices(i, v1,v2,v3);
    
    // -------MANIPULATE----------
    // -------Uncomment to use
    //Rotate(v1, v2, v3, rx, ry, rz);
    //Scale(v1, v2, v3, scaleCoeff);
    //Translate(v1, v2, v3, tx, ty, tz);
    //Shear(v1,v2,v3,sx,sy); // To keep things simple, only shear x,y
    
    glBegin(GL_POINTS);
      glVertex2i((int)v1[0],(int)v1[1]);
      glVertex2i((int)v2[0],(int)v2[1]);
      glVertex2i((int)v3[0],(int)v3[1]);
      
      // ---------LINE DRAWING-----------
      // -------Uncomment to use
      //DoDDA(v1,v2,v3);        // Simple DDA lines
      DoBresenham(v1,v2,v3);    // Standard bresenham/midpoint lines
      //DoWuLine(v1,v2,v3);     // Wu patterns
      //DoGS(v1,v2,v3);         // Gupta-Sproul + anti-aliasing
      //DoAAL(v1,v2,v3,1);      // Anti-aliased lines, hard to see
      //DoEFLA(v1,v2,v3);       // An Extremely Fast Line Algorithm (found on the web)
    glEnd();
  }
  rx += 0.5;
  ry += 0.5;
  rz += 0.5;
  
  glFlush(); // Output everything
  
  /*
  GLdouble mdl[16];
  double camera_org[3];
  glGetDoublev(GL_MODELVIEW_MATRIX, mdl);
  camera_org[0] = -(mdl[0] * mdl[12] + mdl[1] * mdl[13] + mdl[2] * mdl[14]);
  camera_org[1] = -(mdl[4] * mdl[12] + mdl[5] * mdl[13] + mdl[6] * mdl[14]);
  camera_org[2] = -(mdl[8] * mdl[12] + mdl[9] * mdl[13] + mdl[10] * mdl[14]);
  printf("Camera x,y,z: %f,%f,%f \n", camera_org[0], camera_org[1], camera_org[2]);
  */
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
  
  glEnable(GL_BLEND);
  glBlendFunc(GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA);
  glClearColor(0.0,0.0,0.0,0.0); // Set the bg colour
  
  glutDisplayFunc(myDisplay);// Callback function
  glutIdleFunc(myDisplay); // Display this while idling
  
  /* ---- To Set Interactive ---------
   * Comment out glutIdleFunc() and uncomment glutKeyboardFunc(),
   * glutMouseFunc(), glutMotionFunc() to switch to an interactive 'mode'.
   */
  //glutKeyboardFunc(myKB);
  //glutMouseFunc(myMouse);
  //glutMotionFunc(myMotion);
  
  glutMainLoop();// Display everything and wait
}