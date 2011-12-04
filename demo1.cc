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
#include "Color.h"

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
float z_buffer[640][480] = {0};
//vector< vector<float> > z_buffer;

Color lx; // Result colour
Color Lx; // Light colour
Color Ax; // Ambient colour
Color Dx; // Diffuse colour
Color Sx; // Specular colour
float    Ka; // Ambient coefficient (intensity)
float    Kd; // Diffuse coefficient
float    Ks; // Specular coefficient
float    Att; // Attenuation coefficient
int      small_n; // Shine/roughness (specular hardness)
Vector3f light; // Light vector
Vector3f view; // View vector

Color ambient;

#define swap(a,b)       {a^=b; b^=a; a^=b;}
#define absolute(i,j,k) ( (i-j)*(k = ( (i-j)<0 ? -1 : 1)))
#define max(a,b)        ( (a<b) ? b : a)
#define max3(a,b,c)     ( max( (a), max((b),(c)) ) )
#define min(a,b)        ( (a>b) ? b : a)
#define min3(a,b,c)     ( min( (a), min((b),©) ) )

// Vector3f stuff
Vector3f Cross(Vector3f& v1, Vector3f& v2)
{
  float x = v2[1] * v1[2] - v2[2] * v1[1];
  float y = v2[2] * v1[0] - v2[0] * v1[2];
  float z = v2[0] * v1[1] - v2[1] * v1[0];

  return Vector3f(x, y, z);	
}

float dot(Vector3f& v1, Vector3f& v2)
{
  float dot = (v1[0]*v2[0]) + (v1[1]*v2[1]) + (v1[2]*v2[2]);
  return dot;
}

void normalise(Vector3f& v)
{
  const float Length = sqrtf(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
  float x = v[0]/Length;
  float y = v[1]/Length;
  float z = v[2]/Length;
  v[0] = x;
  v[1] = y;
  v[2] = z;
}

float length(Vector3f& v)
{
  float len = sqrtf( (v[0]*v[0]) + (v[1]*v[1]) + (v[2]*v[2]) );
  return len;
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

void SetPixel(int x, int y, Color &col) {
  glColor4f(col.R, col.G, col.B, 1.0);
  //glColor4f(0.4,0.8,0.4,1.0);
  glVertex2i(x,y);
}

void putpixel(int x, int y, float alpha) {
  glColor4f(1.0,1.0,1.0,alpha);
  glVertex2i(x,y);
}

/**********************************/
/*         Rasterization          */
/**********************************/

void populateNormals()
{
  for(int i=0; i<trig.vNum(); i++)
  {
    Vector3f v = Vector3f(0.0,0.0,0.0);
    vertex_normals.push_back(v);
  }
  
  for(int i=0; i<trig.trigNum(); i++)
  {
    Vector3f v = Vector3f(0.0,0.0,0.0);
    face_normals.push_back(v);
  }
}

void computeNormals()
{
  int num_trigs = trig.trigNum();
  int num_vertices = trig.vNum();
  int *num_faces_for_vertex = new int[num_vertices];
  
  Vector3f v1,v2,v3;
  
  // Calculate the face normals
  for (int i=0; i<num_trigs-1; i++)
  {
    trig.getTriangleVertices(i, v1,v2,v3);
    //Vector3f face_normal = face_normals[i];
    
    Vector3f v1tov2 = Vector3f( (v2[0]-v1[0]), (v2[1]-v1[1]), (v2[2]-v1[2]) );
    Vector3f v2tov3 = Vector3f( (v3[0]-v2[0]), (v3[1]-v2[1]), (v3[2]-v2[2]) );
    
    Vector3f face_normal = Cross(v1tov2, v2tov3); // should be our face normal
    //face_normal[0] = v1tov2[0];
    //face_normal[1] = v1tov2[1];
    //face_normal[2] = v1tov2[2];
    
    normalise(face_normal); //Normalise the vector????
    face_normals[i] = face_normal;
  }
  
  // Calculate vertex normals (averaging the face normals)
  for (int i=0; i<num_trigs-1; i++)
  {
    Vector3f f_norm = face_normals[i];
    Triangle triggy = trig.getTriangle(i);
    
    for (int j=0; j<3; j++)
    {
      int v_index = triggy._vertex[j];
      num_faces_for_vertex[ v_index ] += 1;
      Vector3f v_norm = vertex_normals[ v_index ];
      
      v_norm[0] += f_norm[0];
      v_norm[1] += f_norm[1];
      v_norm[2] += f_norm[2];
      
      vertex_normals[ v_index ] = v_norm; 
    }
  }
  
  for (int i=0; i<num_vertices; i++)
  {
    int num_faces = num_faces_for_vertex[i];
    Vector3f v = vertex_normals[i];
    
    float x = v[0]/num_faces;
    float y = v[1]/num_faces;
    float z = v[2]/num_faces;
    
    v[0] = x;
    v[1] = y;
    v[2] = z;
    
    vertex_normals[i] = v;
  }
  /*
  Vector3f v = vertex_normals[0];
  printf("Vertex 0\n");
  printVector3f(v);
  printf("\n");
  */
  // Because we're preparing the array after run time
  delete [] num_faces_for_vertex;
}

void drawSpan(const Span &span, int y, int trignum)
{
	Vector3f v1,v2,v3;
	trig.getTriangleVertices(trignum,v1,v2,v3);
	
	//Triangle triggy = trig.getTriangle(trignum);
	
	int X1 = (int)v1[0];
	int X2 = (int)v2[0];
	int X3 = (int)v3[0];
	int Y1 = (int)v1[1];
	int Y2 = (int)v2[1];
	int Y3 = (int)v3[1];
	
	float detA = (float)( (X1*Y2)-(X1*Y3)-(X2*Y1)+(X2*Y3)+(X3*Y1)-(X3*Y2) );
	
	int v_ind1 = trig.getTriangleVertexIndex1(trignum);
	int v_ind2 = trig.getTriangleVertexIndex2(trignum);
	int v_ind3 = trig.getTriangleVertexIndex3(trignum);
	Vector3f v1_n = vertex_normals[v_ind1];
	Vector3f v2_n = vertex_normals[v_ind2];
	Vector3f v3_n = vertex_normals[v_ind3];
	
	//printf("indices: %i, %i, %i\n", v_ind1, v_ind2, v_ind3);
	
	for (int x = span.X1; x<span.X2; x++) {
    // Get the barycentric coords
    float detA1 = (float)( (x*Y2)-(x*Y3)-(X2*y)+(X2*Y3)+(X3*y)-(X3*Y2) );
    float detA2 = (float)( (X1*y)-(X1*Y3)-(x*Y1)+(x*Y3)+(X3*Y1)-(X3*y) );
    float detA3 = (float)( (X1*Y2)-(X1*y)-(X2*Y1)+(X2*y)+(x*Y1)-(x*Y2) );
  	
  	float u = detA1 / detA;
  	float v = detA2 / detA;
  	float w = detA3 / detA;
  	
  	float pz = u*v1[2] + v*v2[2] + w*v3[2];
  	//if (pz < 0.0) continue; // Don't draw anything at the back..
  	//if (y<0 || x<0) {continue;}
  	
  	if (pz < z_buffer[y+108][x+190]) {
      continue;
  	} else {
      z_buffer[y+108][x+190] = pz;
  	}
  	
  	Vector3f point_norm = Vector3f(u*v1_n[0] + v*v2_n[0] + w*v3_n[0], u*v1_n[1] + v*v2_n[1] + w*v3_n[1], u*v1_n[2] + v*v2_n[2] + w*v3_n[2]);
  	
  	Vector3f L_vect = Vector3f(x-light[0], y-light[1], pz-light[2]);
  	normalise(L_vect);
  	Vector3f V_vect = Vector3f(x-view[0], y-view[1], pz-view[2]); // pz should never be larger than view[2]=200
  	normalise(V_vect);
  	//normalise(point_norm);
  	
  	//r = i - (2 * n * dot(i, n))
  	float reflec = dot(V_vect, point_norm);
  	Vector3f ref2 = Vector3f(V_vect[0] - reflec*2*point_norm[0], V_vect[1] - reflec*2*point_norm[1], V_vect[2] - reflec*2*point_norm[2]);
  	normalise(ref2);
  	
  	float nDotL = dot(point_norm, L_vect);
  	float vDotR = dot(V_vect, ref2);
  	vDotR = pow(vDotR, small_n);
  	//vDotR = fabs(vDotR);
  	//nDotL = fabs(nDotL);
  	
  	float red = Ax.R*Ka + Att*Lx.R*( Kd*Dx.R*nDotL + Ks*Sx.R*vDotR );
  	float green = Ax.G*Ka + Att*Lx.G*( Kd*Dx.G*nDotL + Ks*Sx.G*vDotR );
  	float blue = Ax.B*Ka + Att*Lx.B*( Kd*Dx.B*nDotL + Ks*Sx.B*vDotR );
  	
  	if ((u + v + w) > 1.01) {
      printf("What the fuck? Lambdas: %f, %f, %f \n", u, v, w);
      printf("Added: %f\n", u+v+w);
      printf("For triangle: %i\n", trignum);
      exit(1);
  	}
   	
   	
   	
   	Color col = Color(red, green, blue, 1.0f);
   	//printf("RGB: %f, %f, %f\n", red,green,blue);
   	//printf("nDotl, vDotR: %f, %f\n", nDotL, vDotR);
   	//exit(1);
   	
    SetPixel(x, y, col);
	}
}

void drawSpansBetweenEdges(const Edge &e1, const Edge &e2, int trignum)
{
	// calculate difference between the y coordinates
	// of the first edge and return if 0
	float e1ydiff = (float)(e1.Y2 - e1.Y1);
	if(e1ydiff == 0.0f)
		return;

	// calculate difference between the y coordinates
	// of the second edge and return if 0
	float e2ydiff = (float)(e2.Y2 - e2.Y1);
	if(e2ydiff == 0.0f)
		return;

	// calculate differences between the x coordinates
	// and colors of the points of the edges
	float e1xdiff = (float)(e1.X2 - e1.X1);
	float e2xdiff = (float)(e2.X2 - e2.X1);
	Color e1colordiff = (e1.Color2 - e1.Color1);
	Color e2colordiff = (e2.Color2 - e2.Color1);

	// calculate factors to use for interpolation
	// with the edges and the step values to increase
	// them by after drawing each span
	float factor1 = (float)(e2.Y1 - e1.Y1) / e1ydiff;
	float factorStep1 = 1.0f / e1ydiff;
	float factor2 = 0.0f;
	float factorStep2 = 1.0f / e2ydiff;

	// loop through the lines between the edges and draw spans
	for(int y = e2.Y1; y < e2.Y2; y++) {
		// create and draw span
		Span span(e1.Color1 + (e1colordiff * factor1),
		          e1.X1 + (int)(e1xdiff * factor1),
		          e2.Color1 + (e2colordiff * factor2),
		          e2.X1 + (int)(e2xdiff * factor2));
		drawSpan(span, y, trignum);

		// increase factors
		factor1 += factorStep1;
		factor2 += factorStep2;
	}
}

void drawTriangle(const Color &color1, float x1, float y1,
                  const Color &color2, float x2, float y2,
                  const Color &color3, float x3, float y3, int trignum)
{
	// create edges for the triangle
	Edge edges[3] = {
		Edge(color1, (int)x1, (int)y1, color2, (int)x2, (int)y2),
		Edge(color2, (int)x2, (int)y2, color3, (int)x3, (int)y3),
		Edge(color3, (int)x3, (int)y3, color1, (int)x1, (int)y1)
	};

	int maxLength = 0;
	int longEdge = 0;

	// find edge with the greatest length in the y axis
	for(int i = 0; i < 3; i++) {
		int length = edges[i].Y2 - edges[i].Y1;
		if(length > maxLength) {
			maxLength = length;
			longEdge = i;
		}
	}

	int shortEdge1 = (longEdge + 1) % 3;
	int shortEdge2 = (longEdge + 2) % 3;

	// draw spans between edges; the long edge can be drawn
	// with the shorter edges to draw the full triangle
	drawSpansBetweenEdges(edges[longEdge], edges[shortEdge1], trignum);
	drawSpansBetweenEdges(edges[longEdge], edges[shortEdge2], trignum);
}

/**********************************/
/*        Line algorithms         */
/**********************************/
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
/*          Line Drawing          */
/**********************************/

void DoBresenham(Vector3f v1, Vector3f v2, Vector3f v3) {
  Bresenham(v1[0], v1[1], v2[0], v2[1]);
  Bresenham(v1[0], v1[1], v3[0], v3[1]);
  Bresenham(v2[0], v2[1], v3[0], v3[1]);
}

/**********************************/
/*            Display             */
/**********************************/

void myDisplay()
{
  glClear(GL_COLOR_BUFFER_BIT); // Clear OpenGL Window
  
  int trignum = trig.trigNum();
  Vector3f v1, v2, v3;
  
  glColor4f(1.0,1.0,1.0,1.0); // Needed for every algorithm except anti-aliasing
  Color col1 = Color(1.0f,0.0f,0.0f,1.0f);
  Color col2 = Color(0.0f,1.0f,0.0f,1.0f);
  Color col3 = Color(0.0f,0.0f,1.0f,1.0f);
  
  for (int i = 0; i < trignum-1; i++) {
    trig.getTriangleVertices(i, v1,v2,v3);
    
    //Rotate(v1, v2, v3, rx, ry, rz);
    
    glBegin(GL_POINTS);
      drawTriangle(col1, v1[0], v1[1],
                   col2, v2[0], v2[1],
                   col3, v3[0], v3[1],
                   i);
    glEnd();
  }
  
  glFlush(); // Output everything
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
  
  // set our parameters
  Lx = Color(1.0f,1.0f,1.0f,1.0f); // Light colour is white
  Ax = Color(0.6f,0.8f,0.4f,1.0f); // Ambient colour is greeny
  Dx = Color(1.0f,1.0f,1.0f,1.0f); // Diffuse set
  Sx = Color(1.0f,1.0f,1.0f,1.0f); // Specular set to nothing
  Ka = 1.0f; // Ambient coefficient (intensity)
  Kd = 0.8f; // Diffuse coefficient
  Ks = 0.1f; // Specular coefficient
  Att = 0.6f; // Attenuation coefficient
  small_n = 6; // Shine/roughness
  
  //Vector3f big_N; // Surface normal of pixel
  Vector3f ll = Vector3f(0.0f,0.0f,200.0f);
  light = ll; // Vector from light to pixel
  //Vector3f reflection; // Reflection vector from pixel
  Vector3f vv = Vector3f(0.0f,0.0f,200.0f);
  view = vv; // Vector from pixel to eye (or viceversa) (0.0,0.0,1.0)
  
  // We can precomputer the ambient colour
  ambient = Color(Ax.R*Ka*Dx.R, Ax.G*Ka*Dx.G, Ax.B*Ka*Dx.B);
  
  populateNormals();
  computeNormals();
  
  glutDisplayFunc(myDisplay);// Callback function
  //glutIdleFunc(myDisplay); // Display this while idling
  
  glutMainLoop();// Display everything and wait
}