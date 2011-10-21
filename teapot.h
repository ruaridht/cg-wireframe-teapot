#ifndef TEAPOT_TEAPOT_H
#define TEAPOT_TEAPOT_H

#include "headers.h"

using namespace std;

// Global forward declarations
struct Options {
    Options() { quiet = verbose = false;
                imageFile = ""; }
    bool quiet, verbose;
    string imageFile;
};

// Global constants
#define TEAPOT_VERSION "0.1"
#define WIDTH "640"
#define HEIGHT "480"
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

/* These aren't used.
inline float fmax(float f1,float f2, float f3) {
	float f = f1;

	if (f < f2) f = f2;
	if (f < f3) f = f3;

	return f;
}

inline float fmin(float f1,float f2, float f3) {
	float f = f1;

	if (f > f2) f = f2;
	if (f > f3) f = f3;

	return f;
}
*/

#endif // TEAPOT_TEAPOT_H