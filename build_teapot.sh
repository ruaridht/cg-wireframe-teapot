#!/bin/bash

# Build script, Ruaridh Thomson

CC="g++"
OS_TYPE="Darwin"
LDFLAGS="-lGL -lglut -lGLU"

if [ $OS_TYPE == "Darwin" ]
	then
  	LDFLAGS="-framework GLUT -framework OpenGL"
fi

$CC -o teapot teapot.cpp parser.cpp triangle.cpp wireframe.cpp $LDFLAGS

#./teapot scenes/MIT_teapot.obj