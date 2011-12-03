#!/bin/bash

# Build script, Ruaridh Thomson

CC="g++"
OS_TYPE="Darwin"
LDFLAGS="-lGL -lglut -lGLU"

if [ $OS_TYPE == "Darwin" ]
	then
  	LDFLAGS="-framework GLUT -framework OpenGL"
fi

$CC -o demo demo1.cc Color.cpp $LDFLAGS

./demo scenes/MIT_teapot.obj