# Makefile

CPP = g++
LDFLAGS = -lGLU -lGL -lX11 -framework GLUT -framework OpenGL
OBJS = main.o shapes.o viewing.o menu.o
HEADERS = shapes.h viewing.h menu.h

INCLUDE = -I/opt/local/include
LIBDIR  = -L/opt/local/lib
CFLAGS = -Wall -g $(INCLUDE) $(LIBDIR) -framework GLUT -framework OpenGL

# Compile the program.

all: t03m01

t03m01: $(OBJS)
	$(CPP) $(CFLAGS) -o t03m01 $(OBJS) $(LDFLAGS)

main.o: main.cpp $(HEADERS)
	$(CPP) $(CFLAGS) -c main.cpp

shapes.o: shapes.cpp $(HEADERS)
	$(CPP) $(CFLAGS) -c shapes.cpp

viewing.o: viewing.cpp $(HEADERS)
	$(CPP) $(CFLAGS) -c viewing.cpp

menu.o: menu.cpp $(HEADERS)
	$(CPP) $(CFLAGS) -c menu.cpp

clean:
	rm -f t03m01 $(OBJS)