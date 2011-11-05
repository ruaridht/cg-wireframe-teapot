README

Computer Graphics Practical 1
Ruaridh Thomson
s0786036

This readme is slightly long, sorry.

To run:
The program will compile and run using the following command
./build.sh

----------Before running-----------
Some code will need to be commented/uncommented to demonstrate the different algorithms.
Apologies for requiring code to be commented/uncommented, though viewing various techniques proved to be much easier by commenting/uncommenting than passing arguments to the program (due to the various combinations).

By default the program is automatic (no mouse/keyboard input) as this better demonstrates the algorithms.  The interactive mode only uses the bresenham/midpoint algorithm.

-----------Capabilities------------
The program can:
- Scale the object.
- Rotate along the x,y,z axis.
- Translate along the x,y,z axis.
- Shear along the x,y axis (for simplicity).

Draw lines via:
- DDA Algorithm
- Bresenhams/midpoint
- Xiaolin Wu's patterns
- Gupta-Sproull
- AAL (bresenham variant for anti-aliasing)
- EFLA (Extremely Fast Line Algorithm -  http://www.edepot.com/algorithm.html)

-------Switching algorithms--------
The program needs to be re-compiled and run (using ./build.sh) when changes are made.  All code that concerns running the program can be found at the bottom of the file demo.cc, below the comment "For Marker".

There are two methods at the bottom of the file, main() and myDisplay().
An interactive mode can be achieved by following 'To Set Interactive' in main().

Otherwise, from myDisplay() (just above main) different line algorithms can be used by commenting/uncommenting them (under "LINE DRAWING").

Also from myDisplay(), rotation, scale, translation, shear can be used by uncommenting the desired manipulation.

If desired the rotation speed, translation, scale and shear factors can be changed.

---------Interactive mode --------
To control while interactive:
w,s - rotate z
a,d - rotate y

Left mouse click and hold and move to rotate the object.

---------Notes on line algorithms--------
Anti-aliasing(AA) draws thicker lines and is better viewed while the object is scaled up.

DDA performs as expected.
Bresenham/midpoint also performs as expected.
Wu lines are thick, likely expecting anti-aliasing.  Though I could make it draw well while anti-aliasing.
Gupta-Sproul anti-aliases to an extent.  (The anti-aliasing may not be completely accurate)
AAL should draw and anti-alias as expected.
EFLA draws as expected.

-----------Final comments--------------
Everything should work with little complication.  Translating along the z-axis does not produce any clear result (due to the camera?).

The program has been tested to run on DICE, please email me if there are any problems.
s0786036@sms.ed.ac.uk.