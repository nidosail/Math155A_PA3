#pragma once

// 
// MyInitial.h   ---  Header file for Initial.cpp.
// 
//   Sets up and renders the initial (alphabet letter)
//   for the Math 155A project.
//
//   Comes supplied with the code for Sam Buss's "S".
//

// These objects take care of generating and loading VAO's, VBO's and EBO's,
//    They are used for the ellipsoids and cylinders in the initial.
//    They are extern'ed here to let InitialProj.cpp use them too,
//    but MyInitial.cpp is responsible for maintaining them.
extern GlGeomSphere unitSphere;
extern GlGeomCylinder unitCylinder;

// These variables control the animation.
// MyInitial.cpp has the primary responsibility of maintaining them.
// But, they are extern'd here to let the main program change them 
//    in response to keystrokes.
extern double animateIncrement;
extern double currentTime;
extern bool spinMode;
extern bool singleStep;


//
// Function Prototypes
//
void MySetupInitialGeometries();   // Called once, before rendering begins.
void MyRemeshGeometries();         // Called when mesh changes, must update initial's goemetries.

void MyRenderInitial();


