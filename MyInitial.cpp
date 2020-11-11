//
//  MyInitial.cpp   - Version 0.4 - October 24, 2020
//
//   Sets up and renders the initial (alphabet letter)
//   for the Math 155A project.
//
//   Comes supplied with some code for a cylinders, a torus,
//      and a revolving ellipsoid.
//
//  THIS FILE IS WHAT YOU WILL MODIFY FOR PROJECT #3.
//  IT WILL ALSO BE USED FOR PROJECTS #4 and #5.
//


// Use the static library (so glew32.dll is not needed):
#define GLEW_STATIC
#include <GL/glew.h> 
#include <GLFW/glfw3.h>

#include "LinearR3.h"		// Adjust paths as needed.
#include "LinearR4.h"		
#include "GlGeomSphere.h"
#include "GlGeomCylinder.h"
#include "GlGeomTorus.h"
#include <cmath>

// Enable standard input and output via printf(), etc.
// Put this include *after* the includes for glew and GLFW!
#include <stdio.h>

#include "MyInitial.h"
#include "InitialProj.h"
//#include "SurfaceProj.h"      // Replace previous line with this line for Project #4

// ***********
// Animation Controls.  
//     Optional: you may customize these for Project #3.
// ***********

// YOU MAY WISH TO RE-DO THIS FOR YOUR CUSTOM ANIMATION.  
double animateIncrement = 0.01;   // Make bigger to speed up animation, smaller to slow it down.
double currentTime = 0.0;         // Current "time" for the animation.
double maxTime = 2.0;             // Time cycles back to 0 after reaching maxTime.

// These two variables control whether running or paused.
//  IT IS OPTIONAL TO INCLUDE THIS FUNCTIONALITY IN YOUR PROGRAMMING ASSIGNMENT
bool spinMode = true;
bool singleStep = false;


// These objects take care of generating and loading VAO's, VBO's and EBO's,
//    rendering ellipsoids and cylinders
// IF YOU ADDED EXTRA TORII, THEY SHOULD BE DECLARED HERE
GlGeomSphere unitSphere;
GlGeomCylinder unitCylinder;
GlGeomTorus torus1(6, 6, 0.1f);   // Set default mesh resolutions and inner radius
                                  // Initialize multiple tori if they have different inder radii.

// **********************
// This sets up a sphere and a cylinder and a torus needed for the "Initial" (the 3-D alphabet letter)
//  This routine is called only once, for the first initialization.
//  IF YOU ADD EXTRA TORII, THEY NEED TO HANDLED HERE.
// **********************
void MySetupInitialGeometries() {
    // Initialize the sphere and cylinder, and give them the vertPos, vertNormal, and vertTexCoords locations.
    MyRemeshGeometries();
    unitSphere.InitializeAttribLocations(vPos_loc,vNormal_loc,vTexcoords_loc);
    unitCylinder.InitializeAttribLocations(vPos_loc, vNormal_loc, vTexcoords_loc);
    torus1.InitializeAttribLocations(vPos_loc, vNormal_loc, vTexcoords_loc);

    check_for_opengl_errors();
}

// *********************
// This is called when geometric shapes are initialized.
// And is called again whenever the mesh resolution changes.
// IF YOU ADD EXTRA TORII, THEY NEED TO BE HANDLED HERE
// ********************
void MyRemeshGeometries() {
    unitSphere.Remesh(meshRes, meshRes);              // Number of slices and stacks both set to meshRes
    unitCylinder.Remesh(meshRes, meshRes, meshRes);   // Number of slices, stacks and rings all set to meshRes
    torus1.Remesh(meshRes, meshRes);                  // Number of rings and number of sides per ring.
}

// *************************************
// Render the initial (3D alphabet letter)
// THIS CODE IS THE CORE PART TO RE_WRITE FOR YOUR 155A PROJECT  ****************************
// ************
void MyRenderInitial() {
    // Compute the "currentTime" for the animation.
    //    As initially implemented, CurrentTime goes from 0.0 to 1.0, and then back to 0.0
    //    THIS IS SPECIFIC TO THE ANIMATION IN THE DEMO.
    //    FOR PROJECT 3 YOU MAY DO SOMETHING DIFFERENT, FOR INSTANCE, SIMILAR TO WHAT 
	//       THE SOLAR SYSTEM PROJECT DID.
    //  
    if (spinMode) {
        currentTime += animateIncrement;
        if (currentTime >= maxTime) {
            spinMode = false;  // stop animation
        }
        if (singleStep) {
            spinMode = false;       // If in single step mode, turn off future animation
        }
    }

    // Render the letter "X" (Sam's initial) with two cylinders,
    //    Plus a revolving ellipsoid.
   
    // Setup the main Model view matrix
    LinearMapR4 mat1 = viewMatrix;              // Base off of viewMatrix
    mat1.Mult_glTranslate(-2.5, 2.0, -2.5);     // Center of the letter

    // Make the X partgreen-ish (YOU ARE ENCOURAGED TO ALTER COLORS)
    glVertexAttrib3f(vColor_loc, currentTime / maxTime, currentTime/maxTime, currentTime / maxTime);
    
    // Draw middle span of the N
    LinearMapR4 mat2 = mat1;                   // Modelview matrix for cylinders & sphere
    mat2.Mult_glTranslate(0.0, 0.0, 0.0);         // Translate slightly towards the viewer  
    mat2.Mult_glRotate(1/currentTime, 0.0, 1.0, 0.0); // revolve around axis
    mat2.Mult_glRotate(PIsixths, 0.0, 0.0, 1.0); // Rotate -30 degrees
    mat2.Mult_glScale(0.4, 2.2, 0.2);          // Scale the cylinder, to thiner, flater and taller 
    mat2.DumpByColumns(matEntries);
    glUniformMatrix4fv(modelviewMatLocation, 1, false, matEntries);
    unitCylinder.Render();

    // Draw the left cylinder of the N
    mat2 = mat1;                               // Back to the main Modelview matrix
    mat2.Mult_glRotate(1/currentTime, 0.0, 1.0, 0.0); //revolve with middle span, decelerating
    mat2.Mult_glTranslate(-3.2 + (currentTime/1), 0.0, 0.0);        // Translate slightly away from the viewer
    mat2.Mult_glScale(0.4, 2.0, 0.2);          // Scale the cylinder, to thiner, flater and taller 
    mat2.DumpByColumns(matEntries);
    glUniformMatrix4fv(modelviewMatLocation, 1, false, matEntries);
    unitCylinder.Render();

    // Draw the right cylinder of the N
    mat2 = mat1;                               // Back to the main Modelview matrix
    mat2.Mult_glRotate(1 / currentTime, 0.0, 1.0, 0.0); //revolve with middle span
    mat2.Mult_glTranslate(3.2 - (currentTime / 1), 0.0, 0.0);        // Translate in and 
    mat2.Mult_glScale(0.4, 2.0, 0.2);          // Scale the cylinder, to thiner, flater and taller 
    mat2.DumpByColumns(matEntries);
    glUniformMatrix4fv(modelviewMatLocation, 1, false, matEntries);
    unitCylinder.Render();

    
    // Render top torus
    glVertexAttrib3f(vColor_loc, 0.8f, 0.0f, 0.0f);  // Red color (slightly darkened)
    mat2 = mat1;                              // Back to the main Modelview matrix
    mat2.Mult_glRotate(PIsixths, 0.0, 0.0, 1.0);
    mat2.Mult_glTranslate(0.0, -2.5 + pow(3.0, currentTime), 0.0);
    mat2.Mult_glScale(3 - pow(2.0, currentTime - .3));                   // Uniform scaling
    mat2.DumpByColumns(matEntries);
    glUniformMatrix4fv(modelviewMatLocation, 1, false, matEntries);
    if (currentTime <= maxTime - .1) {
        torus1.Render();
    }

    // Render bottom torus
    glVertexAttrib3f(vColor_loc, 0.8f, 0.0f, 0.0f);  // Red color (slightly darkened)
    mat2 = mat1;                              // Back to the main Modelview matrix
    mat2.Mult_glRotate(PIsixths, 0.0, 0.0, 1.0);
    mat2.Mult_glTranslate(0.0, 2.5 - pow(3.0, currentTime), 0.0);
    mat2.Mult_glScale(3 - pow(2.0, currentTime - .3));                   // Uniform scaling
    mat2.DumpByColumns(matEntries);
    glUniformMatrix4fv(modelviewMatLocation, 1, false, matEntries);
    if (currentTime <= maxTime - .1) {
        torus1.Render();
    }

    // Render base torus 1
    glVertexAttrib3f(vColor_loc, 0.5f, 0.5f, 1.0f);  // blueish color (slightly lightened)
    mat2 = mat1;                              // Back to the main Modelview matrix
    mat2.Mult_glRotate(1 / currentTime, 0.0, 1.0, 0.0);
    mat2.Mult_glTranslate(-3.2 + (currentTime / 1), -3.4 +  3*(2 / pow(2.0, currentTime)), 0.0);
    mat2.Mult_glScale(0.4, 2.0, 0.2);                   // Uniform scaling
    mat2.DumpByColumns(matEntries);
    glUniformMatrix4fv(modelviewMatLocation, 1, false, matEntries);
    torus1.Render();

    // Render base torus 2
    glVertexAttrib3f(vColor_loc, 0.5f, 0.5f, 1.0f);  // blueish color (slightly lightened)
    mat2 = mat1;                              // Back to the main Modelview matrix
    mat2.Mult_glRotate(1 / currentTime, 0.0, 1.0, 0.0);
    mat2.Mult_glTranslate(3.2 - (currentTime / 1), 3.4 - 3 * (2 / pow(2.0, currentTime)), 0.0);
    mat2.Mult_glScale(0.4, 2.0, 0.2);                   // Uniform scaling
    mat2.DumpByColumns(matEntries);
    glUniformMatrix4fv(modelviewMatLocation, 1, false, matEntries);
    torus1.Render();

    // render i
    glVertexAttrib3f(vColor_loc, 0.4*(currentTime - 1), 0.3 * (currentTime - 1), 0.8*(currentTime - 1));  // purpleish?
    LinearMapR4 imat = mat1;                               // Back to the main Modelview matrix
    imat.Mult_glTranslate(13.0 - currentTime*5, -1.0, 0.0);        // Translate in and 
    imat.Mult_glScale(0.4, 0.7, 0.2);          // Scale the cylinder, to thiner, flater and taller 
    imat.DumpByColumns(matEntries);
    glUniformMatrix4fv(modelviewMatLocation, 1, false, matEntries);
    if (currentTime >= 1.0) {
        unitCylinder.Render();
    }
    
    // render dot on i
    glVertexAttrib3f(vColor_loc, 0.4*(pow(2.0, currentTime) - 3), 0.3 * (pow(2.0, currentTime) - 3), 0.8*(pow(2.0, currentTime) - 3));  // purpleish?
    mat2 = mat1;
    mat2.Mult_glTranslate(13.0 - pow(3.16, currentTime), 0.5, 0.0);
    mat2.Mult_glScale(0.4, 0.4, 0.2);
    mat2.DumpByColumns(matEntries);
    glUniformMatrix4fv(modelviewMatLocation, 1, false, matEntries);
    if (currentTime >= 1.5) {
        unitSphere.Render();
    }

    /*
    // create c location and scaling, all that is needed is rotation
    glVertexAttrib3f(vColor_loc, 0.2f, 0.0f, 0.8f);  // purpleish?
    LinearMapR4 clocmat = mat1;
    clocmat.Mult_glTranslate(4.0, -1.0, 0.0);
    clocmat.Mult_glScale(0.2, 0.6, 0.2);

    // render center of c
    clocmat.DumpByColumns(matEntries);
    glUniformMatrix4fv(modelviewMatLocation, 1, false, matEntries);
    if (currentTime >= 1.5) {
        unitCylinder.Render();
    }

    //create top part of c
    glVertexAttrib3f(vColor_loc, 0.2f, 0.0f, 0.8f);  // purpleish?
    LinearMapR4 ctop = clocmat;
    ctop.Mult_glRotate(PIhalves, 0.0, 0.0, 1.0);
    ctop.Mult_glTranslate(0.2, 0.2, 0.0);
    ctop.Mult_glScale(0.0, .66, 0.0);
    ctop.DumpByColumns(matEntries);
    glUniformMatrix4fv(modelviewMatLocation, 1, false, matEntries);
    if (currentTime >= 1.5) {
        unitCylinder.Render();
    }
 
    // create bottom part of c
    glVertexAttrib3f(vColor_loc, 0.2f, 0.0f, 0.8f);  // purpleish?
    LinearMapR4 cbot = clocmat;
    cbot.Mult_glTranslate(0.2, -0.2, 0.0);
    cbot.Mult_glScale(0.0, .66, 0.0);
    cbot.Mult_glRotate(PIhalves, 0.0, 0.0, 1.0);
    cbot.DumpByColumns(matEntries);
    glUniformMatrix4fv(modelviewMatLocation, 1, false, matEntries);
    if (currentTime >= 1.5) {
        unitCylinder.Render();
    }

    




    /*
    // Render the revolving ellipsoid
    glVertexAttrib3f(vColor_loc, 0.6f, 0.4f, 1.0f);  // Blue/Magenta-ish color
    mat2 = mat1;                              // Back to the main Modelview matrix
    mat2.Mult_glRotate(currentTime*PI2, 0.0, 1.0, 0.0);   // PI2 is 2*pi (defined in MathMisc.h)
    mat2.Mult_glTranslate(0.0, 0.0, 1.8);     // Pull towards the viewer (for revolution)
    mat2.Mult_glScale(0.5, 0.2, 1.0);         // Nonuniform scaling  
    mat2.DumpByColumns(matEntries);
    glUniformMatrix4fv(modelviewMatLocation, 1, false, matEntries);
    unitSphere.Render();
    */
}








