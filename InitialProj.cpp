/*
 * InitialProj.cpp - Version 0.4 - October 24, 2020
 *
 * Starting code for Math 155A, Project #3,
 * 
 * Author: Sam Buss
 *
 * Software accompanying POSSIBLE SECOND EDITION TO the book
 *		3D Computer Graphics: A Mathematical Introduction with OpenGL,
 *		by S. Buss, Cambridge University Press, 2003.
 *
 * Software is "as-is" and carries no warranty.  It may be used without
 *   restriction, but if you modify it, please change the filenames to
 *   prevent confusion between different versions.
 * Bug reports: Sam Buss, sbuss@ucsd.edu.
 * Web page: http://math.ucsd.edu/~sbuss/MathCG2
 */

// These libraries are needed to link the program.
// First five are usually provided by the system.
#pragma comment(lib,"opengl32.lib")
#pragma comment(lib,"glu32.lib")
#pragma comment(lib,"glfw3.lib")
#pragma comment(lib,"glew32s.lib")
#pragma comment(lib,"glew32.lib")


// Use the static library (so glew32.dll is not needed):
#define GLEW_STATIC

#include <GL/glew.h> 
#include <GLFW/glfw3.h>

#include "LinearR3.h"		// Adjust path as needed.
#include "LinearR4.h"		// Adjust path as needed.
#include "GlGeomSphere.h"
#include "GlGeomCylinder.h"
#include "GlGeomTorus.h"
#include "GlShaderMgr.h"

// Enable standard input and output via printf(), etc.
// Put this include *after* the includes for glew and GLFW!
#include <stdio.h>

#include "InitialProj.h"
#include "MyInitial.h"



// ********************
// Animation controls and state infornation
// ********************

// These variables control the view direction.
//    The arrow keys are used to change these values.
double viewAzimuth = 0.25;	// Angle of view up/down (in radians)
double viewDirection = 0.0; // Rotation of view around y-axis (in radians)
double deltaAngle = 0.01;	// Change in view angle for each up/down/left/right arrow key press
double extraViewDistance = 0.0;  // Added to distance of viewer from center of the scene
LinearMapR4 viewMatrix;		// The current view matrix, based on viewAzimuth and viewDirection.

// These two variables control how triangles are rendered.
bool wireframeMode = false;	// Equals true for polygon GL_FILL mode. False for polygon GL_LINE mode.
bool cullBackFaces = true; // Equals true to cull backfaces. Equals false to not cull backfaces. 

// The next variable controls the resoluton of the meshes for cylinders and spheres.
int meshRes=10;             // Resolution of the meshes (slices, stacks, and rings all equal)

// this boolean determines whether to reverse the spin or not
bool spinRev = false;

// ************************
// General data helping with setting up VAO (Vertex Array Objects)
//    and Vertex Buffer Objects.
// ***********************

// We create one shader program: consisting of a vertex shader and a fragment shader
unsigned int shaderProgram1;
const unsigned int vPos_loc = 0;     // Corresponds to "location = 0" in the verter shader definitions
const unsigned int vColor_loc = 1;   // Corresponds to "location = 1" in the verter shader definitions
const unsigned int vNormal_loc = UINT_MAX;          // UINT_MAX means no normal vector used by the vertex shader
const unsigned int vTexcoords_loc = UINT_MAX;       // UINT_MAX means no texture coordinates used by the vertex shader
unsigned int projMatLocation;						// Location of the projectionMatrix in the shader program.
unsigned int modelviewMatLocation;					// Location of the modelviewMatrix in the shader program.

const char* projMatName = "projectionMatrix";		// Name of the uniform variable projectionMatrix
const char* modelviewMatName = "modelviewMatrix";	// Name of the uniform variable modelviewMatrix

//  The Projection matrix: Controls the "camera view/field-of-view" transformation
//     Generally is the same for all objects in the scene.
LinearMapR4 theProjectionMatrix;		//  The Projection matrix: Controls the "camera/view" transformation

// A ModelView matrix controls the placement of a particular object in 3-space.
//     It is generally different for each object.
// The array matEntries holds the matrix values as floats to be loaded into the shader program. 
float matEntries[16];		// Holds 16 floats (since cannot load doubles into a shader that uses floats)

// *****************************
// These variables set the dimensions of the perspective region we wish to view.
// They are used to help form the projection matrix and the view matrix
// All rendered objects lie in the rectangular prism centered on the z-axis
//     equal to (-Xmax,Xmax) x (-Ymax,Ymax) x (Zmin,Zmax)
// Be sure to leave some slack in these values, to allow for rotations, etc.
// The model/view matrix can be used to move objects to this position
// THESE VALUES MAY NEED AD-HOC ADJUSTMENT TO GET THE SCENE TO BE VISIBLE.
const double Xmax = 7.0;                // Control x dimensions of viewable scene
const double Ymax = 5.0;                // Control y dimensions of viewable scene
const double Zmin = -9.0, Zmax = 9.0;   // Control z dimensions of the viewable scene
// zNear equals the distance from the camera to the z = Zmax plane
const double zNear = 20.0;              // Make this value larger or smaller to affect field of view.


// ************************
// General data helping with setting up VAO (Vertex Array Objects)
//    and Vertex Buffer Objects.
// ***********************

const int NumObjects = 1;
const int iFloor = 0;

unsigned int myVBO[NumObjects];  // a Vertex Buffer Object holds an array of data
unsigned int myVAO[NumObjects];  // a Vertex Array Object - holds info about an array of vertex data;
unsigned int myEBO[NumObjects];  // an Element Buffer Object - holds indices of vertices for glDrawElements

// *************************
// mySetupGeometries defines the scene data, especially vertex  positions and colors.
//    - It also loads all the data into the VAO's (Vertex Array Objects) and
//      into the VBO's (Vertex Buffer Objects).
// This routine is only called once to initialize the data.
// *************************
void mySetupGeometries() {
 
    // Allocate the needed Vertex Array Objects (VAO's) and Vertex Buffer Objects (VBO')
    glGenVertexArrays(NumObjects, &myVAO[0]);
    glGenBuffers(NumObjects, &myVBO[0]);
    glGenBuffers(NumObjects, &myEBO[0]);

    mySetupFloorGeometry();
    MySetupInitialGeometries();

    mySetViewMatrix();

    check_for_opengl_errors();   // Really a great idea to check for errors -- esp. good for debugging!
}

// **********************
// This sets up the data structures needed for the floor
//   The "baby" version is fixed rectangle.
//   The real version is can be remeshed on the fly.
// **********************
void mySetupFloorGeometry() {
    check_for_opengl_errors();
    // First Geometry: Forms the rectangular floor 
    //   as a triangle fan, centered at the center.
    float floorVerts[] = {
        // Positions			
        5.0, 0.0, 5.0,      // Front right vertex - Vertex 0
        -5.0, 0.0, 5.0,     // Front left vertex - Vertex 1
        5.0, 0.0, -5.0,     // Back right vertex - Vertex 2
        -5.0, 0.0, -5.0,    // Back left vertex - Vertex 3
        0.0, 0.0, 0.0,      // Center point, Vertex 4
    };
    // The indices of the vertices in the order needed for the triangle fan.
    //   These are called "elements".
    unsigned int floorVertexIndices[]{
        4, 1, 0, 2, 3, 1
    };

    glBindVertexArray(myVAO[iFloor]);
    glBindBuffer(GL_ARRAY_BUFFER, myVBO[iFloor]);
    glBufferData(GL_ARRAY_BUFFER, sizeof(floorVerts), floorVerts, GL_STATIC_DRAW);
    glVertexAttribPointer(vPos_loc, 3, GL_FLOAT, GL_FALSE, 0, (void*)0);	// Store vertices in the VBO
    glEnableVertexAttribArray(vPos_loc);									// Enable the stored vertices
 
                                                                            // Bind the Element Buffer object to the VAO, and load its data from the array.
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, myEBO[iFloor]);
    check_for_opengl_errors();
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(floorVertexIndices), floorVertexIndices, GL_STATIC_DRAW);
}

void mySetViewMatrix() {
    // Set the view matrix. Sets view distance, and view direction.
    // The final translation is done because the ground plane lies in the xz-plane,
    //    se the center of the scene is about 3 or 4 units above the origin.
    // YOU MAY NEED TO ADJUST THE FINAL TRANSLATION.
    viewMatrix.Set_glTranslate(0.0, 0.0, -(Zmax + zNear)+extraViewDistance);      // Translate to be in front of the camera
    viewMatrix.Mult_glRotate(viewAzimuth, 1.0, 0.0, 0.0);	    // Rotate viewAzimuth radians around x-axis
    viewMatrix.Mult_glRotate(-viewDirection, 0.0, 1.0, 0.0);    // Rotate -viewDirection radians around y-axis
    viewMatrix.Mult_glTranslate(0.0, -3.5, 0.0);                // Translate the scene down the y-axis so the center is near the origin.
}

// *************************************
// Main routine for rendering the scene
// myRenderScene() is called every time the scene needs to be redrawn.
// mySetupGeometries() has already created the vertex and buffer objects
//    and the model view matrices.
// setup_shaders() has already created the shader programs.
// *************************************
void myRenderScene() {
   
    // Clear the rendering window
    static const float black[] = { 0.0f, 0.0f, 0.0f, 0.0f };
    const float clearDepth = 1.0f;
    glClearBufferfv(GL_COLOR, 0, black);
    glClearBufferfv(GL_DEPTH, 0, &clearDepth);	// Must pass in a *pointer* to the depth
 
    MyRenderFloor();
    MyRenderExtraCylinders();
    MyRenderInitial();

    check_for_opengl_errors();   // Really a great idea to check for errors -- esp. good for debugging!
}

// ***********************************
// Draw the 10x10 ground plane in the xz-plane as a two triangle triangle strip.
// *********
void MyRenderFloor() {
    viewMatrix.DumpByColumns(matEntries);
    glUniformMatrix4fv(modelviewMatLocation, 1, false, matEntries);
    glBindVertexArray(myVAO[iFloor]);
    glVertexAttrib3f(vColor_loc, 1.0f, 0.4f, 0.4f);			// A light red color (R, G, B values).
    glDrawElements(GL_TRIANGLE_FAN, 6, GL_UNSIGNED_INT, 0);    // Use the EBO for this Triangle Fan (since glDrawElements)
}

// ***********************************
// Render two extra cylinders in front right section
// **********
void MyRenderExtraCylinders() {
    LinearMapR4 mat1 = viewMatrix;
    mat1.Mult_glTranslate(2.5, 0.0, 2.5);   // Center of front right quadrant of ground plane
    LinearMapR4 mat2 = mat1;
    mat2.Mult_glTranslate(0.0, 0.3, 0.0);   // Position cylinder       
    mat2.Mult_glScale(2.75, 0.2, 2.75);     // It is a very flat cylinder
    mat2.DumpByColumns(matEntries);
    glUniformMatrix4fv(modelviewMatLocation, 1, false, matEntries);
    glVertexAttrib3f(vColor_loc, 0.9f, 0.2f, 0.7f);     // Make the cylinder magenta-ish
    unitCylinder.Render();
    mat2 = mat1;                             // Base off of mat1 again
    mat2.Mult_glTranslate(0.0, 1.6, 0.0);   // Position cylinder barely above cylinder   
    mat2.Mult_glScale(0.8, 1.0, 0.8);
    mat2.DumpByColumns(matEntries);
    glUniformMatrix4fv(modelviewMatLocation, 1, false, matEntries);
    glVertexAttrib3f(vColor_loc, 0.7f, 0.9f, 0.0f);     // Make the cylinder a yellow-ish
    unitCylinder.Render();
}

void my_setup_SceneData() {

	mySetupGeometries();

    GlShaderMgr::LoadShaderSource("InitialProj.glsl");
    shaderProgram1 = GlShaderMgr::CompileAndLinkAll();
	
	// Get the locations of the projection and model view matrices in the shader programs.
	projMatLocation = glGetUniformLocation(shaderProgram1, projMatName);
	modelviewMatLocation = glGetUniformLocation(shaderProgram1, modelviewMatName);

	check_for_opengl_errors();   // Really a great idea to check for errors -- esp. good for debugging!
}

// *******************************************************
// Process all key press events.
// This routine is called each time a key is pressed or released.
// *******************************************************
void key_callback(GLFWwindow* window, int key, int scancode, int action, int mods) {
    static const double Pi = 3.1415926535f;
    if (action == GLFW_RELEASE) {
        return;			// Ignore key up (key release) events
    }
    switch (key) {
    case GLFW_KEY_ESCAPE:
        glfwSetWindowShouldClose(window, true);
        return;
    case 'G':
            spinMode = true; // make the animation happen when you press g
            currentTime = 0.0;
            return;
    case 'R':
        if (singleStep) {			// If ending single step mode
            singleStep = false;
            spinMode = true;		// Restart animation
        }
        else {
            spinMode = !spinMode;	// Toggle animation on and off.
        }
        return;
    case 'S':
        singleStep = true;
        spinMode = true;
        return;
    case 'W':		// Toggle wireframe mode
        if (wireframeMode) {
            wireframeMode = false;
            glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
        }
        else {
            wireframeMode = true;
            glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
        }
        return;
    case 'C':		// Toggle backface culling
        cullBackFaces = !cullBackFaces;     // Negate truth value of cullBackFaces
        if (cullBackFaces) {
            glEnable(GL_CULL_FACE);
        }
        else {
            glDisable(GL_CULL_FACE);
        }
        return;
    case 'M':
        if (mods & GLFW_MOD_SHIFT) {            
            meshRes = meshRes<79 ? meshRes+1 : 80;  // Uppercase 'M'
        }
        else {
            meshRes = meshRes>4 ? meshRes-1 : 3;    // Lowercase 'm'
        }
        MyRemeshGeometries();
        return;
    case 'F':
        if (mods & GLFW_MOD_SHIFT) {                // If upper case 'F'
            animateIncrement *= sqrt(2.0);			// Double the animation time step after two key presses
        }
        else {                                      // Else lose case 'f',
            animateIncrement *= sqrt(0.5);			// Halve the animation time step after two key presses
        }
        return;
    case GLFW_KEY_UP:
        viewAzimuth = Min(viewAzimuth + 0.01, PIhalves - 0.05);
        break;
    case GLFW_KEY_DOWN:
        viewAzimuth = Max(viewAzimuth - 0.01, -PIhalves + 0.05);
        break;
    case GLFW_KEY_RIGHT:
        viewDirection += 0.01;
        if (viewDirection > PI) {
            viewDirection -= PI2;
        }
        break;
    case GLFW_KEY_LEFT:
        viewDirection -= 0.01;
        if (viewDirection < -PI) {
            viewDirection += PI2;
        }
        break;
    case GLFW_KEY_HOME:
        extraViewDistance = Min(4.0, extraViewDistance + 0.2);
        break;
    case GLFW_KEY_END:
        extraViewDistance = Max(-8.0, extraViewDistance - 0.2);
        break;
    }
    mySetViewMatrix();
}


// *************************************************
// This function is called with the graphics window is first created,
//    and again whenever it is resized.
// The Projection View Matrix is typically set here.
//    But this program does not use any transformations or matrices.
// *************************************************
void window_size_callback(GLFWwindow* window, int width, int height) {
	// Define the portion of the window used for OpenGL rendering.
	glViewport(0, 0, width, height);

	// Setup the projection matrix as a perspective view.
	// The complication is that the aspect ratio of the window may not match the
	//		aspect ratio of the scene we want to view.
	double w = (width == 0) ? 1.0 : (double)width;
	double h = (height == 0) ? 1.0 : (double)height;
	double windowXmax, windowYmax;
    double aspectFactor = w * Ymax / (h * Xmax);   // == (w/h)/(Xmax/Ymax), ratio of aspect ratios
	if (aspectFactor>1) {
		windowXmax = Xmax * aspectFactor;
		windowYmax = Ymax;
	}
	else {
		windowYmax = Ymax / aspectFactor;
		windowXmax = Xmax;
	}

	// Using the max & min values for x & y & z that should be visible in the window,
	//		we set up the orthographic projection.
    double zFar = zNear + Zmax - Zmin + 8.0f;
	theProjectionMatrix.Set_glFrustum(-windowXmax, windowXmax, -windowYmax, windowYmax, zNear, zFar);

	if (glIsProgram(shaderProgram1)) {
		glUseProgram(shaderProgram1);
		theProjectionMatrix.DumpByColumns(matEntries);
		glUniformMatrix4fv(projMatLocation, 1, false, matEntries);
	}
	check_for_opengl_errors();   // Really a great idea to check for errors -- esp. good for debugging!
}

void my_setup_OpenGL() {
	
	glEnable(GL_DEPTH_TEST);	// Enable depth buffering
	glDepthFunc(GL_LEQUAL);		// Useful for multipass shaders

	// Set polygon drawing mode for front and back of each polygon
    glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);

    glEnable(GL_CULL_FACE);

	check_for_opengl_errors();   // Really a great idea to check for errors -- esp. good for debugging!
}

void error_callback(int error, const char* description)
{
	// Print error
	fputs(description, stderr);
}

void setup_callbacks(GLFWwindow* window) {
	// Set callback function for resizing the window
	glfwSetFramebufferSizeCallback(window, window_size_callback);

	// Set callback for key up/down/repeat events
	glfwSetKeyCallback(window, key_callback);

	// Set callbacks for mouse movement (cursor position) and mouse botton up/down events.
	// glfwSetCursorPosCallback(window, cursor_pos_callback);
	// glfwSetMouseButtonCallback(window, mouse_button_callback);
}

int main() {
	glfwSetErrorCallback(error_callback);	// Supposed to be called in event of errors. (doesn't work?)
	glfwInit();
#if defined(__APPLE__) || defined(__linux__)
	glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
	glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
	glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
	glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
#endif
	const int initWidth = 800;
	const int initHeight = 600;
	GLFWwindow* window = glfwCreateWindow(initWidth, initHeight, "Letter (Initial) Demo", NULL, NULL);
	if (window == NULL) {
		printf("Failed to create GLFW window!\n");
		return -1;
	}
	glfwMakeContextCurrent(window);

	if (GLEW_OK != glewInit()) {
		printf("Failed to initialize GLEW!.\n");
		return -1;
	}

	// Print info of GPU and supported OpenGL version
	printf("Renderer: %s\n", glGetString(GL_RENDERER));
	printf("OpenGL version supported %s\n", glGetString(GL_VERSION));
#ifdef GL_SHADING_LANGUAGE_VERSION
	printf("Supported GLSL version is %s.\n", (char *)glGetString(GL_SHADING_LANGUAGE_VERSION));
#endif
    printf("Using GLEW version %s.\n", glewGetString(GLEW_VERSION));

	printf("------------------------------\n");
    printf("Use Arrow Keys to adjust the view angle.\n   ");
    printf("Use Home/End to adjust view distance.\n");
    printf("Press 'r' or 'R' (Run) to toggle(off and on) running the animation.\n");
	printf("Press 's' or 'S' to single - step animation.\n");
    printf("Press 'c' or 'C' (cull) to toggle whether back faces are culled.\n");
    printf("Press 'w' or 'W' (wireframe) to toggle whether wireframe or fill mode.\n");
    printf("Press 'M' (mesh) to increase the mesh resolution.\n");
    printf("Press 'm' (mesh) to decrease the mesh resolution.\n");
    printf("Press 'F'(faster) or 'f' (slower) to speed up or slow down the animation.\n");
    printf("Press ESCAPE to exit.\n");
	
    setup_callbacks(window);
   
	// Initialize OpenGL, the scene and the shaders
    my_setup_OpenGL();
	my_setup_SceneData();
 	window_size_callback(window, initWidth, initHeight);

    // Loop while program is not terminated.
	while (!glfwWindowShouldClose(window)) {
	
		myRenderScene();				// Render into the current buffer
		glfwSwapBuffers(window);		// Displays what was just rendered (using double buffering).

		// Poll events (key presses, mouse events)
		glfwWaitEventsTimeout(1.0/60.0);	    // Use this to animate at 60 frames/sec (timing is NOT reliable)
		// glfwWaitEvents();					// Or, Use this instead if no animation.
		// glfwPollEvents();					// Use this version when animating as fast as possible
	}

	glfwTerminate();
	return 0;
}

// If an error is found, it could have been caused by any command since the
//   previous call to check_for_opengl_errors()
// To find what generated the error, you can try adding more calls to
//   check_for_opengl_errors().
char errNames[8][36] = {
	"Unknown OpenGL error",
	"GL_INVALID_ENUM", "GL_INVALID_VALUE", "GL_INVALID_OPERATION",
	"GL_INVALID_FRAMEBUFFER_OPERATION", "GL_OUT_OF_MEMORY",
	"GL_STACK_UNDERFLOW", "GL_STACK_OVERFLOW" };
bool check_for_opengl_errors() {
	int numErrors = 0;
	GLenum err;
	while ((err = glGetError()) != GL_NO_ERROR) {
		numErrors++;
		int errNum = 0;
		switch (err) {
		case GL_INVALID_ENUM:
			errNum = 1;
			break;
		case GL_INVALID_VALUE:
			errNum = 2;
			break;
		case GL_INVALID_OPERATION:
			errNum = 3;
			break;
		case GL_INVALID_FRAMEBUFFER_OPERATION:
			errNum = 4;
			break;
		case GL_OUT_OF_MEMORY:
			errNum = 5;
			break;
		case GL_STACK_UNDERFLOW:
			errNum = 6;
			break;
		case GL_STACK_OVERFLOW:
			errNum = 7;
			break;
		}
		printf("OpenGL ERROR: %s.\n", errNames[errNum]);
	}
    if (numErrors != 0) {
        assert("OpenGL error: Check console window for error messages!" && false);
    }
	return (numErrors != 0);
}