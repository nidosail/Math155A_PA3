/*
* GlGeomBezier.h - Version 0.3 - June 5, 2020
*
* C++ class for rendering Bezier patches in Modern OpenGL.
*   A GlGeomBezier object encapsulates a VAO, VBO, and VEO,
*   which can be used to render a sphere.
*   The number of slices and stacks can be varied.
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

#ifndef GLGEOM_BEZIER_H
#define GLGEOM_BEZIER_H

#include <assert.h>

#include "GlGeomBase.h"

// GlGeomBezier
//     Generates vertices, normals, and texture coordinates for Bezier patches.
//     Can be remeshed dynamically.
//     Stores and renders multiple Bezier patches at once.
//     Each patch must have the same uOrder and vOrder, and the same resolutions
// Supports:
//    (1) Allocating and loading a VAO, VBO, and EBO
//    (2) Rendering the Bezier patches with OpenGL.
// How to use:
//     * First call the constructor GlGeomBezier() and possibly ReMesh()
//             to set the numbers of mesh resolutions.
//             These numbers can be changed by calling ReMesh().
//     * Then call LoadControlPts() to load all the control points
//             for all the Bezier patches at once.
//     * Then call InitializeAttribLocations() to specify whether to use
//          normals and texture coordinates and to specify
//          the locations in the VBO buffer for the shader program.
//          The also allocates and loads the VAO, VBO and EBO.
//     * Call Render() to render all the Bezier patches
//     * Call RenderPatch(i) to render just the i-th patch.
//     * The routines Render(...) issues the the glDrawElements commands 
//       for the Bezier patch(es) using the VAO, VBO and EBO.

class GlGeomBezier : public GlGeomBase
{
public:
    GlGeomBezier() : GlGeomBezier(8, 8) {}
    GlGeomBezier( int uMeshResolution, int vMeshResolution);
    ~GlGeomBezier();

    // Load, all at once, all the control points for all the Bezier patches.
    // uOrder, vOrder: the order of the Bezier patch in u and v directions.
    //     "order" is the same as "degree + 1"
    //     A single patch has uOrder*vOrder many control points.
    // numCoordinates - number of coordinates per control point.
    //     Use "3" for (x,y,z) values in 3-space.  
    //     Use "4" for (x,y,z,w) values, homogeneous representations of points in 3-space.
    //     When "4" is used: it is only for the control points
    //           The VBO still stores (x,y,z) only. (Shader program should use a vec3.)
    // controlPoints: an array of numPatches*uOrder*vOrder*numCoordinates many doubles,
    //     specifying numPatches*uOrder*vOrder many vertices.
    // LoadControlPts makes a copy of the control points in case ReMeshing is needed.
    //     It is possible to pass in a null pointer, and then just the parameters are saved (For future compatibility.)
    void LoadControlPts(int uOrder, int vOrder, int numCoordinates, int numPatches, double* controlPoints);

    // Disable all copy and assignment operators for a GlGeomBezier object.
    //     If you need to pass it to/from a function, use references or pointers
    //     and be sure that there are no implicit copy or assignment operations!
    GlGeomBezier(const GlGeomBezier&) = delete;
    GlGeomBezier& operator=(const GlGeomBezier&) = delete;
    GlGeomBezier(GlGeomBezier&&) = delete;
    GlGeomBezier& operator=(GlGeomBezier&&) = delete;

public:

    // Remesh: re-mesh to change the mesh resolution.
    // Can be called either before or after InitializeAttribLocations(), but it is
    //    more efficient if Remesh() is called first, or if the constructor sets the mesh resolution;
    //    and InitializeAttribLocations() is called afterwards.
    void Remesh(int uMeshResolution, int vMeshResolution);

    // Allocate the VAO, VBO, and EBO.
    // Set up info about the Vertex Attribute Locations
    // LoadControlPts(...) must be called before InitializeAttribLocations.
    // This must be called before Render() or RenderPatch() is first called.
    // First parameter is the location for the vertex position vector in the shader program.
    // Second parameter is the location for the vertex normal vector in the shader program.
    // Third parameter is the location for the vertex 2D texture coordinates in the shader program.
    // The second and third parameters are optional (use UINT_MAX to omit).
    void InitializeAttribLocations(
        unsigned int pos_loc, unsigned int normal_loc = UINT_MAX, unsigned int texcoords_loc = UINT_MAX);

    // Render the Bezier patches.  Must call InitializeAttribLocations first.
    void Render();                      // Render all the patches
    void RenderPatch(int i);            // Render the i-th patch only

    // GetNumElementsMax() returns the maximum number of elements in the EBO
    // GetNumElementsRender() returns the actual number of elements in the EBO (for rendering)
    // GetNumVerticesTexCoords() returns the number of vertices when there are texture coordinates
    // GetNumVerticesNoTexCoords() returns the number of vertices when no texture coordinates
    int GetNumElementsMax() const { return numPatches * 6 * uMeshRes * vMeshRes; }
    int GetNumElementsRender() const { return 3*GetFirstTriInPatch(numPatches); }
    int GetNumVerticesTexCoords() const { return (numPatches * (uMeshRes + 1) * (vMeshRes + 1)); }
    int GetNumVerticesNoTexCoords() const { return GetNumVerticesTexCoords(); }

    int GetuMeshRes() const { return uMeshRes; }
    int GetvMeshRes() const { return vMeshRes; }
    int GetNumPatches() const { return numPatches; }

    // CalcVboAndEbo- return all VBO vertex information, and EBO elements for GL_TRIANGLES drawing.
    // See GlGeomBase.h for additional information
    void CalcVboAndEbo(float* VBOdataBuffer, unsigned int* EBOdataBuffer,
                        int vertPosOffset, int vertNormalOffset, int vertTexCoordsOffset, 
                        unsigned int stride);

protected:
    void BezierMultiEval(
        double* cntlPts, int stride, 
        int order, int alphaNumerator, int alphaDenominator, 
        double* destVal, double* destDeriv = 0);

private:
    int uMeshRes;   // Mesh resolution in the "u" direction
    int vMeshRes;   // Mesh resolution in the "v" direction

    int uOrder = 0;             // Bezier order in u direction (order = degree+1)
    int vOrder = 0;             // Bezier order in v direction (order = degree+1)
    int numPatches = 0;         // Number of patches loaded
    int numCoordinates = 0;     // Should equal 3, or 4 (E.g.., ordinary coordinates or homogenenous coordinates)
    double* controlPts = 0;     // Has numPatches*uOrder*vOrder*numCoordinates values (if non-zero)

    int* firstTriInPatch = 0;   // firstTriInPatch[i] = index of first triangle in the i-th patch

    int GetFirstTriInPatch(int i) const { assert(i <= numPatches); return firstTriInPatch[i]; }
    int GetNumTrisInPatch(int i) const;
    bool EqualVerts(int i, int j, float* vboVertPtr, int stride) const;

    bool VboEboLoaded = false;

    void PreRender();

    void CalcCornerNormals(const double* controlPointsPtr, double* retCornerNormals, double* retNearZeroSq);
    int CalcOneCornerNormal(const double* cornerPtr, int uStride, int vStride, double dest[3], double nearZeroSq);
    double MaxAbsEntry(double* cntlPt);

    // Basic vector functions (here just to avoid bringing in a vector package such as LinearMapR3)
    static void VecCopy(const double* source, double* retCopy, int numComponents);
    static void VecDiff(const double* minuhend, const double* subtrahend, double* retDifference, int numComponents);
    static void VecCrossProd(const double* aVec, const double* bVec, double* retCrossProd);
    static void VecMultScalar(const double* multiplicand, double multiplier, double* retProduct, int numComponents);
    static void VecAddScaled(const double* multiplicand, double multiplier, double* retAccum, int numComponents);
    static void VecSetZero(double* retZero, int numComponents);
    static double VecDotProd(const double* vecA, const double* vecB, int numComponents);
    static double VecNormSq(const double* vec, int numComponents);
    static bool VecIsZero(const double* vec, int numComponents);
};

// Constructor
inline GlGeomBezier::GlGeomBezier(int uMeshResolution, int vMeshResolution)
{
    // Set initial mesh resolutions
    uMeshRes = uMeshResolution;
    vMeshRes = vMeshResolution;
}

// Destructor
inline GlGeomBezier::~GlGeomBezier()
{
    delete[] controlPts;
}

inline int GlGeomBezier::GetNumTrisInPatch(int i) const 
{
    assert(i < numPatches);
    return  firstTriInPatch[i + 1] - firstTriInPatch[i];
}


#endif  // GLGEOM_BEZIER_H