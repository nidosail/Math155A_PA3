/*
* GlGeomBezier.cpp - Version 0.3 - June 5, 2020
*
* C++ class for rendering bezier patches in Modern OpenGL.
*   A GlGeomBezier object encapsulates a VAO, a VBO, and an EBO,
*   which can be used to render a set of Bezier patches.
*   The u- and v-mesh resolution can be varied.
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

// Use the static library (so glew32.dll is not needed):
#define GLEW_STATIC
#include <GL/glew.h> 
#include <GLFW/glfw3.h>

#include "GlGeomBezier.h"
#include "MathMisc.h"
#include "assert.h"

// ************************************
// LoadControlPts
//      Load, all at once, all the control points for all the Bezier patches.
//      See header file for information on the parameters
void GlGeomBezier::LoadControlPts(int uOrder, int vOrder,
                                    int numCoordinates, int numPatches,
                                    double* controlPoints)
{
    assert(uOrder > 0 && vOrder > 0);
    assert(numCoordinates == 3 || numCoordinates <= 4);           // Points in R^3 or homogeneous representations
    int numVals = uOrder * vOrder * numCoordinates * numPatches;  // Number of floats for control points
    if (controlPoints == 0) {
        delete[] controlPts;      // If no control points are passed in
        controlPts = 0;
        delete[] firstTriInPatch;
        firstTriInPatch = 0;
        numVals = 0;            // Disable copying of values below
    }
    else if (controlPts != 0) {   // If control points already allocated
        int oldNumVals = this->uOrder * this->vOrder * this->numCoordinates * this->numPatches;
        if (numVals != oldNumVals) {
            delete[] controlPts; // Delete to re-allocate (potentially inefficient)
            controlPts = new double[numVals];
        }
        if (numPatches != this->numPatches) {
            delete[] firstTriInPatch;
            firstTriInPatch = new int[numPatches + 1];
        }
    }
    else {
        assert(controlPts == 0);
        assert(firstTriInPatch == 0);
        controlPts = new double[numVals];
        firstTriInPatch = new int[numPatches + 1];
    }
    this->uOrder = uOrder;
    this->vOrder = vOrder;
    this->numCoordinates = numCoordinates;
    this->numPatches = numPatches;

    for (int i = 0; i < numVals; i++) {
        *(controlPts + i) = *(controlPoints + i);
    }
}

void GlGeomBezier::Remesh(int uMeshResolution, int vMeshResolution)
{
    if (uMeshResolution == uMeshRes && vMeshResolution == vMeshRes) {
        return;
    }
    uMeshRes = uMeshResolution;
    vMeshRes = vMeshResolution;
    VboEboLoaded = false;
}

void GlGeomBezier::CalcVboAndEbo(float* VBOdataBuffer, unsigned int* EBOdataBuffer,
    int vertPosOffset, int vertNormalOffset, int vertTexCoordsOffset, unsigned int stride)
{
    assert(vertPosOffset >= 0 && stride > 0);
    assert(uOrder > 1 && vOrder > 1 && numPatches > 0 && controlPts != 0);
    bool calcNormals = (vertNormalOffset >= 0);       // Should normals be calculated?
    bool calcTexCoords = (vertTexCoordsOffset >= 0);  // Should texture coordinates be calculated?

    int numCntlPtsEntriesJ = uOrder * (vMeshRes+1) * numCoordinates;   // Size for J-slices control points
    double* sliceCntlPtsJ = new double[numCntlPtsEntriesJ];
    int numCntlPtsEntriesI = (uMeshRes+1) * vOrder * numCoordinates;   // Size for I-slices control points
    double* sliceCntlPtsI = (calcNormals ? new double[numCntlPtsEntriesI] : (double*)0);
    int numValuesInPatch = uOrder * vOrder * numCoordinates;
    for (int patchNum = 0; patchNum < numPatches; patchNum++) {
        // Calculate vertices for patch number patchNum
        double* patchPtr = controlPts + patchNum * numValuesInPatch;            // Pointer to patch's control points
        double nearZeroSq;              // Threshhold at which a "normal" is considered be null

        // First, calculate the control points for the "horizontal" paths (with j constant)
        for (int j = 0; j <= vMeshRes; j++) {
            for (int i = 0; i < uOrder; i++) {
                // Calculate i-th control point for j-th horizontal slice
                BezierMultiEval(patchPtr + i * numCoordinates, uOrder * numCoordinates,
                    vOrder, j, vMeshRes,
                    sliceCntlPtsJ + (j * uOrder + i) * numCoordinates);
            }
        }
        // Second, if needed for the computation of normal vectors,
        //    calculate the control points for the "vertical" paths (with i constant)
        if (calcNormals) {
            for (int i = 0; i <= uMeshRes; i++) {
                for (int j = 0; j < vOrder; j++) {
                    // Calculate j-th control point for i-th horizontal slice
                    BezierMultiEval(patchPtr + j * uOrder * numCoordinates, numCoordinates,
                        uOrder, i, uMeshRes,
                        sliceCntlPtsI + (i * vOrder + j) * numCoordinates);
                }
            }
        }
        // Third, if needed, calculate normals at corners, including certain degenerate cases.
        float* vboSubbufferPtr = VBOdataBuffer + stride * (patchNum * (uMeshRes + 1) * (vMeshRes + 1));
        if (calcNormals) {
            // Precompute the normals for the four corner vertices
            //    with special coding to handle degenerate cases
            // Then copy them into the VBO buffer
            double cornerNormals[4][3];     // Precomputed normals at the four corners
            CalcCornerNormals(patchPtr, &cornerNormals[0][0], &nearZeroSq);
            for (int i = 0; i <= 1; i++) {
                for (int j = 0; j <= 1; j ++) {
                    float* vboPtr = vboSubbufferPtr + stride * (j * vMeshRes * (uMeshRes + 1) + i * uMeshRes);
                    float* normalPtr = vboPtr + vertNormalOffset;
                    int cornerIdx = 2*j + i;
                    *(normalPtr++) = (float)cornerNormals[cornerIdx][0];
                    *(normalPtr++) = (float)cornerNormals[cornerIdx][1];
                    *normalPtr = (float)cornerNormals[cornerIdx][2];
                }
            }
        }

        // Fourth, generate coordinates and normals for each mesh vertex
        for (int i = 0; i <= uMeshRes; i++) {
            // Handle the (i,*) vertices
            for (int j = 0; j <= vMeshRes; j++) {
                // Handle vertex (i,j)
                float* vboPtr = vboSubbufferPtr + stride * (j * (uMeshRes+1) + i);
                double tempPositionJ[4];            // Max numCoordinates is 4
                double tempDerivativeJ[4];          // Partial wrt i (ie., j is fixed.
                BezierMultiEval(sliceCntlPtsJ + (j * uOrder * numCoordinates), numCoordinates,
                    uOrder, i, uMeshRes,
                    tempPositionJ, tempDerivativeJ);
                if (numCoordinates == 4) {
                    // Divide by w component to convert to R^3 position
                    VecMultScalar(tempPositionJ, 1.0 / tempPositionJ[3], tempPositionJ, 3);
                }
                // Copy the three coordinates (do not use homogeneous representation in VBO)
                for (int t = 0; t < 3; t++) {
                    *(vboPtr + vertPosOffset + t) = (float)tempPositionJ[t];
                }
                if (calcNormals) {
                    float* normalPtr = vboPtr + vertNormalOffset;
                    double normalNormSq;
                    if ((i == 0 || i == uMeshRes) && (j == 0 || j == vMeshRes)) {
                        // At a corner, do nothing
                    }
                    else {
                        // Not at a corner
                        double tempPositionI[4];        // Max numCoordinates is 4
                        double tempDerivativeI[4];      // Partial wrt j (ie., i is fixed).
                        BezierMultiEval(sliceCntlPtsI + i * vOrder * numCoordinates, numCoordinates,
                            vOrder, j, vMeshRes, tempPositionI, tempDerivativeI);
                        // tempPositionI is equal to tempPositionJ (except not divided by w, since not used)
                        // tempDerivativeI is the partial wrt j (ie., i is fixed).
                        // tempDerivativeJ is the partial wrt i (ie., j is fixed.
                        if (numCoordinates == 4) {
                            // (x/w)' = (x'*w -x*w')/w^2 = (x' -(x/w)*w')/w
                            // (x/w) is equal to tempPositionI[0] and to tempPositionJ[0]. Similarly for (y/w) and (z/w)
                            // x', y', z', w' are in tempDerivativeI[] and tempDerivativeJ[], respectively for wrt j, i
                            // We skip dividing by w since the result is normalized anyway.
                            VecAddScaled(tempPositionJ, -tempDerivativeI[3], tempDerivativeI, 3);
                            VecAddScaled(tempPositionJ, -tempDerivativeJ[3], tempDerivativeJ, 3);
                        }
                        double normalVec[3];
                        VecCrossProd(tempDerivativeJ, tempDerivativeI, normalVec);
                        // Check if the calculated normal was (close to) zero
                        normalNormSq = VecNormSq(normalVec, 3);
                        if (normalNormSq < nearZeroSq) {
                            if (i == 0 || i == uMeshRes) {
                                if (VecNormSq(tempDerivativeI, 3) < nearZeroSq) {
                                    // This calculation valid if the edge is degenerate. (not checked)
                                    // It is also valid if the if second partial w.r.t. v is also zero.
                                    // Go to neighboring row, ii=1 or ii = uOrder-2. Find partial along control polygon
                                    int ii = (i == 0) ? 1 : uOrder - 2;
                                    double tempPositionI2[4];        // Values from neighboring row of origial control point
                                    double tempDerivativeI2[4];      // Partial wrt j (ie., ii is fixed).
                                    BezierMultiEval(patchPtr + ii * numCoordinates, uOrder * numCoordinates,
                                        vOrder, j, vMeshRes,
                                        tempPositionI2, tempDerivativeI2);
                                    if (numCoordinates == 4) {
                                        VecMultScalar(tempPositionI2, 1.0 / tempPositionI2[3], tempPositionI2, 3);
                                        VecAddScaled(tempPositionI2, -tempDerivativeI2[3], tempDerivativeI2, 3);
                                    }
                                    VecCrossProd(tempDerivativeJ, tempDerivativeI2, normalVec);
                                }
                            }
                            else if (j == 0 || j == vMeshRes) {
                                if (VecNormSq(tempDerivativeJ, 3) < nearZeroSq) {
                                    // This calculation valid if the edge is degenerate. (not checked)
                                    // It is also valid if the if second partial w.r.t. v is also zero.
                                    // Go to neighboring column, j=1 or j = vOrder-2. Find partial along control polygon
                                    int jj = (j == 0) ? 1 : vOrder - 2;
                                    double tempPositionJ2[4];        // Values from neighboring column
                                    double tempDerivativeJ2[4];      // Partial wrt i (ie., jj is fixed).
                                    BezierMultiEval(patchPtr + jj * uOrder * numCoordinates, numCoordinates,
                                        uOrder, i, uMeshRes,
                                        tempPositionJ2, tempDerivativeJ2);
                                    if (numCoordinates == 4) {
                                        VecMultScalar(tempPositionJ2, 1.0 / tempPositionJ2[3], tempPositionJ2, 3);
                                        VecAddScaled(tempPositionJ2, -tempDerivativeJ2[3], tempDerivativeJ2, 3);
                                    }
                                    VecCrossProd(tempDerivativeJ2, tempDerivativeI, normalVec);
                                }
                            }
                            normalNormSq = VecNormSq(normalVec, 3);
                        }
                        // End of special edge calculation for degenerate edge
                        double normSqInv = 1.0 / sqrt(normalNormSq);
                        *(normalPtr++) = (float)(normalVec[0] * normSqInv);
                        *(normalPtr++) = (float)(normalVec[1] * normSqInv);
                        *normalPtr = (float)(normalVec[2] * normSqInv);
                    }
                }
                // Calculate texture coordinate
                if (calcTexCoords) {
                    float* tcPtr = vboPtr + vertTexCoordsOffset;
                    *tcPtr = (float)i / (float)uMeshRes;
                    *(tcPtr + 1) = (float)j / (float)vMeshRes;
                }
            }
        }
    }
    delete[] sliceCntlPtsI;
    delete[] sliceCntlPtsJ;
 
    // Set up the element array. Somewhat wastefully, each triangle has its own
    // entries in the EBO, so as to fit the framework used by GlGeomBase.
    firstTriInPatch[0] = 0;
    int idx = 0;
    unsigned int* eboPtr = EBOdataBuffer;
    for (int patchNum = 0; patchNum < numPatches; patchNum++) {
        for (int j = 0; j < vMeshRes; j++) {
            for (int i = 0; i < uMeshRes; i++) {
                if (!(EqualVerts(idx, idx + 1, VBOdataBuffer+vertPosOffset, stride) 
                        || EqualVerts(idx, idx + (uMeshRes + 1), VBOdataBuffer + vertPosOffset, stride))) {
                    *(eboPtr++) = idx;
                    *(eboPtr++) = idx + 1;
                    *(eboPtr++) = idx + (uMeshRes + 1);
                }
                if (!(EqualVerts(idx + (uMeshRes + 1), idx + (uMeshRes + 2), VBOdataBuffer + vertPosOffset, stride)
                        || EqualVerts(idx + (uMeshRes + 2), idx + 1, VBOdataBuffer + vertPosOffset, stride))) {
                    *(eboPtr++) = idx + (uMeshRes + 1);
                    *(eboPtr++) = idx + 1;
                    *(eboPtr++) = idx + (uMeshRes + 2);
                }
                idx++;
            }
            idx++;
        }
        idx += uMeshRes+1;
        firstTriInPatch[patchNum + 1] = (int)((eboPtr - EBOdataBuffer) / 3);
    }

}

// retCornerNormals is a pointer to an array where the
//    four normal vectors are returned.
void GlGeomBezier::CalcCornerNormals(const double* controlPointsPtr, 
                                     double* retCornerNormals, double* retNearZeroSq)
{
    // Fix a threshhold for deeming a normal to be zero (up to roundoff error)
    //   Any vector with NormSq < nearZeroSq will be deemed to be non-zero.
    double allNormSq = VecNormSq(controlPointsPtr, uOrder * vOrder * numCoordinates);
    double nearZeroSq = 1.0e-25 * allNormSq / (double)(uOrder * vOrder);
    *retNearZeroSq = nearZeroSq;

    // Find the normals at the four corners
    for (int i = 0; i <= 1; i++) {
        // i is the value of u
        int uStride = (1 - 2 * i) * numCoordinates;
        for (int j = 0; j <= 1; j++) {
            // j is the value of v
            int vStride = (1 - 2 * j) * uOrder * numCoordinates;
            int cornerIdx = i * (uOrder - 1) + j * uOrder * (vOrder - 1);
            const double* cornerPtr = controlPointsPtr + cornerIdx * numCoordinates;
            double* dest = retCornerNormals + (2 * j + i) * 3;
            CalcOneCornerNormal(cornerPtr, uStride, vStride, dest, nearZeroSq);
            if (i + j == 1) {
                VecMultScalar(dest, -1.0, dest, 3);  // Flip sign to correct for crossproduct direction
            }
        }
    }
}

// Attempt to calculate the normal at a corner vertex using a crossproduct.
//   cornerIdx is the index to the corner vertex.  
//   uStride and vStride give displacements to the adjacent vertices.
// Answer is returned in dest[3].
// Handles both ordinary coordinates and homogeneous coordinates
//   and adjacent vertices which are points at infinity
// Returns 1 if successful with crossproduct method.
// Returns 2 if successful with special second order methods.
// Returns 0 if unsuccessful.
int GlGeomBezier::CalcOneCornerNormal(const double* cornerPtr, int uStride, int vStride, 
                                      double dest[3], double nearZeroSq)
{
    if (numCoordinates == 3) {
        double P10minusP00[3];
        double P01minusP00[3];
        VecDiff(cornerPtr + uStride, cornerPtr, P10minusP00, 3);
        VecDiff(cornerPtr + vStride, cornerPtr, P01minusP00, 3);
        double crossProduct[3];
        VecCrossProd(P10minusP00, P01minusP00, crossProduct);
        double cpMagSq = VecNormSq(crossProduct, 3);
        if (cpMagSq > nearZeroSq) {
            // Use the crossproduct as the normal vector (after normalization)
            VecMultScalar(crossProduct, 1.0 / sqrt(cpMagSq), &dest[0], 3);
            return 1;   
        }
        // Use the "special" method
        // First: Compute partialUV, equal to 
        //             (partialU)(partialV)(0,0) / (uOrder*vOrder).
        double P11minusP01[3];
        VecDiff(cornerPtr + uStride + vStride, cornerPtr + vStride,
                P11minusP01, 3);
        double partialUV[3]; // (partial U)(partial V) / (uOrder*vOrder)
        VecDiff(P11minusP01, P10minusP00, partialUV, 3);
        VecMultScalar(partialUV, (uOrder - 1) * (vOrder - 1), partialUV, 3);
        // Second, check if either of the two partials partialU and partialV are zero
        // If so, crossproduct nonzero one with partialUV to obtain the normal
        double partialU[3];
        double partialV[3];
        VecMultScalar(P10minusP00, uOrder - 1, partialU, 3);
        VecMultScalar(P01minusP00, vOrder - 1, partialV, 3);
        double partialUmagSq = VecNormSq(partialU, 3);
        double partialVmagSq = VecNormSq(partialV, 3);
        if (partialUmagSq < nearZeroSq) {
            VecCrossProd(partialUV, partialV, crossProduct);
        }
        else if (partialVmagSq < nearZeroSq) {
            VecCrossProd(partialU, partialUV, crossProduct);
        }
        else {
            // Both partials are non-zero
            // Use crossproduct with a directional second derivative
            double partialUdotPartialV = VecDotProd(partialU, partialV, 3);
            double alpha = 1.0, beta = 1.0;
            if (partialUmagSq < partialVmagSq) {
                beta = -alpha * partialUdotPartialV / partialVmagSq;
            }
            else {
                alpha = -beta * partialUdotPartialV / partialUmagSq;
            }
            double partialUU[3] = { 0.0, 0.0, 0.0 };  // Value zero in case k_u==2
            double partialVV[3] = { 0.0, 0.0, 0.0 };  // Value zero in case k_v==2
            if (uOrder > 2) {
                // Set partialUU, second derivative w.r.t. u.
                double P20minusP10[3];
                VecDiff(cornerPtr + 2 * uStride, cornerPtr + uStride, P20minusP10, 3);
                VecDiff(P20minusP10, P10minusP00, partialUU, 3);
                VecMultScalar(partialUU, (uOrder - 1) * (uOrder - 2), partialUU, 3);
            }
            if (vOrder > 2) {
                // Set partialVV, second derivative w.r.t. v.
                double P02minusP01[3];
                VecDiff(cornerPtr + 2 * vStride, cornerPtr + vStride, P02minusP01, 3);
                VecDiff(P02minusP01, P01minusP00, partialVV, 3);
                VecMultScalar(partialVV, (vOrder - 1) * (vOrder - 2), partialVV, 3);
            }
            // Calculate second derivative in the direction D = alpha*partialU + beta*partialV.
            double partialDD[3];
            VecMultScalar(partialUU, alpha * alpha, partialDD, 3);
            VecAddScaled(partialVV, beta * beta, partialDD, 3);
            VecAddScaled(partialUV, 2 * alpha * beta, partialDD, 3);
            if (partialUmagSq < partialVmagSq) {
                VecCrossProd(partialDD, partialV, crossProduct);
            }
            else {
                VecCrossProd(partialU, partialDD, crossProduct);
            }
        }
        cpMagSq = VecNormSq(crossProduct, 3); 
        // If normal is still zero, leave as is.
        double cpMagInv = (cpMagSq > nearZeroSq) ? 1.0 / sqrt(cpMagSq) : 1.0;
        VecMultScalar(crossProduct, cpMagInv, &dest[0], 3);
        return 2;
    }
    else {
        // numCoordinates==4 - homogeneous coordinates.
        const double* P00 = cornerPtr;
        // partialUhg, partialVhg - partial derivatives w.r.t. u and v -- for homogeneous coordinates
        double partialUhg[4];
        double partialVhg[4];
        double P10minusP00hg[4];
        double P01minusP00hg[4];
        VecDiff(cornerPtr + uStride, cornerPtr, P10minusP00hg, 4);
        VecMultScalar(P10minusP00hg, uOrder - 1, partialUhg, 4);
        VecDiff(cornerPtr + vStride, cornerPtr, P01minusP00hg, 4);
        VecMultScalar(P01minusP00hg, vOrder - 1, partialVhg, 4);
        // partialU and partialV are partials wrt u, v in R^3 
        double partialU[3];                               // Partial w.r.t. u --- in R^3
        double partialV[3];                               // Partial w.r.t. v --- in R^3
        double P00wInv = 1.0 / (*(P00 + 3));              // 1.0/P00.w
        double P00wInvSq = P00wInv* P00wInv;              // (1.0/P00.w)^2
        double WuOverWsq = *(partialUhg + 3) * P00wInvSq; // w_u/P00.w^2  (w_u =  partial w w.r.t u)
        double WvOverWsq = *(partialVhg + 3) * P00wInvSq; // w_v/P00.w^2  (w_v =  partial w w.r.t v)
        VecMultScalar(partialUhg, P00wInv, partialU, 3);
        VecAddScaled(P00, -WuOverWsq, partialU, 3);
        VecMultScalar(partialVhg, P00wInv, partialV, 3);
        VecAddScaled(P00, -WvOverWsq, partialV, 3);
        // Try (partial U) crossproduct (partial V)
        double crossProduct[3];
        VecCrossProd(partialU, partialV, crossProduct);
        double cpMagSq = VecNormSq(crossProduct, 3);
        if (cpMagSq > nearZeroSq) {
            // Use the crossproduct as the normal vector (after normalization)
            VecMultScalar(crossProduct, 1.0 / sqrt(cpMagSq), &dest[0], 3);
            return 1;
        }
        // Use the "special" method
        // First: Compute partialUVhg, equal to (partialU)(partialV)(0,0), for homogeneous coordinates
        //        and partialUV, equal to (partialU)(partialV)(0,0), for R^3 coordinates
        double partialUVhg[4]; // (partial U)(partial V) -- using homogeneous coordinates
        VecDiff(cornerPtr + uStride + vStride, cornerPtr + vStride, partialUVhg, 4);
        VecMultScalar(partialUVhg, (uOrder - 1)* (vOrder - 1), partialUVhg, 4);
        VecAddScaled(partialUhg, -(vOrder - 1), partialUVhg, 4);
        double partialUV[3];     // (partial U)(partial V) -- for R^3
        VecMultScalar(partialUVhg, P00wInv, partialUV, 3);
        VecAddScaled(partialUhg, -WvOverWsq, partialUV, 3);
        VecAddScaled(partialVhg, -WuOverWsq, partialUV, 3);
        VecAddScaled(P00, 2 * (*(P00 + 3)) * WuOverWsq * WvOverWsq - partialUVhg[3] * P00wInvSq, partialUV, 3);
        double partialUmagSq = VecNormSq(partialU, 3);
        double partialVmagSq = VecNormSq(partialV, 3);
        if (partialVmagSq < nearZeroSq) {
            VecCrossProd(partialU, partialUV, crossProduct);
        }
        else if (partialUmagSq < nearZeroSq) {
            VecCrossProd(partialUV, partialV, crossProduct);
        }
        else {
            // Both partials are non-zero
            // Use crossproduct with a directional second derivative
            double partialUdotPartialV = VecDotProd(partialU, partialV, 3);
            double alpha = 1.0, beta = 1.0;
            if (partialUmagSq < partialVmagSq) {
                beta = -alpha * partialUdotPartialV / partialVmagSq;
            }
            else {
                alpha = -beta * partialUdotPartialV / partialUmagSq;
            }
            double partialUUhg[4] = { 0.0, 0.0, 0.0, 0.0 };  // Value zero in case k_u==2
            double partialVVhg[4] = { 0.0, 0.0, 0.0, 0.0 };  // Value zero in case k_u==2
            if (uOrder > 2) {
                // Set partialUUhg, second derivative w.r.t. u. in homogeneous coordinates
                double P20minusP10hg[4];
                VecDiff(cornerPtr + 2 * uStride, cornerPtr + uStride, P20minusP10hg, 4);
                VecDiff(P20minusP10hg, P10minusP00hg, partialUUhg, 4);
                VecMultScalar(partialUUhg, (uOrder - 1) * (uOrder - 2), partialUUhg, 4);
            }
            if (vOrder > 2) {
                // Set partialVV, second derivative w.r.t. v., in homogeneous coordinates
                double P02minusP01hg[4];
                VecDiff(cornerPtr + 2 * vStride, cornerPtr + vStride, P02minusP01hg, 4);
                VecDiff(P02minusP01hg, P01minusP00hg, partialVVhg, 4);
                VecMultScalar(partialVVhg, (vOrder - 1) * (vOrder - 2), partialVVhg, 4);
            }
            double partialUU[3];    // Second derivative in R^3 coordinates
            double partialVV[3];    // Second derivative in R^3 coordinates
            VecMultScalar(partialUUhg, P00wInv, partialUU, 3);
            VecAddScaled(partialUhg, -2.0 * partialUhg[3] * P00wInvSq, partialUU, 3);
            VecAddScaled(P00, (2.0 * partialUhg[3] * partialUhg[3] * P00wInv - partialUUhg[3])* P00wInvSq, partialUU, 3);
            VecMultScalar(partialVVhg, P00wInv, partialVV, 3);
            VecAddScaled(partialVhg, -2.0 * partialVhg[3] * P00wInvSq, partialVV, 3);
            VecAddScaled(P00, (2.0 * partialVhg[3] * partialVhg[3] * P00wInv - partialVVhg[3])* P00wInvSq, partialVV, 3);
            // Calculate second derivative in the direction D = alpha*partialU + beta*partialV.
            double partialDD[3];
            VecMultScalar(partialUU, alpha * alpha, partialDD, 3);
            VecAddScaled(partialVV, beta * beta, partialDD, 3);
            VecAddScaled(partialUV, 2 * alpha * beta, partialDD, 3);
            if (partialUmagSq < partialVmagSq) {
                VecCrossProd(partialDD, partialV, crossProduct);
            }
            else {
                VecCrossProd(partialU, partialDD, crossProduct);
            }
        }
        cpMagSq = VecNormSq(crossProduct, 3);
        VecMultScalar(crossProduct, 1.0 / sqrt(cpMagSq), &dest[0], 3);
        return 2;
    }
}

double GlGeomBezier::MaxAbsEntry(double* cntlPt)
{
    double ret = 0.0;
    for (int i = 0; i < 3; i++) {
        double xyz = *(cntlPt + i);
        if (numCoordinates == 4) {
            xyz /= *(cntlPt + 3);
        }
        // ret = Max(ret, xyz, -xyz);
        ret = (ret < xyz) ? xyz : (ret < -xyz) ? -xyz : ret;
    }
    return ret;
}

bool GlGeomBezier::EqualVerts(int i, int j, float* vboVertPtr, int stride) const
{
    float* vertI = vboVertPtr + stride * i;
    float* vertJ = vboVertPtr + stride * j;
    return (*vertI == *vertJ && *(vertI + 1) == *(vertJ + 1) && *(vertI + 2) == *(vertJ + 2));
}


void GlGeomBezier::InitializeAttribLocations(
    unsigned int pos_loc, unsigned int normal_loc, unsigned int texcoords_loc)
{
    // The call to GlGeomBase::InitializeAttribLocations will further call
    //   GlGeomBezier::CalcVboAndEbo()

    GlGeomBase::InitializeAttribLocations(pos_loc, normal_loc, texcoords_loc);
    VboEboLoaded = true;
}

// Main computation for values and derivatives of Bezier curves of arbitrary degree
// Parameters:
//     cntlPts - pointer to the control points.
//               Each control point consists numCoordinates many floats
//     stride - stride between control points. Measured in floats.
//     order - The order (degree plus 1) of the Bezier curve.
//     alphaNumerator & alphaDenominator: divide to get value where curve is evaluated
//             They are specified separately to avoid round off error causing 
//             different results if a curve is traversed in both directions 
//             (i.e., with control points reversed)
//     destVal- Pointer to where the results are returned (numCoordinates many values)
//     destDeriv - Pointer to where the derivative is returned. 
//                 If null, no derivative returned
void GlGeomBezier::BezierMultiEval(
    double* cntlPts, int stride,
    int order, int alphaNumerator, int alphaDenominator,
    double* destVal, double* destDeriv)
{
    assert(alphaDenominator > 0);
    assert(order > 1);
    bool reverseLerpDirection = (alphaNumerator > alphaDenominator / 2);
    double alpha;
    if (reverseLerpDirection) {
        // Reverse the traversal order
        cntlPts += stride * (order - 1);
        stride = -stride;
        alpha = (double)(alphaDenominator-alphaNumerator) / (double)alphaDenominator;
    }
    else {
        alpha = (double)alphaNumerator / (double)alphaDenominator;
    }
    double beta = 1.0 - alpha;
    double* workSpace = new double[(order * (order + 1)) >> 1];
    for (int i = 0; i < numCoordinates; i++) {
        // Handle i-th coordinate
        // First, copy control points into the workSpace array
        for (int j = 0; j < order; j++) {
            *(workSpace + j) = *(cntlPts + j * stride + i);
        }
        // De Casteljau algorithm
        double* fromPtr = workSpace;
        double* toPtr = workSpace + order;
        for (int j = 0; j < order - 1; j++) {
            for (int k = 0; k < order - j - 1; k++) {
                *toPtr = beta * (*fromPtr) + alpha * (*(fromPtr + 1));
                toPtr++;
                fromPtr++;
            }
            fromPtr++;
        }
        // Copy answers to correct locations
        double* ansPtr = workSpace + (((order * (order + 1)) >> 1) - 1);
        *(destVal + i) = *ansPtr;
        if (destDeriv != 0) {
            double deriv = (*(ansPtr - 1) - *(ansPtr - 2));
            if (reverseLerpDirection) {
                deriv = -deriv;
            }
            *(destDeriv + i) = deriv;
        }
    }
   
    delete[] workSpace;
}

// **********************************************
// These routines do the rendering.
// If the patches' VAO, VBO, EBO need to be loaded, it does this first.
// **********************************************

void GlGeomBezier::PreRender()
{
    GlGeomBase::PreRender();

    if (!VboEboLoaded) {
        ReInitializeAttribLocations();
    }
}

// Render entire torus as triangles
void GlGeomBezier::Render()
{
    PreRender();
    GlGeomBase::Render();
}

// Render the i-th patch as triangles
void GlGeomBezier::RenderPatch(int i)
{
    PreRender();

    int patchEboIndex = 3*GetFirstTriInPatch(i);
    GlGeomBase::RenderEBO(GL_TRIANGLES, 3*GetNumTrisInPatch(i), patchEboIndex);
}

// **********************************************
// Simple vector functions (without classes)
// **********************************************
void GlGeomBezier::VecCopy(const double* source, double* retCopy, int numComponents) {
    for (int i = 0; i < numComponents; i++) {
        *(retCopy++) = *(source++);
    }
}

void GlGeomBezier::VecDiff(const double* minuhend, const double* subtrahend, double* retDifference, int numComponents) {
    for (int i = 0; i < numComponents; i++) {
        *(retDifference++) = *(minuhend++) - *(subtrahend++);
    }
}

double GlGeomBezier::VecDotProd(const double* vecA, const double* vecB, int numComponents)
{
    double ret = 0.0;
    for (int i = 0; i < numComponents; i++) {
        ret += (*(vecA++)) * (*(vecB++));
    }
    return ret;
}

void GlGeomBezier::VecCrossProd(const double* aVec, const double* bVec, double* retCrossProd) {
    *retCrossProd = (*(aVec + 1)) * (*(bVec + 2)) - (*(aVec + 2)) * (*(bVec + 1));
    *(retCrossProd + 1) = (*(aVec + 2)) * (*bVec) - (*(aVec)) * (*(bVec + 2));
    *(retCrossProd + 2) = (*aVec) * (*(bVec + 1)) - (*(aVec + 1)) * (*bVec);
}

double GlGeomBezier::VecNormSq(const double* vec, int numComponents) {
    double ret = 0.0;
    for (int i = 0; i < numComponents; i++) {
        double u = *(vec++);
        ret += u * u;
    }
    return ret;
}

void GlGeomBezier::VecMultScalar(const double* multiplicand, double multiplier, double* retProduct, int numComponents) {
    for (int i = 0; i < numComponents; i++) {
        *(retProduct++) = *(multiplicand++) * multiplier;
    }
}

// Add (multiplier)*(multiplicand vector) to (retAccum vector)
void GlGeomBezier::VecAddScaled(const double* multiplicand, double multiplier, double* retAccum, int numComponents) {
    for (int i = 0; i < numComponents; i++) {
        *(retAccum++) += *(multiplicand++) * multiplier;
    }
}

// Set equal to the zero vector
void GlGeomBezier::VecSetZero(double* retZero, int numComponents) {
    for (int i = 0; i < numComponents; i++) {
        *(retZero++) = 0.0;
    }
}

bool GlGeomBezier::VecIsZero(const double* vec, int numComponents)
{
    for (int i = 0; i < numComponents; i++) {
        if (*(vec++) != 0.0) {
            return false;
        }
    }
    return true;
}
