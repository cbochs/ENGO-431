/*
 * The purpose of this header is to provide the design matrices to the user that
 * are required when attempting complete an affine and/or similarity transformation
 */

#pragma once

#include "LeastSquares.h"
#include "Point.h"

/** similarity_A
 * Computes an N-by-4 design matrix for a linear similarity transformation
 * 
 * @param A			  - the desired output design matrix
 * @param coordinates - "from" coordinates for determinng the design matrix
 */
void similarity_A(Matrix &A, PointList &coordinates);

/** similarity_N
 * Computes a 4-by-4 normal matrix for a linear similarity transformation
 * 
 * N = A_trans * A
 * 
 * @param N			  - the desired output normal matrix
 * @param coordinates - "from" coordinates for determining the normal matrix
 */
void similarity_N(Matrix &N, PointList &coordinates);

/** similarity_u
 * Computes a 4-by-1 normal vector for a linear similarity transformation
 * 
 * u = A_trans * w
 * w = l - f(x)
 * 
 * @param u			   - the desired output normal vector
 * @param coordinates  - "from" coordinates for determining the normal vector
 * @param observations - "to" coordinates acting as the misclosure vector with f(x) = 0
 */
void similarity_u(Matrix &u, PointList &coordinates, PointList &observations);

/** affine_A
 * Computes an N-by-6 design matrix for a linear affine transformation
 * 
 * @param A			  - the desired output design matrix
 * @param coordinates - "from" coordinates for determining the design matrix
 */
void affine_A(Matrix &A, PointList &coordinates);

/** affine_N
 * Computes a 6-by-6 normal matrix for a linear affine transformation
 * 
 * N = A_trans * A
 * 
 * @param N			  - the desired output normal matrix
 * @param coordinates - "from" coordinates for determining the normal matrix
 */
void affine_N(Matrix &N, PointList &coordinates);

/** affine_u
 * Computes a 6-by-1 normal vector for a linear affine transformation
 *
 * u = A_trans * w
 * w = l - f(x)
 * 
 * @param u			   - the desired output normal vector
 * @param coordinates  - "from" coordinates for determining the normal vector
 * @param observations - "to" coordinates acting as the misclosure vector with f(x) = 0
 */
void affine_u(Matrix &u, PointList &coordinates, PointList &observations);
