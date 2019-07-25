#pragma once

#include <string>
#include "Matrix.h"

/** delta
 * Computes the delta vector in a Least-Squares adjustment
 * 
 * delta = x_hat - x
 * delta = -(A_trans * P * A)^-1 * A_trans * P * w
 * delta = -(A_trans * A)^-1 * A_trans * w
 * delta = -N^-1 * u
 * 
 * @param A		- design matrix for the adjustment (constant)
 * @param P		- weight matrix for the adjustment
 * @param w		- misclosure vector for the adjustment
 * 
 * @return - the desired output delta vector
 */
Matrix delta(Matrix &A, const Matrix &P, const Matrix &w);
Matrix delta(Matrix &A, const Matrix &w);

/** belowTolerance
* checks a delta vector against a given tolerance vector and returns whether
* or not the delta is completely under the tolerance values
*
* @param delta		 - an n-by-1 delta vector to be checked
* @param tolerances - an n-by-1 tolerances vector to be compared to
*
* @return			 - true if ALL the delta values are below their corresponding tolerance
*/
bool belowTolerances(const Matrix &delta, const Matrix &tolerances);

/** misclosure
 * Computes the misclosure vector in a Least-Squares adjustment
 * 
 * w = f(x) - l
 * 
 * @param obs - the observations of the adjustment
 * @param est - the estimated f(x) values
 * 
 * @return    - the desired output misclosure vector
 */
Matrix misclosure(const Matrix &obs, const Matrix &est);

/** residuals
 * Computes the residual vector in a Least-Squares Adjustment
 * 
 * v = l_hat - l
 * v = A * delta + w
 * 
 *
 * @param A		- design matrix for the adjustment (constant)
 * @param delta - delta vector fot the adjustment
 * @param W		- misclosure vector for the adjustment
 * 
 * @return		- the desired output residual vector
 */
Matrix residuals(Matrix &A, const Matrix &delta, const Matrix &w);

/** aposteriori
 * Computes the aposteriori of the Least-Squares adjustment
 * 
 * apost = (v_trans * P * v) / d.o.f
 * apost = (v_trans * v) / d.o.f
 * 
 * v = (A * delta + w || l_hat - l)
 * 
 * @param v - the residual vector of the adjustment (constant)
 * @param P - the weight matrix of the adjuesmtment
 * 
 * @return  - the aposteriori value of the LS adjustment
 */
double aposteriori(Matrix &v, const Matrix &P, double dof);
double aposteriori(Matrix &v, double dof);

/** cofactorMatrix
 * Computes the cofactor matrix of the Least-Squares adjustment
 * 
 * Qx = (A_trans * P * A)^-1
 * Qx = (A_trans * A)^-1
 * 
 * @param A	 - the design matrix for the adjustment
 * @param P  - the weight matrix for the adjustment
 * 
 * @return   - the desired output cofactor matrix
 */
Matrix cofactorMatrix(Matrix &A, const Matrix &P);
Matrix cofactorMatrix(Matrix &A);

/** unknownCovariance
 * Computes the estimated covariance matrix of the Least-Squares adjustment
 * 
 * Cx = aposteriori * Qx
 * 
 * @param aposteriori - apost. for the adjustment
 * @Param Qx		  - the cofactor matrix for the adjustment
 * 
 * @return			  - the desired output estimated covariance matrix
 */
Matrix unknownCovariance(Matrix &A, const Matrix &P, double aposteriori);
Matrix unknownCovariance(Matrix &A, double aposteriori);

/** observeCovariance
 * Computes the corrected covariance matrix of the Least-Squares adjustment
 * 
 * Cl = A * Cx * A_trans
 * 
 * @param A  - the design matrix for the adjustment
 * @param Cx - the estimated covariance matrix for the adjustment
 * 
 * @return   - the desired output corrected covariance matrix
 */
Matrix observeCovariance(Matrix &A, const Matrix &P);
Matrix observeCovariance(Matrix &A);

/** residualCovariance
 * Computes the residual covaraince matrix of the Least-Squares adjustment
 * 
 * Cv = C - Cl
 * 
 * @param C  - the observervations covariance matrix for the adjustment
 * @param Cl - the corrected covariance matrix for the adjustment
 * 
 * @return	 - the desired output residual covariance matrix
 */
Matrix residualCovariance(const Matrix &C, Matrix &A, const Matrix &P);
Matrix residualCovariance(Matrix &A);

/** correlation
 * Computes the correlation matrix of the unknowns in a Least-Squares adjustment
 * 
 * @param A  - the design matrix for the adjustment
 * 
 * @return   - the desired output correlation matrix
 */
Matrix correlation(Matrix &A);

/** eye
 * Computes an n-by-n iidentity matrix
 * 
 * eye = diag(1, 1, ..., 1, 1)
 *  
 * @param size - the size of the matrix
 * 
 * @return - the desired output identity matrix
 */
Matrix eye(int size);