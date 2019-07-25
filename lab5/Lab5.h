/* The purpose of this header is to provide space for all functions/methods that do
 * not have a solid role within the classes established already. These methods will later
 * become implemented, but for the purpose of this lab they will reside in Lab5.h
 */

#pragma once

#include<vector>

#include "Point.h"
#include "Resection.h"

const double pi = acos(-1.0);

/** printStatistics
* prints the residuals (misclosure), correlation, covariance, and redundancy of the resection
*
* @param absolute - the resection object
*/
void printStatistics(Resection &resection);

/** printParameters
* prints the parameters of the resection
*
* @param absolute - the resection object
*/
void printParameters(Resection &resection);

/** isTesting
* a simple method to determine whether the use is inputting test data
* or an entire set of RTK data
*
* @return - true: using test data; false: not using test data
*/
bool isTesting();

/** rad2deg
* at least I didn't have to make a is_even function this time. Yes, this
* is converts a given radian angle to its decimal degree equivilant
*
* @param rad - the input angle in radians
*
* @return    - an output angle in degrees
*/
double rad2deg(double rad);

/** deg2rad
* do i even need to say what this means? As per mentioned above, this function
* does the opposite: converting a decimal degree angle to radians
*
* @param deg - the input angle in degrees
*
* @return	  - an output angle in radians
*/
double deg2rad(double deg);