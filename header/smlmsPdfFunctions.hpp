/* ######################################################################
* File Name: smlmsPdfFunctions.hpp
* Project: ermine
* Version: 19.02
* Creation Date: 31.03.2017
* Created By Sebastian Malkusch
* <malkusch@chemie.uni-frankfurt.de>
* Goethe University of Frankfurt
* Physical and Theoretical Chemistry
* Single Molecule Biophysics
###################################################################### */


#ifndef SMLMSPDFFUNCTIONS_HPP
#define SMLMSPDFFUNCTIONS_HPP
namespace SMLMS{
/* jump Distance Distribution */
double judiPdf(double *val, double *para);
double judiSuperPosPdf(double *val, double *para);
double expectedDiffCoeff(double dist, double dt, double sigma);
double expectedDistance(double diffCoeff, double dt, double sigma);

/* Normal Distribution */
//double normalPdf(double *val, double *para);
//double normalSuperPos(double *val, double *para);

/* dwell time analysis function */
double monoExpDecFunc(double *val, double *para);
}/*SMLMS*/
#endif /* SmlmsPdfFunctions_hpp */

