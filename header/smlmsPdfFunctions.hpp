/* ######################################################################
* File Name: smlmsPdfFunctions.hpp
* Project: SMLMS
* Version: 16.03
* Creation Date: 31.03.2017
* Created By Sebastian Malkusch
* <malkusch@chemie.uni-frankfurt.de>
* Goethe University of Frankfurt
* Physical and Theoretical Chemistry
* Single Molecule Biophysics
###################################################################### */


#ifndef SmlmsPdfFunctions_hpp
#define SmlmsPdfFunctions_hpp
namespace SMLMS{
/* jump Distance Distribution */
double judiPdf(double *val, double *para);
double judiSuperPosPdf(double *val, double *para);
double expectedDiffCoeff(double dist, double dt, double sigma);
double expectedDistance(double diffCoeff, double dt, double sigma);

/* Normal Distribution */
//double normalPdf(double *val, double *para);
//double normalSuperPos(double *val, double *para);

}/*SMLMS*/
#endif /* SmlmsPdfFunctions_hpp */

