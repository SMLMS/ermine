/* ######################################################################
* File Name:
* Project: 
* Version:
* Creation Date:
* Created By Sebastian Malkusch
* <malkusch@chemie.uni-frankfurt.de>
* Goethe University of Frankfurt
* Physical and Theoretical Chemistry
* Single Molecule Biophysics
###################################################################### */


#include <iostream>
#include <cmath>
#include "header/smlmsPdfFunctions.hpp"

namespace SMLMS{

double judiPdf(double *r, double *para){
        double y=0;
        y = (para[0]*2*r[0])/(4*para[1]*para[2]+2*std::pow(para[3],2))*
                std::exp((-std::pow(r[0],2))/(4*para[1]*para[2]+2*pow(para[3],2)));
        return y;
}

double judiSuperPosPdf(double *r, double *para){
	double y=0;
	int states = (int) para[0];
	for (int i=0; i<states; i++){
		y += judiPdf(&r[0], &para[1+(i*4)]);
	}
	return y;
}
}/* SMLMS */
