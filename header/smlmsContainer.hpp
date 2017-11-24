/* ######################################################################
* File Name: container.hpp
* Project: SMLMS
* Version:16.02
* Creation Date:03.03.2016
* Created By Sebastian Malkusch
* <malkusch@chemie.uni-frankfurt.de>
* Goethe University of Frankfurt
* Physical and Theoretical Chemistry
* Single Molecule Biophysics
###################################################################### */


#ifndef Smlms_hpp
#define Smlms_hpp
#include <string>
#include <vector>

namespace SMLMS{
	struct Molecule{
		int trace;
		int frame;
		double x;
		double y;
		double z;
		int state;
		double intensity;
		double precision;
	};

	struct ROI{
		double minX;
		double maxX;
		double minY;
		double maxY;
		double minFrame;
		double maxFrame;
		double minIntensity;
		double maxIntensity;
	};

	struct Jump{
		int trace;
		double jumpDistance;
		int state;
	};
	
	struct Settings{
    		double pxl;
    		double dt;
    		double prec;
	};

	struct HmmStatistics{
		unsigned states;
		unsigned symbols;
		double logLikelihood;
		unsigned dof;
		double bic;
		double aic;
	};

	struct ModelState{
		unsigned fix;
		double value;
		double min;
		double max;
	};
}/* SMLMS */
#endif /* Smlms_hpp */

