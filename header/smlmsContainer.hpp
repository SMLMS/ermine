/* ######################################################################
* File Name: smlmsContainer.hpp
* Project: ermine
* Version: 19.02
* Creation Date:03.03.2016
* Created By Sebastian Malkusch
* <malkusch@chemie.uni-frankfurt.de>
* Goethe University of Frankfurt
* Physical and Theoretical Chemistry
* Single Molecule Biophysics
###################################################################### */


#ifndef SMLMSCONTAINER_HPP
#define SMLMSCONTAINER_HPP

namespace SMLMS{
	struct Molecule{
		unsigned trace;
		unsigned frame;
		double x;
		double y;
		double z;
		unsigned state;
		double intensity;
		double precision;
	};

	struct TraceStatistics{
		unsigned index;
		unsigned begin;
		unsigned end;
		unsigned length;
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
		unsigned trace;
		double jumpDistance;
		unsigned state;
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
		double value;
		unsigned fix;
		double min;
		double max;
	};
}/* SMLMS */
#endif /* Smlms_hpp */

