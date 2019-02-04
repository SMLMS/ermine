/* ######################################################################
* File Name: smlmsRandom.hpp
* Project: ermine
* Version: 19.02
* Creation Date: 2016
* Created By Sebastian Malkusch
* <malkusch@chemie.uni-frankfurt.de>
* Goethe University of Frankfurt
* Physical and Theoretical Chemistry
* Single Molecule Biophysics
###################################################################### */


#ifndef SMLMSRANDOM_HPP
#define SMLMSRANDOM_HPP

#include <boost/random/mersenne_twister.hpp>

namespace SMLMS{

class SMLMSRandom{
	private:		
		boost::random::mt19937 _seed;
	public:
		SMLMSRandom();
		void updateSeed();
		double generateRandomDouble(double min, double max);
};/* SMLMSRandom */


}/* SMLMS */

#endif /* smlmsRandom_hpp */

