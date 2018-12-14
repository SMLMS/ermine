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


#ifndef smlmsRandom_hpp
#define smlmsRandom_hpp

#include <boost/random/mersenne_twister.hpp>
#include <boost/mpi.hpp>
#include <boost/mpi/collectives.hpp>

namespace SMLMS{

class SMLMSRandom{
	private:
		friend class boost::serialization::access;
		template<class Archive>
		void serialize(Archive & ar, const unsigned int version){
			ar & _seed;
		}
		boost::random::mt19937 _seed;
	public:
		SMLMSRandom();
		double generateRandomDouble(double min, double max);
};/* SMLMSRandom */


}/* SMLMS */

#endif /* smlmsRandom_hpp */

