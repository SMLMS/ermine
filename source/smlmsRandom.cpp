/* ######################################################################
* File Name: smlmsRandom.cpp
* Project: ermine
* Version: 19.02
* Creation Date: 2016
* Created By Sebastian Malkusch
* <malkusch@chemie.uni-frankfurt.de>
* Goethe University of Frankfurt
* Physical and Theoretical Chemistry
* Single Molecule Biophysics
###################################################################### */


#include <ctime>
#include "header/smlmsRandom.hpp"
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_real_distribution.hpp>
#include <boost/random/variate_generator.hpp>

namespace SMLMS{

SMLMSRandom::SMLMSRandom(){
	updateSeed();
}

void SMLMSRandom::updateSeed(){
	boost::random::mt19937 seed(std::time(0));
	_seed = seed;
}

double SMLMSRandom::generateRandomDouble(double min, double max){
	boost::random::uniform_real_distribution<> dist(min,max);
  	boost::variate_generator<boost::mt19937&, boost::random::uniform_real_distribution<> > random(_seed,dist);
  	return random();
}

}/* SMLMS */
