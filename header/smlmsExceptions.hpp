/* ######################################################################
* File Name: ermineExceptions.hpp
* Project: ermine
* Version: 19.02
* Creation Date: 02.12.2016
* Created By Sebastian Malkusch
* <malkusch@chemie.uni-frankfurt.de>
* Goethe University of Frankfurt
* Physical and Theoretical Chemistry
* Single Molecule Biophysics
###################################################################### */


#ifndef SMLMSEXCEPTIONS_HPP
#define SMLMSEXCEPTIONS_HPP
#include <string>

namespace SMLMS{

enum {
	SMLMS_SUCCESS = 0,
	SMLMS_FAILURE = -1
};


class SmlmsError{
	private:
		std::string _errorMessage;
	public:
		SmlmsError(std::string);
		std::string what();
};/* SmlmsError */

} /* SMLMS */
#endif /* SmlmsExceptions_hpp */

