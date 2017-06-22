/* ######################################################################
* File Name: ermineExceptions.hpp
* Project: ermine
* Version: 16.12
* Creation Date: 02.12.2016
* Created By Sebastian Malkusch
* <malkusch@chemie.uni-frankfurt.de>
* Goethe University of Frankfurt
* Physical and Theoretical Chemistry
* Single Molecule Biophysics
###################################################################### */


#ifndef SmlmsExceptions_hpp
#define SmlmsExceptions_hpp
#include <string>

namespace SMLMS{

class SmlmsError{
	private:
		std::string _errorMessage;
	public:
		SmlmsError(std::string);
		std::string what();
};/* SmlmsError */

} /* SMLMS */
#endif /* SmlmsExceptions_hpp */
