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


#ifndef ermineExceptions_hpp
#define ermineExceptions_hpp
#include <string>

namespace SMLMS{

class ErmineParserError{
	private:
		std::string _errorMessage;
	public:
		ErmineParserError(std::string);
		std::string what();
}; /* ErmineParserError */


class SMLMSFolderError{
	private:
		std::string _errorMessage;
	public:
		SMLMSFolderError(std::string);
		std::string what();
}; /* SMLMSFolderError */

} /* SMLMS */
#endif /* ermineExceptions_hpp */

