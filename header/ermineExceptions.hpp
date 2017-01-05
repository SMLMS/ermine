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

class NoFileName{
	private:
		std::string _errorMessage;
	public:
		NoFileName();
		std::string returnError();

}; /* NoFileName */


class WrongFileName{
	private:
		std::string _fileName;
		std::string _errorMessage;
	public:
		WrongFileName(std::string &);
		std::string returnError();
}; /* WrongFileName */

class NoAlgorithm{
	private:
		std::string _errorMessage;
	public:
		NoAlgorithm();
		std::string returnError();
}; /*  NoAlgorithm */

class WrongAlgorithm{
	private:
		std::string _algorithmName;
		std::string _errorMessage;
	public:
		WrongAlgorithm(std::string &);
		std::string returnError();
}; /* WrongAlgorithm */


} /* SMLMS */
#endif /* ermineExceptions_hpp */

