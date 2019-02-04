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


#ifndef ERMINEEXCEPTIONS_HPP
#define ERMINEEXCEPTIONS_HPP
#include <string>

namespace SMLMS{

class ErmineError{
	private:
		std::string _errorMessage;
	public:
		ErmineError(std::string);
		std::string what();
};/* ErmineError */

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

class ErmineFileNameError{
	private:
		std::string _errorMessage;
	public:
		ErmineFileNameError(std::string);
		std::string what();
}; /* SMLMSFileNameError*/

class SMLMSMicroscopeError{
	private:
		std::string _errorMessage;
	public:
		SMLMSMicroscopeError(std::string);
		std::string what();
}; /* SMLMSMicroscopeError */

class SMLMSMoleculesError{
	private:
		std::string _errorMessage;
	public:
		SMLMSMoleculesError(std::string);
		std::string what();
}; /* SMLMSMoleculesError */

class ErmineJudiError{
	private:
		std::string _errorMessage;
	public:
		ErmineJudiError(std::string);
		std::string what();
}; /* ErmineJudiExceptions */

} /* SMLMS */
#endif /* ermineExceptions_hpp */

