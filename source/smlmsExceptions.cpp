/* ######################################################################
* File Name: ermineExceptions.cpp
* Project: ermine
* Version: 19.02
* Creation Date: 02.12.2016
* Created By Sebastian Malkusch
* <malkusch@chemie.uni-frankfurt.de>
* Goethe University of Frankfurt
* Physical and Theoretical Chemistry
* Single Molecule Biophysics
###################################################################### */

#include <iostream>
#include <string>
#include <sstream>
#include "header/smlmsExceptions.hpp"
namespace SMLMS{

SmlmsError::SmlmsError(std::string errorMessage){
	_errorMessage = "SMLMS Error: ";
	_errorMessage.append(errorMessage);
	_errorMessage.append("\ntype './ermine --help (-h)' for help!\n");
}

std::string SmlmsError::what(){
	return _errorMessage;
} 

}/* SMLMS */
