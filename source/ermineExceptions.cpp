/* ######################################################################
* File Name: ermineExceptions.cpp
* Project: ermine
* Version: 16.12
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
#include "header/ermineExceptions.hpp"
namespace SMLMS{
ErmineError::ErmineError(std::string errorMessage){
	_errorMessage = "Ermine Parser Error: ";
	_errorMessage.append(errorMessage);
	_errorMessage.append("\ntype './ermine --help (-h)' for help!\n");
}

std::string ErmineError::what(){
	return _errorMessage;
} 

ErmineParserError::ErmineParserError(std::string errorMessage){
	_errorMessage = "Ermine Parser Error: ";
	_errorMessage.append(errorMessage);
	_errorMessage.append("\ntype './ermine --help (-h)' for help!\n");
}

std::string ErmineParserError::what(){
	return _errorMessage;
}

SMLMSFolderError::SMLMSFolderError(std::string errorMessage){
	_errorMessage = "SMLMS Folder Error: ";
	_errorMessage.append(errorMessage);
	_errorMessage.append("\n");
}

std::string SMLMSFolderError::what(){
	return _errorMessage;
}

ErmineFileNameError::ErmineFileNameError(std::string errorMessage){
	_errorMessage = "Ermine File Name Error: ";
	_errorMessage.append(errorMessage);
	_errorMessage.append("\n");
}

std::string ErmineFileNameError::what(){
	return _errorMessage;
}



}/* SMLMS */
