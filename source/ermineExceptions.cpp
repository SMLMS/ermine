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
#include "header/ermineExceptions.hpp"
namespace SMLMS{
ErmineError::ErmineError(std::string errorMessage){
	_errorMessage = "Oops, an ermine parser error occurred: ";
	_errorMessage.append(errorMessage);
	_errorMessage.append("\ntype './ermine --help (-h)' for help!\n");
}

std::string ErmineError::what(){
	return _errorMessage;
} 

ErmineParserError::ErmineParserError(std::string errorMessage){
	_errorMessage = "Oops, an ermine parser error occurred: ";
	_errorMessage.append(errorMessage);
	_errorMessage.append("\ntype './ermine --help (-h)' for help!\n");
}

std::string ErmineParserError::what(){
	return _errorMessage;
}

SMLMSFolderError::SMLMSFolderError(std::string errorMessage){
	_errorMessage = "Oops, an ermine folder error occurred: ";
	_errorMessage.append(errorMessage);
	_errorMessage.append("\n");
}

std::string SMLMSFolderError::what(){
	return _errorMessage;
}

ErmineFileNameError::ErmineFileNameError(std::string errorMessage){
	_errorMessage = "Oops, an ermine file name error occurred: ";
	_errorMessage.append(errorMessage);
	_errorMessage.append("\n");
}

std::string ErmineFileNameError::what(){
	return _errorMessage;
}

SMLMSMicroscopeError::SMLMSMicroscopeError(std::string errorMessage){
	_errorMessage = "Oops an ermine microscope error occurred: ";
	_errorMessage.append(errorMessage);
	_errorMessage.append("\n");
}

std::string SMLMSMicroscopeError::what(){
	return _errorMessage;
}

SMLMSMoleculesError::SMLMSMoleculesError(std::string errorMessage){
	_errorMessage = "Oops, an ermine molecule list error occurred: ";
	_errorMessage.append(errorMessage);
	_errorMessage.append("\n");
}

std::string SMLMSMoleculesError::what(){
	return _errorMessage;
}

ErmineJudiError::ErmineJudiError(std::string errorMessage){
	_errorMessage = "Oops, an ermine jump distance distribution error occurred:";
	_errorMessage.append(errorMessage);
	_errorMessage.append("\n");
}

std::string  ErmineJudiError::what(){
	return _errorMessage;
}
}/* SMLMS */
