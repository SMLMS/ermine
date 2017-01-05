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

SMLMS::NoFileName::NoFileName(){
	std::stringstream message;
	message<<"Error: No filename was given to ermine"<<std::endl<<"type --help(-h) for help."<<std::endl;
	_errorMessage = message.str();
}

std::string SMLMS::NoFileName::returnError(){
	return _errorMessage;
}

SMLMS::WrongFileName::WrongFileName(std::string &name){
	_fileName=name;
	std::stringstream message;
	message<<"Error: Wrong filename \nermine could not open "<<_fileName<<"\nmake sure it does exist!"<<std::endl;
	_errorMessage = message.str();
}

std::string SMLMS::WrongFileName::returnError(){
	return _errorMessage;
}

SMLMS::NoAlgorithm::NoAlgorithm(){
	std::stringstream message;
	message<<"Error: No algorithm was chosen by user"<<std::endl<<"type --help (-h) for help."<<std::endl;
	_errorMessage = message.str();
}

std::string SMLMS::NoAlgorithm::returnError(){
	return _errorMessage;
}

SMLMS::WrongAlgorithm::WrongAlgorithm(std::string &name){
	_algorithmName=name;
	std::stringstream message;
	message<<_algorithmName<<" is not a valid algorithm for ermine!."<<std::endl<<"type --help (-h) for help."<<std::endl;
	_errorMessage = message.str();
}

std::string SMLMS::WrongAlgorithm::returnError(){
	return _errorMessage;
}
