/* ######################################################################
* File Name: smlmsMicroscope.cpp
* Project: ermine
* Version: 19.02
* Creation Date:11.03.2016
* Created By Sebastian Malkusch
* <malkusch@chemie.uni-frankfurt.de>
* Goethe University of Frankfurt
* Physical and Theoretical Chemistry
* Single Molecule Biophysics
###################################################################### */

#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <iomanip>
#include "header/smlmsMicroscope.hpp"
#include "header/smlmsExceptions.hpp"

namespace SMLMS{
/* Constructor */
Microscope::Microscope(){
	//std::cout<<"Microscope Constructor called."<<std::endl;
}

Microscope::Microscope(double initPxlSize, double initIntTime, double initLocPrec){
	setPxlSize(initPxlSize);
	setIntTime(initIntTime);
	setLocPrec(initLocPrec);
}
/* Destructor */
Microscope::~Microscope(){
	//std::cout<<"Microscope removed from Heap!"<<std::endl;
}

/* Copy Constructor*/
Microscope::Microscope(const Microscope &obj){
	//std::cout<<"Microscope Copy Constructor Called."<<std::endl;
	setPxlSize(obj._pxlSize);
	setIntTime(obj._intTime);
	setLocPrec(obj._locPrec);
}
/* Elementary Functions */

void Microscope::setPxlSize(double initPxlSize){
	_pxlSize = initPxlSize;
}

double Microscope::pxlSize(){
	return _pxlSize;
}

double Microscope::pxlSize()const{
	return _pxlSize;
}

void Microscope::setIntTime(double initIntTime){
	_intTime = initIntTime;
}

double Microscope::intTime(){
	return _intTime;
}

double Microscope::intTime()const{
	return _intTime;
}

void Microscope::setLocPrec(double initLocPrec){
	_locPrec = initLocPrec;
}

double Microscope::locPrec(){
	return _locPrec;
}

double Microscope::locPrec()const{
	return _locPrec;
}

/* Special Functions */
void Microscope::clearMicroscope(){
	setPxlSize(0);
	setIntTime(0);
	setLocPrec(0);
}

void Microscope::loadMicroscope(const std::string &name){
	clearMicroscope();
	std::string line;
	std::vector<double> input(3);
	std::fstream inFile(name.data());
	if (inFile.is_open()){
		int i=0;
		while(std::getline(inFile, line)){
			if( !line.find("#")) continue;
			std::stringstream lineContent(line);
			lineContent>>input.at(i);
			i+=1;
		}
		inFile.close();
		setPxlSize(input.at(0));
		setIntTime(input.at(1));
		setLocPrec(input.at(2));	
	}
	else{
		std::stringstream errorMessage;
		errorMessage<<name<<" could not be opened!"<<std::endl;
		SMLMS::SmlmsError error(errorMessage.str());
		throw error;	
	}	
}

void Microscope::saveMicroscope(const std::string &folderName){
	std::string name = folderName;
	name.append("/microscope.txt");
	std::ofstream outFile;
	outFile.open(name.data());
	outFile<<std::scientific;
	outFile<<std::setprecision(6);
	/* header line */
	outFile<<"# Ermine Microscope File"<<std::endl;
	outFile<<"# pxl size[nm]"<<std::endl;
	outFile<<"# integration time [s]"<<std::endl;
	outFile<<"# localization precision [nm]"<<std::endl;
	/* data */
	outFile<<pxlSize()<<std::endl;
	outFile<<intTime()<<std::endl;
	outFile<<locPrec()<<std::endl;
	/* close */
	outFile.close();
}

}/* SMLMS */
