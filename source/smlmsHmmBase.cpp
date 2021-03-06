/* ######################################################################
* File Name: smlmsHmmBase.cpp
* Project: ermine
* Version: 19.02
* Creation Date: 23.03.2016
* Created By Sebastian Malkusch
* <malkusch@chemie.uni-frankfurt.de>
* Goethe University of Frankfurt
* Physical and Theoretical Chemistry
* Single Molecule Biophysics
###################################################################### */


#include <iostream>
#include <vector>
#include <string>
#include <stdexcept>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <math.h>
#include "header/smlmsExceptions.hpp"
#include "header/smlmsHmmBase.hpp"
#include "header/smlmsMatrix.hpp"
#include "header/smlmsPhysModBrownLatDiff.hpp"

namespace SMLMS{
/* constructor */
HMMBase::HMMBase (){
	//std::cout<<"HMM constructor called."<<std::endl;
	_logLikelihood = 0.0;
	_dof = 0;
	_bic = 0.0;
	_aic = 0.0;
}

HMMBase::HMMBase(unsigned states, unsigned symbols){
	//std::cout<<"HMM base constructor called"<<std::endl;
	//set everytging
	setStateNumber(states);
	setSymbolNumber(symbols);
	initHMM();
	checkHMM();
	calcDof();
	_logLikelihood = 0.0;
	_bic = 0.0;
	_aic = 0.0;
}

/* destructor */
HMMBase::~HMMBase (){
	//std::cout<<"HMMBase removed from Heap!"<<std::endl;
}

/* copy constructor */
HMMBase::HMMBase(const HMMBase &obj){
	//std::cout<<"HMM copy constructor called."<<std::endl;
	_stateNumber = obj._stateNumber;
	_symbolNumber  = obj._symbolNumber;
	_symbolInterval = obj._symbolInterval;
	_equiPDF = obj._equiPDF;
	_transPDF = obj._transPDF;
	_obsPDF = obj._obsPDF;
	_equiCDF = obj._equiCDF;
	_transCDF = obj._transCDF;
	_obsCDF = obj._obsCDF;
	_obsAlphabet = obj._obsAlphabet;
	_logLikelihood = obj._logLikelihood;
	_dof = obj._dof;
	_bic = obj._bic;
	_aic = obj._aic;
	_stopCrit = obj._stopCrit;
	_maxIt = obj._maxIt;
}

/* elementary functions */
void HMMBase::setStateNumber(unsigned numberOfStates){
	_stateNumber = numberOfStates;
}

unsigned HMMBase::stateNumber(){
	return _stateNumber;
}

void HMMBase::setSymbolNumber(unsigned numberOfSymbols){
	_symbolNumber = numberOfSymbols;
}

unsigned HMMBase::symbolNumber(){
	return _symbolNumber;
}

void HMMBase::setSymbolInterval(double data){
	_symbolInterval = data;
}

double HMMBase::symbolInterval(){
	return _symbolInterval;
}

void HMMBase::setMinValue(double val){
	_minValue = val;
}

double HMMBase::minValue(){
	return _minValue;
}

void HMMBase::setMaxValue(double val){
	_maxValue = val;
}

double HMMBase::maxValue(){
	return _maxValue;
}

void HMMBase::setEquiPDF(SMLMS::Matrix matrix){
	_equiPDF = matrix;
}

SMLMS::Matrix HMMBase::equiPDF(){
	return _equiPDF;
}

void HMMBase::setTransPDF(SMLMS::Matrix matrix){
	_transPDF = matrix;
}

SMLMS::Matrix HMMBase::transPDF(){
	return _transPDF;
}

void HMMBase::setObsPDF(SMLMS::Matrix matrix){
	_obsPDF = matrix;
}

SMLMS::Matrix HMMBase::obsPDF(){
	return _obsPDF;
}

void HMMBase::setObsAlphabet(std::vector<double> data){
	_obsAlphabet = data;
}

std::vector<double> HMMBase::obsAlphabet(){
	return _obsAlphabet;
}

void HMMBase::setLogLikelihood(double data){
	_logLikelihood = data;
}

double HMMBase::logLikelihood(){
	return _logLikelihood;
}

void HMMBase::setDof(unsigned data){
	_dof = data;
}

unsigned HMMBase::dof(){
	return _dof;
}

void HMMBase::setBic(double data){
	_bic = data;
}

double HMMBase::bic(){
	return _bic;
}

void HMMBase::setAic(double data){
	_aic = data;
}

double HMMBase::aic(){
	return _aic;
}

void HMMBase::setStopCrit(double val){
	_stopCrit = val;
}

double HMMBase::stopCrit(){
	return _stopCrit;
}

void HMMBase::setMaxIt(int val){
	_maxIt = val;
}

int HMMBase::maxIt(){
	return _maxIt;
}

/* read/write functions */
void HMMBase::readHMM(std::string const &name){
	unsigned n=0, row=0;
	std::string line;
	std::fstream inFile(name.data());
	if (inFile.is_open()){
		while(std::getline(inFile, line)){
			/* skip comments */
			if(!line.find("#")) continue;
			std::stringstream lineContent(line); 
			/* read state number */
			if (n==0){
				lineContent>>_stateNumber;
				initEqui();
				initTrans();
			}
			/* read para number */
			if (n==1){
				lineContent>>_symbolNumber;
				initObs();
				initObsAlphabet();
			}
			/* read initial probability: pi */
			if (n==2){
				for (unsigned i=0; i<_stateNumber; i++){
					lineContent>>_equiPDF(0,i);
				}
			}
			/* read transition probability: A */
			if(n>2 && n<(3+_stateNumber)){
				row=n-3;
				for (unsigned column=0; column<_stateNumber; column++){
					lineContent>>_transPDF(row, column);
				}
			}
			/* read observation alphabet */
			if(n==3+_stateNumber){
				for(unsigned i=0; i<_symbolNumber; i++) lineContent>>_obsAlphabet.at(i);
				extractParasFromObsAlphabet();
				extractSymbolInterval();
				extractMinValue();
				extractMaxValue();
			}
			/* read observation parameter Matrix: B */
			if(n>(3+_stateNumber) && n<(4+_stateNumber+_stateNumber)){
				row = n-_stateNumber-4;
				for (unsigned column=0; column<_symbolNumber; column++){
					lineContent>>_obsPDF(row, column);
				}
			}
			/* read statistics */
			if(n==4+2*_stateNumber)lineContent>>_logLikelihood;
			if(n==5+2*_stateNumber)lineContent>>_dof;
			if(n==6+2*_stateNumber)lineContent>>_aic;
			if(n==7+2*_stateNumber)lineContent>>_bic;
			/* increase iteration */
			n +=1;
		}
		/* close */
		inFile.close();
	}
	else{
		std::stringstream errorMessage;
		errorMessage<<"could not read hmm from file!"<<std::endl;
		SMLMS::SmlmsError error(errorMessage.str());
		throw error;
	}	

}

void HMMBase::writeHMM(const std::string &folderName){
	std::string name = folderName;
	name.append("/hmm.txt");
	std::ofstream outFile(name.data());
	if (outFile.is_open()){
		outFile<<std::scientific;
		outFile<<std::setprecision(6);
		/* header line */
		outFile<<"# Hidden Markov Model"<<std::endl;
		outFile<<"# ermine format"<<std::endl;
		/* State Number */
		outFile<<"# number of states: N"<<std::endl;
		outFile<<stateNumber()<<std::endl;
		/* para Number */
		outFile<<"# number of symbols in observation matrix: "<<std::endl;
		outFile<<symbolNumber()<<std::endl;
		/* initial probability: pi */
		outFile<<"# initial probability vector: pi"<<std::endl;
		for (unsigned i=0; i<_stateNumber; i++){
			outFile<<_equiPDF.at(0,i)<<"\t";
		}
		outFile<<std::endl;
		/* transition array: A */
		outFile<<"# transition probability matrix: A"<<std::endl;
		for (unsigned row=0; row<_stateNumber; row++){
			for (unsigned column=0; column<_stateNumber; column++){
				outFile<<_transPDF.at(row, column)<<"\t";
			}
			outFile<<std::endl;
		}
		/* observation alphabet */
		outFile<<"# observation alphabet:"<<std::endl;
		for (unsigned i=0; i<_symbolNumber; i++)outFile<<_obsAlphabet.at(i)<<"\t";
		outFile<<std::endl;
		/* observation parameter array: B */
		outFile<<"# observation probability matrix B"<<std::endl;
		for (unsigned row=0; row<_stateNumber; row++){
			for (unsigned column=0; column<_symbolNumber; column++) outFile<<_obsPDF.at(row,column)<<"\t";
			outFile<<std::endl;
		}
		/* print statistics */
		outFile<<"# LogLikelihood:"<<std::endl<<_logLikelihood<<std::endl;
		outFile<<"# DOF:"<<std::endl<<_dof<<std::endl;
		outFile<<"# AIC:"<<std::endl<<_aic<<std::endl;
		outFile<<"# BIC:"<<std::endl<<_bic<<std::endl;
		/* close */
		outFile<<"# go ermine!"<<std::endl;
		outFile.close();
		std::cout<<"\nwrote hmm to "<<name<<std::endl;
	}
	else{
		std::stringstream errorMessage;
		errorMessage<<"could not write hmm to "<<name<<std::endl;
		SMLMS::SmlmsError error(errorMessage.str());
		throw error;
	}
}

/* clear functions */
void HMMBase::clearHMM(){
	_stateNumber=0;
	_symbolNumber=0;
	_equiPDF.clearMatrix();
	_transPDF.clearMatrix();
	_obsPDF.clearMatrix();
	_equiCDF.clearMatrix();
	_transCDF.clearMatrix();
	_obsCDF.clearMatrix();
	_logLikelihood = 0.0;
	_dof = 0;
	_bic = 0.0;
	_aic = 0.0;
}

/* special functions */
/* proof functions */
void HMMBase::checkStateNumber(){
	if(_stateNumber<1){
		std::stringstream errorMessage;
		errorMessage<<"Cannot setup a HmmBase instance without a definite state number!"<<std::endl;
		SMLMS::SmlmsError error(errorMessage.str());
		throw error;
	}
}

void HMMBase::checkSymbolNumber(){
	if(_symbolNumber<1){
		std::stringstream errorMessage;
		errorMessage<<"Cannot setup a hmm without a definite symbol number!"<<std::endl;
		SMLMS::SmlmsError error(errorMessage.str());
		throw error;
	}
}

void HMMBase::checkEqui(){
	if ((_stateNumber != _equiPDF.numberOfColumns()) || (_equiPDF.numberOfRows()!=1)){
		std::stringstream errorMessage;
		errorMessage<<"equilibrum matrix mismatch: Matrix needs to be: "<<_stateNumber<<" x "<<1<<std::endl;
		SMLMS::SmlmsError error(errorMessage.str());
		throw error;
	}
}
	

void HMMBase::checkTrans(){
	if ((_stateNumber != _transPDF.numberOfColumns()) ||(_stateNumber != _transPDF.numberOfRows())){
		std::stringstream errorMessage;
		errorMessage<<"transition matrix mismatch: Matrix needs to be: "<<_stateNumber<<" x "<<_stateNumber<<std::endl;
		SMLMS::SmlmsError error(errorMessage.str());
		throw error;
	}
}

void HMMBase::checkObs(){
	if ((_symbolNumber != _obsPDF.numberOfColumns()) ||(_stateNumber != _obsPDF.numberOfRows())){
		std::stringstream errorMessage;
		errorMessage<<"observation matrix mismatch: Matrix needs to be: "<<_stateNumber<<" x "<<_symbolNumber<<std::endl;
		SMLMS::SmlmsError error(errorMessage.str());
		throw error;
	}
}

void HMMBase::checkAlphabet(){
	if (_obsAlphabet.size() != _symbolNumber){
		std::stringstream errorMessage;
		errorMessage<<"observation alphabet error: no observation alphabet given!"<<_stateNumber<<" x "<<_symbolNumber<<std::endl;
		SMLMS::SmlmsError error(errorMessage.str());
		throw error;
	}
}

void HMMBase::checkProbDen(SMLMS::Matrix & mat){
	double prob;
	for (unsigned i=0; i<mat.numberOfRows(); i++){
		prob = 0;
		for (unsigned j=0; j<mat.numberOfColumns(); j++){
			prob += mat.at(i,j);
		}
		if ((prob < 0.99) || (prob > 1.01)){
			std::cout<<"probabilty desnity warning: PDF integral is not normalized to unity."<<std::endl;
		}
	}
}

void HMMBase::checkHMM(){
	//checkFolderName();
	checkStateNumber();
	checkSymbolNumber();
	checkEqui();
	checkTrans();
	checkObs();
	checkProbDen(_equiPDF);
	checkProbDen(_transPDF);
	checkProbDen(_obsPDF);
}

void HMMBase::checkPhysicalModelCompatibility(SMLMS::PhysicalModelBLD &physMod){
	int errorCounter=0;
	if (physMod.stateNumber() != _stateNumber) errorCounter += 1;
	if (physMod.incNumber() != _symbolNumber) errorCounter += 1;
	if (physMod.minValue() != _minValue) errorCounter += 1;
	if (physMod.maxValue() != _maxValue) errorCounter += 1;
	if (errorCounter > 0){
		std::stringstream errorMessage;
		errorMessage<<"HMM base instance:"<<std::endl;
		errorMessage<<"PhysicalModel (brownian lateral diffusion) and hmmBase are not compatible."<<std::endl;
		SMLMS::SmlmsError error(errorMessage.str());
		throw error;
	}
}

/* init functions */
void HMMBase::initObsAlphabet(){
	/* check stuff */
	checkSymbolNumber();
	std::vector <double> tempVec(_symbolNumber, 0.0);
	_obsAlphabet = tempVec;
}

void HMMBase::initEqui(){
	initEquiPDF();
	initEquiCDF();
}

void HMMBase::initEquiPDF(){
	_equiPDF.clearMatrix();
	SMLMS::Matrix tempMat(1, _stateNumber);
	_equiPDF = tempMat;
	for (unsigned i=0; i<_stateNumber; i++){
		_equiPDF(0,i) = 1.0/_stateNumber;
	}
}

void HMMBase::initEquiCDF(){
	_equiCDF.clearMatrix();
	SMLMS::Matrix tempMat(1, _stateNumber);
	_equiCDF = tempMat;
	calcCDF(_equiPDF, _equiCDF);
}

void HMMBase::initTrans(){
	initTransPDF();
	initTransCDF();
}

void HMMBase::initTransPDF(){
	_transPDF.clearMatrix();
	SMLMS::Matrix tempMat(_stateNumber, _stateNumber);
	_transPDF = tempMat;
	for (unsigned j=0; j<_stateNumber; j++){
		for (unsigned i=0; i<_stateNumber; i++){
			_transPDF(j,i)=1.0/_stateNumber;
		}
	}
}

void HMMBase::initTransCDF(){
	_transCDF.clearMatrix();
	SMLMS::Matrix tempMat(_stateNumber, _stateNumber);
	_transCDF = tempMat;
	calcCDF(_transPDF, _transCDF);
}

void HMMBase::initObs(){
	initObsPDF();
	initObsCDF();
}

void HMMBase::initObsPDF(){
	unsigned i,j;
	_obsPDF.clearMatrix();
	SMLMS::Matrix tempMat(_stateNumber, _symbolNumber);
	_obsPDF = tempMat;
	for (j=0; j<_stateNumber; j++){
		for (i=0; i<_symbolNumber; i++){
			_obsPDF(j,i) = 1.0/_symbolNumber;
		}
	}
}

void HMMBase::initObsCDF(){
	_obsCDF.clearMatrix();
	SMLMS::Matrix tempMat(_stateNumber, _symbolNumber);
	_obsCDF = tempMat;
	calcCDF(_obsPDF, _obsCDF);
}

void HMMBase::initHMM(){
	initEqui();
	initTrans();
	initObs();
	initObsAlphabet();
}

/* init from file */
void HMMBase::initLoadedHMM(){
	/* init */
	initEquiCDF();
	initTransCDF();
	initObsCDF();
	/* calc hmm members */
	extractParasFromObsAlphabet();
	/* normalize */
	normalizeHMM();
	/* calc cdf */
	calcCDF(_equiPDF, _equiCDF);
	calcCDF(_transPDF, _transCDF);
	calcCDF(_obsPDF, _obsCDF);
}

void HMMBase::initFromPhysMod(SMLMS::PhysicalModelBLD &physMod){
	/* adjust state Number */
	setStateNumber(physMod.stateNumber());
	/* adjust Symbol Number */
	setSymbolNumber(physMod.incNumber());
	/* adjust alphabet */
	setObsAlphabet(physMod.alphabet());
	/* initialize pdf container */
	initEquiPDF();
	initTransPDF();
	initObsPDF();
	/* calc equilibrium pdf */
	physMod.updatePi(_equiPDF);
	/* get obs pdf from physical model */
	physMod.normalizeModel();
	_obsPDF = physMod.fitMatrix();
	/* initialize hmm containers */
	initLoadedHMM();
}

/* normalization */
void HMMBase::normalizePDF(SMLMS::Matrix& mat){
	unsigned i,j;
	double area;
	//integrate
	for (i=0; i<mat.numberOfRows(); i++){
		area = 0.0;
		for (j=0; j<mat.numberOfColumns(); j++){
			area += mat.at(i,j); 
		}
		// normalize
		if (area>0.0){
			for (j=0; j<mat.numberOfColumns(); j++){
				mat(i,j)/=area;
			}
		}
		else{
			std::stringstream errorMessage;
			errorMessage<<"PDF error: Integrated Area = 0!"<<std::endl;
			SMLMS::SmlmsError error(errorMessage.str());
			throw error;
		}
	}
}

void HMMBase::normalizeHMM(){
	HMMBase::normalizePDF(_equiPDF);
	HMMBase::normalizePDF(_transPDF);
	HMMBase::normalizePDF(_obsPDF);
	HMMBase::calcCDF(_equiPDF, _equiCDF);
	HMMBase::calcCDF(_transPDF, _transCDF);
	HMMBase::calcCDF(_obsPDF, _obsCDF);
}

/* calc functionc */	
void HMMBase::calcSymbolNumber(){
	_symbolNumber = int(_maxValue-_minValue)/int(_symbolInterval);
	checkSymbolNumber();
}

void HMMBase::calcSymbolInterval(){	
	_symbolInterval = double(_maxValue-_minValue)/double(_symbolNumber);
}

void HMMBase::extractMinValue(){
	_minValue = _obsAlphabet.at(0)-_symbolInterval;
}

void HMMBase::extractMaxValue(){
	_maxValue = _obsAlphabet.at(_symbolNumber-1);
}

void HMMBase::extractSymbolInterval(){
	_symbolInterval = double(_obsAlphabet.at(1)-_obsAlphabet.at(0));
}

void HMMBase::calcObsAlphabetFromParas(){
	checkAlphabet();
	for (unsigned i=0; i<_symbolNumber; i++){
		_obsAlphabet.at(i) = _minValue+((i+1)*_symbolInterval);
	}
}

void HMMBase::extractParasFromObsAlphabet(){
	checkAlphabet();
	extractSymbolInterval();
	extractMinValue();
	extractMaxValue();
}

void HMMBase::calcCDF(const SMLMS::Matrix& pdf, SMLMS::Matrix& cdf){
	double value;
	for (unsigned i=0; i<pdf.numberOfRows(); i++){
		value = 0.0;
		for (unsigned j=0; j<pdf.numberOfColumns(); j++){
			value += pdf.at(i,j);
			cdf.at(i,j,value);
		}
	}
}

void HMMBase::calcDof(){
	_dof = 0;
	/* add equilibrium degrees of freedom */
	_dof += _stateNumber-1;
	/* add transition Matrix degrees of freedom based upon Sriraman et al. J. Phys Chem. 2005 */
	_dof += (_stateNumber+2)*(_stateNumber-1)/2;
	/* add observation degrees of freedom */
	_dof += (_symbolNumber-1)*_stateNumber;
}

void HMMBase::calcDof(SMLMS::PhysicalModelBLD &model){
	_dof = 0;
	SMLMS::Matrix paraMat = model.paraMat();
	/* add equilibrium degrees of freedom */
	for (unsigned i=0; i<_stateNumber; i++){
		if(!paraMat.at(i,1))_dof += 1;
	}
	if (_dof)_dof -= 1;
	/* add transition Matrix degrees of freedom based upon Sriraman et al. J. Phys Chem. 2005 */
	_dof += (_stateNumber+2)*(_stateNumber-1)/2;
	/* add observation degrees of freedom */
	for (unsigned i=0; i<_stateNumber; i++){
		if (!paraMat.at(i,5))_dof += 1;
	}
}

void HMMBase::calcBic(unsigned n){
	/* Bayesian information criterion (BIC) */
	_bic = (_dof*log(n))-(2.0*_logLikelihood);
}

void HMMBase::calcAic(unsigned n){
	/* second order Akaike information criterion (AIC) */
	_aic = (2.0*_dof)-(2.0*_logLikelihood);
	_aic += (2.0*_dof*(_dof + 1))/(n - _dof - 1);
}

void HMMBase::calcModelSelection(unsigned n){
	calcBic(n);
	calcAic(n);
}

/* print functions */
void HMMBase::printStateNumber(){
	std::cout<<std::endl<<"Number of States: "<<_stateNumber<<std::endl;
}

void HMMBase::printMinValue(){
	std::cout<<std::endl<<"Minimum Observation Value: "<<_minValue<<std::endl;
}

void HMMBase::printMaxValue(){
	std::cout<<std::endl<<"Maximum Observation Value: "<<_maxValue<<std::endl;
}

void HMMBase::printSymbolNumber(){
	std::cout<<std::endl<<"Number of Symbols: "<<_symbolNumber<<std::endl;
}

void HMMBase::printSymbolInterval(){
	std::cout<<std::endl<<"Distance between Symbols: "<<_symbolInterval<<std::endl;
}

void HMMBase::printEquiPDF(){
	unsigned i;
	std::cout<<std::endl<<"Equilibrium Matrix (pi): "<<std::endl;
	for (i=0; i<_stateNumber; i++) std::cout<<_equiPDF.at(0,i)<<"\t";
	std::cout<<std::endl;
}

void HMMBase::printTransPDF(){
	unsigned i,j;
	std::cout<<std::endl<<"Transition Matrix (A): "<<std::endl;
	for (i=0; i<_stateNumber; i++){
		for(j=0; j<_stateNumber; j++) std::cout<<_transPDF.at(i,j)<<"\t";
		std::cout<<std::endl;
	}
	std::cout<<std::endl;
}

void HMMBase::printObsPDF(){
	unsigned i,j;
	std::cout<<std::endl<<"Observation Matrix (B):"<<std::endl;
	for (i=0; i<_stateNumber; i++){
		for(j=0; j<_symbolNumber; j++) std::cout<<_obsPDF.at(i,j)<<"\t";
		std::cout<<std::endl;
	}
	std::cout<<std::endl;
}

void HMMBase::printEquiCDF(){
	unsigned i;
	std::cout<<std::endl<<"Cumulative Equilibrium Matrix (pi): "<<std::endl;
	for (i=0; i<_stateNumber; i++) std::cout<<_equiCDF.at(0,i)<<"\t";
	std::cout<<std::endl;
}

void HMMBase::printTransCDF(){
	unsigned i,j;
	std::cout<<std::endl<<"Cumulative Transition Matrix (A): "<<std::endl;
	for (i=0; i<_stateNumber; i++){
		for(j=0; j<_stateNumber; j++) std::cout<<_transCDF.at(i,j)<<"\t";
		std::cout<<std::endl;
	}
	std::cout<<std::endl;
}

void HMMBase::printObsCDF(){
	unsigned i,j;
	std::cout<<std::endl<<"Cumulative Observation Matrix (B):"<<std::endl;
	for (i=0; i<_stateNumber; i++){
		for(j=0; j<_symbolNumber; j++) std::cout<<_obsCDF.at(i,j)<<"\t";
		std::cout<<std::endl;
	}
	std::cout<<std::endl;
}

void HMMBase::printObsAlphabet(){
	unsigned i;
	std::cout<<std::endl<<"Observation Alphabet: "<<std::endl;
	for (i=0; i<_symbolNumber; i++) std::cout<<_obsAlphabet.at(i)<<"\t";
	std::cout<<std::endl;
}

void HMMBase::printLogLikelihood(){
	std::cout<<std::endl<<"Log Likelihood: "<<_logLikelihood<<std::endl;
}

void HMMBase::printDof(){
	std::cout<<std::endl<<"DOF: "<<_dof<<std::endl;
}

void HMMBase::printBic(){
	std::cout<<std::endl<<"BIC: "<<_bic<<std::endl;
}

void HMMBase::printAic(){
	std::cout<<std::endl<<"AIC: "<<_aic<<std::endl;
}

void HMMBase::printHMM(){
	std::cout<<std::endl<<"HMM:"<<std::endl;
	printStateNumber();
	printSymbolNumber();
	printEquiPDF();
	printTransPDF();
	//printObsPDF();
	//printSymbolInterval();
	//printObsAlphabet();
	printLogLikelihood();
	printDof();
	printAic();
	printBic();
}/* printHMM*/

}/* SMLMS */
