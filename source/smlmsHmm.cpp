/* ######################################################################
* File Name: hmm.cpp
* Project: SMLMS
* Version: 16.03
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
#include <ctime>
#include <math.h>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_real_distribution.hpp>
#include <boost/random/variate_generator.hpp>
#include "header/ermineExceptions.hpp"
#include "header/smlmsHmm.hpp"
#include "header/smlmsMatrix.hpp"
#include "header/smlmsRandom.hpp"

namespace SMLMS{
/* constructor */
HMM::HMM (){
	std::cout<<"HMM constructor called."<<std::endl;
	_minl = exp(-99);
	_invMinl = exp(99);
	_fbDone = false;
	_logLikelihood = 0.0;
}

HMM::HMM(unsigned states, unsigned symbols){
	std::cout<<"HMM constructor called"<<std::endl;
	//set everytging
	setStateNumber(states);
	setSymbolNumber(symbols);
	_minl=exp(-99);
	_invMinl=exp(99);
	initHMM();
	checkHMM();
	_fbDone = false;
	_logLikelihood = 0.0;
}

/* destructor */
HMM::~HMM (){
	std::cout<<"HMM removed from Heap!"<<std::endl;
}

/* copy constructor */
HMM::HMM(const HMM &obj){
	std::cout<<"HMM copy constructor called."<<std::endl;
	_stateNumber = obj._stateNumber;
	_symbolNumber  = obj._symbolNumber;
	_equiPDF = obj._equiPDF;
	_transPDF = obj._transPDF;
	_obsPDF = obj._obsPDF;
	_equiCDF = obj._equiCDF;
	_transCDF = obj._transCDF;
	_obsCDF = obj._obsCDF;
	_obsAlphabet = obj._obsAlphabet;
}

/* elementary functions */
void HMM::setStateNumber(unsigned numberOfStates){
	_stateNumber = numberOfStates;
}

unsigned HMM::stateNumber(){
	return _stateNumber;
}

void HMM::setSymbolNumber(unsigned numberOfSymbols){
	_symbolNumber = numberOfSymbols;
}

unsigned HMM::symbolNumber(){
	return _symbolNumber;
}

void HMM::setEquiPDF(SMLMS::Matrix matrix){
	_equiPDF = matrix;
}

SMLMS::Matrix HMM::equiPDF(){
	return _equiPDF;
}

void HMM::setTransPDF(SMLMS::Matrix matrix){
	_transPDF = matrix;
}

SMLMS::Matrix HMM::transPDF(){
	return _transPDF;
}

void HMM::setObsPDF(SMLMS::Matrix matrix){
	_obsPDF = matrix;
}

SMLMS::Matrix HMM::obsPDF(){
	return _obsPDF;
}

void HMM::setObsAlphabet(std::vector<double> data){
	_obsAlphabet = data;
}

std::vector<double> HMM::obsAlphabet(){
	return _obsAlphabet;
}

SMLMS::Matrix HMM::obsProb(){
	return _obsProb;
}

double HMM::logLikelihood(){
	return _logLikelihood;
}

/* read/write functions */

void HMM::readHMM(std::string const &name){
	int n=0, row=0;
	std::string line;
	std::fstream inFile(name.data());
	if (inFile.is_open()){
		while(std::getline(inFile, line)){
			// skip comments
			if(!line.find("#")) continue;
			std::stringstream lineContent(line); 
			if (n==0){// read state number
				lineContent>>_stateNumber;
			}
			if (n==1){//read symbol number
				lineContent>>_symbolNumber;
				HMM::initHMM();
			}
			// read initial probability: pi
			if (n==2){
				for (int i=0; i<_stateNumber; i++){
					lineContent>>_equiPDF(0,1);
				}
			}
			// read transition probability: A
			if(n>2 && n<(3+_stateNumber)){
				row=n-3;
				for (int column=0; column<_stateNumber; column++){
					lineContent>>_transPDF(row, column);
				}
			}
			// read observation parameter: B
			if(n>(2+_stateNumber) && n<(3+_stateNumber+_stateNumber)){
				row = n-_stateNumber-3;
				for (int column=0; column<_symbolNumber; column++){
					lineContent>>_obsPDF(row, column);
				}
			}
			n +=1;
		}
	}
	else{
		std::stringstream errorMessage;
		errorMessage<<"cold not read hmm from file!"<<std::endl;
		SMLMS::SMLMSHmmError	smlmsHmmError(errorMessage.str());
		throw smlmsHmmError;
	}	

}

void HMM::writeHMM(std::string const &name){
	std::ofstream outFile;
	outFile.open(name.data());
	outFile<<std::scientific;
	outFile<<std::setprecision(6);
	// header line
	outFile<<"# Hidden Markov Model"<<std::endl;
	outFile<<"# ermine format"<<std::endl;
	// State Number
	outFile<<"# number of states: N"<<std::endl;
	outFile<<stateNumber()<<std::endl;
	// Symbol Number
	outFile<<"# number of observation symbols: S"<<std::endl;
	outFile<<symbolNumber()<<std::endl;
	// initial probability: pi
	outFile<<"# initial probability vector: pi"<<std::endl;
	for (int i=0; i<_stateNumber; i++){
		outFile<<_equiPDF.at(0,i)<<"\t";
	}
	outFile<<std::endl;
	// transition array: A
	outFile<<"# transition probability matrix: A"<<std::endl;
	for (int row=0; row<_stateNumber; row++){
		for (int column=0; column<_stateNumber; column++){
			outFile<<_transPDF.at(row, column)<<"\t";
		}
		outFile<<std::endl;
	}
	// observation parameter array: B
	outFile<<"# observation probability matrix: B"<<std::endl;
	for (int row=0; row<_stateNumber; row++){
		for (int column=0; column<_symbolNumber; column++){
			outFile<<_obsPDF.at(row,column)<<"\t";
		}
		outFile<<std::endl;
	}
	// close
	outFile<<"# go ermine!"<<std::endl;
	outFile.close();

}

void HMM::printHMM(){
	int i,j;
	//std::stringstream hmmOutput;
	std::cout<<"HMM:\n"
	<<"number of states: "<<_stateNumber<<std::endl
	<<"number of symbols: "<<_symbolNumber<<std::endl
	<<"equilibrium matrix (pi):"<<std::endl;
	for (i=0; i<_stateNumber;i++){
		std::cout<<_equiPDF.at(0,i)<<"\t";
	}
	std::cout<<std::endl;
	std::cout<<"transition Matrix (A):"<<std::endl;
	for (i=0; i<_stateNumber; i++){
		for (j=0; j<_stateNumber; j++){
			std::cout<<_transPDF.at(i,j)<<"\t";
		}
		std::cout<<std::endl;
	}
	std::cout<<"observation Matrix (B):"<<std::endl;
	for (i=0; i<_stateNumber; i++){
		for(j=0; j<_symbolNumber; j++){
			std::cout<<_obsPDF.at(i,j)<<"\t";
		}
		std::cout<<std::endl;
	}
	std::cout<<"observation Alphabet:"<<std::endl;
	for (i=0; i<_symbolNumber; i++ ){
		std::cout<<_obsAlphabet.at(i)<<"\t";
	}
	std::cout<<std::endl;
	std::cout<<"logLikelihood: "<<_logLikelihood<<std::endl;
	std::cout<<std::endl;
}/* printHMM*/


/* special functions */
void HMM::clearHMM(){
	_stateNumber=0;
	_symbolNumber=0;
	_equiPDF.clearMatrix();
	_transPDF.clearMatrix();
	_obsPDF.clearMatrix();
	_equiCDF.clearMatrix();
	_transCDF.clearMatrix();
	_obsCDF.clearMatrix();
	
}

void HMM::checkStateNumber(){
	if(_stateNumber<1){
		std::stringstream errorMessage;
		errorMessage<<"Cannot setup a hmm without a definite state number!"<<std::endl;
		SMLMS::SMLMSHmmError smlmsHmmError(errorMessage.str());
		throw smlmsHmmError;
	}
}

void HMM::checkSymbolNumber(){
	if(_symbolNumber<1){
		std::stringstream errorMessage;
		errorMessage<<"Cannot setup a hmm without a definite symbol number!"<<std::endl;
		SMLMS::SMLMSHmmError smlmsHmmError(errorMessage.str());
		throw smlmsHmmError;
	}
}

void HMM::checkEqui(){
	if ((_stateNumber != _equiPDF.numberOfColumns()) || (_equiPDF.numberOfRows()!=1)){
		std::stringstream errorMessage;
		errorMessage<<"equilibrum Matrix Mismatch: Matrix needs to be: "<<_stateNumber<<" x "<<1<<std::endl;
		SMLMS::SMLMSHmmError smlmsHmmError(errorMessage.str());
		throw smlmsHmmError;
	}
}
	

void HMM::checkTrans(){
	if ((_stateNumber != _transPDF.numberOfColumns()) ||(_stateNumber != _transPDF.numberOfRows())){
		std::stringstream errorMessage;
		errorMessage<<"transition Matrix Mismatch: Matrix needs to be: "<<_stateNumber<<" x "<<_stateNumber<<std::endl;
		SMLMS::SMLMSHmmError smlmsHmmError(errorMessage.str());
		throw smlmsHmmError;
	}
}

void HMM::checkObs(){
	if ((_symbolNumber != _obsPDF.numberOfColumns()) ||(_stateNumber != _obsPDF.numberOfRows())){
		std::stringstream errorMessage;
		errorMessage<<"observation Matrix Mismatch: Matrix needs to be: "<<_stateNumber<<" x "<<_symbolNumber<<std::endl;
		SMLMS::SMLMSHmmError smlmsHmmError(errorMessage.str());
		throw smlmsHmmError;
	}
	if (_obsAlphabet.size() != _symbolNumber){
		std::stringstream errorMessage;
		errorMessage<<"observation Alphabet Error: no Observation Alphabet given!"<<_stateNumber<<" x "<<_symbolNumber<<std::endl;
		SMLMS::SMLMSHmmError smlmsHmmError(errorMessage.str());
		throw smlmsHmmError;
	}
}

void HMM::checkProbDen(SMLMS::Matrix & mat){
	double prob;
	for (int i=0; i<mat.numberOfRows(); i++){
		prob = 0;
		for (int j=0; j<mat.numberOfColumns(); j++){
			prob += mat.at(i,j);
		}
		if ((prob < 0.99) || (prob > 1.01)){
			std::cout<<"probabilty desnity Error: pdf integral is not normalized to unity."<<std::endl;
		}
	}
}

void HMM::checkHMM(){
	checkStateNumber();
	checkSymbolNumber();
	checkEqui();
	checkTrans();
	checkObs();
	checkProbDen(_equiPDF);
	checkProbDen(_transPDF);
	checkProbDen(_obsPDF);
}

void HMM::initEqui(){
	_equiPDF.clearMatrix();
	_equiCDF.clearMatrix();
	SMLMS::Matrix tempMat(1, _stateNumber);
	_equiPDF = tempMat;
	_equiCDF = tempMat;
	for (int i=0; i<_stateNumber; i++){
		_equiPDF(0,i) = 1.0/_stateNumber;
	}
	HMM::calcCDF(_equiPDF, _equiCDF);
}

void HMM::initTrans(){
	_transPDF.clearMatrix();
	_transCDF.clearMatrix();
	SMLMS::Matrix tempMat(_stateNumber, _stateNumber);
	_transPDF = tempMat;
	_transCDF = tempMat;
	for (int j=0; j<_stateNumber; j++){
		for (int i=0; i<_stateNumber; i++){
			_transPDF(j,i)=1.0/_stateNumber;
		}
	}
	HMM::calcCDF(_transPDF, _transCDF);
}

void HMM::initObs(){
	int i,j;
	_obsPDF.clearMatrix();
	_obsCDF.clearMatrix();
	_obsAlphabet.clear();
	SMLMS::Matrix tempMat(_stateNumber, _symbolNumber);
	_obsPDF = tempMat;
	_obsCDF = tempMat;
	for (j=0; j<_stateNumber; j++){
		for (i=0; i<_symbolNumber; i++){
			_obsPDF(j,i) = 1.0/_symbolNumber;
		}
	}
	std::vector<double> tempVec(_symbolNumber);
	for (i=0; i<tempVec.size(); i++){
		tempVec.at(i)=double(i);
	}
	_obsAlphabet = tempVec;
	HMM::calcCDF(_obsPDF, _obsCDF);
}

void HMM::initHMM(){
	initEqui();
	initTrans();
	initObs();
}

void HMM::initAlpha(int obsNumber){
	int i,t;;
	_alpha.clearMatrix();
	SMLMS::Matrix tempMat(_stateNumber, obsNumber);
	_alpha = tempMat;
	for (i=0; i<_stateNumber; i++){
		for (t=0; t<obsNumber; t++){
			_alpha.at(i,t,0.0);
		}
	}
}

void HMM::initNorm(int obsNumber){
	int i;
	_norm.clear();
	std::vector<double> tempVec(obsNumber);
	for(i=0; i<obsNumber; i++) tempVec.at(i)=0.0;
	_norm = tempVec;
}

void HMM::initNormAlpha(int obsNumber){
	int i,t;
	_normAlpha.clearMatrix();
	SMLMS::Matrix tempMat(_stateNumber, obsNumber);
	for(i=0; i<_stateNumber; i++){
		for(t=0; t<obsNumber; t++){
			tempMat.at(i,t,0.0);
		}
	}
	_normAlpha = tempMat;
}

void HMM::initScaleAlpha(int obsNumber){
	int i;
	_scaleAlpha.clear();
	std::vector<int> tempVec(obsNumber);
	_scaleAlpha = tempVec;
	for (i=0; i<obsNumber; i++){
		_scaleAlpha.at(i)=0;
	}
}

void HMM::initBeta(int obsNumber){
	int i,t;
	_beta.clearMatrix();
	SMLMS::Matrix tempMat(_stateNumber, obsNumber);
	_beta = tempMat;
	for (i=0; i<_stateNumber; i++){
		for (t=0; t<obsNumber; t++){
			_beta.at(i,t,0.0);
		}
	}
}

void HMM::initNormBeta(int obsNumber){
	int i,t;
	_normBeta.clearMatrix();
	SMLMS::Matrix tempMat(_stateNumber, obsNumber);
	for(i=0; i<_stateNumber; i++){
		for(t=0; t<obsNumber; t++){
			tempMat.at(i,t,0.0);
		}
	}
	_normBeta = tempMat;
}

void HMM::initScaleBeta(int obsNumber){
	int i;
	_scaleBeta.clear();
	std::vector<int> tempVec(obsNumber);
	_scaleBeta = tempVec;
	for (i=0; i<obsNumber; i++){
		_scaleBeta.at(i)=0;
	}
}

void HMM::initObsProb(int obsNumber){
	int i,j;
	_obsProb.clearMatrix();
	SMLMS::Matrix tempMat(_stateNumber, obsNumber);
	_obsProb = tempMat;
	for (i=0; i<_stateNumber; i++){
		for (j=0; j<obsNumber; j++){
			_obsProb.at(i,j,0.0);
		}
	}
}

void HMM::initXi(){
	int i,j;
	_xi.clearMatrix();
	SMLMS::Matrix tempMat(_stateNumber, _stateNumber);
	_xi = tempMat;
	for (i=0; i<_stateNumber; i++){
		for (j=0; j<_stateNumber; j++){
			_xi.at(i,j,0.0);
		}
	}
}

void HMM::initXiVector(int obsNumber){
	int i,j,t;
	_Xi.clear();
	SMLMS::Matrix tempMat(_stateNumber, _stateNumber);
	for (i=0; i<_stateNumber; i++){
		for (j=0; j<_stateNumber; j++){
			tempMat.at(i,j,0.0);
		}
	}
	for (t=0; t<obsNumber; t++){
		_Xi.push_back(tempMat);
	}
}

void HMM::initGamma(int obsNumber){
	int i,t;
	_gamma.clearMatrix();
	SMLMS::Matrix tempMat(_stateNumber, obsNumber);
	_gamma = tempMat;
	for (i=0; i<_stateNumber; i++){
		for (t=0; t<obsNumber; t++){
			_gamma.at(i,t,0.0);
		}
	}
}

void HMM::initAnalysis(int obsNumber){
	initAlpha(obsNumber);
	initNorm(obsNumber);
	initNormAlpha(obsNumber);
	initScaleAlpha(obsNumber);
	initBeta(obsNumber);
	initNormBeta(obsNumber);
	initScaleBeta(obsNumber);
	initObsProb(obsNumber);
	initXi();
	initXiVector(obsNumber);
	initGamma(obsNumber);
}

/* normalization */
void HMM::normalizePDF(SMLMS::Matrix& mat){
	int i,j;
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
			errorMessage<<"PDF error: Area is 0!"<<std::endl;
			SMLMS::SMLMSHmmError smlmsHmmError(errorMessage.str());
			throw smlmsHmmError;
		}
	}
}

void HMM::normalizeHMM(){
	HMM::normalizePDF(_equiPDF);
	HMM::normalizePDF(_transPDF);
	HMM::normalizePDF(_obsPDF);
	HMM::calcCDF(_equiPDF, _equiCDF);
	HMM::calcCDF(_transPDF, _transCDF);
	HMM::calcCDF(_obsPDF, _obsCDF);
}
/* calc functionc */
void HMM::calcObsAlphabet(double interval){
	for (int i=0; i<_symbolNumber; i++){
		_obsAlphabet.at(i) = (i+1.0)*interval;
	}
}

void HMM::calcCDF(SMLMS::Matrix& pdf, SMLMS::Matrix& cdf){
	double value;
	for (int i=0; i<pdf.numberOfRows(); i++){
		value = 0.0;
		for (int j=0; j<pdf.numberOfColumns(); j++){
			value += pdf.at(i,j);
			cdf.at(i,j,value);
		}
	}
}

/* print functions */
void HMM::printAlpha(){
	int i,t, obsNumber;
	obsNumber = _normAlpha.numberOfColumns();
	std::cout<<"Alpha: "<<std::endl;
	for (t=0; t<obsNumber; t++){
		for(i=0; i<_stateNumber;i++){
			std::cout<<_normAlpha.at(i,t)<<"\t";
		}
		std::cout<<std::endl;
	}
	std::cout<<std::endl;
}

void HMM::printBeta(){
	int i,t, obsNumber;
	obsNumber = _normBeta.numberOfColumns();
	std::cout<<"Beta: "<<std::endl;
	for (t=0; t<obsNumber; t++){
		for(i=0; i<_stateNumber;i++){
			std::cout<<_normBeta.at(i,t)<<"\t";
		}
		std::cout<<std::endl;
	}
	std::cout<<std::endl;
}

void HMM::printNorm(){
	int t, obsNumber;
	obsNumber = _norm.size();
	std::cout<<"Norm: "<<std::endl;
	for (t=0; t<obsNumber; t++){
		std::cout<<_norm.at(t)<<std::endl;
	}
	std::cout<<std::endl;
}

void HMM::printXi(){
	int i,t, obsNumber;
	obsNumber = _xi.numberOfColumns();
	std::cout<<"xi: "<<std::endl;
	for (t=0; t<obsNumber; t++){
		for(i=0; i<_stateNumber; i++){
			std::cout<<_xi.at(i,t)<<"\t";
		}
		std::cout<<std::endl;
	}
	std::cout<<std::endl;
}

void HMM::printGamma(){
	int i,t, obsNumber;
	obsNumber = _gamma.numberOfColumns();
	std::cout<<"gamma: "<<std::endl;
	for (t=0; t<obsNumber; t++){
		for(i=0; i<_stateNumber; i++){
			std::cout<<_gamma.at(i,t)<<"\t";
		}
		std::cout<<std::endl;
	}
	std::cout<<std::endl;
}

/* help functions*/
int HMM::findMatch(SMLMS::Matrix& cdf, int state, double event){
	int eventPdfPosition=0;
	int i;
	if (event<=cdf.at(state,0)){
		eventPdfPosition=0;
	}
	else{
		for (i=1; i<cdf.numberOfColumns(); i++){
			if ((event>cdf.at(state,i-1)) && (event<=cdf.at(state,i))){
				eventPdfPosition=i;
			}
		}
	}
	if (event>cdf.at(state, cdf.numberOfColumns()-1)){
		throw std::out_of_range("Out of range error in probability density function while calculating event probability: event exceeds 1.");
	}
	return eventPdfPosition;
}

int HMM::obsPDFMatch(double event){
	int eventPdfPosition=0;
	int i;
	if (event<=_obsAlphabet.at(0)){
		eventPdfPosition=0;
	}
	else{
		for (i=1; i<_obsAlphabet.size(); i++){
			if ((event>_obsAlphabet.at(i-1)) && (event<=_obsAlphabet.at(i))){
				eventPdfPosition=i;
			}
		}
	}
	if (event>_obsAlphabet.at(_obsAlphabet.size()-1)){
		throw std::out_of_range("Out of range error in probability density function while calculating event probability: event exceeds 1.");
	}
	return eventPdfPosition;
}

bool HMM::obsSymbolMatch(int i, double event){
	if (i==0){
		if (event<=_obsAlphabet.at(0))	return true;
	}
	else{
		if ((event>_obsAlphabet.at(i-1)) && (event<=_obsAlphabet.at(i))) return true;
	}
	return false;
}

void HMM::calcXi2(std::vector<double>& observations){
	int obsNumber = observations.size();
	int t,i,j;
	double denominator;
	SMLMS::Matrix tempMat(_stateNumber, _stateNumber);
	for(t=1;t<obsNumber; t++){
		denominator = 0;
		for(i=0; i<_stateNumber; i++){
			for(j=0; j<_stateNumber; j++){
				denominator += tempMat(i,j)=_normAlpha.at(i, t-1) * _transPDF(i,j) * _obsPDF(j,obsPDFMatch(observations.at(t)))*_normBeta.at(j,t);
			}
		}
		for (i=0; i<_stateNumber; i++){
			for (j=0; j<_stateNumber; j++){
				tempMat(i,j) /= denominator;
			}
		}
		_Xi.at(t) = tempMat;
	}
}

void HMM::calcXi(){
	int obsNumber = _alpha.numberOfColumns();
	int t,i,j,k, scaleDiff;
	double enumerator, denominator, xiTempElement;
	for(i=0; i<_stateNumber; i++){
		for(j=0; j<_stateNumber; j++){
			xiTempElement = 0.0;
			enumerator = 0.0;
			for (t=1; t<obsNumber; t++){
				//enum(i,j) for timepoint t
				/*
				if (_scaleAlpha.at(t-1)!=_scaleAlpha.at(t)){
					scaleDiff = (_scaleAlpha.at(t)-_scaleAlpha.at(t-1));
					
				}
				else scaleDiff=0;
				//std::cout<<scaleDiff<<std::endl;
				
				enumerator += (_alpha.at(i,t-1)*pow(_invMinl,scaleDiff))* _transPDF.at(i,j)*_obsProb.at(j,t)*_beta.at(j,t);
				*/
				enumerator = _normAlpha.at(i,t-1)* _transPDF.at(i,j)*_obsProb.at(j,t)*_normBeta.at(j,t);
				denominator = 0.0;
				//for(k=0; k<_stateNumber; k++) denominator += _alpha(k,t)*_beta(k,t);
				for(k=0; k<_stateNumber; k++) denominator += _normAlpha(k,t)*_normBeta(k,t);
				xiTempElement += enumerator/denominator;
				//xiTempElement += enumerator;
			}
			_xi(i,j)=xiTempElement;
		}
	}
	/*
	// print
	std::cout<<"xi: "<<std::endl;
	for(i=0; i<_stateNumber; i++){
		for(j=0; j<_stateNumber; j++){
			std::cout<<_xi.at(i,j)<<"\t";
		}
		std::cout<<std::endl;
	}
	std::cout<<std::endl;
	*/
}

void HMM::calcGamma2(std::vector<double>& observations){
	int i,j,t,scaleDiff;
	int obsNumber = observations.size();
	double enumerator, denominator, tempGamma;
	for(t=1; t<obsNumber; t++){
		denominator = 0.0;
		//for(i=0; i<_stateNumber; i++) denominator += _normAlpha(i,t)*_normBeta(i,t);
		for(i=0; i<_stateNumber; i++) denominator += _gamma(i,t)= _normAlpha(i,t)*_normBeta(i,t);
		for(i=0; i<_stateNumber; i++){
			_gamma(i,t) /= denominator;
			/*
			enumerator = 0.0;
			for(j=0; j<_stateNumber; j++){
				enumerator += _normAlpha.at(i,t-1)* _transPDF.at(i,j)*_obsPDF(j,obsPDFMatch(observations.at(t)))*_normBeta.at(j,t);
			}
			_gamma(i,t) = enumerator/denominator;
			*/
		}
	}
}

void HMM::calcGamma(){
	int i,j,t,scaleDiff;
	int obsNumber = _alpha.numberOfColumns();
	double enumerator, denominator, tempGamma;
	for(t=1; t<obsNumber; t++){
		denominator = 0.0;
		//for(i=0; i<_stateNumber; i++) denominator += _alpha(i,t)*_beta(i,t);
		for(i=0; i<_stateNumber; i++) denominator += _normAlpha(i,t)*_normBeta(i,t);
		//for(i=0; i<_stateNumber; i++) denominator += _gamma(i,t)= _normAlpha(i,t)*_normBeta(i,t);
		/*
		std::cout<<"likelihood: "<<denominator<<std::endl;
		if (_scaleAlpha.at(t-1)!=_scaleAlpha.at(t)){
			scaleDiff = (_scaleAlpha.at(t)-_scaleAlpha.at(t-1));
		}
		else scaleDiff=0.0;
		*/
		for(i=0; i<_stateNumber; i++){
			//_gamma(i,t) /= denominator;
			
			enumerator = 0.0;
			for(j=0; j<_stateNumber; j++){
				//enumerator += (_alpha.at(i,t-1)*pow(_invMinl,scaleDiff))* _transPDF.at(i,j)*_obsProb.at(j,t)*_beta.at(j,t);
				enumerator += _normAlpha.at(i,t-1)* _transPDF.at(i,j)*_obsProb.at(j,t)*_normBeta.at(j,t);
			}
			_gamma(i,t) = enumerator/denominator;
			
		}
	}
	/*
	std::cout<<"minl: "<<_minl<<"\t"<<"minl^-1"<<_invMinl<<std::endl;
	std::cout<<"Gamma: "<<std::endl;
	for (t=0; t<obsNumber; t++){
		for (i=0;i<_stateNumber; i++){
			std::cout<<_gamma.at(i,t)<<"\t"<<_scaleAlpha.at(t)<<"\t";
		}
		std::cout<<std::endl;
	}
	std::cout<<std::endl;
	*/
}

void HMM::reestimateEquiPDF2(){
	int i;
	for(i=0; i<_stateNumber; i++) _equiPDF(0,i) = _gamma(i,1);
}

void HMM::reestimateEquiPDF(){
	int i,j, scaleDiff;
	double enumerator, denominator, tempPi;
	denominator = 0;
	//std::cout<<"new pi: "<<std::endl;
	for(i=0; i<_stateNumber; i++) denominator += _alpha(i,1)*_beta(i,1);
	if (_scaleAlpha.at(0)!=_scaleAlpha.at(1)) scaleDiff = (_scaleAlpha.at(1)-_scaleAlpha.at(0));
	else scaleDiff = 0;
	//std::cout<<"scale :"<<scaleDiff<<"\tfactor: "<<pow(_invMinl, scaleDiff)<<std::endl;
	for(i=0; i<_stateNumber; i++){
		tempPi = 0;
		for(j=0; j<_stateNumber; j++){
			enumerator = (_alpha.at(i,0)*pow(_invMinl,scaleDiff))*_transPDF.at(i,j)*_obsProb.at(j,1)*_beta.at(j,1);
			tempPi += enumerator/denominator; 
		}
		//std::cout<<tempPi<<"\t";
		_equiPDF(0,i)=tempPi;
	}
	//std::cout<<std::endl;
}

void HMM::reestimateTransPDF2(){
	int i,j,t;
	int obsNumber = _Xi.size();
	double enumerator, denominator;
	//SMLMS::Matrix tempMatrix(_stateNumber, _stateNumber);
	for (i=0; i<_stateNumber; i++){
		for(j=0;j<_stateNumber; j++){
			enumerator=0;
			denominator=0;
			for(t=1; t<obsNumber; t++){
				enumerator += (_Xi.at(t)).at(i,j);
				denominator += _gamma.at(i,t);
			}
			_transPDF(i,j) = enumerator/denominator;
		}
	}
}

void HMM::reestimateTransPDF(){
	int i,j,t;
	int obsNumber=_alpha.numberOfColumns();
	double scale;
	for(i=0; i<_stateNumber; i++){
		scale = 0.0;
		for(t=1; t<obsNumber; t++) scale += _gamma.at(i,t);
		for(j=0; j<_stateNumber; j++) _transPDF(i,j)=_xi.at(i,j)/scale;
	}
	/*
	std::cout<<std::endl;
	std::cout<<"new Trans PDF: "<<std::endl;
	for(i=0; i<_stateNumber; i++){
		for(j=0; j<_stateNumber; j++){
			std::cout<<tempMat.at(i,j)<<"\t";
		}
		std::cout<<std::endl;
	}
	std::cout<<std::endl;
	*/
}

void HMM::reestimateObsPDF2(std::vector<double>& observations){
	int i,j,t;
	int obsNumber = observations.size();
	double enumerator, denominator;
	//SMLMS::Matrix tempMatrix(_stateNumber, _stateNumber);
	for (i=0; i<_stateNumber; i++){
		for(j=0;j<_symbolNumber; j++){
			enumerator=0;
			denominator=0;
			for(t=1; t<obsNumber; t++){
				if (obsSymbolMatch(j, observations.at(t))) enumerator += _gamma.at(i,t);
				denominator += _gamma.at(i,t);
			}
			_obsPDF(i,j) = enumerator/denominator;
		}
	}
}

/* simulate hmm */
void HMM::simulate(std::vector<double>& observations, std::vector<int>& states){
	if (observations.size() != states.size()){
		std::stringstream errorMessage;
		errorMessage<<"dimension error: observation vector and states vector must have the same dimensions."<<std::endl;
		SMLMS::SMLMSHmmError smlmsHmmError(errorMessage.str());
		throw smlmsHmmError;
	}	
	HMM::simulateStates(states);
	HMM::simulateObservations(observations, states);
}

void HMM::simulateStates(std::vector<int>& states){
	double event = 0.0;
	int i, preState;
	// simulate initial state by equiCDF
	event = _randGen.generateRandomDouble(0.0, 1.0);
	states.at(0) = findMatch(_equiCDF, 0, event); 
	for (i=1; i<states.size(); i++){
		event = _randGen.generateRandomDouble(0.0, 1.0);
		preState = states.at(i-1);
		states.at(i) = findMatch(_transCDF, preState, event);
	}
}

void HMM::simulateObservations(std::vector<double>& observations, std::vector<int>& states){
	double event = 0.0;
	int eventPDFPos, i;
	for(i=0; i<observations.size(); i++){
		event = _randGen.generateRandomDouble(0.0, 1.0);
		eventPDFPos = findMatch(_obsCDF, states.at(i), event);
		observations.at(i)=_obsAlphabet.at(eventPDFPos);
	}
}

/* HMM Estimate likelihood */
void HMM::calcObsProb(std::vector<double>& observations){
	int state, obs, obsNumber, obsPDFMatch=0;
	obsNumber = observations.size();
	for (state=0; state<_stateNumber; state++){
		for(obs=0; obs<obsNumber; obs++){
			// Muss definitiev noch verallgemeinert werden
			obsPDFMatch = int(observations.at(obs))-1;
			_obsProb(state,obs)=_obsPDF.at(state, obsPDFMatch);
		}
	}
}

void HMM::forwardAlgorithm(std::vector<double>& observations){
	int i,j,obsNumber,t, underflow;
	double tempNorm;
	obsNumber = observations.size();
	// forward Procedure
        // alpha_t(i) = P(O_1 O_2 ... O_t, q_t = S_i | hmm)
      	
        // Initialize alpha
        // alpha_0(i) = pi(i)*obs(i)
        //_scaleAlpha.at(0)=0;
	underflow =0;
	tempNorm=0.0;
	for (i=0; i<_stateNumber; i++){
		_alpha(i,0)=_equiPDF(0,i) * _obsProb(i,0);
		_normAlpha(i,0)=_equiPDF(0,i) * _obsProb(i,0);
		tempNorm += _normAlpha(i,0);
		if (_alpha.at(i,0)<_minl) underflow = 1;
	}
	_norm.at(0)=1/tempNorm;
	// normalizing
	for(i=0; i<_stateNumber; i++){
		_normAlpha(i,0) *= _norm.at(0);
	}
	// scaling
	if (underflow){
		++_scaleAlpha.at(0);
		for(i=0; i<_stateNumber; i++){
			_alpha(i,0)*=_invMinl;
			++_scaleAlpha.at(0);
		}
	}

	// Induction of alpha
	// alpha_t(j) = sum(i) [alpha_t-1(i) * trans(ij) * obs_j(t)]
	for (t=1; t<obsNumber; t++){
		underflow = 0;
		_scaleAlpha.at(t)=_scaleAlpha.at(t-1);
		tempNorm = 0.0;
		for(j=0; j<_stateNumber; j++){
			for(i=0; i<_stateNumber; i++){
				_alpha(j,t) += _alpha.at(i,t-1) * _transPDF.at(i,j) * _obsProb(j,t);
				_normAlpha(j,t) += _normAlpha.at(i,t-1) * _transPDF.at(i,j) * _obsProb(j,t);
			}
			tempNorm += _normAlpha(j,t);
			if (_alpha.at(j,t)<_minl) underflow=1;
		}
		_norm.at(t) = 1/tempNorm;
		// normalization
		for (j=0; j<_stateNumber; j++){
			_normAlpha(j,t) *= _norm.at(t);
		}
		// scaling to prevent undeflow
		if (underflow){
			++_scaleAlpha.at(t);
			for(j=0; j<_stateNumber; j++){
				_alpha(j,t)*=_invMinl;
			}
		}
	}
}

void HMM::forwardAlgorithm2(std::vector<double>& observations){
	int i,j,t,obsNumber;
	double tempNorm = 0.0;
	obsNumber = observations.size();
	// forward Procedure
        // alpha_t(i) = P(O_1 O_2 ... O_t, q_t = S_i | hmm)
      	
        // Initialize alpha
        // alpha_0(i) = pi(i)*obs(i)
        //_scaleAlpha.at(0)=0;
	for (i=0; i<_stateNumber; i++){
		_normAlpha(i,0)=_equiPDF(0,i) * _obsPDF(i,obsPDFMatch(observations.at(0)));
		tempNorm += _normAlpha(i,0);
	}
	_norm.at(0)=1/tempNorm;
	// normalizing
	for(i=0; i<_stateNumber; i++){
		_normAlpha(i,0) *= _norm.at(0);
	}
	// Induction of alpha
	// alpha_t(j) = sum(i) [alpha_t-1(i) * trans(ij) * obs_j(t)]
	for (t=1; t<obsNumber; t++){
		tempNorm = 0.0;
		for(j=0; j<_stateNumber; j++){
			_normAlpha(j,t)=0;
			for(i=0; i<_stateNumber; i++){
				_normAlpha(j,t) += _normAlpha.at(i,t-1) * _transPDF.at(i,j) * _obsPDF(j,obsPDFMatch(observations.at(t)));
			}
			tempNorm += _normAlpha(j,t);
		}
		_norm.at(t) = 1/tempNorm;
		// normalization
		for (j=0; j<_stateNumber; j++){
			_normAlpha(j,t) *= _norm.at(t);
		}
	}
}

void HMM::backwardAlgorithm(std::vector<double>& observations){
	int i,j,t,obsNumber,underflow;
	obsNumber = observations.size();
	// Backward Procedure
        // beta_t(i) = P(O_t+1 O_t+2 ... O_T | q_t = S_i , hmm)

        // Initialize beta
       	//_scaleBeta.at(obsNumber-1)=0.0;
        for (i=0; i<_stateNumber; i++){
        	_beta(i,obsNumber-1)=1.0;
		// normalization
        	_normBeta(i,obsNumber-1)=1.0*_norm.at(obsNumber-1);
	}
	
	// Induction of beta
	// beta_t-1(i) = sum(i) [trans(ij) * obs_j(t)]
	for (t=obsNumber-1; t>0; t--){
		_scaleBeta.at(t-1) = _scaleBeta.at(t);
		underflow = 0;
		for (i=0; i<_stateNumber; i++){
			for (j=0; j<_stateNumber; j++){
				_beta(i,t-1) += _transPDF.at(i,j) * _obsProb.at(j,t) * _beta.at(j,t);
				_normBeta(i,t-1) += _transPDF.at(i,j) * _obsProb.at(j,t) * _normBeta.at(j,t);
			}
			// normalization
			_normBeta(i,t-1) *= _norm.at(t-1); 
			if (_beta.at(i,t-1)<_minl) underflow=1;
		}
		//scaling
		if (underflow){
			++_scaleBeta.at(t-1);
			for(i=0; i<_stateNumber; i++){
				_beta(i, t-1)*=_invMinl;
			}
		}
	}
}

void HMM::backwardAlgorithm2(std::vector<double>& observations){
	int i,j,t,obsNumber;
	obsNumber = observations.size();
	// Backward Procedure
        // beta_t(i) = P(O_t+1 O_t+2 ... O_T | q_t = S_i , hmm)

        // Initialize beta
       	//_scaleBeta.at(obsNumber-1)=0.0;
        for (i=0; i<_stateNumber; i++){
        	_beta(i,obsNumber-1)=1.0;
		// normalization
        	_normBeta(i,obsNumber-1)=1.0*_norm.at(obsNumber-1);
	}
	
	// Induction of beta
	// beta_t-1(i) = sum(i) [trans(ij) * obs_j(t)]
	for (t=obsNumber-1; t>0; t--){
		for (i=0; i<_stateNumber; i++){
			_normBeta(i,t-1) = 0.0;
			for (j=0; j<_stateNumber; j++){
				_normBeta(i,t-1) += _transPDF.at(i,j) * _obsPDF(j, obsPDFMatch(observations.at(t))) * _normBeta.at(j,t);
			}
			// normalization
			_normBeta(i,t-1) *= _norm.at(t-1); 
		}
	}
}

void HMM::forwardBackward(std::vector<double>& observations){
	int obsNumber = observations.size();
	initAnalysis(obsNumber);
	calcObsProb(observations);
	/*
	for (int i =0; i<obsNumber; i++){
		std::cout<<observations.at(i)<<"\t";
		for(int j=0; j<_stateNumber; j++){
			std::cout<<_obsProb.at(j,i)<<"\t";
		}
		std::cout<<std::endl;
	}
	*/
	forwardAlgorithm(observations);
	/*
	std::cout<<"alpha: "<<std::endl;
	for (int t=0; t<obsNumber;t++){
		for (int i=0; i<_stateNumber; i++){
			std::cout<<_alpha.at(i,t)<<" * minl^ ("<<_scaleAlpha.at(t)<<")\t";
		}
		std::cout<<std::endl;
	}
	*/
	backwardAlgorithm(observations);
	/*
	std::cout<<"\nbeta: "<<std::endl;
	for (int t=0; t<obsNumber;t++){
		for (int i=0; i<_stateNumber; i++){
			std::cout<<_beta.at(i,t)<<" * minl^ ("<<_scaleBeta.at(t)<<")\t";
		}
		std::cout<<std::endl;
	}
	std::cout<<std::endl;
	*/
	_fbDone = true;
}

void HMM::forwardBackward2(std::vector<double>& observations){
	int obsNumber = observations.size();
	initAnalysis(obsNumber);
	forwardAlgorithm2(observations);
	backwardAlgorithm2(observations);
	_fbDone = true;
}

void HMM::estimateLikelihood2(std::vector<double>& observations){
	int t;
	int obsNumber = observations.size();	
	double llc =0.0;
	for(t=0; t<obsNumber; t++){
		llc += log(_norm.at(t));
	}
	_logLikelihood = -1*llc;
	//std::cout<<"Log Likelihood: "<<_logLikelihood<<std::endl;
}

void HMM::estimateLikelihood(std::vector<double>& observations){
	int i,t;
	int obsNumber = observations.size();
	double tempL=0.0;
	double logLikelihood = 0.0;
	for (int t=0; t<obsNumber; t++){
		tempL = 0.0;
		for(int i=0; i<_stateNumber; i++){
			tempL += _alpha.at(i,t)*_beta.at(i,t);
		}
		logLikelihood += log(tempL)+((_scaleAlpha.at(t)+_scaleBeta.at(t))*log(_minl));
	}
	logLikelihood /= obsNumber;
	//_logLikelihood = logLikelihood;
	std::cout<<"log Likelihood nrc++: "<<_logLikelihood<<std::endl;
	double llc =0.0;
	for(t=0; t<obsNumber; t++){
		llc += log(_norm.at(t));
	}
	_logLikelihood = -1*llc;
	std::cout<<"log Likelihood c: "<<-1*llc<<std::endl;
	//std::cout<<"log Likelihood: "<<logLikelihood<<std::endl;
	/*	
	std::cout<<"likelihood: "<<std::endl;
	double prob=0.0;
	double test =.0;
	for (int t=0; t<obsNumber;t++){
		test = 0.0;
		for (int i=0; i<_stateNumber; i++){
			test += _alpha.at(i,t) * _beta.at(i,t);
			std::cout<<_alpha.at(i,t) * _beta.at(i,t)<<"\t";
		}
		std::cout<<test<<std::endl;
		prob += test;
	}
	std::cout<<std::endl;
	for (int t=0; t<obsNumber; t++){
		std::cout<<_scaleAlpha.at(t)<<"\t"<<_scaleBeta.at(t)<<std::endl;
	}
	*/
}

void HMM::baumWelch2(std::vector<double>& observations){
	int i,j,t, obsNumber;
	double tempSum;
	obsNumber = observations.size();
	//estimateLikelihood(observations);
	for (t=0; t<1000; t++){
		//std::cout<<"\niteration: "<<t+1<<"\t obs Number: "<<obsNumber<<std::endl<<std::endl;
		checkHMM();
		forwardBackward2(observations);
		//printNorm();
		calcXi2(observations);
		//printXi();
		calcGamma2(observations);
		//printGamma();
		reestimateEquiPDF2();
		reestimateTransPDF2();
		reestimateObsPDF2(observations);
		normalizeHMM();
		estimateLikelihood2(observations);
		//printHMM();
		_fbDone=false;
	}
	printHMM();
}

void HMM::baumWelch(std::vector<double>& observations){
	int i,j,t, obsNumber;
	double tempSum;
	obsNumber = observations.size();
	//estimateLikelihood(observations);
	for (t=0; t<10; t++){
		std::cout<<"\niteration: "<<t+1<<"\t obs Number: "<<obsNumber<<std::endl<<std::endl;
		checkHMM();
		forwardBackward(observations);
		//printNorm();
		calcXi();
		//printXi();
		calcGamma();
		//printGamma();
		reestimateEquiPDF();
		reestimateTransPDF();
		reestimateObsPDF2(observations);
		estimateLikelihood(observations);
		//printHMM();
		_fbDone = false;
		}
	printHMM();
		/*
		std::cout<<"log Likelihood: "<<logLikelihood()<<std::endl;
		std::cout<<"new Pi: "<<std::endl;
		for(i=0; i<_stateNumber; i++) std::cout<<_equiPDF.at(0,i)<<"\t";
		std::cout<<std::endl<<"new A: "<<std::endl;
		for(i=0; i<_stateNumber; i++){
			for(j=0; j<_stateNumber; j++){
				std::cout<<_transPDF.at(i,j)<<"\t";
			}
			std::cout<<std::endl;
		}
		*/
	/*
	std::cout<<"norm alpha:"<<std::endl;
	for(t=0; t<obsNumber; t++){
		tempSum=0;
		for(i=0; i<_stateNumber; i++){
			tempSum += _normAlpha.at(i,t);
			std::cout<<_normAlpha.at(i,t)<<"\t";
		}
		std::cout<<tempSum<<std::endl;
	}
	std::cout<<"norm beta:"<<std::endl;
	for (t=0; t<obsNumber; t++){
		tempSum = 0;
		for(i=0; i<_stateNumber; i++){
			tempSum += _normBeta.at(i,t);
			std::cout<<_normBeta.at(i,t)<<"\t";
		}
		std::cout<<tempSum<<std::endl;
	}
	*/
}
}/* SMLMS */

