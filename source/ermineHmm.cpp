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
#include "header/ermineHmm.hpp"
#include "header/smlmsMatrix.hpp"
#include "header/ermineJudi.hpp"
//#include "header/probabilityCalculus.hpp"

namespace SMLMS{
/* constructor */
HMM::HMM (){
	std::cout<<"HMM constructor called."<<std::endl;
}
/*
HMM::HMM(unsigned states, double size, double min, double max){
	std::cout<<"HMM constructor called"<<std::endl;
	//set everytging
	setStateNumber(states);
	setMinValue(min);
	setMaxValue(max);
	setBinSize(size);
	//checkHMM();
	//initHMM() states size min max
	
}*/

/* destructor */
HMM::~HMM (){
	std::cout<<"HMM removed from Heap!"<<std::endl;
}

/* copy constructor */
HMM::HMM(const HMM &obj){
	std::cout<<"HMM copy constructor called."<<std::endl;
	_stateNumber = obj._stateNumber;
	_minValue = obj._minValue;
	_maxValue = obj._maxValue;
	_binSize = obj._binSize;
	_stateVector = obj._stateVector;
	_transMatrix = obj._transMatrix;
	_obsAlphabet = obj._obsAlphabet;
	_obsMatrix = obj._obsMatrix;
	_obsParaMatrix = obj._obsParaMatrix;
}

/* elementary functions */
void HMM::setStateNumber(unsigned numberOfStates){
	_stateNumber = numberOfStates;
}

unsigned HMM::stateNumber(){
	return _stateNumber;
}

void HMM::setMinValue(double min){
	_minValue = min;
}

double HMM::minValue(void){
	return _minValue;
}

void HMM::setMaxValue(double max){
	_maxValue = max;
}

double HMM::maxValue(void){
	return _maxValue;
}

void HMM::setBinSize(double size){
	_binSize = size;
}

double HMM::binSize(void){
	return _binSize;
}

void HMM::setStateVector(std::vector<double> const &states){
	_stateVector = states;
}

std::vector<double> HMM::stateVector(){
	return _stateVector;
}

void HMM::setStateVectorCdf(std::vector<double> const&states){
	_stateVectorCdf = states;
}

std::vector<double> HMM::stateVectorCdf(void){
	return _stateVectorCdf;
}

void HMM::setTransMatrix(SMLMS::Matrix const &trans){
	_transMatrix.clearMatrix();
	_transMatrix = trans;
}

SMLMS::Matrix HMM::transMatrix(){
	return _transMatrix;
}

void HMM::setTransMatrixCdf(SMLMS::Matrix const &trans){
	_transMatrixCdf = trans;
}

SMLMS::Matrix HMM::transMatrixCdf(void){
	return _transMatrixCdf;
}

void HMM::setObsMatrix(SMLMS::Matrix const &obs){
	_obsMatrix.clearMatrix();
	_obsMatrix = obs;
}

SMLMS::Matrix HMM::obsMatrix(){
	return _obsMatrix;
}

void HMM::setObsMatrixCdf(SMLMS::Matrix const &obs){
	_obsMatrixCdf = obs;
}

SMLMS::Matrix HMM::obsMatrixCdf(void){
	return _obsMatrixCdf;
}

void HMM::setObsParaMatrix(SMLMS::Matrix const &para){
	_obsParaMatrix.clearMatrix();
	_obsParaMatrix = para;
}

SMLMS::Matrix HMM::obsParaMatrix(){
	return _obsParaMatrix;
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
			if(n==0){/* read state number */
				lineContent>>_stateNumber;
				HMM::initHMM();
			}
			// read initial probability: pi
			if(n==1){
				for (int i=0; i<_stateNumber; i++){
					lineContent>>_stateVector.at(i);
				}
			}
			// read transition probability: A
			if(n>1 && n<(2+_stateNumber)){
				row=n-2;
				for (int column=0; column<_stateNumber; column++){
					lineContent>>_transMatrix(row, column);
				}
			}
			// read observation parameter: B
			if(n>(1+_stateNumber) && n<(2+_stateNumber+_stateNumber)){
				row = n-_stateNumber-2;
				for (int column=0; column<4; column++){
					std::cout<<column<<std::endl;
					lineContent>>_obsParaMatrix(row, column);
				}
			}
			n +=1;
		}
	}	

}

void HMM::writeHMM(std::string const &name){
	std::ofstream outFile;
	outFile.open(name.data());
	outFile<<std::scientific;
	outFile<<std::setprecision(6);
	/* header line */
	outFile<<"# Hidden Markov Model"<<std::endl;
	outFile<<"# ermine format"<<std::endl;
	/* State Number */
	outFile<<"# number of states: N"<<std::endl;
	outFile<<stateNumber()<<std::endl;
	/* initial probability: pi*/
	outFile<<"# initial probability vector: pi"<<std::endl;
	for (int i=0; i<_stateNumber; i++){
		outFile<<_stateVector.at(i)<<"\t";
	}
	outFile<<std::endl;
	/* transition array: A*/
	outFile<<"# transition probability matrix: A"<<std::endl;
	for (int row=0; row<_stateNumber; row++){
		for (int column=0; column<_stateNumber; column++){
			outFile<<_transMatrix(row, column)<<"\t";
		}
		outFile<<std::endl;
	}
	/* observation parameter array: B*/
	outFile<<"# observation probability parameter matrix: B"<<std::endl;
	outFile<<"# A\tD[nm^2/s]\tdt[s]\tsigma[nm]"<<std::endl;
	for (int row=0; row<_stateNumber; row++){
		for (int column=0; column<4; column++){
			outFile<<_obsParaMatrix(row,column)<<"\t";
		}
		outFile<<std::endl;
	}
	/* close */
	outFile<<"# go ermine!"<<std::endl;
	outFile.close();

}

/* special functions */
void HMM::clearHMM(){
	_stateNumber=0;
	_minValue=0;
	_maxValue=0;
	_binSize=0;
	_stateVector.clear();
	_transMatrix.clearMatrix();
	_obsAlphabet.clear();
	_obsMatrix.clearMatrix();
	_obsParaMatrix.clearMatrix();
	_stateVectorCdf.clear();
	_transMatrixCdf.clearMatrix();
	_obsMatrixCdf.clearMatrix();
}

void HMM::checkAlphabet(void){
	if (_obsAlphabet.empty()){
		throw std::out_of_range(" Cannot work with an HMM without an observation alphabet.");
	}
}

void HMM::checkStateNumber(){
	if(_stateNumber<1){
		throw std::out_of_range(" Cannot work with an HMM without a state.");
	}
}

void HMM::checkHMM(){
	checkAlphabet();
	checkStateNumber();

}

void HMM::initStateVector(){
	double stateProb;
	checkStateNumber();
	stateProb = 1.0/_stateNumber;
	_stateVector.clear();
	for (unsigned i=0; i<_stateNumber; i++){
		_stateVector.push_back(stateProb);
	}
}

void HMM::initTransMatrix(){
	checkStateNumber();
	_transMatrix.setNumberOfRows(_stateNumber);
	_transMatrix.setNumberOfColumns(_stateNumber);
	_transMatrix.calcNumberOfElements();
	_transMatrix.initMatrix();	
}	

void HMM::initObsAlphabet(double size, double min, double max){
	setBinSize(size);
	setMinValue(min);
	setMaxValue(max);
	if(_maxValue==0 || _binSize==0){
		throw std::out_of_range(" Cannot work with an HMM without an observation alphabet.");
	}
	_obsAlphabet.clear();
	double binNumber = (_maxValue -_minValue)/_binSize;
	for (int i=0; i<=binNumber; i++){
		_obsAlphabet.push_back(_binSize*i);
	}
}

void HMM::initObsParaMatrix(){
	checkStateNumber();
	_obsParaMatrix.setNumberOfRows(_stateNumber);
	_obsParaMatrix.setNumberOfColumns(4);/* temporär für judi */
	_obsParaMatrix.calcNumberOfElements();
	_obsParaMatrix.initMatrix();
}

void HMM::initObsMatrix(){
	checkHMM();
	_obsMatrix.setNumberOfColumns(_obsAlphabet.size());
	_obsMatrix.setNumberOfRows(_stateNumber);
	_obsMatrix.calcNumberOfElements();
	_obsMatrix.initMatrix();
}

void HMM::initStateVectorCdf(){
	double stateProb;
	checkStateNumber();
	stateProb = 1.0/_stateNumber;
	_stateVectorCdf.clear();
	for (unsigned i=0; i<_stateNumber; i++){
		_stateVectorCdf.push_back(stateProb);
	}
}

void HMM::initTransMatrixCdf(){
	checkStateNumber();
	_transMatrixCdf.setNumberOfRows(_stateNumber);
	_transMatrixCdf.setNumberOfColumns(_stateNumber);
	_transMatrixCdf.calcNumberOfElements();
	_transMatrixCdf.initMatrix();	
}	

void HMM::initObsMatrixCdf(){
	checkHMM();
	_obsMatrixCdf.setNumberOfColumns(_obsAlphabet.size());
	_obsMatrixCdf.setNumberOfRows(_stateNumber);
	_obsMatrixCdf.calcNumberOfElements();
	_obsMatrixCdf.initMatrix();
}


void HMM::initHMM(){
	checkHMM();
	initStateVector();
	initTransMatrix();
	initObsParaMatrix();
	initObsMatrix();
	initStateVectorCdf();
	initTransMatrixCdf();
	initObsMatrixCdf();
}

void HMM::calcObsMatrix(void){
        std::vector<double> tempObsPara(5);
        std::vector<double> tempPdf(_obsAlphabet.size());
        for (int i=0; i<_stateNumber; i++){
                tempObsPara=HMM::getSingleStateObsPara(i); 
                HMM::calcPdfFromPara(tempPdf, tempObsPara); 
                for (int j=0; j<tempPdf.size(); j++){
                        HMM::_obsMatrix(i,j)=tempPdf.at(j);
		}
	}
}

void HMM::calcObsMatrixCdf(void){
	if(HMM::_obsMatrixCdf.numberOfElements()<=0){
		throw std::out_of_range("obsMatrixCdf needs to be initialized before calculated!");
	}
	for (int state=0; state<_stateNumber; state++){			
		_obsMatrixCdf(state, 0)=_obsMatrix(state,0)*_binSize;
		for (int i=1; i<_obsMatrixCdf.numberOfColumns(); i++){
				_obsMatrixCdf(state, i)=_obsMatrixCdf(state,i-1)+(_obsMatrix(state,i)*_binSize);
		}
	}

}

void HMM::calcStateVectorCdf(void){
	if (HMM::_stateVector.size()!=_stateNumber){
		throw std::out_of_range("state Vector Cdf needs to be initialized before calculated!");
	}
	_stateVectorCdf.at(0) = _stateVector.at(0);
	for (int i=1; i<_stateNumber; i++){
		_stateVectorCdf.at(i) = _stateVectorCdf.at(i-1)+_stateVector.at(i); 
	}
}

void HMM::calcTransMatrixCdf(void){
	if(HMM::_transMatrixCdf.numberOfElements()!=_stateNumber*_stateNumber){ 
		throw std::out_of_range("transMatrixCdf needs to be initialized before calculated!");
	}
	for (int state=0; state<_stateNumber; state++){
		_transMatrixCdf(state, 0)=_transMatrix(state,0);
		for (int i=1; i<_stateNumber; i++){
			_transMatrixCdf(state,i)=_transMatrixCdf(state,i-1)+_transMatrix(state,i);
		}
	}
}

void HMM::calcCdf(void){
	HMM::calcStateVectorCdf();
	HMM::calcTransMatrixCdf();
	HMM::calcObsMatrixCdf();
}

/* core functions */
void HMM::initialize(SMLMS::JumpDistanceList &jumpData, std::string &name){
	checkHMM();
	std::cout<<_binSize<<std::endl;
	std::cout<<_minValue<<std::endl;
	std::cout<<_maxValue<<std::endl;
	std::cout<<_obsParaMatrix.numberOfColumns()<<std::endl;
	std::cout<<_obsParaMatrix.numberOfRows()<<std::endl;
	std::cout<<_obsParaMatrix.numberOfElements()<<std::endl;
	//SMLMS::BrownianLatDiff pdfCalc(_binSize, _minValue, _maxValue, _obsParaMatrix);
	//std::vector<double> jumpList;
	//jumpData.getAllJumps(jumpList);
	//pdfCalc.calcPdfByHist(jumpList);
	//pdfCalc.fitPdf(name);
}

void HMM::simulate(unsigned stepNumber, unsigned traceNumber, SMLMS::JumpDistanceList &simJudi){
	checkHMM();
	calcObsMatrix();
	normalizeHMM();

	//simStateSequence(stepNumber, traceNumber, simJudi);
	//for i in trace:
		//setTraceNumber
		//calc first state
		//for j in stepNumber
			//setTraceNumber
			//calc init state based on prior
	//simObsSequence:
	//for i in jumps
		//getStateNumber
		//calc obs based on state number
	//SMLMS::BrownianLatDiff pdfCalc(_binSize, _minValue, _maxValue, _obsParaMatrix);
	
}

}/* SMLMS */

