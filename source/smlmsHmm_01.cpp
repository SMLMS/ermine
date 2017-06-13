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
#include "header/smlmsHmm_01.hpp"
#include "header/smlmsMatrix.hpp"
#include "header/smlmsRandom.hpp"

namespace SMLMS{
/* constructor */
HMM::HMM (){
	std::cout<<"HMM constructor called."<<std::endl;
	_fbDone = false;
	_logLikelihood = 0.0;
}

HMM::HMM(unsigned states, unsigned symbols, unsigned obs){
	std::cout<<"HMM constructor called"<<std::endl;
	//set everytging
	setStateNumber(states);
	setSymbolNumber(symbols);
	setObsNumber(obs);
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
	//_obsSeq = obj._obsSeq;
	//_stateSeq = obj._stateSeq;
	//_obsProb = obj._obsProb;
	_obsNumber = obj._obsNumber;
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
	_normLH = obj._normLH;
	_alpha = obj._alpha;
	_beta = obj._beta;
	_xi = obj._xi;
	_gamma = obj._gamma;
	_logLikelihood = obj._logLikelihood;
	_fbDone = obj._fbDone;
}

/* elementary functions */

/*
void HMM::setObsSeq(std::vector<double> obs){
	_obsSeq = obs;
}

std::vector<double> HMM::obsSeq(){
	return _obsSeq;
}

void HMM::setStateSeq(std::vector<int> states){
	_stateSeq = states;
}

std::vector<int> HMM::stateSeq(){
	return _stateSeq;
}
*/

void HMM::setObsNumber(unsigned numberOfObs){
	_obsNumber = numberOfObs;
}

unsigned HMM::obsNumber(){
	return _obsNumber;
}

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

void HMM::setSymbolInterval(double data){
	_symbolInterval = data;
}

double HMM::symbolInterval(){
	return _symbolInterval;
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

std::vector<double> HMM::normLH(){
	return _normLH;
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
		SMLMS::ErmineError	smlmsHmmError(errorMessage.str());
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

/* clear functions */
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

/* special functions */
/* proof functions */
void HMM::checkObsNumber(){
	if(_obsNumber<1){
		std::stringstream errorMessage;
		errorMessage<<"Cannot setup a hmm without a definite number of observations!"<<std::endl;
		SMLMS::ErmineError smlmsErmineError(errorMessage.str());
		throw smlmsErmineError;
	}
}

void HMM::checkObsLength(unsigned seqSize){
	if(_obsNumber != seqSize){
		std::stringstream errorMessage;
		errorMessage<<"Observation sequence does not match pre-set sequence length!"<<std::endl;
		SMLMS::ErmineError smlmsErmineError(errorMessage.str());
		throw smlmsErmineError;
	}
}

void HMM::checkStateNumber(){
	if(_stateNumber<1){
		std::stringstream errorMessage;
		errorMessage<<"Cannot setup a hmm without a definite state number!"<<std::endl;
		SMLMS::ErmineError smlmsHmmError(errorMessage.str());
		throw smlmsHmmError;
	}
}

void HMM::checkSymbolNumber(){
	if(_symbolNumber<1){
		std::stringstream errorMessage;
		errorMessage<<"Cannot setup a hmm without a definite symbol number!"<<std::endl;
		SMLMS::ErmineError smlmsHmmError(errorMessage.str());
		throw smlmsHmmError;
	}
}

void HMM::checkEqui(){
	if ((_stateNumber != _equiPDF.numberOfColumns()) || (_equiPDF.numberOfRows()!=1)){
		std::stringstream errorMessage;
		errorMessage<<"equilibrum Matrix Mismatch: Matrix needs to be: "<<_stateNumber<<" x "<<1<<std::endl;
		SMLMS::ErmineError smlmsHmmError(errorMessage.str());
		throw smlmsHmmError;
	}
}
	

void HMM::checkTrans(){
	if ((_stateNumber != _transPDF.numberOfColumns()) ||(_stateNumber != _transPDF.numberOfRows())){
		std::stringstream errorMessage;
		errorMessage<<"transition Matrix Mismatch: Matrix needs to be: "<<_stateNumber<<" x "<<_stateNumber<<std::endl;
		SMLMS::ErmineError smlmsHmmError(errorMessage.str());
		throw smlmsHmmError;
	}
}

void HMM::checkObs(){
	if ((_symbolNumber != _obsPDF.numberOfColumns()) ||(_stateNumber != _obsPDF.numberOfRows())){
		std::stringstream errorMessage;
		errorMessage<<"observation Matrix Mismatch: Matrix needs to be: "<<_stateNumber<<" x "<<_symbolNumber<<std::endl;
		SMLMS::ErmineError smlmsHmmError(errorMessage.str());
		throw smlmsHmmError;
	}
	if (_obsAlphabet.size() != _symbolNumber){
		std::stringstream errorMessage;
		errorMessage<<"observation Alphabet Error: no Observation Alphabet given!"<<_stateNumber<<" x "<<_symbolNumber<<std::endl;
		SMLMS::ErmineError smlmsHmmError(errorMessage.str());
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
			std::cout<<"probabilty desnity warning: pdf integral is not normalized to unity."<<std::endl;
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

/* init functions */
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
	for (i=0; i<_symbolNumber; i++){
		tempVec.at(i)=double(i);
	}
	_obsAlphabet = tempVec;
	HMM::calcCDF(_obsPDF, _obsCDF);
}

/*
void HMM::initObsProb(){
	_obsProb.clearMatrix();
	SMLMS::Matrix tempMat(_stateNumber, _obsNumber);
	_obsProb = tempMat;
	resetObsProb();
}

void HMM::resetObsProb(){
	int i,t;
	for (int i=0; i<_stateNumber; i++){
		for (int t=0; t<_stateNumber; t++){
			_obsProb(i,t)=0.0;
		}
	}
}
*/

void HMM::initHMM(){
	initEqui();
	initTrans();
	initObs();
	//initObsProb();
}

void HMM::initAlpha(){
	_alpha.clearMatrix();
	SMLMS::Matrix tempMat(_stateNumber, _obsNumber);
	_alpha = tempMat;
	resetAlpha();
}

void HMM::resetAlpha(){
	int i,t;
	for (i=0; i<_stateNumber; i++){
		for (t=0; t<_obsNumber; t++){
			_alpha.at(i,t,0.0);
		}
	}
}

void HMM::initNormLH(){
	_normLH.clear();
	std::vector<double> tempVec(_obsNumber);
	_normLH = tempVec;
	resetNormLH();
}

void HMM::resetNormLH(){
	int t;
	for (t=0; t<_obsNumber; t++){
		_normLH.at(t)=0.0;
	}
}

void HMM::initBeta(){
	_beta.clearMatrix();
	SMLMS::Matrix tempMat(_stateNumber, _obsNumber);
	_beta = tempMat;
	resetBeta();
}

void HMM::resetBeta(){
	int i,t;
	for (i=0; i<_stateNumber; i++){
		for (t=0; t<_obsNumber; t++){
			_beta.at(i,t,0.0);
		}
	}
}

void HMM::initXi(){
	int i,j,t;
	_xi.clear();
	SMLMS::Matrix tempMat(_stateNumber, _stateNumber);
	for (i=0; i<_stateNumber; i++){
		for (j=0; j<_stateNumber; j++){
			tempMat.at(i,j,0.0);
		}
	}
	for (t=0; t<_obsNumber; t++){
		_xi.push_back(tempMat);
	}
}

void HMM::resetXi(){
	int i,j,t;
	SMLMS::Matrix tempMat(_stateNumber, _stateNumber);
	for (i=0; i<_stateNumber; i++){
		for (j=0; j<_stateNumber; j++){
			tempMat.at(i,j,0.0);
		}
	}
	for (t=0; t<_obsNumber; t++){
		_xi.at(t)=tempMat;
	}
}

void HMM::initGamma(){
	_gamma.clearMatrix();
	SMLMS::Matrix tempMat(_stateNumber, _obsNumber);
	_gamma = tempMat;
	resetGamma();
}

void HMM::resetGamma(){
	int i,t;
	for (i=0; i<_stateNumber; i++){
		for (t=0; t<_obsNumber; t++){
			_gamma.at(i,t,0.0);
		}
	}
}

void HMM::initAnalysis(){
	initAlpha();
	initNormLH();
	initBeta();
	initXi();
	initGamma();
}

void HMM::resetAnalysis(){
	//resetObsProb();
	resetAlpha();
	resetNormLH();
	resetBeta();
	resetXi();
	resetGamma();
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
			errorMessage<<"PDF error: Integrated Area = 0!"<<std::endl;
			SMLMS::ErmineError smlmsHmmError(errorMessage.str());
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

void HMM::calcCDF(const SMLMS::Matrix& pdf, SMLMS::Matrix& cdf){
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

void HMM::printObsNumber(){
	std::cout<<std::endl<<"Number of Observations: "<<_obsNumber<<std::endl;
}

void HMM::printStateNumber(){
	std::cout<<std::endl<<"Number of States: "<<_stateNumber<<std::endl;
}

void HMM::printSymbolNumber(){
	std::cout<<std::endl<<"Number of Symbols: "<<_symbolNumber<<std::endl;
}

/*
void HMM::printObsSeq(){
	int t;
	std::cout<<std::endl<<"Observation sequence: "<<std::endl;
	for (t=0; t<_obsNumber; t++) std::cout<<_obsSeq.at(t)<<"\t";
	std::cout<<std::endl;
}

void HMM::printObsProb(){
	int i,t;
	std::cout<<std::endl<<"Observation Sequence Probability by State: "<<std::endl;
	for (t=0; t<_obsNumber; t++){
		for (i=0; i<_stateNumber; i++) std::cout<<_obsProb.at(i,t)<<"\t";
		std::cout<<std::endl;
	}
	std::cout<<std::endl;
}

void HMM::printStateSeq(){
	int t;
	std::cout<<std::endl<<"State sequence: "<<std::endl;
	for (t=0; t<_obsNumber; t++) std::cout<<_stateSeq.at(t)<<"\t";
	std::cout<<std::endl;
}
*/

void HMM::printEquiPDF(){
	int i;
	std::cout<<std::endl<<"Equilibrium Matrix (pi): "<<std::endl;
	for (i=0; i<_stateNumber; i++) std::cout<<_equiPDF.at(0,i)<<"\t";
	std::cout<<std::endl;
}

void HMM::printTransPDF(){
	int i,j;
	std::cout<<std::endl<<"Transition Matrix (A): "<<std::endl;
	for (i=0; i<_stateNumber; i++){
		for(j=0; j<_stateNumber; j++) std::cout<<_transPDF.at(i,j)<<"\t";
		std::cout<<std::endl;
	}
	std::cout<<std::endl;
}

void HMM::printObsPDF(){
	int i,j;
	std::cout<<std::endl<<"Observation Matrix (B):"<<std::endl;
	for (i=0; i<_stateNumber; i++){
		for(j=0; j<_symbolNumber; j++) std::cout<<_obsPDF.at(i,j)<<"\t";
		std::cout<<std::endl;
	}
	std::cout<<std::endl;
}

void HMM::printEquiCDF(){
	int i;
	std::cout<<std::endl<<"Cumulative Equilibrium Matrix (pi): "<<std::endl;
	for (i=0; i<_stateNumber; i++) std::cout<<_equiCDF.at(0,i)<<"\t";
	std::cout<<std::endl;
}

void HMM::printTransCDF(){
	int i,j;
	std::cout<<std::endl<<"Cumulative Transition Matrix (A): "<<std::endl;
	for (i=0; i<_stateNumber; i++){
		for(j=0; j<_stateNumber; j++) std::cout<<_transCDF.at(i,j)<<"\t";
		std::cout<<std::endl;
	}
	std::cout<<std::endl;
}

void HMM::printObsCDF(){
	int i,j;
	std::cout<<std::endl<<"Cumulative Observation Matrix (B):"<<std::endl;
	for (i=0; i<_stateNumber; i++){
		for(j=0; j<_symbolNumber; j++) std::cout<<_obsCDF.at(i,j)<<"\t";
		std::cout<<std::endl;
	}
	std::cout<<std::endl;
}

void HMM::printObsAlphabet(){
	int i;
	std::cout<<std::endl<<"Observation Alphabet: "<<std::endl;
	for (i=0; i<_symbolNumber; i++) std::cout<<_obsAlphabet.at(i)<<"\t";
	std::cout<<std::endl;
}

void HMM::printAlpha(){
	int i,t;
	std::cout<<std::endl<<"Alpha: "<<std::endl;
	for (t=0; t<_obsNumber; t++){
		for(i=0; i<_stateNumber;i++) std::cout<<_alpha.at(i,t)<<"\t";
		std::cout<<std::endl;
	}
	std::cout<<std::endl;
}

void HMM::printBeta(){
	int i,t;
	std::cout<<std::endl<<"Beta: "<<std::endl;
	for (t=0; t<_obsNumber; t++){
		for(i=0; i<_stateNumber;i++) std::cout<<_beta.at(i,t)<<"\t";
		std::cout<<std::endl;
	}
	std::cout<<std::endl;
}

void HMM::printNormLH(){
	int t;
	std::cout<<std::endl<<"Normalization Factor: "<<std::endl;
	for (t=0; t<_obsNumber; t++) std::cout<<_normLH.at(t)<<std::endl;
	std::cout<<std::endl;
}

void HMM::printXi(){
	int i,j,t;
	std::cout<<std::endl<<"Xi: "<<std::endl;
	for (t=0; t<_obsNumber; t++){
		for (i=0; i<_stateNumber; i++){
			for (j=0; j<_stateNumber; j++) std::cout<<_xi.at(t).at(i,j)<<"\t";
		}
		std::cout<<std::endl;
	}
	std::cout<<std::endl;
}

void HMM::printGamma(){
	int i,t;
	std::cout<<std::endl<<"Gamma: "<<std::endl;
	for (t=0; t<_obsNumber; t++){
		for(i=0; i<_stateNumber; i++) std::cout<<_gamma.at(i,t)<<"\t";
		std::cout<<std::endl;
	}
	std::cout<<std::endl;
}

void HMM::printLogLikelihood(){
	std::cout<<std::endl<<"Log Likelihood: "<<_logLikelihood<<std::endl;
}

void HMM::printHMM(){
	std::cout<<std::endl<<"HMM:"<<std::endl;
	printObsNumber();
	printStateNumber();
	printSymbolNumber();
	printEquiPDF();
	printTransPDF();
	printObsPDF();
	printObsAlphabet();
	printLogLikelihood();
}/* printHMM*/

/* help functions*/
int HMM::findMatch(const SMLMS::Matrix& cdf, int state, double event){
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
		std::stringstream errorMessage;
		errorMessage<<"Out of range error in probability density function while claculating event probability: "<<event<<"exceeds number of columns = "<<cdf.numberOfColumns()<<std::endl;
		SMLMS::ErmineError smlmsHmmError(errorMessage.str());
		throw smlmsHmmError;
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
		std::stringstream errorMessage;
		errorMessage<<"Out of range error in probability density function while claculating event probability: "<<event<<"exceeds number of columns = "<<_obsAlphabet.size()<<std::endl;
		SMLMS::ErmineError smlmsHmmError(errorMessage.str());
		throw smlmsHmmError;
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

void HMM::calcXi(const std::vector<double> &obsSeq){
	int t,i,j;
	double denominator;
	SMLMS::Matrix tempMat(_stateNumber, _stateNumber);
	for (t=1; t<_obsNumber; t++){
		denominator = 0;
		for (i=0; i<_stateNumber; i++){
			for (j=0; j<_stateNumber; j++){
				denominator += tempMat(i,j)=_alpha.at(i, t-1) * _transPDF(i,j) * _obsPDF(j,obsPDFMatch(obsSeq.at(t)))*_beta.at(j,t);
			}
		}
		for (i=0; i<_stateNumber; i++){
			for (j=0; j<_stateNumber; j++){
				tempMat(i,j) /= denominator;
			}
		}
		_xi.at(t) = tempMat;
	}
}

void HMM::calcGamma(){
	int i,j,t,scaleDiff;
	double enumerator, denominator, tempGamma;
	for (t=1; t<_obsNumber; t++){
		denominator = 0.0;
		for(i=0; i<_stateNumber; i++) denominator += _gamma(i,t)= _alpha(i,t)*_beta(i,t);
		for(i=0; i<_stateNumber; i++) _gamma(i,t) /= denominator;
	}
}

void HMM::reestimateEquiPDF(){
	int i;
	for(i=0; i<_stateNumber; i++) _equiPDF(0,i) = _gamma(i,1);
}

void HMM::reestimateTransPDF(){
	int i,j,t;
	double enumerator, denominator;
	for (i=0; i<_stateNumber; i++){
		for (j=0;j<_stateNumber; j++){
			enumerator=0;
			denominator=0;
			for (t=1; t<_obsNumber; t++){
				enumerator += (_xi.at(t)).at(i,j);
				denominator += _gamma.at(i,t);
			}
			_transPDF(i,j) = enumerator/denominator;
		}
	}
}

void HMM::reestimateObsPDF(const std::vector<double> &obsSeq){
	int i,j,t;
	double enumerator, denominator;
	for (i=0; i<_stateNumber; i++){
		for (j=0;j<_symbolNumber; j++){
			enumerator=0;
			denominator=0;
			for (t=1; t<_obsNumber; t++){
				if (obsSymbolMatch(j, obsSeq.at(t))) enumerator += _gamma.at(i,t);
				denominator += _gamma.at(i,t);
			}
			_obsPDF(i,j) = enumerator/denominator;
		}
	}
}
/* core functions*/
/* core fundtions simulation */
void HMM::simulate(std::vector<double> &obsSeq, std::vector<int> &stateSeq){
	if (obsSeq.size() != stateSeq.size()){
	std::stringstream errorMessage;
	errorMessage<<"dimension error: observation vector and states vector must have the same dimensions."<<std::endl;
	SMLMS::ErmineError smlmsHmmError(errorMessage.str());
	throw smlmsHmmError;
	}	
	HMM::simulateStates(stateSeq);
	HMM::simulateObservations(obsSeq, stateSeq);
}

void HMM::simulateStates(std::vector<int> &stateSeq){
	double event = 0.0;
	int t, preState;
	checkObsLength(stateSeq.size());
	// simulate initial state by equiCDF
	event = _randGen.generateRandomDouble(0.0, 1.0);
	stateSeq.at(0) = findMatch(_equiCDF, 0, event); 
	for (t=1; t<_obsNumber; t++){
		event = _randGen.generateRandomDouble(0.0, 1.0);
		preState = stateSeq.at(t-1);
		stateSeq.at(t) = findMatch(_transCDF, preState, event);
	}
}

void HMM::simulateObservations(std::vector<double> &obsSeq, const std::vector<int> &stateSeq){
	double event = 0.0;
	int eventPDFPos, t;
	checkObsLength(obsSeq.size());
	checkObsLength(stateSeq.size());
	for(t=0; t<_obsNumber; t++){
		event = _randGen.generateRandomDouble(0.0, 1.0);
		eventPDFPos = findMatch(_obsCDF, stateSeq.at(t), event);
		obsSeq.at(t)=_obsAlphabet.at(eventPDFPos);
	}
}

/* core functions: HMM Estimate likelihood */
/*
void HMM::calcObsProb(){
	int state, obs, obsPDFMatch=0;
	for (state=0; state<_stateNumber; state++){
		for(obs=0; obs<_obsNumber; obs++){
			// Muss definitiev noch verallgemeinert werden
			obsPDFMatch = int(_obsSeq.at(obs))-1;
			_obsProb(state,obs)=_obsPDF.at(state, obsPDFMatch);
		}
	}
}
*/

void HMM::forwardAlgorithm(const std::vector<double> &obsSeq){
	int i,j,t;
	double tempNorm = 0.0;
	checkObsLength(obsSeq.size());
	// forward Procedure
        // alpha_t(i) = P(O_1 O_2 ... O_t, q_t = S_i | hmm)
      	
        // Initialize alpha
        // alpha_0(i) = pi(i)*obs(i)
	for (i=0; i<_stateNumber; i++){
		_alpha(i,0)=_equiPDF(0,i) * _obsPDF(i,obsPDFMatch(obsSeq.at(0)));
		tempNorm += _alpha(i,0);
	}
	_normLH.at(0)=1/tempNorm;
	// normalizing
	for (i=0; i<_stateNumber; i++){
		_alpha(i,0) *= _normLH.at(0);
	}
	// Induction of alpha
	// alpha_t(j) = sum(i) [alpha_t-1(i) * trans(ij) * obs_j(t)]
	for (t=1; t<_obsNumber; t++){
		tempNorm = 0.0;
		for (j=0; j<_stateNumber; j++){
			_alpha(j,t)=0;
			for (i=0; i<_stateNumber; i++){
				_alpha(j,t) += _alpha.at(i,t-1) * _transPDF.at(i,j) * _obsPDF(j,obsPDFMatch(obsSeq.at(t)));
			}
			tempNorm += _alpha(j,t);
		}
		_normLH.at(t) = 1/tempNorm;
		// normalization
		for (j=0; j<_stateNumber; j++){
			_alpha(j,t) *= _normLH.at(t);
		}
	}
}

void HMM::backwardAlgorithm(const std::vector<double> &obsSeq){
	int i,j,t;
	checkObsLength(obsSeq.size());
	// Backward Procedure
        // beta_t(i) = P(O_t+1 O_t+2 ... O_T | q_t = S_i , hmm)

        // Initialize beta
       	//_scaleBeta.at(obsNumber-1)=0.0;
        for (i=0; i<_stateNumber; i++){
        	_beta(i,_obsNumber-1)=1.0;
		// normalization
        	_beta(i,_obsNumber-1)=1.0*_normLH.at(_obsNumber-1);
	}
	
	// Induction of beta
	// beta_t-1(i) = sum(i) [trans(ij) * obs_j(t)]
	for (t=_obsNumber-1; t>0; t--){
		for (i=0; i<_stateNumber; i++){
			_beta(i,t-1) = 0.0;
			for (j=0; j<_stateNumber; j++){
				_beta(i,t-1) += _transPDF.at(i,j) * _obsPDF(j, obsPDFMatch(obsSeq.at(t))) * _beta.at(j,t);
			}
			// normalization
			_beta(i,t-1) *= _normLH.at(t-1); 
		}
	}
}



void HMM::forwardBackward(const std::vector<double> &obsSeq){
	initAnalysis();
	forwardAlgorithm(obsSeq);
	backwardAlgorithm(obsSeq);
	_fbDone = true;
}

void HMM::estimateLikelihood(){
	int t;
	double llc =0.0;
	for (t=0; t<_obsNumber; t++){
		llc += log(_normLH.at(t));
	}
	_logLikelihood = -1*llc;
}

void HMM::baumWelch(const std::vector<double> &obsSeq){
	//int i,j,t;
	//double tempSum;
	checkHMM();
	forwardBackward(obsSeq);
	calcXi(obsSeq);
	calcGamma();
	reestimateEquiPDF();
	reestimateTransPDF();
	reestimateObsPDF(obsSeq);
	normalizeHMM();
	estimateLikelihood();
	_fbDone=false;
}

void HMM::train(const std::vector<double> &obsSeq, double stopCrit, int maxIt){
	int n=0;
	double ll=0.0;
	double lldiff = 0.0;
	bool crit=true;
	while (crit){
		ll=_logLikelihood;
		baumWelch(obsSeq);
		n +=1;
		lldiff = fabs(ll-_logLikelihood);
		std::cout<<n<<"\t"<<_logLikelihood<<"\t"<<ll<<"\t"<<lldiff<<std::endl;
		if(lldiff <= stopCrit){
			std::cout<<std::endl<<"Training ended by reaching stop criterium."<<std::endl;
			crit = false;
		}
		else if (n>=maxIt){
			std::cout<<std::endl<<"Training ended by reaching maximum number iterations."<<std::endl;
			crit = false;
		}
	}
}
}/* SMLMS */

