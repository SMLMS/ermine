/* ######################################################################
* File Name: smlmsHmmUnique.cpp
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
#include <sstream>
#include <ctime>
#include <algorithm>
#include <math.h>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_real_distribution.hpp>
#include <boost/random/variate_generator.hpp>
#include "header/smlmsExceptions.hpp"
#include "header/smlmsHmmBase.hpp"
#include "header/smlmsHmmUnique.hpp"
#include "header/smlmsMatrix.hpp"
#include "header/smlmsRandom.hpp"

namespace SMLMS{
/* constructor */
HMMUnique::HMMUnique():SMLMS::HMMBase(){
	//std::cout<<"HMMUnique constructor called."<<std::endl;
	_fbDone = false;
}

HMMUnique::HMMUnique(unsigned states, unsigned symbols):SMLMS::HMMBase(states, symbols){
	//std::cout<<"HMMUnique constructor called"<<std::endl;
	//set everytging
	setObsNumber(0);
	_fbDone = false;
}

HMMUnique::HMMUnique(unsigned states, unsigned symbols, unsigned obs):SMLMS::HMMBase(states, symbols){
	//std::cout<<"HMMUnique constructor called"<<std::endl;
	//set everytging
	setObsNumber(obs);
	_fbDone = false;
}

/* copy constructor */
HMMUnique::HMMUnique(const HMMUnique &obj){
	//std::cout<<"HMMUnique copy constructor called."<<std::endl;
	_randGen = obj._randGen;
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
	_transPDFNumer = obj._transPDFNumer;
	_obsPDFNumer = obj._obsPDFNumer;
	_pdfDenominator = obj._pdfDenominator;
	_logLikelihood = obj._logLikelihood;
	_bic = obj._bic;
	_aic = obj._aic;
	_fbDone = obj._fbDone;
}

/* elementary functions */
void HMMUnique::setObsNumber(unsigned numberOfObs){
	_obsNumber = numberOfObs;
}

unsigned HMMUnique::obsNumber(){
	return _obsNumber;
}

std::vector<double> HMMUnique::normLH(){
	return _normLH;
}

SMLMS::Matrix HMMUnique::transPDFNumer(){
	return _transPDFNumer;
}

SMLMS::Matrix HMMUnique::obsPDFNumer(){
	return _obsPDFNumer;
}

SMLMS::Matrix HMMUnique::pdfDenominator(){
	return _pdfDenominator;
}

/* special functions */
/* proof functions */
void HMMUnique::checkObsNumber(){
	if(_obsNumber<1){
		std::stringstream errorMessage;
		errorMessage<<"Cannot setup a hmm without a definite number of observations!"<<std::endl;
		SMLMS::SmlmsError error(errorMessage.str());
		throw error;
	}
}

void HMMUnique::checkObsLength(unsigned seqSize){
	if(_obsNumber != seqSize){
		std::stringstream errorMessage;
		errorMessage<<"Observation sequence does not match pre-set sequence length!"<<std::endl;
		SMLMS::SmlmsError error(errorMessage.str());
		throw error;
	}
}

/* init functions */
void HMMUnique::initNormLH(){
	_normLH.clear();
	std::vector<double> tempVec(_obsNumber);
	_normLH = tempVec;
	resetNormLH();
}

void HMMUnique::resetNormLH(){
	unsigned t;
	for (t=0; t<_obsNumber; t++){
		_normLH.at(t)=0.0;
	}
}

void HMMUnique::initAlpha(){
	_alpha.clearMatrix();
	SMLMS::Matrix tempMat(_stateNumber, _obsNumber);
	_alpha = tempMat;
	resetAlpha();
}

void HMMUnique::resetAlpha(){
	unsigned i,t;
	for (i=0; i<_stateNumber; i++){
		for (t=0; t<_obsNumber; t++){
			_alpha.at(i,t,0.0);
		}
	}
}


void HMMUnique::initBeta(){
	_beta.clearMatrix();
	SMLMS::Matrix tempMat(_stateNumber, _obsNumber);
	_beta = tempMat;
	resetBeta();
}

void HMMUnique::resetBeta(){
	unsigned i,t;
	for (i=0; i<_stateNumber; i++){
		for (t=0; t<_obsNumber; t++){
			_beta.at(i,t,0.0);
		}
	}
}

void HMMUnique::initXi(){
	unsigned i,j,t;
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

void HMMUnique::resetXi(){
	unsigned i,j,t;
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

void HMMUnique::initGamma(){
	_gamma.clearMatrix();
	SMLMS::Matrix tempMat(_stateNumber, _obsNumber);
	_gamma = tempMat;
	resetGamma();
}

void HMMUnique::resetGamma(){
	unsigned i,t;
	for (i=0; i<_stateNumber; i++){
		for (t=0; t<_obsNumber; t++){
			_gamma.at(i,t,0.0);
		}
	}
}

void HMMUnique::initTransPDFNumer(){
	_transPDFNumer.clearMatrix();
	SMLMS::Matrix tempMat(_stateNumber, _stateNumber);
	_transPDFNumer = tempMat;
	resetTransPDFNumer();
}

void HMMUnique::resetTransPDFNumer(){
	unsigned i,j;
	for (i=0; i<_stateNumber; i++){
		for (j=0; j<_stateNumber; j++){
			_transPDFNumer.at(i,j,0.0);
		}
	}
}

void HMMUnique::initObsPDFNumer(){
	_obsPDFNumer.clearMatrix();
	SMLMS::Matrix tempMat(_stateNumber, _symbolNumber);
	_obsPDFNumer = tempMat;
	resetObsPDFNumer();
}

void HMMUnique::resetObsPDFNumer(){
	unsigned i,j;
	for (i=0; i<_stateNumber; i++){
		for (j=0; j<_symbolNumber; j++){
			_obsPDFNumer.at(i,j,0.0);
		}
	}
}

void HMMUnique::initPDFDenominator(){
	_pdfDenominator.clearMatrix();
	SMLMS::Matrix tempMat(1, _stateNumber);
	_pdfDenominator = tempMat;
	resetPDFDenominator();
}

void HMMUnique::resetPDFDenominator(){
	unsigned i;
	for (i=0; i<_stateNumber; i++) _pdfDenominator.at(0,i,0.0);
}

void HMMUnique::initAnalysis(){
	initAlpha();
	initNormLH();
	initBeta();
	initXi();
	initGamma();
	initTransPDFNumer();
	initObsPDFNumer();
	initPDFDenominator();
}

void HMMUnique::resetAnalysis(){
	resetAlpha();
	resetNormLH();
	resetBeta();
	resetXi();
	resetGamma();
	resetTransPDFNumer();
	resetObsPDFNumer();
	resetPDFDenominator();
}

/* print functions */
void HMMUnique::printObsNumber(){
	std::cout<<std::endl<<"Number of observations: "<<_obsNumber<<std::endl;
}

void HMMUnique::printNormLH(){
	unsigned t;
	std::cout<<std::endl<<"Normalization factor: "<<std::endl;
	for (t=0; t<_obsNumber; t++) std::cout<<_normLH.at(t)<<std::endl;
	std::cout<<std::endl;
}

void HMMUnique::printAlpha(){
	unsigned i,t;
	std::cout<<std::endl<<"Alpha: "<<std::endl;
	for (t=0; t<_obsNumber; t++){
		for(i=0; i<_stateNumber;i++) std::cout<<_alpha.at(i,t)<<"\t";
		std::cout<<std::endl;
	}
	std::cout<<std::endl;
}

void HMMUnique::printBeta(){
	unsigned i,t;
	std::cout<<std::endl<<"Beta: "<<std::endl;
	for (t=0; t<_obsNumber; t++){
		for(i=0; i<_stateNumber;i++) std::cout<<_beta.at(i,t)<<"\t";
		std::cout<<std::endl;
	}
	std::cout<<std::endl;
}

void HMMUnique::printXi(){
	unsigned i,j,t;
	std::cout<<std::endl<<"Xi: "<<std::endl;
	for (t=0; t<_obsNumber; t++){
		for (i=0; i<_stateNumber; i++){
			for (j=0; j<_stateNumber; j++) std::cout<<_xi.at(t).at(i,j)<<"\t";
		}
		std::cout<<std::endl;
	}
	std::cout<<std::endl;
}

void HMMUnique::printGamma(){
	unsigned i,t;
	std::cout<<std::endl<<"Gamma: "<<std::endl;
	for (t=0; t<_obsNumber; t++){
		for(i=0; i<_stateNumber; i++) std::cout<<_gamma.at(i,t)<<"\t";
		std::cout<<std::endl;
	}
	std::cout<<std::endl;
}

void HMMUnique::printTransPDFNumer(){
	unsigned i,j;
	std::cout<<std::endl<<"Transition matrix numer:"<<std::endl;
	for (i=0; i<_stateNumber; i++){
		for (j=0; j<_stateNumber; j++) std::cout<<_transPDFNumer.at(i,j)<<"\t";
		std::cout<<std::endl;
	}
	std::cout<<std::endl;
}

void HMMUnique::printObsPDFNumer(){
	unsigned i,j;
	std::cout<<std::endl<<"Observation matrix numer:"<<std::endl;
	for (i=0; i<_stateNumber; i++){
		for (j=0; j<_symbolNumber; j++) std::cout<<_obsPDFNumer.at(i,j)<<"\t";
		std::cout<<std::endl;
	}
	std::cout<<std::endl;
}

void HMMUnique::printPDFDenominator(){
	unsigned i;
	std::cout<<std::endl<<"Transition/Observation matrix denominator:"<<std::endl;
	for (i=0; i<_stateNumber; i++) std::cout<<_pdfDenominator.at(0,i)<<"\t";
	std::cout<<std::endl;
}
/* help functions*/
int HMMUnique::findMatch(const SMLMS::Matrix& cdf, int state, double event){
	int eventPdfPosition=0;
	unsigned i;
	// binary search
	std::vector<double> cdfState(cdf.numberOfColumns());
	for (i=0; i<cdf.numberOfColumns(); i++) cdfState.at(i) = cdf.at(state,i);
	auto lb = std::lower_bound(cdfState.begin(), cdfState.end(), event);
	eventPdfPosition = lb - cdfState.begin();
	return eventPdfPosition;
}

int HMMUnique::obsPDFMatch(double event){
	int eventPdfPosition=0;
	// binary search
	auto lb = std::lower_bound(_obsAlphabet.begin(), _obsAlphabet.end(), event);
	eventPdfPosition = lb - _obsAlphabet.begin();
	return eventPdfPosition;
}

bool HMMUnique::obsSymbolMatch(int i, double event){
	if (i==0){
		if (event<=_obsAlphabet.at(0))	return true;
	}
	else{
		if ((event>_obsAlphabet.at(i-1)) && (event<=_obsAlphabet.at(i))) return true;
	}
	return false;
}

void HMMUnique::calcXi(const std::vector<double> &obsSeq){
	unsigned t,i,j;
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

void HMMUnique::calcGamma(){
	unsigned i,t;
	double denominator;
	for (t=1; t<_obsNumber; t++){
		denominator = 0.0;
		for(i=0; i<_stateNumber; i++) denominator += _gamma(i,t)= _alpha(i,t)*_beta(i,t);
		for(i=0; i<_stateNumber; i++) _gamma(i,t) /= denominator;
	}
}

void HMMUnique::reestimateEquiPDF(){
	unsigned i;
	for(i=0; i<_stateNumber; i++) _equiPDF(0,i) = _gamma(i,1);
}

void HMMUnique::estimateTransPDFNumer(){
	unsigned i,j,t;
	resetTransPDFNumer();
	for (i=0; i<_stateNumber; i++){
		for (j=0;j<_stateNumber; j++){
			for (t=1; t<_obsNumber; t++)_transPDFNumer(i,j) += (_xi.at(t)).at(i,j);
		}
	}
}

void HMMUnique::reestimateTransPDF(){
	unsigned i,j;
	estimateTransPDFNumer();
	estimatePDFDenominator();
	for (i=0; i<_stateNumber; i++){
		for (j=0;j<_stateNumber; j++) _transPDF(i,j) = _transPDFNumer.at(i,j)/_pdfDenominator.at(0,i);
	}
}

void HMMUnique::estimateObsPDFNumer(const std::vector<double> &obsSeq){
	unsigned i,j,t;
	resetObsPDFNumer();
	for (i=0; i<_stateNumber; i++){
		for (j=0;j<_symbolNumber; j++){
			for (t=1; t<_obsNumber; t++){
				if (obsSymbolMatch(j, obsSeq.at(t))) _obsPDFNumer(i,j) += _gamma.at(i,t);
			}
		}
	}
}

void HMMUnique::estimatePDFDenominator(){
	unsigned i,t;
	resetPDFDenominator();
	for (i=0; i<_stateNumber; i++){
		for (t=0; t<_obsNumber; t++) _pdfDenominator(0,i) += _gamma.at(i,t);
	}
}

void HMMUnique::reestimateObsPDF(const std::vector<double> &obsSeq){
	unsigned i,j;
	estimateObsPDFNumer(obsSeq);
	estimatePDFDenominator();
	for (i=0; i<_stateNumber; i++){
		for (j=0;j<_symbolNumber; j++) _obsPDF(i,j) = _obsPDFNumer.at(i,j)/_pdfDenominator.at(0,i);
	}
}

/* core functions*/
/* core fundtions simulation */
void HMMUnique::simulate(std::vector<double> &obsSeq, std::vector<int> &stateSeq){
	if (obsSeq.size() != stateSeq.size()){
		std::stringstream errorMessage;
		errorMessage<<"Dimension error: Observation vector and states vector must have the same dimensions."<<std::endl;
		SMLMS::SmlmsError error(errorMessage.str());
		throw error;
	}	
	HMMUnique::simulateStates(stateSeq);
	HMMUnique::simulateObservations(obsSeq, stateSeq);
}

void HMMUnique::simulateStates(std::vector<int> &stateSeq){
	double event = 0.0;
	unsigned t;
	int preState;
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

void HMMUnique::simulateObservations(std::vector<double> &obsSeq, const std::vector<int> &stateSeq){
	double event = 0.0;
	int eventPDFPos;
	unsigned t;
	checkObsLength(obsSeq.size());
	checkObsLength(stateSeq.size());
	for(t=0; t<_obsNumber; t++){
		event = _randGen.generateRandomDouble(0.0, 1.0);
		eventPDFPos = findMatch(_obsCDF, stateSeq.at(t), event);
		obsSeq.at(t)=_obsAlphabet.at(eventPDFPos);
	}
}

/* core functions: HMM Estimate likelihood */

void HMMUnique::forwardAlgorithm(const std::vector<double> &obsSeq){
	unsigned i,j,t;
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

void HMMUnique::backwardAlgorithm(const std::vector<double> &obsSeq){
	unsigned i,j,t;
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



void HMMUnique::forwardBackward(const std::vector<double> &obsSeq){
	initAnalysis();
	forwardAlgorithm(obsSeq);
	backwardAlgorithm(obsSeq);
	_fbDone = true;
}

void HMMUnique::estimateLikelihood(){
	unsigned t;
	double llc =0.0;
	for (t=0; t<_obsNumber; t++){
		llc += log(_normLH.at(t));
	}
	_logLikelihood = -1*llc;
	calcModelSelection(_obsNumber);
}

void HMMUnique::baumWelch(const std::vector<double> &obsSeq){
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

void HMMUnique::train(const std::vector<double> &obsSeq){
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
		if(lldiff <= _stopCrit){
			std::cout<<std::endl<<"Training ended by reaching stop criterium."<<std::endl;
			crit = false;
		}
		else if (n>=_maxIt){
			std::cout<<std::endl<<"Training ended by reaching maximum number iterations."<<std::endl;
			crit = false;
		}
	}
}

void HMMUnique::viterbi(std::vector<int> &stateSeq, const std::vector<double> &obsSeq){
	unsigned i,j,t;
	//proof seqLength
	// Initialization
	double tempNorm = 0.0;
	double tempProb =0.0;
	double maxVal=0.0;
	int maxInd=0;
	SMLMS::Matrix delta(_stateNumber, _obsNumber);
	SMLMS::Matrix psi(_stateNumber, _obsNumber);
	//initAnalysis();
	//forwardAlgorithm(obsSeq);		
	//Induction
	for (i=0; i<_stateNumber; i++){
		//delta(i,0)=_alpha(i,0);
		delta(i,0)=_equiPDF(0,i) * _obsPDF(i,obsPDFMatch(obsSeq.at(0)));
		tempNorm += delta(i,0);
		psi(i,0)=0.0;
	}
	for (i=0; i<_stateNumber; i++){
		delta(i,0) /=tempNorm;
		tempProb = delta(i,0);
		if (maxVal<tempProb)psi(i,0)=i;
	}
	//Recursion
	for (t=1; t<_obsNumber; t++){
		tempNorm = 0.0;
		for (j=0; j<_stateNumber; j++){
			tempProb = 0.0;
			//maxVal = _alpha.at(0,t);
			//maxVal= delta.at(0,t-1) * _transPDF.at(0,j) * _obsPDF(j,obsPDFMatch(obsSeq.at(t)));
			maxVal = 0.0;
			maxInd = 0;
			for (i=0; i<_stateNumber; i++){
				//tempProb = _alpha.at(j,t);
				tempProb=delta.at(i,t-1) * _transPDF.at(i,j) *_obsPDF(j,obsPDFMatch(obsSeq.at(t)));
				//tempProb = delta.at(j,t-1) * _transPDF.at(i,j) * _obsPDF(i,obsPDFMatch(obsSeq.at(t)));
				if(maxVal < tempProb){
					maxVal = tempProb;
					maxInd = i;
				}
			}
			delta(j,t) = maxVal;
			tempNorm += maxVal;
			psi(j,t) = maxInd;
		}
		for (j=0; j<_stateNumber;j++) delta(j,t) /= tempNorm;
	}
	//Backtracking
	//maxVal = delta(0,_obsNumber-1);
	maxVal = 0.0;
	maxInd = 0;
	for (i=0; i<_stateNumber; i++){
		if (maxVal < delta(i, _obsNumber-1)){
			maxVal = delta(i, _obsNumber-1);
			maxInd = i;
		}
	}
	stateSeq.at(_obsNumber-1)=maxInd;
	for (t=_obsNumber-1; t>1; t--) stateSeq.at(t-1)=psi.at(stateSeq.at(t),t);
}
		
}/* SMLMS */

