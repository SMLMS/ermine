/* ######################################################################
* File Name: smlmsHmmSequence.cpp
* Project: SMLMS
* Version: 17.02
* Creation Date: 23.02.2017
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
#include <boost/progress.hpp>
#include "header/smlmsContainer.hpp"
#include "header/smlmsExceptions.hpp"
#include "header/smlmsHmmBase.hpp"
#include "header/smlmsHmmUnique.hpp"
#include "header/smlmsHmmSequence.hpp"
#include "header/smlmsMatrix.hpp"
#include "header/smlmsJudi.hpp"
#include "header/smlmsPhysModBrownLatDiff.hpp"

namespace SMLMS{
/* constructor */
HMMSequence::HMMSequence():SMLMS::HMMBase(){
	std::cout<<"HMM constructor called."<<std::endl;
	setTraceNumber(0);
	_seqLogLikelihood = 0.0;
	_fbSeqDone = false;
}

HMMSequence::HMMSequence(unsigned states, unsigned symbols):SMLMS::HMMBase(states, symbols){
	std::cout<<"HMM constructor called"<<std::endl;
	//set everytging
	setTraceNumber(0);
	_seqLogLikelihood = 0.0;
	_fbSeqDone = false;
}

HMMSequence::HMMSequence(unsigned states, unsigned symbols, unsigned trace):SMLMS::HMMBase(states, symbols){
	std::cout<<"HMM constructor called"<<std::endl;
	//set everytging
	setTraceNumber(trace);
	_seqLogLikelihood = 0.0;
	_fbSeqDone = false;
}

/* copy constructor */
HMMSequence::HMMSequence(const HMMSequence &obj){
	std::cout<<"HMM copy constructor called."<<std::endl;
	_traceNumber = obj._traceNumber;
	_stateNumber = obj._stateNumber;
	_modelAdjustInd = obj._modelAdjustInd;
	_seqLogLikelihood = obj._seqLogLikelihood;
	_symbolNumber  = obj._symbolNumber;
	_symbolInterval = obj._symbolInterval;
	_equiPDF = obj._equiPDF;
	_transPDF = obj._transPDF;
	_obsPDF = obj._obsPDF;
	_equiCDF = obj._equiCDF;
	_transCDF = obj._transCDF;
	_obsCDF = obj._obsCDF;
	_obsAlphabet = obj._obsAlphabet;
	_equiPDFNumer = obj._equiPDFNumer;
	_transPDFNumer = obj._transPDFNumer;
	_obsPDFNumer = obj._obsPDFNumer;
	_pdfDenominator = obj._pdfDenominator;
	_logLikelihood = obj._logLikelihood;
	_fbSeqDone = obj._fbSeqDone;
}

/* elementary functions */
void HMMSequence::setTraceNumber(unsigned val){
	_traceNumber = val;
}

unsigned HMMSequence::traceNumber(){
	return _traceNumber;
}

void HMMSequence::setModelAdjustInd(bool val){
	_modelAdjustInd = val;
}

bool HMMSequence::modelAdjustInd(){
	return _modelAdjustInd;
}

double HMMSequence::seqLogLikelihood(){
	return _seqLogLikelihood;
}

/* proof functions */
void HMMSequence::checkTraceNumber(){
	if(_traceNumber<1){
		std::stringstream errorMessage;
		errorMessage<<"Cannot setup a hmm sequence without a definite number of traces!"<<std::endl;
		SMLMS::SmlmsError error(errorMessage.str());
		throw error;
	}
}

void HMMSequence::checkSimulationDimension(const SMLMS::JumpDistanceList &judi, unsigned obsVal){
	unsigned seqLength = obsVal * _traceNumber;
	unsigned judiLength = judi.getNumberOfJumps();
	if (judiLength != seqLength){
		std::stringstream errorMessage;
		errorMessage<<"Cannot simulate an hmm sequence: Dimension Mismatch!"<<std::endl;
		SMLMS::SmlmsError error(errorMessage.str());
		throw error;
	}
}

/* init functions */
void HMMSequence::initSeqLogLikelihood(){
	_seqLogLikelihood = 0.0;
}

void HMMSequence::initEquiPDFNumer(){
	checkStateNumber();
	_equiPDFNumer.clearMatrix();
	SMLMS::Matrix tempMat(1, _stateNumber);
	_equiPDFNumer = tempMat;
}

void HMMSequence::initTransPDFNumer(){
	checkStateNumber();
	_transPDFNumer.clearMatrix();
	SMLMS::Matrix tempMat(_stateNumber, _stateNumber);
	_transPDFNumer = tempMat;
}


void HMMSequence::initObsPDFNumer(){
	checkStateNumber();
	checkSymbolNumber();
	_obsPDFNumer.clearMatrix();
	SMLMS::Matrix tempMat(_stateNumber, _symbolNumber);
	_obsPDFNumer = tempMat;
}

void HMMSequence::initPDFDenominator(){
	checkStateNumber();
	_pdfDenominator.clearMatrix();
	SMLMS::Matrix tempMat(1, _stateNumber);
	_pdfDenominator = tempMat;
}

void HMMSequence::initTrainingSequences(){
	initSeqLogLikelihood();
	initEquiPDFNumer();
	initTransPDFNumer();
	initObsPDFNumer();
	initPDFDenominator();
	_fbSeqDone = false;
}

/* clear functions */
void HMMSequence::resetEquiPDFNumer(){
	for (int i=0; i<_stateNumber; i++)_equiPDFNumer.at(0,i,0.0);
}

void HMMSequence::resetTransPDFNumer(){
	for(int i=0; i<_stateNumber; i++){
		for(int j=0; j<_stateNumber; j++){
			_transPDFNumer.at(i,j,0.0);
		}
	}
}

void HMMSequence::resetObsPDFNumer(){
	for(int i=0; i<_stateNumber; i++){
		for(int j=0; j<_symbolNumber; j++){
			_obsPDFNumer.at(i,j,0.0);
		}
	}
}

void HMMSequence::resetPDFDenominator(){
	for (int i=0; i<_stateNumber; i++)_pdfDenominator.at(0,i,0.0);
}

void HMMSequence::resetTrainingSequences(){
	resetEquiPDFNumer();
	resetTransPDFNumer();
	resetObsPDFNumer();
	resetPDFDenominator();
	initSeqLogLikelihood();
	_fbSeqDone = false;
}

/* print functions */
void HMMSequence::printTraceNumber(){
	std::cout<<"Trace number: "<<_traceNumber<<std::endl<<std::endl;
}

void HMMSequence::printSeqLogLikelihood(){
	std::cout<<"Sequence LogLikelihood: "<<_seqLogLikelihood<<std::endl<<std::endl;
}

void HMMSequence::printEquiPDFNumer(){
	std::cout<<"Equilibrium Matrix Numer:"<<std::endl;
	for (int i=0; i<_stateNumber; i++){
		std::cout<<_equiPDFNumer.at(0,i)<<"\t";
	}
	std::cout<<std::endl<<std::endl;
}

void HMMSequence::printTransPDFNumer(){
	std::cout<<"Transition Matrix Numer:"<<std::endl;
	for(int i=0; i<_stateNumber; i++){
		for(int j=0; j<_stateNumber; j++){
			std::cout<<_transPDFNumer.at(i,j)<<"\t";
		}
		std::cout<<std::endl;
	}
	std::cout<<std::endl;
}

void HMMSequence::printObsPDFNumer(){
	std::cout<<"Observation Matrix Numer:"<<std::endl;
	for(int i=0; i<_stateNumber; i++){
		for(int j=0; j<_symbolNumber; j++){
			std::cout<<_obsPDFNumer.at(i,j)<<"\t";
		}
		std::cout<<std::endl;
	}
	std::cout<<std::endl;
}

void HMMSequence::printPDFDenominator(){
	std::cout<<"Transition/Observation Matrix Denominator:"<<std::endl;
	for(int i=0; i<_stateNumber; i++){
		std::cout<<_pdfDenominator.at(0,i)<<"\t";
	}
	std::cout<<std::endl<<std::endl;
}

/* help functions */
void HMMSequence::initHelpUniqueHMM(SMLMS::HMMUnique& hmm, unsigned obsNumber){
	hmm.setStateNumber(_stateNumber);
	hmm.setSymbolNumber(_symbolNumber);
	hmm.setSymbolInterval(_symbolInterval);
	hmm.setEquiPDF(_equiPDF);
	hmm.setTransPDF(_transPDF);
	hmm.setObsPDF(_obsPDF);
	hmm.setObsAlphabet(_obsAlphabet);
	hmm.setObsNumber(obsNumber);
	hmm.normalizeHMM();
	hmm.checkObsNumber();
	hmm.checkHMM();
}

void HMMSequence::increaseSimulationSequence(SMLMS::JumpDistanceList &judi, unsigned trace, const std::vector<double> &obs, const std::vector<int> &states){
	SMLMS::Jump jump;
	int obsNumber = obs.size(); 
	for (int i=0; i<obsNumber; i++){
		jump.trace=trace;
		jump.jumpDistance=obs.at(i);
		jump.state=states.at(i);
		judi.addJumpToEnd(jump);
	}
}

// model adjustment

void HMMSequence::reestimateEquiPDF(){
	int i;
	for (i=0; i<_stateNumber; i++){
		_equiPDF(0,i) = _equiPDFNumer.at(0,i)/_traceNumber; 
	}
}

void HMMSequence::reestimateTransPDF(){
	int i,j;
	for (i=0; i<_stateNumber; i++){
		for (j=0; j<_stateNumber; j++) _transPDF(i,j)=_transPDFNumer.at(i,j)/_pdfDenominator(0,i);
	}
}

void HMMSequence::reestimateObsPDF(){
	int i,j;
	for (i=0; i<_stateNumber; i++){
		for (j=0; j<_symbolNumber; j++) _obsPDF(i,j)=_obsPDFNumer.at(i,j)/_pdfDenominator(0,i);
	}
}

void HMMSequence::reestimateHMM(){
	reestimateEquiPDF();
	reestimateTransPDF();
	reestimateObsPDF();
}

bool HMMSequence::benchmarkTrainingsResult(double llResult, int itStep){
	if(llResult<_stopCrit){
		std::cout<<std::endl<<"Reached stop criteium after "<<itStep<<" iterations."<<std::endl;
		return false;
	}
	if(itStep>_maxIt){
		std::cout<<std::endl<<"Stop criterium not reached after maximal number of iterations."<<std::endl;
		return false;
	}
	return true;
}

/* core functions */
void HMMSequence::simulateSequence(unsigned obsNumber, SMLMS::JumpDistanceList &judi){
	/* check stuff */
	checkStateNumber();
	/* run simulation */
	SMLMS::JumpDistanceList tempJudi;
	std::vector<int> tempStateList(obsNumber);
	std::vector<double> tempObsList(obsNumber);
	SMLMS::HMMUnique tempHMM(_stateNumber, _symbolNumber);
	initHelpUniqueHMM(tempHMM, obsNumber);
	//simulate
	for (int i=0; i<_traceNumber; i++){
		tempHMM.simulate(tempObsList, tempStateList);
		increaseSimulationSequence(tempJudi, i+1, tempObsList, tempStateList);
	}
	tempJudi.calcTraceNumber();
	judi = tempJudi; 
}

void HMMSequence::estimateSeqLikelihood(const SMLMS::JumpDistanceList &judi){
	// test length
	checkStateNumber();
	judi.checkTraceNumber();
	// init
	initTrainingSequences();
	std::vector<double> obs;
	SMLMS::HMMUnique tempHmm(_stateNumber, _symbolNumber);
	unsigned obsNumber = 0;
	unsigned cummulativeObsNumber = 0;
	for (int i=1; i<judi.traceNumber()+1; i++){
		obs = judi.getTraceJumps(i);
		obsNumber = obs.size();
		cummulativeObsNumber += obsNumber;
		initHelpUniqueHMM(tempHmm, obsNumber);
		tempHmm.forwardBackward(obs);
		tempHmm.estimateLikelihood();
		_seqLogLikelihood += tempHmm.logLikelihood();
	}
	_logLikelihood = _seqLogLikelihood;
	estimateSeqBic(cummulativeObsNumber);
	estimateSeqAic();
}


void HMMSequence::estimateSeqBic(unsigned obsNumber){
	int p = (_stateNumber * (_stateNumber -1)) + (_stateNumber-1) + (_stateNumber * (_symbolNumber - 1));
	_bic = (-2.0 * _logLikelihood) + p * log(obsNumber);
}

void HMMSequence::estimateSeqAic(){
	int p = (_stateNumber * (_stateNumber -1)) + (_stateNumber-1) + (_stateNumber * (_symbolNumber - 1));
	_aic = (-2.0 * _logLikelihood) + p * 2.0;
}


void HMMSequence::baumWelchSequence(const SMLMS::JumpDistanceList &judi){
	// test length
	checkStateNumber();
	judi.checkTraceNumber();
	// init
	resetTrainingSequences();
	std::vector<double> obs;
	SMLMS::HMMUnique tempHmm(_stateNumber, _symbolNumber);
	unsigned obsNumber = 0;
	for (int i=1; i<judi.traceNumber()+1; i++){
		obs = judi.getTraceJumps(i);
		obsNumber = obs.size();
		initHelpUniqueHMM(tempHmm, obsNumber);
		tempHmm.baumWelch(obs);
		tempHmm.estimateLikelihood();
		_equiPDFNumer += tempHmm.equiPDF();
		_transPDFNumer += tempHmm.transPDFNumer();
		_obsPDFNumer += tempHmm.obsPDFNumer();
		_pdfDenominator += tempHmm.pdfDenominator();
		_seqLogLikelihood += tempHmm.logLikelihood();
	}
}

void HMMSequence::trainSequence(const SMLMS::JumpDistanceList &judi){
	int it = 0;
	double llResult = 0.0;
	double lastLl = 0.0;
	// set trace Number
	setTraceNumber(judi.traceNumber());
	// test length
	checkStateNumber();
	checkTraceNumber();
	judi.checkTraceNumber();
	// init
	std::cout<<std::endl<<"The ermine is training:"<<std::endl;
	boost::progress_display show_progress(_maxIt+1);	
	initTrainingSequences();
	estimateSeqLikelihood(judi);
	llResult = 0-_logLikelihood;
	lastLl = _logLikelihood;
	while (benchmarkTrainingsResult(llResult, it)){
		baumWelchSequence(judi);
		reestimateHMM();
		estimateSeqLikelihood(judi);
		llResult = _logLikelihood - lastLl;
		lastLl = _logLikelihood;
		++show_progress;
		it +=1;
	}
	std::cout<<std::endl;
}

void HMMSequence::trainPhysModSequence(const SMLMS::JumpDistanceList &judi, SMLMS::PhysicalModelBLD &physMod){
	int it = 0;
	double llResult = 0.0;
	double lastLl = 0.0;
	// set trace Number
	setTraceNumber(judi.traceNumber());
	// test length
	checkStateNumber();
	checkTraceNumber();
	judi.checkTraceNumber();
	// init
	std::cout<<std::endl<<"The ermine is training:"<<std::endl;
	boost::progress_display show_progress(_maxIt+1);	
	initTrainingSequences();
	estimateSeqLikelihood(judi);
	llResult = 0-_logLikelihood;
	lastLl = _logLikelihood;
	while (benchmarkTrainingsResult(llResult, it)){
		baumWelchSequence(judi);
		reestimateHMM();
		/* refit model */
		physMod.baumWelchFit(_equiPDF, _obsPDF);
		normalizeHMM();
		estimateSeqLikelihood(judi);
		llResult = _logLikelihood - lastLl;
		lastLl = _logLikelihood;
		++show_progress;
		it +=1;
	}
	std::cout<<std::endl;
}

void HMMSequence::estimateStateSequence(SMLMS::JumpDistanceList& judi){
	std::vector<int> states;
	std::vector<double> obs;
	std::cout<<std::endl<<"The ermine is estimating states:"<<std::endl;
	// set trace number
	setTraceNumber(judi.traceNumber());
	boost::progress_display show_progress(_traceNumber+1);	
	SMLMS::HMMUnique tempHmm(_stateNumber, _symbolNumber);
	unsigned obsNumber = 0;
	for (int i=1; i<_traceNumber+1; i++){
		// get jump of trace
		obs.clear();
		obs = judi.getTraceJumps(i);
		obsNumber = obs.size();
		// get states of trace
		states.clear();
		states = judi.getTraceStates(i);
		// calc viterbi
		initHelpUniqueHMM(tempHmm, obsNumber);
		tempHmm.viterbi(states, obs);
		// set states of trace
		judi.setTraceStates(i, states);
		// update progress
		++show_progress;
	}
}

} /* SMLMS */
