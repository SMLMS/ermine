/* ######################################################################
* File Name: ermineParser.cpp
* Project: SMLMS
* Version: 1611
* Creation Date: 04.11.2016
* Created By Sebastian Malkusch
* <malkusch@chemie.uni-frankfurt.de>
* Goethe University of Frankfurt
* Physical and Theoretical Chemistry
* Single Molecule Biophysics
###################################################################### */

#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <sstream>
#include <unistd.h> /* for getopt */
#include <boost/program_options.hpp>
#include "header/ermineParser.hpp"
#include "header/ermineExceptions.hpp"

namespace po=boost::program_options;
namespace SMLMS{

// constructors
// default constructor
ErmineParser::ErmineParser(){
	// Initialize default parser arguments
	std::map<std::string, int> alphabet;
	//batch, mol2judi, initialize, simulate, estimate, train, viterbi
	alphabet.insert(std::make_pair("batch",0));
	alphabet.insert(std::make_pair("mol2judi",0));
	alphabet.insert(std::make_pair("initialize",0));
	alphabet.insert(std::make_pair("simulate",0));
	alphabet.insert(std::make_pair("estimate",0));
	alphabet.insert(std::make_pair("train",0));
	alphabet.insert(std::make_pair("path",0));
	setAlgorithmAlphabet(alphabet);
	setAlgorithmArgument("train");
	setStopCritArgument(0.01);	
	setFileNameArgument("");
	setFolderNameArgument("");
	setJumpIntervalArgument(10);
	setMinDistArgument(10);
	setMaxDistArgument(500);
	setTimeIntervalArgument(0.0);
	setTraceLengthArgument(20);
	setParticleArgument(1000);
}// ErmineParser()

// Destructor
ErmineParser::~ErmineParser(){
	std::cout<< "ErmineParser removed from heap." <<std::endl;
}// ~ErmineParser

// Assessor functions for ErmineParser arguments
void ErmineParser::setAlgorithmAlphabet(std::map<std::string, int> parserArg){
	_algorithmAlphabet = parserArg;
}

std::map<std::string, int> ErmineParser::algorithmAlphabet(){
	return _algorithmAlphabet;
}

void ErmineParser::setAlgorithmArgument(std::string parserArg){
	_algorithmArgument = parserArg;
}
std::string ErmineParser::algorithmArgument(){
	return _algorithmArgument;
}

void ErmineParser::setStopCritArgument(double parserArg){
	_stopCritArgument = parserArg;
}
double ErmineParser::stopCritArgument(){
	return _stopCritArgument;
}

void ErmineParser::setFileNameArgument(std::string parserArg){
	_fileNameArgument = parserArg;
}
std::string ErmineParser::fileNameArgument(){
	return _fileNameArgument;
}

void ErmineParser::setFolderNameArgument(std::string parserArg){
	_folderNameArgument = parserArg;
}

std::string ErmineParser::folderNameArgument(){
	return _folderNameArgument;
}

void ErmineParser::setJumpIntervalArgument(int parserArg){
	_jumpIntervalArgument = parserArg;
}
int ErmineParser::jumpIntervalArgument(){
	return _jumpIntervalArgument;
}
void ErmineParser::setMinDistArgument(int parserArg){
	_minDistArgument = parserArg;
}
int ErmineParser::minDistArgument(){
	return _minDistArgument;
}

void ErmineParser::setMaxDistArgument(int parserArg){
	_maxDistArgument = parserArg;
}
int ErmineParser::maxDistArgument(){
	return _maxDistArgument;
}

void ErmineParser::setTimeIntervalArgument(double parserArg){
	_timeIntervalArgument = parserArg;
}

double ErmineParser::timeIntervalArgument(){
	return _timeIntervalArgument;
}

void ErmineParser::setDurationArgument(double parserArg){
	_durationArgument = parserArg;
}

double ErmineParser::durationArgument(){
	return _durationArgument;
}

void ErmineParser::setTraceLengthArgument(int parserArg){
	_traceLengthArgument = parserArg;
}

int ErmineParser::traceLengthArgument(){
	return _traceLengthArgument;
}

void ErmineParser::setParticleArgument(int parserArg){
	_particleArgument=parserArg;
}

int ErmineParser::particleArgument(){
	return _particleArgument;
}

void ErmineParser::printHelp(){
	std::string helpHeader, helpLine, line;
	helpHeader = "\nEstimate Reaction-rates by Markov-based Investigation of Nanoscopy Experiments (ermine):\n";
	helpLine = "\n type: ./ermine -h (to view help)\n\n";
	line       = "\n-----------------------------------------------------------------------------------------\n";

	std::cout<<line<< helpHeader<<helpLine<<line <<std::endl;
}

void ErmineParser::parseArguments(po::variables_map &vm){
	// parse filename
	if(vm.count("file")<1){
		SMLMS::NoFileName noFileNameError;
		throw noFileNameError;
	}
	else{
		setFileNameArgument(vm["file"].as<std::string>());
	}

	// parse algorithm
	if(vm.count("algorithm")<1){
		SMLMS::NoAlgorithm noAlgorithmError;
		throw noAlgorithmError;
	}

	setAlgorithmArgument(vm["algorithm"].as<std::string>());

	if(proofAlgorithmArgument()){
		SMLMS::WrongAlgorithm wrongAlgorithmError(_algorithmArgument);
		throw wrongAlgorithmError;
	}
}	

int ErmineParser::proofAlgorithmArgument(){
	if(_algorithmAlphabet.find(_algorithmArgument) != _algorithmAlphabet.end()){
		std::cout<<0<<std::endl;
		return 0;	
	}	
	else{
		std::cout<<1<<std::endl;
		return 1;
	}
}

void ErmineParser::calcTraceLength(){
	_traceLengthArgument = int(_durationArgument/_timeIntervalArgument);
}
}//SMLMS
