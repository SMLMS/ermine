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
#include <sstream>
#include <unistd.h> /* for getopt */
#include <boost/program_options.hpp>
#include "header/ermineParser.hpp"

namespace po=boost::program_options;
namespace SMLMS{

// constructors
// default constructor
ErmineParser::ErmineParser(){
	// Initialize parser indicators
	setAlgorithmIndicator(0);
	setStopCritIndicator(0);
	setFilenameIndicator(0);
	setJumpIntervalIndicator(0);
	setMinDistIndicator(0);
	setMaxDistIndicator(0);
	setLengthIndicator(0);
	setParticleIndicator(0);
	// Initialize default parser arguments
	setAlgorithmArgument("train");
	setStopCritArgument(0.01);	
	setFileNameArgument("");
	setFolderNameArgument("");
	setJumpIntervalArgument(10);
	setMinDistArgument(10);
	setMaxDistArgument(500);
	setLengthArgument(20);
	setParticleArgument(1000);
}// ErmineParser()

// Destructor
ErmineParser::~ErmineParser(){
	std::cout<< "ErmineParser removed from heap." <<std::endl;
}// ~ErmineParser

// Assessor functions for ErmienParser indicators
void ErmineParser::setAlgorithmIndicator(int parserInd){
	_algorithmIndicator = parserInd;
}

int ErmineParser::algorithmIndicator(){
	return _algorithmIndicator;
}

void ErmineParser::setStopCritIndicator(int parserInd){
	_stopCritIndicator = parserInd;
}
int ErmineParser::stopCritIndicator(){
	return _stopCritIndicator;
}

void ErmineParser::setFilenameIndicator(int parserInd){
	_filenameIndicator = parserInd;
}
int ErmineParser::filenameIndicator(){
	return _filenameIndicator;
}

void ErmineParser::setJumpIntervalIndicator(int parserInd){
	_jumpIntervalIndicator = parserInd;
}
int ErmineParser::jumpIntervalIndicator(){
	return _jumpIntervalIndicator;
}

void ErmineParser::setMinDistIndicator(int parserInd){
	_minDistIndicator = parserInd;
}
int ErmineParser::minDistIndicator(){
	return _minDistIndicator;
}

void ErmineParser::setMaxDistIndicator(int parserInd){
	_maxDistIndicator = parserInd;
}
int ErmineParser::maxDistIndicator(){
	return _maxDistIndicator;
}

void ErmineParser::setLengthIndicator(int parserInd){
	_lengthIndicator= parserInd;
}
int ErmineParser::lengthIndicator(){
	return _lengthIndicator;
}

void ErmineParser::setParticleIndicator(int parserInd){
	_particleIndicator=parserInd;
}

int ErmineParser::particleIndicator(){
	return _particleIndicator;
}

// Assessor functions for ErmineParser arguments
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

void ErmineParser::setLengthArgument(int parserArg){
	_lengthArgument = parserArg;
}

int ErmineParser::lengthArgument(){
	return _lengthArgument;
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
	helpLine = "\n type\n.\ermine -h\nto view help\n";
	line       = "\n-----------------------------------------------------------------------------------------\n";

	std::cout<<line<< helpHeader<<helpLine<<line <<std::endl;
}

void ErmineParser::parseArguments(){
	std::cout<<"so"<<std::endl;
}
}//SMLMS
