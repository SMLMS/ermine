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
#include "header/smlmsFolder.hpp"

namespace po=boost::program_options;
namespace SMLMS{

// constructors
// default constructor
ErmineParser::ErmineParser(){
	// Initialize default parser arguments
	std::map<std::string, int> alphabet;
	//batch, mol2judi, initialize, simulate, likelihood, train, viterbi
	alphabet.insert(std::make_pair("batch",0));
	alphabet.insert(std::make_pair("mol2judi",0));
	alphabet.insert(std::make_pair("initialize",0));
	alphabet.insert(std::make_pair("simulate",0));
	alphabet.insert(std::make_pair("likelihood",0));
	alphabet.insert(std::make_pair("train",0));
	alphabet.insert(std::make_pair("path",0));
	setAlgorithmAlphabet(alphabet);
	setAlgorithmArgument("train");
	setStopCritArgument(0.01);	
	setFileNameArgument("");
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

void ErmineParser::setFolderArgument(SMLMS::SMLMSFolder parserArg){
	_folderArgument = parserArg;
}

SMLMS::SMLMSFolder ErmineParser::folderArgument(){
	return _folderArgument;
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

void ErmineParser::printAlgorithmHelp(){
	std::string line;
	line       = "\n-----------------------------------------------------------------------------------------\n";
	
	std::stringstream message;
	message<<"possible algorithms in ermine are:"<<std::endl
	<<"algorithm\t\tdescription\t\t\t\tessential parameters"<<std::endl
	<<"batch:\t\tmerge several .trc data sets\t\t\t(-f, -a)"<<std::endl
	<<"mol2judi:\tcalculate judi from .trc file.\t\t\t(to be annonced)"<<std::endl
	<<"initialize:\treturns an initial guess for a hmm.\t\t(to be announced)"<<std::endl
	<<"simulate:\tcalculates a mchmm simulation.\t\t\t(to be announced)"<<std::endl
	<<"likelihood:\tcalculates the likelihood of a given model\t(to be announced)"<<std::endl
	<<"train:\t\ttrains a hmm on a given data set by Baum-Welch.\t(to be announced)"<<std::endl
	<<"path:\t\testimates the most likely path by Viterbi.\t(to be announced)"<<std::endl;

	std::cout<<line<<message.str()<<std::endl;
}

void ErmineParser::parseArguments(po::variables_map &vm){
	std::string parameter;
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
	
	// parse stopCrit
	if (vm.count("stopCrit")>0){
		if (vm["stopCrit"].as<double>()<0.0) {
			parameter = "stopCrit";
			SMLMS::WrongDataType wrongDataTypeError(parameter);
			throw wrongDataTypeError;
		}
		setStopCritArgument(vm["stopCrit"].as<double>());
	}
	// parse jumpInterval
	if (vm.count("jumpInterval")>0){
		
		if (vm["jumpInterval"].as<int>()<0){
			parameter = "jumpInterval";
			SMLMS::WrongDataType wrongDataTypeError(parameter);
			throw wrongDataTypeError;
		}
		setJumpIntervalArgument(vm["jumpInterval"].as<int>());
	}
	// parse minDist
	if (vm.count("minDist")>0){
		if (vm["minDist"].as<int>()<0.0){
			parameter = "minDist";
			SMLMS::WrongDataType wrongDataTypeError(parameter);
			throw wrongDataTypeError;
		}
		setMinDistArgument(vm["minDist"].as<int>());
	}
	// parse maxdist
	if (vm.count("maxDist")>0){
		if(vm["maxDist"].as<int>()<0){
			parameter = "maxDist";
			SMLMS::WrongDataType wrongDataTypeError(parameter);
			throw wrongDataTypeError;
		}
		setMaxDistArgument(vm["maxDist"].as<int>());
	}
	// parse time
	if (vm.count("time")>0){
		if (vm["time"].as<double>()<0){
			parameter = "time";
			SMLMS::WrongDataType wrongDataTypeError(parameter);
			throw wrongDataTypeError;
		}
		setTimeIntervalArgument(vm["time"].as<double>());
	}
	// parse duration
	if (vm.count("duration")>0){
		if (vm["duration"].as<double>()<0){
			parameter = "duration";
			SMLMS::WrongDataType wrongDataTypeError(parameter);
			throw wrongDataTypeError;
		}
		setDurationArgument(vm["duration"].as<double>());
	}
	// parse particles
	if (vm.count("particles")>0){
		if (vm["particles"].as<int>()<0){
			parameter = "particles";
			SMLMS::WrongDataType wrongDataTypeError(parameter);
			throw wrongDataTypeError;
		}
		setParticleArgument(vm["particles"].as<int>());
	}
	// calculate _traceLength
	calcTraceLength();
	// create result Folder
	_folderArgument.extractFolderName(_fileNameArgument, _algorithmArgument);
	
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

void ErmineParser::writeErmineParser(){
	std::cout<<"writing parser"<<std::endl;
}

void ErmineParser::printArguments(){
	std::cout<<"filename: "<<_fileNameArgument<<std::endl;
	std::cout<<"algorithm: "<<_algorithmArgument<<std::endl;
	std::cout<<"folder: "<<_folderArgument.folderName()<<std::endl;
	std::cout<<"stopCrit: "<<_stopCritArgument<<std::endl;
	std::cout<<"jump interval: "<<_jumpIntervalArgument<<std::endl;
	std::cout<<"min Dist: "<<_minDistArgument<<std::endl;
	std::cout<<"max Dist: "<<_maxDistArgument<<std::endl;
	std::cout<<"time: "<<_timeIntervalArgument<<std::endl;
	std::cout<<"duration: "<<_durationArgument<<std::endl;
	std::cout<<"trace length: "<<_traceLengthArgument<<std::endl;
	std::cout<<"particles: "<<_particleArgument<<std::endl;
}
}//SMLMS
