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
#include <cmath>
#include <map>
#include <sstream>
#include <fstream>
#include <unistd.h> /* for getopt */
#include <sys/stat.h>
#include <boost/program_options.hpp>
#include "header/ermineParser.hpp"
#include "header/ermineExceptions.hpp"
#include "header/smlmsFolder.hpp"

namespace po=boost::program_options;
namespace SMLMS{

// constructors
// default constructor
ErmineParser::ErmineParser(){
	std::cout<<"Ermine Parser created"<<std::endl;
	// Initialize default parser arguments
	std::map<std::string, int> alphabet;
	//batch, mol2judi, initialize, simulate, likelihood, train, viterbi
	alphabet.insert(std::make_pair("batch",0));
	alphabet.insert(std::make_pair("mol2judi",0));
	alphabet.insert(std::make_pair("initPhysMod",0));
	alphabet.insert(std::make_pair("fitPhysMod",0));
	alphabet.insert(std::make_pair("initHMM",0));
	alphabet.insert(std::make_pair("simulate",0));
	alphabet.insert(std::make_pair("evaluate",0));
	alphabet.insert(std::make_pair("train",0));
	alphabet.insert(std::make_pair("bestPath",0));
	alphabet.insert(std::make_pair("dwellTime",0));
	setAlgorithmAlphabet(alphabet);
	setAlgorithmArgument("train");
	setStopCritArgument(0.01);
	setMaxItArgument(100);
	setFileNameArgument("");
	setJumpIntervalArgument(10.0);
	setMinDistArgument(10.0);
	setMaxDistArgument(500.0);
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

void ErmineParser::setMaxItArgument(int itArg){
	_maxItArgument = itArg;
}

int ErmineParser::maxItArgument(){
	return _maxItArgument;
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

void ErmineParser::setFolderNameArgument(std::string parserArg){
	_folderNameArgument = parserArg;
}

std::string ErmineParser::folderNameArgument(){
	return _folderNameArgument;
}

void ErmineParser::setJumpIntervalArgument(double parserArg){
	_jumpIntervalArgument = parserArg;
}
double ErmineParser::jumpIntervalArgument(){
	return _jumpIntervalArgument;
}
void ErmineParser::setMinDistArgument(double parserArg){
	_minDistArgument = parserArg;
}
double ErmineParser::minDistArgument(){
	return _minDistArgument;
}

void ErmineParser::setMaxDistArgument(double parserArg){
	_maxDistArgument = parserArg;
}

double ErmineParser::maxDistArgument(){
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
	<<"batch:\t\tmerge several .trc data sets\t\t\t(-a, -f)"<<std::endl
	<<"mol2judi:\tcalculate judi from .trc file.\t\t\t(-a, -f)"<<std::endl
	<<"initPhysMod:\treturns an initial .mod file.\t\t\t(-a, -f)"<<std::endl
	<<"fitPhysMod:\tfits a physical model to a given judi.\t\t(to be announced)"<<std::endl
	<<"initHMM:\treturns an initial guess for a hmm.\t\t(-a -f -j --midDist --maxDist)"<<std::endl
	<<"simulate:\tcalculates a mchmm simulation.\t\t\t(to be announced)"<<std::endl
	<<"evaluate:\tevaluates how well a given model fits a distinct data set\t(to be announced)"<<std::endl
	<<"train:\t\ttrains a hmm on a given data set by Baum-Welch.\t(to be announced)"<<std::endl
	<<"bestPath:\t\testimates the most likely path of hidden states by Viterbi.\t(to be announced)"<<std::endl
	<<"dwellTime:\t\tretimates the model transition rates form an optimized path."<<std::endl;

	std::cout<<line<<message.str()<<std::endl;
}

void ErmineParser::proofFilename(po::variables_map &vm){
	if(vm.count("file")<1){
		std::stringstream errorMessage;
		errorMessage<<"No filename was given!"<<std::endl;;
		SMLMS::ErmineParserError ermineParserError(errorMessage.str());
		throw ermineParserError;
	}
	struct stat buffer;
	std::string name=vm["file"].as<std::string>();
	if ((stat (name.c_str(), &buffer))){
		std::stringstream errorMessage;
		errorMessage<<vm["file"].as<std::string>()<<" does not exist!"<<std::endl;
		SMLMS::ErmineParserError ermineParserError(errorMessage.str());
		throw ermineParserError;
	}
}

void ErmineParser::proofAlgorithm(po::variables_map &vm){
	if(vm.count("algorithm")<1){
		std::stringstream errorMessage;
		errorMessage<<"No algorithm was given!"<<std::endl;
		SMLMS::ErmineParserError ermineParserError(errorMessage.str());
		throw ermineParserError;
	}
}
	
void ErmineParser::proofAlgorithmArgument(){
	if(_algorithmAlphabet.find(_algorithmArgument) == _algorithmAlphabet.end()){
		std::stringstream errorMessage;
		errorMessage<<_algorithmArgument<<" is not a valid algorithm!"<<std::endl;
		SMLMS::ErmineParserError ermineParserError(errorMessage.str());
		throw ermineParserError;	
	}	
}

void ErmineParser::proofStopCrit(po::variables_map &vm){
	if ((vm.count("stopCrit")>0) && (vm["stopCrit"].as<double>()<0.0)){
		std::stringstream errorMessage;
		errorMessage<<"stopCrit needs to be a positive floating point number!"<<std::endl;
		SMLMS::ErmineParserError ermineParserError(errorMessage.str());
		throw ermineParserError;
	}
}

void ErmineParser::proofMaxIt(po::variables_map &vm){
	if ((vm.count("maxIt")>0) && (vm["maxIt"].as<int>()<0)){
		std::stringstream errorMessage;
		errorMessage<<"maxIt needs to be a positive integer!"<<std::endl;
		SMLMS::ErmineParserError ermineParserError(errorMessage.str());
		throw ermineParserError;
	}
}

void ErmineParser::proofJumpInterval(po::variables_map &vm){
	if ((vm.count("jumpInterval")>0) && (vm["jumpInterval"].as<double>()<0)){
		std::stringstream errorMessage;
		errorMessage<<"jumpInterval needs to be a positive integer!"<<std::endl;
		SMLMS::ErmineParserError ermineParserError(errorMessage.str());
		throw ermineParserError;
		}
}

void ErmineParser::proofMinDist(po::variables_map &vm){
	if ((vm.count("minDist")>0) && (vm["minDist"].as<double>()<0.0)){
		std::stringstream errorMessage;
		errorMessage<<"minDist needs to be a positive integer!"<<std::endl;
		SMLMS::ErmineParserError ermineParserError(errorMessage.str());
		throw ermineParserError;
		}
}

void ErmineParser::proofMaxDist(po::variables_map &vm){
	if ((vm.count("maxDist")>0) && (vm["maxDist"].as<double>()<0)){
		std::stringstream errorMessage;
		errorMessage<<"maxDist needs to be a positive integer!"<<std::endl;
		SMLMS::ErmineParserError ermineParserError(errorMessage.str());
		throw ermineParserError;
		}
}

void ErmineParser::proofTime(po::variables_map &vm){
	if ((vm.count("time")>0) && (vm["time"].as<double>()<0)){
		std::stringstream errorMessage;
		errorMessage<<"time needs to be a positive floating point number!"<<std::endl;
		SMLMS::ErmineParserError ermineParserError(errorMessage.str());
		throw ermineParserError;
		}
}

void ErmineParser::proofDuration(po::variables_map &vm){
	if ((vm.count("duration")>0) && (vm["duration"].as<double>()<0)){
		std::stringstream errorMessage;
		errorMessage<<"duration needs to be a positive floating point number!"<<std::endl;
		SMLMS::ErmineParserError ermineParserError(errorMessage.str());
		throw ermineParserError;
		}
}

void ErmineParser::proofTraceLengthRest(double rest){
	if ((rest != 0.0)){
		std::stringstream errorMessage;
		errorMessage<<"Given time and Duration do not match!"<<std::endl;
		SMLMS::ErmineParserError ermineParserError(errorMessage.str());
		throw ermineParserError;
		}

}

void ErmineParser::proofParticles(po::variables_map &vm){
	if ((vm.count("particles")>0) && (vm["particles"].as<int>()<0)){
		std::stringstream errorMessage;
		errorMessage<<"particles needs to be a positive integer!"<<std::endl;
		SMLMS::ErmineParserError ermineParserError(errorMessage.str());
		throw ermineParserError;
		}
}

void ErmineParser::proofJumpIntervalValidity(){
	double interval = (_maxDistArgument - _minDistArgument) / _jumpIntervalArgument;
	double length;
	double rest = std::modf(interval, &length); 
	if ((rest != 0.0)){
		std::stringstream errorMessage;
		errorMessage<<"Given observation interval (min max) and jump distance interval do not match!"<<std::endl;
		SMLMS::ErmineParserError ermineParserError(errorMessage.str());
		throw ermineParserError;
	}
}

void ErmineParser::calcTraceLength(){
	double length = _durationArgument/_timeIntervalArgument;
	double tempLength;
	double rest = std::modf(length, &tempLength);
	proofTraceLengthRest(rest);
	_traceLengthArgument = int(tempLength);
}


void ErmineParser::parseArguments(po::variables_map &vm){
	// parse filename
	proofFilename(vm);
	setFileNameArgument(vm["file"].as<std::string>());
	// parse algorithm
	proofAlgorithm(vm);
	setAlgorithmArgument(vm["algorithm"].as<std::string>());
	proofAlgorithmArgument();
	// parse stopCrit
	proofStopCrit(vm);
	setStopCritArgument(vm["stopCrit"].as<double>());
	// parse maxIt
	proofMaxIt(vm);
	setMaxItArgument(vm["maxIt"].as<int>());
	// parse minDist
	proofMinDist(vm);
	setMinDistArgument(vm["minDist"].as<double>());
	// parse maxDist
	proofMaxDist(vm);
	setMaxDistArgument(vm["maxDist"].as<double>());
	// parse jumpInterval
	proofJumpInterval(vm);
	setJumpIntervalArgument(vm["jumpInterval"].as<double>());
	proofJumpIntervalValidity();
	// parse time
	proofTime(vm);
	setTimeIntervalArgument(vm["time"].as<double>());
	// parse duration
	proofDuration(vm);
	setDurationArgument(vm["duration"].as<double>());
	// parse particles
	proofParticles(vm);
	setParticleArgument(vm["particles"].as<int>());
	// calc trace length
	calcTraceLength();
	// create result Folder
	extractFolderName();
	_folderArgument.setFolderName(_folderNameArgument);
}
	
void ErmineParser::extractFolderName(){
/* Extracts a foldername for the calculated results from the given input filename:
 * it strips the filename suffix, that must be ".txt"
 */
	std::string suffix(".txt");
	std::size_t found = _fileNameArgument.find(suffix);	
	if (found!=std::string::npos){
		_folderNameArgument=_fileNameArgument;
		_folderNameArgument.erase(_fileNameArgument.find(suffix),4);
		_folderNameArgument.append("_");
		_folderNameArgument.append(_algorithmArgument);
	}
	else{
		std::stringstream errorMessage;
		errorMessage<<"Input file needs the suffix '.txt'!"<<std::endl;
		SMLMS::SMLMSFolderError smlmsFolderError(errorMessage.str());
		throw smlmsFolderError;
	}
}


void ErmineParser::makeFolder(){
	/* Proofs weather folder exists.
 	 *if not: Folder will be created.
	*/ 	
	_folderArgument.createFolder();
}


void ErmineParser::writeErmineParser(){
	std::string outFileName;
	outFileName = _folderNameArgument;
	outFileName.append("/parserArgs.txt");
	std::ofstream parserFile;
  	parserFile.open (outFileName.data());
  	parserFile << "# Ermine parser"<<std::endl
	<<"# filename:\n"<<_fileNameArgument<<std::endl
	<<"# algorithm:\n "<<_algorithmArgument<<std::endl
	<<"# folder:\n"<<_folderNameArgument<<std::endl
	<<"# stop criterium:\n"<<_stopCritArgument<<std::endl
	<<"# maximum iterations:\n"<<_maxItArgument<<std::endl
	<<"# jump interval [nm]:\n"<<_jumpIntervalArgument<<std::endl
	<<"# min distance [nm]:\n"<<_minDistArgument<<std::endl
	<<"# max distance [nm]:\n"<<_maxDistArgument<<std::endl
	<<"# time interval [s]:\n"<<_timeIntervalArgument<<std::endl
	<<"# following parameters are only of interest for algorithm 'simulate'!"<<std::endl
	<<"# duration [s]:\n"<<_durationArgument<<std::endl
	<<"# trace length:\n"<<_traceLengthArgument<<std::endl
	<<"# particles:\n"<<_particleArgument<<std::endl;	
  	parserFile.close();
	std::cout<<"Writing parser to: "<<outFileName<<std::endl;
}
	
void ErmineParser::printArguments(){
	std::cout<<std::endl<<"GivenArguments:"<<std::endl;
	std::cout<<"filename: "<<_fileNameArgument<<std::endl;
	std::cout<<"algorithm: "<<_algorithmArgument<<std::endl;
	std::cout<<"folder: "<<_folderNameArgument<<std::endl;
	std::cout<<"stopCrit: "<<_stopCritArgument<<std::endl;
	std::cout<<"maxIt: "<<_maxItArgument<<std::endl;
	std::cout<<"jump interval: "<<_jumpIntervalArgument<<std::endl;
	std::cout<<"min Dist: "<<_minDistArgument<<std::endl;
	std::cout<<"max Dist: "<<_maxDistArgument<<std::endl;
	std::cout<<"time: "<<_timeIntervalArgument<<std::endl;
	std::cout<<"duration: "<<_durationArgument<<std::endl;
	std::cout<<"trace length: "<<_traceLengthArgument<<std::endl;
	std::cout<<"particles: "<<_particleArgument<<std::endl;
	std::cout<<std::endl;
}
}//SMLMS
