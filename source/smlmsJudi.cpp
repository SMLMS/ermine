/* ######################################################################
* File Name: judi.cpp
* Project: SMLMS
* Version: 18.09
* Creation Date: 04.03.2016
* Created By Sebastian Malkusch
* <malkusch@chemie.uni-frankfurt.de>
* Goethe University of Frankfurt
* Physical and Theoretical Chemistry
* Single Molecule Biophysics
###################################################################### */


#include <iostream>
#include <vector>
#include <string>
#include <math.h>
#include <fstream>
#include <sstream>
#include <iomanip>
#include "header/smlmsContainer.hpp"
#include "header/smlmsMolecules.hpp"
#include "header/smlmsJudi.hpp"
#include "header/smlmsExceptions.hpp"

namespace SMLMS{
/* Constructor */
JumpDistanceList::JumpDistanceList(){
	std::cout<<"JumpDistanceListConstructor called"<<std::endl;
};
/* Destructor */
JumpDistanceList::~JumpDistanceList(){
	std::cout<<"JumpDistanceList removed from Heap!"<<std::endl;
};
/* Copy Constructor */
JumpDistanceList::JumpDistanceList(const JumpDistanceList &obj){
	std::cout<<"JumpDistanceList copy constructor called."<<std::endl;
	_jumpDistanceList=obj._jumpDistanceList;
	_traceStatList = obj._traceStatList;
	_traceNumber = obj._traceNumber;
};
/* Operator overload */
SMLMS::Jump& SMLMS::JumpDistanceList::operator()(unsigned index){
	return _jumpDistanceList.at(index);
}

SMLMS::Jump SMLMS::JumpDistanceList::operator()(unsigned index) const{
	return _jumpDistanceList.at(index);
}

/* Elementary Functions */
void JumpDistanceList::setJumpDistanceList(std::vector<Jump> &inputList){
	_jumpDistanceList = inputList;
	calcTraceNumber();
	calcTraceStatList();
}

std::vector<Jump> JumpDistanceList::jumpDistanceList(){
	return _jumpDistanceList;
}

std::vector<Jump> JumpDistanceList::jumpDistanceList()const{
	return _jumpDistanceList;
}

unsigned JumpDistanceList::traceNumber(){
	return _traceNumber;
}

unsigned JumpDistanceList::traceNumber()const{
	return _traceNumber;
}

/* proof functions */
void JumpDistanceList::checkTraceNumber(){
	if(_traceNumber < 1 ){
		std::stringstream errorMessage;
		errorMessage<<"Jump Distance List has a length of 0 traces. Did you forget to calculate the trace length?"<<std::endl;
		SMLMS::SmlmsError error(errorMessage.str());
		throw error;
	}
}

void JumpDistanceList::checkTraceNumber() const{
	if(_traceNumber < 1 ){
		std::stringstream errorMessage;
		errorMessage<<"Jump Distance List has a length of 0 traces. Did you forget to calculate the trace length?"<<std::endl;
		SMLMS::SmlmsError error(errorMessage.str());
		throw error;
	}
}

void JumpDistanceList::checkTraceNumber(unsigned tempTrace){
	if(tempTrace == 0 ){
		std::stringstream errorMessage;
		errorMessage<<"Jump Distance List has no trace of index 0. Please type a valid trace number >0 !"<<std::endl;
		SMLMS::SmlmsError error(errorMessage.str());
		throw error;
	}
	if(tempTrace > _traceNumber ){
		std::stringstream errorMessage;
		errorMessage<<"Trace index exceeds number of valid traces. Please type a valid trace number <"<<_traceNumber<<std::endl;
		SMLMS::SmlmsError error(errorMessage.str());
		throw error;
	}
}

void JumpDistanceList::checkTraceNumber(unsigned tempTrace)const{
	if(tempTrace == 0 ){
		std::stringstream errorMessage;
		errorMessage<<"Jump Distance List has no trace of index 0. Please type a valid trace number >0 !"<<std::endl;
		SMLMS::SmlmsError error(errorMessage.str());
		throw error;
	}
	if(tempTrace > _traceNumber ){
		std::stringstream errorMessage;
		errorMessage<<"Trace index exceeds number of valid traces. Please type a valid trace number <="<<_traceNumber<<std::endl;
		SMLMS::SmlmsError error(errorMessage.str());
		throw error;
	}
}

/* Special Functions */
void JumpDistanceList::calcTraceNumber(){
	unsigned tempTraceNumber, tn, seqLen;
	SMLMS::Jump jump;
	tempTraceNumber = 0;
	seqLen = _jumpDistanceList.size();
	for(unsigned i=0; i<seqLen; i++){
		jump = _jumpDistanceList.at(i);
		tn = jump.trace;
		if (tn>tempTraceNumber)tempTraceNumber=tn;
	}
	_traceNumber = tempTraceNumber;
}

void JumpDistanceList::clearJumpDistanceList(){
	_jumpDistanceList.clear();
	_jumpDistanceList.resize(0);
	clearTraceStatList();
}


SMLMS::Jump JumpDistanceList::calculateJump(Molecule &initMol, Molecule &finalMol){
	//claculate jumpdistance
	double dist = sqrt(
	(initMol.x-finalMol.x)*(initMol.x-finalMol.x)+
	(initMol.y-finalMol.y)*(initMol.y-finalMol.y)
	);
	//create jump
	SMLMS::Jump jump;
	jump.trace=finalMol.trace;
	jump.jumpDistance=dist;
	jump.state=finalMol.state;
	return jump;
}

void JumpDistanceList::addJumpToEnd(Jump &jump){
	//pushback Jump
	_jumpDistanceList.push_back(jump);
}

Jump JumpDistanceList::getJump(int position){
	return _jumpDistanceList.at(position);
}

Jump JumpDistanceList::getJump(int position) const{
	return _jumpDistanceList.at(position);
}

std::vector<double> JumpDistanceList::getTraceJumps(unsigned traceNumber){
	/*
	std::vector<double> traceJumps;
	for(unsigned i=0; i<_jumpDistanceList.size(); i++){
		if (_jumpDistanceList.at(i).trace==traceNumber){
			traceJumps.push_back(_jumpDistanceList.at(i).jumpDistance);
		}
	}
	return traceJumps;
	*/
	checkTraceNumber(traceNumber);
	int i = 0;
	while (_traceStatList.at(i).index != traceNumber){
		++i;
	}
	std::vector<double> traceJumps;
	traceJumps.reserve(_traceStatList.at(i).length);
	for (unsigned j=_traceStatList.at(i).begin; j<_traceStatList.at(i).end +1; j++){
		traceJumps.push_back(_jumpDistanceList.at(j).jumpDistance);
	}
	traceJumps.shrink_to_fit();
	return traceJumps;
}

std::vector<double> JumpDistanceList::getTraceJumps(unsigned traceNumber) const{
	/*
	std::vector<double> traceJumps;
	for(unsigned i=0; i<_jumpDistanceList.size(); i++){
		if (_jumpDistanceList.at(i).trace==traceNumber){
			traceJumps.push_back(_jumpDistanceList.at(i).jumpDistance);
		}
	}
	return traceJumps;
	*/
	checkTraceNumber(traceNumber);
	int i = 0;
	while (_traceStatList.at(i).index != traceNumber){
		++i;
	}
	std::vector<double> traceJumps;
	traceJumps.reserve(_traceStatList.at(i).length);
	for (unsigned j=_traceStatList.at(i).begin; j<_traceStatList.at(i).end +1; j++){
		traceJumps.push_back(_jumpDistanceList.at(j).jumpDistance);
	}
	traceJumps.shrink_to_fit();
	return traceJumps;
}

void JumpDistanceList::setTraceStates(unsigned traceNumber, std::vector<int> &stateTrace){
	unsigned traceLength = stateTrace.size();
	unsigned traceStartPosition=_jumpDistanceList.size();
	unsigned jumpsPerTrace=0;
	for (unsigned i=0; i<_jumpDistanceList.size(); i++){
		if(_jumpDistanceList.at(i).trace==traceNumber){
			jumpsPerTrace+=1;
			if (i<=traceStartPosition)traceStartPosition=i;
		}
	}
	if (traceLength != jumpsPerTrace){
		std::cout<<"Error stateTrace length ("<<traceLength<<") and jumpTrace length("<<jumpsPerTrace<<") do not match!"<<std::endl;
	}	
	else{
		for (unsigned i=0; i<traceLength; i++){
			_jumpDistanceList.at(i+traceStartPosition).state=stateTrace.at(i);
		} 
	}
	
}

std::vector<int> JumpDistanceList::getTraceStates(unsigned traceNumber){
	checkTraceNumber(traceNumber);
	int n = 0;
	while (_traceStatList.at(n).index != traceNumber){
		++n;
	}
	std::vector<int> traceStates(_traceStatList.at(n).length);
	for(unsigned i=0; i<_traceStatList.at(n).length; i++){
		traceStates.at(i)=(_jumpDistanceList.at(i+_traceStatList.at(n).begin).state);
	}
	return traceStates;
}

std::vector<int> JumpDistanceList::getTraceStates(unsigned traceNumber) const{
	/*
	std::vector<int> traceStates;
	for(unsigned i=0; i<_jumpDistanceList.size(); i++){
		if (_jumpDistanceList.at(i).trace==traceNumber){
			traceStates.push_back(_jumpDistanceList.at(i).state);
		}
	}
	return traceStates;
	*/
	checkTraceNumber(traceNumber);
	int n = 0;
	while (_traceStatList.at(n).index != traceNumber){
		++n;
	}
	std::vector<int> traceStates(_traceStatList.at(n).length);
	for(unsigned i=0; i<_traceStatList.at(n).length; i++){
		traceStates.at(i)=(_jumpDistanceList.at(i+_traceStatList.at(n).begin).state);
	}
	return traceStates;
}

int JumpDistanceList::getNumberOfJumps(){
	return _jumpDistanceList.size();
}

int JumpDistanceList::getNumberOfJumps()const{
	return _jumpDistanceList.size();
}

void JumpDistanceList::getAllJumps(std::vector<double> &jumpList){
	int n=getNumberOfJumps();
	SMLMS::Jump tempJump;
	jumpList.clear();
	jumpList.resize(0);
	jumpList.reserve(n);
	for (int i=0; i<n; i++){
		tempJump = getJump(i);
		jumpList.push_back(tempJump.jumpDistance);
	}
	jumpList.shrink_to_fit();
}

std::vector<double> JumpDistanceList::getAllJumpsOfState(unsigned state){
	int n=getNumberOfJumps();
	int j = 0;
	int i;
	SMLMS::Jump tempJump;
	std::vector<double> jumpList;
	for (i=0; i<n; i++){
		tempJump = getJump(i);
		if (tempJump.state==state){
			++j;
		}
	}
	jumpList.reserve(j);
	for (i=0; i<n; i++){
		tempJump = getJump(i);
		if (tempJump.state == state){
			jumpList.push_back(tempJump.jumpDistance);
		}
	}
	jumpList.shrink_to_fit();
	return jumpList;
}

void JumpDistanceList::clearTraceStatList(){
	_traceStatList.clear();
	_traceStatList.resize(0);
}

void JumpDistanceList::calcTraceStatList(){
	checkTraceNumber();
	SMLMS::Jump tempJump;
	SMLMS::TraceStatistics tempTraceStat;
	int tempTraceLength = 1;
	int n=getNumberOfJumps();
	int formerTraceNumber = 0;
	int actualTraceNumber;
	clearTraceStatList();
	_traceStatList.reserve(_traceNumber);
	for (int i=0; i<n; i++){
		tempJump = _jumpDistanceList.at(i);
		actualTraceNumber = tempJump.trace;
		if (formerTraceNumber == 0){
			tempTraceStat.index = tempJump.trace;
			tempTraceStat.begin = i;
			formerTraceNumber = actualTraceNumber;
		}
		else if((formerTraceNumber != 0) and (formerTraceNumber<actualTraceNumber)){
			tempTraceStat.length = tempTraceLength;
			tempTraceStat.end = tempTraceStat.begin + tempTraceLength - 1;
			_traceStatList.push_back(tempTraceStat);
			tempTraceStat.index = tempJump.trace;
			tempTraceStat.begin = i;
			tempTraceLength = 1;
			formerTraceNumber = actualTraceNumber;
		}
		else ++tempTraceLength;
		
	}
	tempTraceStat.length = tempTraceLength;
	tempTraceStat.end = tempTraceStat.begin + tempTraceLength - 1;
	_traceStatList.push_back(tempTraceStat);
	_traceStatList.shrink_to_fit();
}

void JumpDistanceList::reserveMemory(const std::string &name){
	std::string line;
	std::fstream inFile(name.data());
	int lineNumber=0;
	if (inFile.is_open()){
		while(std::getline(inFile, line)){
			lineNumber ++;	
		}
		inFile.close();
		_jumpDistanceList.reserve(lineNumber);
	}
	else{
		std::stringstream errorMessage;
		errorMessage<<"Unable to read molecule list form: "<<name<<std::endl;
		SMLMS::SmlmsError error(errorMessage.str());
		throw error;

	}
}

void JumpDistanceList::readJumpDistanceList(const std::string &name){
	clearJumpDistanceList();
	SMLMS::Jump jump;
	std::string line;
	std::fstream inFile(name.data());
	// allocate momory for _jumpDistanceList
	reserveMemory(name);
	if (inFile.is_open()){
		while(std::getline(inFile, line)){
			// skip comments
			if(!line.find("#")) continue;
			std::stringstream lineContent(line);
			lineContent>>jump.trace>>jump.jumpDistance>>jump.state;
			JumpDistanceList::addJumpToEnd(jump);	
		}
		inFile.close();
		_jumpDistanceList.shrink_to_fit();
		calcTraceNumber();
		calcTraceStatList();
	}
	else{
	std::stringstream errorMessage;
	errorMessage<<"unable to read jumpdistance list from:"<<name<<std::endl;
	SMLMS::SmlmsError error(errorMessage.str());
	throw error;
	}
}

void JumpDistanceList::writeJumpDistanceList(const std::string &folderName){
	std::string name = folderName;
	name.append("/judi.txt");
	SMLMS::Jump jump;
	std::ofstream outFile;
	outFile.open(name.data());
	outFile<<std::scientific;
	outFile<<std::setprecision(6);
	/* write header lines */
	outFile<<"# ERMINE Judi file"<<std::endl;
	outFile<<"# trace \tjumpDistance[nm]\tstate"<<std::endl;
	/* write data */
	for (unsigned i=0; i<_jumpDistanceList.size(); i++){
		jump=getJump(i);
		outFile<<jump.trace<<"\t"<<jump.jumpDistance<<"\t"<<jump.state<<std::endl;
	}
	outFile.close();
}

void JumpDistanceList::calcJumpDistanceList(SMLMS::MoleculeList &molList){
	SMLMS::Jump jump;
	SMLMS::Molecule initMol, finalMol;
	int jumpNumber;
	clearJumpDistanceList();
	finalMol = molList.getMolecule(molList.getNumberOfMolecules());
	jumpNumber = molList.getNumberOfMolecules()-finalMol.trace;
	_jumpDistanceList.reserve(jumpNumber);
	for (int i=1; i<molList.getNumberOfMolecules(); i++){
		if (molList.getMolecule(i-1).trace==molList.getMolecule(i).trace){
			initMol = molList.getMolecule(i-1);
			finalMol = molList.getMolecule(i);
			jump = calculateJump(initMol, finalMol);
			addJumpToEnd(jump);
		}
	}	
	_jumpDistanceList.shrink_to_fit();
	calcTraceNumber();
	calcTraceStatList();
}

void JumpDistanceList::transferStatesToMoleculeList(SMLMS::MoleculeList &molList){
	int i,j,k;
	int start, stop;
	int maxTrace;
	std::vector<int> states;
	maxTrace  = traceNumber();
	for (i=1; i<maxTrace+1; i++){
		if (molList.getTraceIndices(i, start, stop)){
			if((stop-start)==0){
				std::cout<<"warning: ermine cannot determine any state information on trace "<<i<<" of length 1."<<std::endl;
			}
			states = getTraceStates(i);
			if(unsigned(stop-start) != states.size()){
				std::stringstream errorMessage;
				errorMessage<<"error: trace: "<<i<<" from judi does not fit the molecule lists correspondent."<<std::endl;
				SMLMS::SmlmsError error(errorMessage.str());
				throw error;
			}
			else{
				k=0;
				// one molecule entry more than calculated jumps !!
				molList.setMoleculeState(start, states.at(k));
				for(j=start+1; j<stop+1; j++){
					molList.setMoleculeState(j, states.at(k));
					k +=1;
				}
			}
		}
		else{
			std::cout<<"warning: ermine has no information about trace: "<<i<<"."<<std::endl;
		}
	}	
}
}/* namespace */
