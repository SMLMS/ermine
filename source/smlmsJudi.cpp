/* ######################################################################
* File Name: judi.cpp
* Project: SMLMS
* Version: 17.03
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

/* Special Functions */
void JumpDistanceList::calcTraceNumber(){
	unsigned tempTraceNumber, tn, seqLen;
	SMLMS::Jump jump;
	tempTraceNumber = 0;
	seqLen = _jumpDistanceList.size();
	for(int i=0; i<seqLen; i++){
		jump = _jumpDistanceList.at(i);
		tn = jump.trace;
		if (tn>tempTraceNumber)tempTraceNumber=tn;
	}
	_traceNumber = tempTraceNumber;
}

void JumpDistanceList::clearJumpDistanceList(){
	_jumpDistanceList.clear();
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

std::vector<double> JumpDistanceList::getTraceJumps(int traceNumber){
	std::vector<double> traceJumps;
	for(int i=0; i<_jumpDistanceList.size(); i++){
		if (_jumpDistanceList.at(i).trace==traceNumber){
			traceJumps.push_back(_jumpDistanceList.at(i).jumpDistance);
		}
	}
	return traceJumps;	
}

std::vector<double> JumpDistanceList::getTraceJumps(int traceNumber) const{
	std::vector<double> traceJumps;
	for(int i=0; i<_jumpDistanceList.size(); i++){
		if (_jumpDistanceList.at(i).trace==traceNumber){
			traceJumps.push_back(_jumpDistanceList.at(i).jumpDistance);
		}
	}
	return traceJumps;	
}

void JumpDistanceList::setTraceStates(int traceNumber, std::vector<int> &stateTrace){
	int traceLength = stateTrace.size();
	int traceStartPosition=_jumpDistanceList.size();
	int jumpsPerTrace=0;
	for (int i=0; i<_jumpDistanceList.size(); i++){
		if(_jumpDistanceList.at(i).trace==traceNumber){
			jumpsPerTrace+=1;
			if (i<=traceStartPosition)traceStartPosition=i;
		}
	}
	if (traceLength != jumpsPerTrace){
		std::cout<<"Error stateTrace length ("<<traceLength<<") and jumpTrace length("<<jumpsPerTrace<<") do not match!"<<std::endl;
	}	
	else{
		for (int i=0; i<traceLength; i++){
			_jumpDistanceList.at(i+traceStartPosition).state=stateTrace.at(i);
		} 
	}
	
}

std::vector<int> JumpDistanceList::getTraceStates(int traceNumber){
	std::vector<int> traceStates;
	for(int i=0; i<_jumpDistanceList.size(); i++){
		if (_jumpDistanceList.at(i).trace==traceNumber){
			traceStates.push_back(_jumpDistanceList.at(i).state);
		}
	}
	return traceStates;
}

std::vector<int> JumpDistanceList::getTraceStates(int traceNumber) const{
	std::vector<int> traceStates;
	for(int i=0; i<_jumpDistanceList.size(); i++){
		if (_jumpDistanceList.at(i).trace==traceNumber){
			traceStates.push_back(_jumpDistanceList.at(i).state);
		}
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
	for (int i=0; i<n; i++){
		tempJump = getJump(i);
		jumpList.push_back(tempJump.jumpDistance);
	}
}

std::vector<double> JumpDistanceList::getAllJumpsOfState(int state){
	int n=getNumberOfJumps();
	SMLMS::Jump tempJump;
	std::vector<double> jumpList;
	for (int i=0; i<n; i++){
		tempJump = getJump(i);
		if (tempJump.state==state){
			jumpList.push_back(tempJump.jumpDistance);
		}
	}
	return jumpList;
}

void JumpDistanceList::readJumpDistanceList(const std::string &name){
	clearJumpDistanceList();
	SMLMS::Jump jump;
	std::string line;
	std::fstream inFile(name.data());
	if (inFile.is_open()){
		while(std::getline(inFile, line)){
			/* skip comments */
			if(!line.find("#")) continue;
			std::stringstream lineContent(line);
			lineContent>>jump.trace>>jump.jumpDistance>>jump.state;
			JumpDistanceList::addJumpToEnd(jump);	
		}
		inFile.close();
		calcTraceNumber();
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
	for (int i=0; i<_jumpDistanceList.size(); i++){
		jump=getJump(i);
		outFile<<jump.trace<<"\t"<<jump.jumpDistance<<"\t"<<jump.state<<std::endl;
	}
	outFile.close();
}

void JumpDistanceList::calcJumpDistanceList(SMLMS::MoleculeList &molList){
	SMLMS::Jump jump;
	SMLMS::Molecule initMol, finalMol;
	clearJumpDistanceList();
	for (int i=1; i<molList.getNumberOfMolecules(); i++){
		if (molList.getMolecule(i-1).trace==molList.getMolecule(i).trace){
			initMol = molList.getMolecule(i-1);
			finalMol = molList.getMolecule(i);
			jump = calculateJump(initMol, finalMol);
			addJumpToEnd(jump);
		}
	}	
}

void JumpDistanceList::transferStatesToMoleculeList(SMLMS::MoleculeList &molList){
	int i,j,k;
	int start, stop;
	int maxTrace;
	std::vector<int> states;
	maxTrace  = traceNumber();
	for (i=1; i<maxTrace+1; i++){
		if (molList.getTraceIndices(i, start, stop)){
			states = getTraceStates(i);
			if((stop-start) != states.size()){
				std::stringstream errorMessage;
				errorMessage<<"error: trace: "<<i<<" from judi does not fit the molecule lists correspondent."<<std::endl;
				SMLMS::SmlmsError error(errorMessage.str());
				throw error;
			}
			if((stop-start)==0){
				std::cout<<"warning: ermine cannot determine any state information on trace "<<i<<" of length 1."<<std::endl;
			}
			else{
				k=0;
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
	// judi get trace states
	// molList get idices of trace
	// compare traces
	// if true: transfer states
	// else: throw error
}
}/* namespace */
