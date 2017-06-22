/* ######################################################################
* File Name: judi.cpp
* Project: SMLMS
* Version: 16.02
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
#include "header/ermineJudi.hpp"
#include "header/ermineExceptions.hpp"

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

/* Special Functions */

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
};

Jump JumpDistanceList::getJump(int position){
	return _jumpDistanceList.at(position);
};

std::vector<double> JumpDistanceList::getTraceJumps(int traceNumber){
	std::vector<double> traceJumps;
	for(int i=0; i<_jumpDistanceList.size(); i++){
		if (_jumpDistanceList.at(i).trace==traceNumber){
			traceJumps.push_back(_jumpDistanceList.at(i).jumpDistance);
		}
	}
	return traceJumps;	
};

void JumpDistanceList::setTraceState(int traceNumber, std::vector<int> &stateTrace){
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
			_jumpDistanceList.at(i+traceStartPosition).trace=stateTrace.at(i);
		} 
	}
	
};

std::vector<int> JumpDistanceList::getTraceStates(int traceNumber){
	std::vector<int> traceStates;
	for(int i=0; i<_jumpDistanceList.size(); i++){
		if (_jumpDistanceList.at(i).trace==traceNumber){
			traceStates.push_back(_jumpDistanceList.at(i).state);
		}
	}
	return traceStates;
};

int JumpDistanceList::getNumberOfJumps(){
	return _jumpDistanceList.size();
};

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

void JumpDistanceList::readJumpDistanceList(std::string name){
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
	}
	else{
	std::stringstream errorMessage;
	errorMessage<<"unable to read jumpdistance list from:"<<name<<std::endl;
	SMLMS::ErmineJudiError ermineJudiError(errorMessage.str());
	throw ermineJudiError;
	}
}

void JumpDistanceList::writeJumpDistanceList(std::string name){
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

}/* namespace */