/* ######################################################################
* File Name: statement.cpp
* Project: SMLMS
* Version: 1602
* Creation Date: 16.03.2016
* Created By Sebastian Malkusch
* <malkusch@chemie.uni-frankfurt.de>
* Goethe University of Frankfurt
* Physical and Theoretical Chemistry
* Single Molecule Biophysics
###################################################################### */


#include <iostream>
#include <string>
#include "header/ermineStatement.hpp"

namespace SMLMS{
/* constructor */
Statement::Statement(){
	setStart();
	setFinished();
	setBatch();
	setMol2Judi();
	setInitialize();
	setSimulate();
	setEstimate();
	setLearn();
	setViterbi();
}

/* destructor */
Statement::~Statement(){
	std::cout<<_finished<<std::endl;
}

/* copy constructor */
Statement::Statement(const Statement &obj){
	_start=obj._start;
	_finished=obj._finished;
	_batch=obj._batch;
	_mol2judi=obj._mol2judi;
	_initialize=obj._initialize;
	_simulate=obj._simulate;
	_estimate=obj._estimate;
	_learn=obj._learn;
	_viterbi=obj._viterbi;
}

/* elementary functions */
void Statement::setStart(){
	_start.clear();
	_start.append("Ermine by single molecule biophysics:");
}

void Statement::setFinished(){
	_finished.clear();
	_finished.append("The ermine got the job done and is going to sleep.");	
}

void Statement::setBatch(){
	_batch.clear();
	_batch.append("The ermine scurries in batch mode.");
}

void Statement::setMol2Judi(){
	_mol2judi.clear();
	_mol2judi.append("The ermine scurries in mol2judi mode.");
}

void Statement::setInitialize(){
	_initialize.clear();
	_initialize.append("The ermine scurries in initialize mode.");
}

void Statement::setSimulate(){
	_simulate.clear();
	_simulate.append("The ermine scurries in simulate mode.");
}

void Statement::setEstimate(){
	_estimate.clear();
	_estimate.append("The ermine scurries in estimate mode.");
}

void Statement::setLearn(){
	_learn.clear();
	_learn.append("The ermine scurries in learn mode.");
}

void Statement::setViterbi(){
	_viterbi.clear();
	_viterbi.append("The ermine scurries in viterbi mode.");
}

/* print functions */

void Statement::printStart(){
	std::cout<<_start<<std::endl;
}

void Statement::printFinished(){
	std::cout<<_finished<<std::endl;
}

void Statement::printBatch(){
	std::cout<<_batch<<std::endl;
}

void Statement::printMol2Judi(){
	std::cout<<_mol2judi<<std::endl;
}

void Statement::printInitialize(){
	std::cout<<_initialize<<std::endl;
}

void Statement::printSimulate(){
	std::cout<<_simulate<<std::endl;
}

void Statement::printEstimate(){
	std::cout<<_estimate<<std::endl;
}

void Statement::printLearn(){
	std::cout<<_learn<<std::endl;
}

void Statement::printViterbi(){
	std::cout<<_viterbi<<std::endl;
}
}/* SMLMS  */

