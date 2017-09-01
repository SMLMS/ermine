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
	setInitPhysMod();
	setFitPhysMod();
	setInitHMM();
	setSimulate();
	setEvaluate();
	setTrain();
	setBestPath();
	setDwellTime();
	setTransferStates();
	setTidy();
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
	_initPhysMod=obj._initPhysMod;
	_fitPhysMod=obj._fitPhysMod;
	_initHMM=obj._initHMM;
	_simulate=obj._simulate;
	_evaluate=obj._evaluate;
	_train=obj._train;
	_bestPath=obj._bestPath;
	_dwellTime = obj._dwellTime;
	_transferStates = obj._transferStates;
	_tidy = obj._tidy;
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
	_batch.append("The ermine calculates a batch from single judi files.");
}

void Statement::setMol2Judi(){
	_mol2judi.clear();
	_mol2judi.append("The ermine calculates a judi file from a trc file.");
}

void Statement::setInitPhysMod(){
	_initPhysMod.clear();
	_initPhysMod.append("The ermine initializes a physical model.");
}

void Statement::setFitPhysMod(){
	_fitPhysMod.clear();
	_fitPhysMod.append("The ermine fits a physical model to a judi data set.");
}

void Statement::setInitHMM(){
	_initHMM.clear();
	_initHMM.append("The ermine initializes an HMM.");
}

void Statement::setSimulate(){
	_simulate.clear();
	_simulate.append("The ermine simulates a judi data set from an HMM.");
}

void Statement::setEvaluate(){
	_evaluate.clear();
	_evaluate.append("The ermine evaluates the likelihood of a judi data set to be described by an HMM.");
}

void Statement::setTrain(){
	_train.clear();
	_train.append("The ermine trains an HMM with a judi data set.");
}

void Statement::setBestPath(){
	_bestPath.clear();
	_bestPath.append("The ermine evaluates the best state path through a given data set based on an HMM.");
}

void Statement::setDwellTime(){
	_dwellTime.clear();
	_dwellTime.append("The ermine evaluates the dwell times of all states based on the best path through a judi data set");
}

void Statement::setTransferStates(){
	_transferStates.clear();
	_transferStates.append("The ermine transfers all states from a judi to a trc data set.");
}

void Statement::setTidy(){
	_tidy.clear();
	_tidy.append("\nThe ermine is clearing the heap.\n");
}
/* print functions */
void Statement::printStart(){
	std::cout<<std::endl<<_start<<std::endl;
}

void Statement::printFinished(){
	std::cout<<std::endl<<_finished<<std::endl;
}

void Statement::printBatch(){
	std::cout<<std::endl<<_batch<<std::endl;
}

void Statement::printMol2Judi(){
	std::cout<<std::endl<<_mol2judi<<std::endl;
}

void Statement::printInitPhysMod(){
	std::cout<<std::endl<<_initPhysMod<<std::endl;
}

void Statement::printFitPhysMod(){
	std::cout<<std::endl<<_fitPhysMod<<std::endl;
}

void Statement::printInitHMM(){
	std::cout<<std::endl<<_initHMM<<std::endl;
}

void Statement::printSimulate(){
	std::cout<<std::endl<<_simulate<<std::endl;
}

void Statement::printEvaluate(){
	std::cout<<std::endl<<_evaluate<<std::endl;
}

void Statement::printTrain(){
	std::cout<<std::endl<<_train<<std::endl;
}

void Statement::printBestPath(){
	std::cout<<std::endl<<_bestPath<<std::endl;
}

void Statement::printDwellTime(){
	std::cout<<std::endl<<_dwellTime<<std::endl;
}

void Statement::printTransferStates(){
	std::cout<<std::endl<<_transferStates<<std::endl;
}

void Statement::printTidy(){
	std::cout<<std::endl<<_tidy<<std::endl;
}
}/* SMLMS  */

