/* ######################################################################
* File Name: Molecules
* Project: SMLMS
* Version:16.02
* Creation Date:01.03.2016
* Created By Sebastian Malkusch
* <malkusch@chemie.uni-frankfurt.de>
* Goethe University of Frankfurt
* Physical and Theoretical Chemistry
* Single Molecule Biophysics
###################################################################### */


#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <iomanip>
#include "header/smlmsMicroscope.hpp"
#include "header/smlmsMolecules.hpp"
#include "header/smlmsExceptions.hpp"

namespace SMLMS{
/* Constructor */
MoleculeList::MoleculeList(){
	std::cout<<"MoleculeList Constructor called."<<std::endl;
}
/* Destructor */
MoleculeList::~MoleculeList(){
	std::cout<<"MoleculeList removed from Heap!"<<std::endl;
}
/* Copy Constructor */
MoleculeList::MoleculeList(const MoleculeList &obj){
	std::cout<<"Copy Constructor called."<<std::endl;
	_roi = obj._roi;
	_moleculeList = obj._moleculeList;
}

/* Operator overload */
SMLMS::Molecule& SMLMS::MoleculeList::operator()(unsigned index){
	return _moleculeList.at(index);
}

SMLMS::Molecule SMLMS::MoleculeList::operator()(unsigned index) const{
	return _moleculeList.at(index);
}
/* Elementary functions */
void MoleculeList::setRoi(ROI &inputRoi){
	_roi = inputRoi;
}

ROI MoleculeList::roi(){
	return _roi;
}

ROI MoleculeList::roi()const{
	return _roi;
}

void MoleculeList::setMoleculeList(std::vector<Molecule> &inputList){
	_moleculeList = inputList;
}

std::vector<Molecule> MoleculeList::moleculeList(){
	return _moleculeList;
}

/*  special functions */
void MoleculeList::addMoleculeToEnd(Molecule &mol){
	_moleculeList.push_back(mol);
}

Molecule MoleculeList::getMolecule(int position){
	return _moleculeList.at(position);
}

Molecule MoleculeList::getMolecule(int position)const{
	return _moleculeList.at(position);
}

void MoleculeList::addMoleculeList(MoleculeList &molList){
	int maxTrace;
	int molNum;
	Molecule mol;
	maxTrace = _moleculeList.at(_moleculeList.size()-1).trace;
	molNum = molList.getNumberOfMolecules();
	for(int i=0; i<molNum; i++){
		mol = molList.getMolecule(i);
		mol.trace += maxTrace;
		addMoleculeToEnd(mol);
	}	
}

void MoleculeList::deleteMolecule(int position){
	_moleculeList.erase(_moleculeList.begin() + position);
}

int MoleculeList::getNumberOfMolecules(){
	return _moleculeList.size();
}

int MoleculeList::getNumberOfMolecules() const{
	return _moleculeList.size();
}

bool MoleculeList::getTraceIndices(int trace, int &start, int &stop){
	bool exists = false;
	Molecule mol;
	int i;
	int molNumber = getNumberOfMolecules();
	for (i=0; i<molNumber; i++){
		mol = getMolecule(i);
		if (mol.trace == trace){
			start = i;
			exists = true;
			break;
		}
	}
	if (exists){
		for(i=start; i<molNumber; i++){
			mol = getMolecule(i);
			if (mol.trace == trace){
				stop = i;
			}
			else{break;}
		}
	}
	return exists;
}

void MoleculeList::setMoleculeState(int index, int state){
	_moleculeList.at(index).state = state;
}

/* load functions */
void MoleculeList::readMoleculeList(std::string locListName, std::string roiName){
	readLocList(locListName);
	readROI(roiName);
}

void MoleculeList::readTrcList(SMLMS::Microscope &microscope, std::string name){
	readLocList(name);
	for (int i=0; i<getNumberOfMolecules(); i++){
		_moleculeList.at(i).x=_moleculeList.at(i).x*microscope.pxlSize();
		_moleculeList.at(i).y=_moleculeList.at(i).y*microscope.pxlSize();
		_moleculeList.at(i).state=0;
	}
}

void MoleculeList::readLocList(std::string name){
	clearLocList();
	SMLMS::Molecule mol;
	std::string line;
	std::fstream inFile(name.data());
	if (inFile.is_open()){
		while(std::getline(inFile, line)){
			// skip comments //
			if (!line.find("#")) continue;
			std::stringstream lineContent(line);
			lineContent>>mol.trace>>mol.frame>>mol.x>>mol.y>>mol.state>>mol.intensity;
			// precision not implemented yet
			mol.precision = 0.0;
			addMoleculeToEnd(mol);
		}
		inFile.close();
	}
	else{
		std::stringstream errorMessage;
		errorMessage<<"Unable to red molecule list form: "<<name<<std::endl;
		SMLMS::SmlmsError error(errorMessage.str());
		throw error;
	}

}/* readLocList */

void MoleculeList::readROI(std::string name){
	clearROI();
	std::string line;
	std::vector<double> input(8);
	std::fstream inFile(name.data());
	if (inFile.is_open()){
		int i=0;
		while (std::getline(inFile, line)){
			// skip comments
    			if( !line.find("#")) continue;
			std::stringstream lineContent(line);
			lineContent>>input.at(i)>>input.at(i+1);
			i+=2;
		}
		inFile.close();
		_roi.minX=input.at(0);
		_roi.maxX=input.at(1);
		_roi.minY=input.at(2);
		_roi.maxY=input.at(3);
		_roi.minFrame=input.at(4);
		_roi.maxFrame=input.at(5);
		_roi.minIntensity=input.at(6);
		_roi.maxIntensity=input.at(7);
	}
	else{
		std::stringstream errorMessage;
		errorMessage<<"Unable to read roi from: "<<name<<std::endl;
		SMLMS::SmlmsError error(errorMessage.str());
		throw error;
		
	}
}/* readRoi */

/* write functions */

void MoleculeList::writeMoleculeList(std::string locListName, std::string roiName){
	writeLocList(locListName);
	writeROI(roiName);
}/* writeMoleculeList */

void MoleculeList::writeLocList(std::string name){
	SMLMS::Molecule mol;
	std::ofstream outFile;
	outFile.open(name.data());
	outFile<<std::scientific;
	outFile<<std::setprecision(6);
	/* header lines*/
	outFile<<"# ERMINE TRC file"<<std::endl;
	outFile<<"# trace\tframe\tx\ty\tstate\tintensity"<<std::endl;
	for(int i=0; i<_moleculeList.size(); i++){
		mol=getMolecule(i);
		outFile<<mol.trace<<"\t"<<mol.frame<<"\t"<<mol.x<<"\t"<<mol.y<<"\t"<<mol.state<<"\t"<<mol.intensity<<std::endl;	
	}	
	outFile.close();
}/* writeLocList */

void MoleculeList::writeROI(std::string name){
	std::ofstream outFile;
	outFile.open(name.data());
	outFile<<std::scientific;
	outFile<<std::setprecision(6);
	outFile<<"# ERMINE ROI file"<<std::endl;
	outFile<<"# min\tmax"<<std::endl;
	outFile<<"# x"<<std::endl;
	outFile<<_roi.minX<<"\t"<<_roi.maxX<<std::endl;
	outFile<<"# y"<<std::endl;
	outFile<<_roi.minY<<"\t"<<_roi.maxY<<std::endl;
	outFile<<"# frame"<<std::endl;
	outFile<<_roi.minFrame<<"\t"<<_roi.maxFrame<<std::endl;
	outFile<<"# intensity"<<std::endl;
	outFile<<_roi.minIntensity<<"\t"<<_roi.maxIntensity;	
	outFile.close();
}/* writeROI*/

/* clear functions */
void MoleculeList::clearMoleculeList(){
	clearLocList();
	clearROI();
}/* clearMoleculeList */

void MoleculeList::clearLocList(){
	_moleculeList.clear();
}/* clearLocList */

void MoleculeList::clearROI(){
	_roi.minX=0;
	_roi.maxX=0;
	_roi.minY=0;
	_roi.maxY=0;
	_roi.minFrame=0;
	_roi.maxFrame=0;
	_roi.minIntensity=0;
	_roi.maxIntensity=0;
};

/* filter functions */
void MoleculeList::filterMoleculeList(){
	for (int i=_moleculeList.size()-1; i>-1; i--){
		//std::cout<<i<<std::endl;
		Molecule tempMol=getMolecule(i);
		int j=0;
		if((tempMol.x<_roi.minX)
		||(tempMol.x>_roi.maxX)
		||(tempMol.y<_roi.minY)
		||(tempMol.y>_roi.maxY)
		||(tempMol.frame<_roi.minFrame-1)
		||(tempMol.frame>_roi.maxFrame)	
		||(tempMol.intensity<_roi.minIntensity)
		||(tempMol.intensity>_roi.maxIntensity)) j+=1;
		if(j>0){
			//std::cout<<"molecule ("<<i<<") should be deleted"<<std::endl;
			deleteMolecule(i);
		}
	}
}
}/*namespace*/
