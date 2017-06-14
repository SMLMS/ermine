/* ######################################################################
* File Name: filenames.cpp
* Project: SMLMS
* Version:1602
* Creation Date: 16.03.2016
* Created By Sebastian Malkusch
* <malkusch@chemie.uni-frankfurt.de>
* Goethe University of Frankfurt
* Physical and Theoretical Chemistry
* Single Molecule Biophysics
###################################################################### */


#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <iomanip>
#include "header/ermineFilenames.hpp"
#include "header/ermineExceptions.hpp"

namespace SMLMS{

/* constructor */
FileNames::FileNames(){
	std::cout<<"FileName constructor called."<<std::endl;
}

/* destructor */
FileNames::~FileNames(){
	std::cout<<"FileNames removed from Heap!"<<std::endl;
}
/* copy-constructor*/
FileNames::FileNames(const FileNames &obj){
	std::cout<<"FileNames copy-constructor called."<<std::endl;
	setSourceFileName(obj._sourceFileName);
	setFolderName(obj._folderName);
	setMicroscopeName(obj._microscopeName);
	setRoiName(obj._roiName);
	setJudiName(obj._judiName);
	setHmmName(obj._hmmName);
	setModelName(obj._modelName);
	setMolListName(obj._molListName);
	setTrcNames(obj._trcNames);
}
/* elementary functions */
void FileNames::setSourceFileName(std::string name){
	_sourceFileName = name;
}

std::string FileNames::sourceFileName(){
	return _sourceFileName;
}

void FileNames::setFolderName(std::string name){
	_folderName = name;
}

std::string FileNames::folderName(){
	return _folderName;
}

void FileNames::setMicroscopeName(std::string name){
	_microscopeName=name;
}

std::string FileNames::microscopeName(){
	return _microscopeName;
}

void FileNames::setRoiName(std::string name){
	_roiName=name;
}

std::string FileNames::roiName(){
	return _roiName;
}

void FileNames::setJudiName(std::string name){
	_judiName=name;
}

std::string FileNames::judiName(){
	return _judiName;
}

void FileNames::setHmmName(std::string name){
	_hmmName=name;
}

std::string FileNames::hmmName(){
	return _hmmName;
}

void FileNames::setModelName(std::string name){
	_modelName=name;
}

std::string FileNames::modelName(){
	return _modelName;
}

void FileNames::setMolListName(std::string name){
	_molListName=name;
}

std::string FileNames::molListName(){
	return _molListName;
}

void FileNames::setTrcNames(std::vector<std::string> names){
	_trcNames=names;
}

std::vector<std::string> FileNames::trcNames(){
	return _trcNames;
}

/* special functions */
int FileNames::trcNumber(){
	return _trcNames.size();
}

std::string FileNames::getTrcName(int position){
	return _trcNames.at(position);
}

void FileNames::addTrcName(std::string name){
	_trcNames.push_back(name);
}

void FileNames::clearFileNames(){
	_sourceFileName.clear();
	_folderName.clear();
	_microscopeName.clear();
	_roiName.clear();
	_judiName.clear();
	_hmmName.clear();
	_modelName.clear();
	_molListName.clear();
	_trcNames.clear();
}

void FileNames::readNamesFromSourceFile(){
	std::cout<<std::endl<<"reading filenames from: "<<_sourceFileName<<std::endl;
	std::string tempFileName;
	std::string line;
	std::fstream inFile(_sourceFileName.data());
	if(inFile.is_open()){
		while(std::getline(inFile, line)){
			if( !line.find("#")) continue;
			else{
				std::stringstream lineContent(line);
				lineContent>>tempFileName;
				if (tempFileName.find(".mic") !=std::string::npos){
					setMicroscopeName(tempFileName);
					std::cout<<"Microscope file is: "<<_microscopeName<<std::endl;
				}
				else if (tempFileName.find(".roi")!= std::string::npos){
					setRoiName(tempFileName);
					std::cout<<"ROI file is: "<<_roiName<<std::endl;
				}
				else if (tempFileName.find(".jud")!= std::string::npos){
					setJudiName(tempFileName);
					std::cout<<"Judi file is: "<<_judiName<<std::endl;
				}
				else if (tempFileName.find(".hmm")!= std::string::npos){
					setHmmName(tempFileName);
					std::cout<<"HMM file is: "<<_hmmName<<std::endl;
				}
				else if (tempFileName.find(".mod")!= std::string::npos){
					setModelName(tempFileName);
					std::cout<<"Physical Model file is: "<<_modelName<<std::endl;
				}
				else if (tempFileName.find(".mol")!= std::string::npos){
					setMolListName(tempFileName);
					std::cout<<"Mol file is: "<<_molListName<<std::endl;
				}
				else if (tempFileName.find(".trc")!= std::string::npos){
					addTrcName(tempFileName);
				}
				else{
					std::stringstream errorMessage;
					errorMessage<<tempFileName<<" is of an unkown data type."<<std::endl;
					SMLMS::ErmineFileNameError ermineFileNameError(errorMessage.str());
					throw ermineFileNameError;;
				}
			}
		}
	}
	else{
		std::stringstream errorMessage;
		errorMessage<<_sourceFileName<<" could not be opened."<<std::endl;
		SMLMS::ErmineFileNameError ermineFileNameError(errorMessage.str());
		throw ermineFileNameError;
	}
	inFile.close();
	std::cout<<std::endl<<"trc files are:"<<std::endl;
	for (int i=0; i<trcNumber(); i++){
		std::cout<<getTrcName(i)<<std::endl;
	}
}

}/* SMLMS */
