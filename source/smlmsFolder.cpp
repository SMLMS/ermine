/* ######################################################################
* File Name:
* Project: 
* Version:
* Creation Date:
* Created By Sebastian Malkusch
* <malkusch@chemie.uni-frankfurt.de>
* Goethe University of Frankfurt
* Physical and Theoretical Chemistry
* Single Molecule Biophysics
###################################################################### */


#include <iostream>
#include <string>
#include "header/smlmsFolder.hpp"

namespace SMLMS{

// Constructor
SMLMSFolder::SMLMSFolder(){
	std::cout<<"Folder created"<<std::endl;
}

// Destructor
SMLMSFolder::~SMLMSFolder(){
	std::cout<<"SMLMSFolder removed from heap!"<<std::endl;
}

//Assessor functions of class folder
void SMLMSFolder::setFolderName(std::string name){
	_folderName = name;
}

std::string SMLMSFolder::folderName(){
	return _folderName;
}

// specific Class Functions

void SMLMSFolder::extractFolderName(std::string &fileName, std::string &algorithm){
/* Extracts a foldername for the calculated results from the given input filename:
 * it strips the filename suffix, that must be ".txt"
 */
	std::string suffix(".txt");
	if (fileName.find(suffix)){
		std::string name=fileName;
		name.erase(fileName.find(suffix),4);
		name.append("_");
		name.append(algorithm);
		_folderName = name;
	}
	std::cout<<_folderName<<std::endl;

}

void SMLMSFolder::printFolderName(){
/* prints the foldername top the terminal */
	std::cout<<"Results will be written to: "<<_folderName<<std::endl;
}

int SMLMSFolder::checkFolder(){
/* checks wether folder exists */
	if (_folderName.size()<1){
		std::cout<<"no folderName chosen. Use extractFolderName first!"<<std::endl;
		return 1;
	}
	
}

void SMLMSFolder::createFolder(){
/* checks wether folder already exists. If not it creates it. */
	if (checkFolder()){
		std::cout<<"created result Folder:"<<_folderName<<std::endl;
	}
}
} /* SMLMS */
