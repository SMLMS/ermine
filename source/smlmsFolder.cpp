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
#include <sstream>
#include <sys/types.h>
#include <sys/stat.h>
#include "header/smlmsFolder.hpp"
#include "header/ermineExceptions.hpp"

namespace SMLMS{

// Constructor
SMLMSFolder::SMLMSFolder(){
	std::cout<<"SMLMSFolder created"<<std::endl;
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
void SMLMSFolder::printFolderName(){
/* prints the foldername top the terminal */
	std::cout<<"Results will be written to: "<<_folderName<<std::endl;
}

int SMLMSFolder::checkFolder(){
/* checks wether folder exists */
	if (_folderName.size()<1){
		std::stringstream errorMessage;
		errorMessage<<"_folderName is not defined!"<<std::endl;
		SMLMS::SMLMSFolderError smlmsFolderError(errorMessage.str());
		throw smlmsFolderError;
	}

	struct stat info;
	stat(_folderName.c_str(), &info);
	if( info.st_mode & S_IFDIR ){
		// S_ISDIR() doesn't exist on my windows 
		std::cout<<_folderName<<" already exists. Results will be replaced."<<std::endl;
		return 1;
	}
	else{
		std::cout<<_folderName<<" does not exist."<<std::endl;
		return 0;
	}
	
}

void SMLMSFolder::createFolder(){
/* checks wether folder already exists. If not it creates it. */
	if (checkFolder()<1){
		#ifdef __linux__
			std::cout<<
       			mkdir(_folderName.c_str(), 0777); 
   		#else
       			//need to implement for Windows
			std::stringstream errorMessage;
			errorMessage<<"creating a directory is not yet implemented for Windows OS. Shame on the lazy developer!"<<std::endl;
			SMLMS::SMLMSFolderError smlmsFolderError(errorMessage.str());
			throw smlmsFolderError;
   		#endif
		std::cout<<std::endl<<"created result Folder:"<<_folderName<<std::endl;
	}
}
} /* SMLMS */
