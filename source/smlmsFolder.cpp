/* ######################################################################
* File Name: smlmsFolder.cpp
* Project: ermine
* Version: 19.02
* Creation Date: 01.12.2016 
* Created By Sebastian Malkusch
* <malkusch@chemie.uni-frankfurt.de>
* Goethe University of Frankfurt
* Physical and Theoretical Chemistry
* Single Molecule Biophysics
###################################################################### */


#include <iostream>
#include <string>
#include <sstream>
#include <boost/filesystem.hpp>
//#include <sys/types.h>
//#include <sys/stat.h>
#include "header/smlmsFolder.hpp"
#include "header/ermineExceptions.hpp"

namespace fs = boost::filesystem;
namespace SMLMS{

// Constructor
SMLMSFolder::SMLMSFolder(){
	//std::cout<<"SMLMSFolder created"<<std::endl;
}

// Destructor
SMLMSFolder::~SMLMSFolder(){
	//std::cout<<"SMLMSFolder removed from heap!"<<std::endl;
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

	fs::path p (_folderName.c_str());   // p reads clearer than argv[1] in the following code

  	if (exists(p)) {   
		// does p actually exist?
    	if (is_regular_file(p)) {
			// is p a regular file?   
      		std::stringstream errorMessage;
			errorMessage<<p<<" is a file, expected a folder!"<<std::endl;
			SMLMS::SMLMSFolderError smlmsFolderError(errorMessage.str());
			throw smlmsFolderError;
		}
		else if (is_directory(p)){
			// is p a directory?
      		std::cout<<_folderName<<" already exists. Results will be replaced."<<std::endl;
			return 1;
		}
		
    	else{
			std::stringstream errorMessage;
			errorMessage<<p<<" exists, but is neither a regular file nor a directory!"<<std::endl;
			SMLMS::SMLMSFolderError smlmsFolderError(errorMessage.str());
			throw smlmsFolderError;
  		}
	}
	else{
		std::cout<<_folderName<<" does not exist."<<std::endl;
		return 0;	
	}
}

void SMLMSFolder::createFolder(){
	/* checks wether folder already exists. If not it creates it. */
	fs::path p (_folderName.c_str());   // p reads clearer than argv[1] in the following code
	if (checkFolder()<1){
	 	if(create_directory(p)){
			std::cout<<p<<" has been created"<<std::endl;
    	}
		else{
			std::stringstream errorMessage;
			errorMessage<<"Oops, the ermine was not able to create "<<_folderName<<std::endl;
			SMLMS::SMLMSFolderError smlmsFolderError(errorMessage.str());
			throw smlmsFolderError;
		}
	}
}
} /* SMLMS */
