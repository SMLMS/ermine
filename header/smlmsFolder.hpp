/* ######################################################################
* File Name: smlmsFolder.hpp
* Project: ermine
* Version: 19.02
* Creation Date: 2016
* Created By Sebastian Malkusch
* <malkusch@chemie.uni-frankfurt.de>
* Goethe University of Frankfurt
* Physical and Theoretical Chemistry
* Single Molecule Biophysics
###################################################################### */


#ifndef SMLMSFOLDER_HPP
#define SMLMSFOLDER_HPP

#include<string>

namespace SMLMS{

class SMLMSFolder{
	private:
		std::string _folderName;
	public:
		// Constructor
		SMLMSFolder();
		// Destructor
		~SMLMSFolder();
		// Assessor functions of program arguments
		void setFolderName(std::string);
		std::string folderName();
		// specific Class Functions
		void printFolderName();
		int checkFolder();
		void createFolder();
}; /* SMLMSFolder */

} /* SMLMS */
#endif /* smlmsFolder_hpp */

