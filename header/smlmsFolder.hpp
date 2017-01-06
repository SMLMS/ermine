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


#ifndef smlmsFolder_hpp
#define smlmsFolder_hpp

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
		void extractFolderName(std::string &fileName, std::string &algorithm);
		void printFolderName();
		int checkFolder();
		void createFolder();
}; /* SMLMSFolder */

} /* SMLMS */
#endif /* smlmsFolder_hpp */

