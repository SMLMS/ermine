/* ######################################################################
* File Name: filenames.hpp
* Project: SMLMS
* Version: 1602
* Creation Date: 15.03.2016
* Created By Sebastian Malkusch
* <malkusch@chemie.uni-frankfurt.de>
* Goethe University of Frankfurt
* Physical and Theoretical Chemistry
* Single Molecule Biophysics
###################################################################### */


#ifndef ErmineFilenames_hpp
#define ErmineFilenames_hpp
#include <vector>
#include <string>

namespace SMLMS{
class FileNames{
	private:
		std::string _sourceFileName;
		std::string _folderName;
		std::string _microscopeName;
		std::string _roiName;
		std::string _judiName;
		std::string _hmmName;
		std::string _modelName;
		std::string _molListName;
		std::string _archiveName;
		std::vector<std::string> _trcNames;

	public:
		/* constructor */
		FileNames();
		/* destructor */
		~FileNames();
		/* copy-constructor */
		FileNames (const FileNames &);
		/* elementary functions */
		void setSourceFileName(std::string);
		std::string sourceFileName();
		void setFolderName(std::string);
		std::string folderName();
		void setMicroscopeName(std::string);
		std::string microscopeName();
		void setRoiName(std::string);
		std::string roiName();
		void setJudiName(std::string);
		std::string judiName();
		void setHmmName(std::string);
		std::string hmmName();
		void setModelName(std::string);
		std::string modelName();
		void setMolListName(std::string);
		std::string molListName();
		void setTrcNames(std::vector<std::string>);
		std::vector<std::string> trcNames();
		void setArchiveName(std::string);
		std::string archiveName();
		/* special functions */
		int trcNumber();
		std::string getTrcName(int);
		void addTrcName(std::string);
		void clearFileNames();
		void readNamesFromSourceFile();
		/* proof functions */
		int proofModel();
}; /* FileNames */
}/* SMLMS */
#endif /* ErmineFilenames_hpp */

