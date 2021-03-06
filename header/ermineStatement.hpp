/* ######################################################################
* File Name: ermineStatement.hpp
* Project: ermine
* Version: 19.02
* Creation Date:
* Created By Sebastian Malkusch
* <malkusch@chemie.uni-frankfurt.de>
* Goethe University of Frankfurt
* Physical and Theoretical Chemistry
* Single Molecule Biophysics
###################################################################### */


#ifndef ERMINESTATEMENT_HPP
#define ERMINESTATEMENT_HPP
#include <string>

namespace SMLMS{
class Statement{
	private:
		std::string _start;
		std::string _finished;
		std::string _batch;
		std::string _mol2judi;
		std::string _initPhysMod;
		std::string _fitPhysMod;
		std::string _initHMM;
		std::string _simulate;
		std::string _evaluate;
		std::string _train;
		std::string _bestPath;
		std::string _dwellTime;
		std::string _transferStates;
		std::string _wholeCell;
		std::string _archive;
		std::string _extract;
		std::string _tidy;
	public:
		/* constructor */
		Statement();
		/* destructor */
		~Statement();
		/* copy constructor */
		Statement(const Statement &);
		/* elementary functions */
		void setStart();
		void setFinished();
		void setBatch();
		void setMol2Judi();
		void setInitPhysMod();
		void setFitPhysMod();
		void setInitHMM();
		void setSimulate();
		void setEvaluate();
		void setTrain();
		void setBestPath();
		void setDwellTime();
		void setTransferStates();
		void setWholeCell();
		void setArchive();
		void setExtract();
		void setTidy();
		/* print functions */
		void printInit();
		void printStart();
		void printFinished();
		void printBatch();
		void printMol2Judi();
		void printInitPhysMod();
		void printFitPhysMod();
		void printInitHMM();
		void printSimulate();
		void printEvaluate();
		void printTrain();
		void printBestPath();
		void printDwellTime();
		void printTransferStates();
		void printWholeCell();
		void printArchive();
		void printExtract();
		void printTidy();
};/* Satement */
}/* SMLMS */
#endif /* Name_hpp */

