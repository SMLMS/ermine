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


#ifndef Name_hpp
#define Name_hpp
#include <string>

namespace SMLMS{
class Statement{
	private:
		std::string _start;
		std::string _finished;
		std::string _batch;
		std::string _mol2judi;
		std::string _initialize;
		std::string _simulate;
		std::string _estimate;
		std::string _learn;
		std::string _viterbi;
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
		void setInitialize();
		void setSimulate();
		void setEstimate();
		void setLearn();
		void setViterbi();
		void setTidy();
		/* print functions */
		void printInit();
		void printStart();
		void printFinished();
		void printBatch();
		void printMol2Judi();
		void printInitialize();
		void printSimulate();
		void printEstimate();
		void printLearn();
		void printViterbi();
		void printTidy();
};/* Satement */
}/* SMLMS */
#endif /* Name_hpp */

