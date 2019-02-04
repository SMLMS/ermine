/* ######################################################################
* File Name: smlmsJudi.hpp
* Project: ermine
* Version: 19.02
* Creation Date: 03.03.2016
* Created By Sebastian Malkusch
* <malkusch@chemie.uni-frankfurt.de>
* Goethe University of Frankfurt
* Physical and Theoretical Chemistry
* Single Molecule Biophysics
###################################################################### */


#ifndef SMLMSJUDI_HPP
#define SMLMSJUDI_HPP
#include <vector>
#include <string>
#include "header/smlmsMolecules.hpp"
#include "header/smlmsContainer.hpp"


namespace SMLMS{

 
class JumpDistanceList{
	private:
		std::vector<Jump> _jumpDistanceList;
		std::vector<TraceStatistics> _traceStatList;
		unsigned _traceNumber;
	public:
		/* constructor */
		JumpDistanceList();
		/* destructor */
		~JumpDistanceList();
		/* copy-constructor */
		JumpDistanceList(const JumpDistanceList &);		
		/* operator overload */
		Jump& operator()(unsigned);
		Jump operator()(unsigned) const;
		/* elementary functions */
		void setJumpDistanceList(std::vector<Jump> &);
		std::vector<Jump> jumpDistanceList();
		std::vector<Jump> jumpDistanceList() const;
		unsigned traceNumber();
		unsigned traceNumber() const;
		/* proof functions */
		void checkTraceNumber();
		void checkTraceNumber()const;
		void checkTraceNumber(unsigned tempTrace);
		void checkTraceNumber(unsigned tempTrace)const;
		/* special functions*/
		void calcTraceNumber();
		Jump calculateJump(Molecule &initMol, Molecule &finalMol);
		void addJumpToEnd(Jump &);
		Jump getJump(int);
		Jump getJump(int) const;
		std::vector<double> getTraceJumps(unsigned traceNumber);
		std::vector<double> getTraceJumps(unsigned traceNumber) const;
		void setTraceStates(unsigned traceNumber, std::vector<int> &stateTrace);
		std::vector<int> getTraceStates(unsigned TraceNumber);
		std::vector<int> getTraceStates(unsigned TraceNumber) const;
		int getNumberOfJumps();
		int getNumberOfJumps() const;
		void getAllJumps(std::vector<double> &);
		std::vector<double> getAllJumpsOfState(unsigned state);
		void clearTraceStatList();
		void calcTraceStatList();
		void clearJumpDistanceList();
		void reserveMemory(const std::string& fileName);
		void readJumpDistanceList(const std::string& fileName);
		void writeJumpDistanceList(const std::string &folderName);
		void calcJumpDistanceList(MoleculeList &);
		void transferStatesToMoleculeList(MoleculeList &);
};/* JumpDistanceList  */
}/* SMLMS */
#endif /* Judi_hpp */

