/* ######################################################################
* File Name: judi.hpp
* Project: SMLMS
* Version: 17.03
* Creation Date: 03.03.2016
* Created By Sebastian Malkusch
* <malkusch@chemie.uni-frankfurt.de>
* Goethe University of Frankfurt
* Physical and Theoretical Chemistry
* Single Molecule Biophysics
###################################################################### */


#ifndef Judi_hpp
#define Judi_hpp
#include <vector>
#include <string>
#include "header/smlmsMolecules.hpp"
#include "header/smlmsContainer.hpp"


namespace SMLMS{

 
class JumpDistanceList{
	private:
		std::vector<Jump> _jumpDistanceList;
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
		/* special functions*/
		void calcTraceNumber();
		Jump calculateJump(Molecule &initMol, Molecule &finalMol);
		void addJumpToEnd(Jump &);
		Jump getJump(int);
		Jump getJump(int) const;
		std::vector<double> getTraceJumps(int traceNumber);
		std::vector<double> getTraceJumps(int traceNumber) const;
		void setTraceStates(int traceNumber, std::vector<int> &stateTrace);
		std::vector<int> getTraceStates(int TraceNumber);
		std::vector<int> getTraceStates(int TraceNumber) const;
		int getNumberOfJumps();
		int getNumberOfJumps() const;
		void getAllJumps(std::vector<double> &);
		std::vector<double> getAllJumpsOfState(int state);
		void clearJumpDistanceList();
		void readJumpDistanceList(std::string const&);
		void writeJumpDistanceList(std::string &);
		void calcJumpDistanceList(MoleculeList &);
		//void transferStatesToMoleculeList(MoleculeList &);
};/* JumpDistanceList  */
}/* SMLMS */
#endif /* Judi_hpp */

