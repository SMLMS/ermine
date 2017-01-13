/* ######################################################################
* File Name: judi.hpp
* Project: SMLMS
* Version: 1602
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
		/* special functions*/
		Jump calculateJump(Molecule &initMol, Molecule &finalMol);
		void addJumpToEnd(Jump &);
		Jump getJump(int);
		std::vector<double> getTraceJumps(int traceNumber);
		void setTraceState(int traceNumber, std::vector<int> &stateTrace);
		std::vector<int> getTraceStates(int TraceNumber);
		int getNumberOfJumps();
		void getAllJumps(std::vector<double> &);
		std::vector<double> getAllJumpsOfState(int state);
		void clearJumpDistanceList();
		void readJumpDistanceList(std::string);
		void writeJumpDistanceList(std::string);
		void calcJumpDistanceList(MoleculeList &);
		//void transferStatesToMoleculeList(MoleculeList &);
};/* JumpDistanceList  */
}/* SMLMS */
#endif /* Judi_hpp */

