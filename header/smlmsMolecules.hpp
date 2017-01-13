/* ######################################################################
* File Name: Molecules
* Project: SMLMS
* Version:16.02
* Creation Date:29.02.2016
* Created By Sebastian Malkusch
* <malkusch@chemie.uni-frankfurt.de>
* Goethe University of Frankfurt
* Physical and Theoretical Chemistry
* Single Molecule Biophysics
###################################################################### */


#ifndef Molecules_hpp
#define Molecules_hpp
#include <vector>
#include <string>
#include "header/smlmsContainer.hpp"
#include "header/smlmsMicroscope.hpp"

namespace SMLMS{

class MoleculeList{
	private:
		ROI _roi;
		std::vector<Molecule> _moleculeList;
	public:
	/*Constructor */
		MoleculeList();
	/* Copy Constructor */
		MoleculeList(const MoleculeList &);
	/* Destructor */
		~MoleculeList();
	/* operator overload */
	Molecule& operator()(unsigned);
	Molecule operator()(unsigned) const;
	/* elementary functions */
		void setRoi(ROI&);
		ROI roi();
		void setMoleculeList(std::vector<Molecule>&);
		std::vector<Molecule> moleculeList();
	/* special functions */
		void addMoleculeToEnd(Molecule&);
		Molecule getMolecule(int);
		void addMoleculeList(MoleculeList &);
		void deleteMolecule(int);
		int getNumberOfMolecules();
		//load functions
		void readMoleculeList(std::string locsName, std::string roiName);
		void readTrcList(SMLMS::Microscope&, std::string);
		void readLocList(std::string);
		void readROI(std::string);
		//write functions
		void writeMoleculeList(std::string moleculeListName, std::string roiName);
		void writeROI(std::string);
		void writeLocList(std::string);
		//clear functions
		void clearMoleculeList();
		void clearLocList();
		void clearROI();
		//sort functions
		void filterMoleculeList();
};/*MoleculeList*/
}/* SMLMS */
#endif /* Molecules_hpp */

