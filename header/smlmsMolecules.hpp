/* ######################################################################
* File Name: smlmsMolecules.hpp
* Project: ermine
* Version: 19.02
* Creation Date:29.02.2016
* Created By Sebastian Malkusch
* <malkusch@chemie.uni-frankfurt.de>
* Goethe University of Frankfurt
* Physical and Theoretical Chemistry
* Single Molecule Biophysics
###################################################################### */


#ifndef SMLMSMOLECULES_HPP
#define SMLMSMOLECULES_HPP
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
		ROI roi()const;
		void setMoleculeList(std::vector<Molecule>&);
		std::vector<Molecule> moleculeList();
	/* special functions */
		void addMoleculeToEnd(Molecule&);
		Molecule getMolecule(int);
		Molecule getMolecule(int) const;
		void addMoleculeList(MoleculeList &);
		void deleteMolecule(int);
		int getNumberOfMolecules();
		int getNumberOfMolecules() const;
		bool getTraceIndices(unsigned trace, int &start, int &stop);
		void setMoleculeState(int index, int state);
		//load functions
		void readMoleculeList(const std::string &locsName,const std::string &roiName);
		void readTrcList(SMLMS::Microscope&, const std::string&);
		void readLocList(const std::string&);
		void readROI(const std::string&);
		//write functions
		void writeMoleculeList(const std::string &folderName);
		void writeROI(const std::string &folderName);
		void writeLocList(const std::string &folderName);
		//memory management  functions
		void reserveMemory(const std::string &fileName);
		void clearMoleculeList();
		void clearLocList();
		void clearROI();
		//sort functions
		void filterMoleculeList();
};/*MoleculeList*/
}/* SMLMS */
#endif /* Molecules_hpp */

