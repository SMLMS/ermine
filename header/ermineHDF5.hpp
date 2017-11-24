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


#ifndef ERMINEHDF5_hpp
#define ERMINEHDF5_hpp

#include <string>
#include "header/smlmsContainer.hpp"
#include "header/smlmsMicroscope.hpp"
#include "header/smlmsMolecules.hpp"
#include "header/smlmsJudi.hpp"
#include "header/smlmsHmmSequence.hpp"
#include "header/smlmsPhysModBrownLatDiff.hpp"
#include "H5Cpp.h"

namespace SMLMS{

 
class HDF5{
	private:
		/* file Name */
		H5std_string _fileName;
		/* group names */
		H5std_string _gnData;
		H5std_string _gnHmm;
		H5std_string _gnModel;
		H5std_string _gnSettings;
		/* data set names */
		H5std_string _dsHmmStatistics;
		H5std_string _dsJudi;
		H5std_string _dsMicroscope;
		H5std_string _dsMol;
		H5std_string _dsAlphabet;
		H5std_string _dsObsMat;
		H5std_string _dsRoi;
		H5std_string _dsTransMat;
		H5std_string _dsWeightMat;
		H5std_string _dsModelWeight;
		H5std_string _dsModelDiff;
		/* HDF5 attributes */
		H5::H5File* _file;
		H5::Group* _group;
		H5::DataSet* _data;
		H5::DataSpace* _space;
		H5::CompType* _compType;
	public:
		/* constructor */
		HDF5();
		/* destructor */
		~HDF5();
		/* copy-constructor */
		HDF5(const HDF5 &);
		/* elementary functions */
		void setFileName(std::string name);
		std::string fileName(void);
		/* read functions */
		int openFileForReading(void); 
		int readHmm();
		int readJudi();
		int readMicroscopope();
		int readMol();
		int readPhysMod();
		int readRoi();
		/* write functions */
		int writeMicroscopeData(SMLMS::Microscope &);
		int writeRoiData(SMLMS::ROI &);
		int writeMolData(SMLMS::MoleculeList &);
		int writeJudiData(SMLMS::JumpDistanceList &);
		int writeStatisticData(SMLMS::HMMSequence &);
		int writeWeightMatData(SMLMS::HMMSequence &);
		int writeTransMatData(SMLMS::HMMSequence &);
		int writeObsMatData(SMLMS::HMMSequence &);
		int writeAlphabetData(SMLMS::HMMSequence &);
		int writePhysModData(SMLMS::PhysicalModelBLD &);
		int writeHmm();
		/* create file */
		int createFile(void);
		/* create group functions */
		int createGroup(H5std_string&);
		int createDataGroup(void);
		int createHmmGroup(void);
		int createModelGroup(void);
		int createSettingsGroup(void);
		/* create data set functions */
		int createMicroscopeData(void);
		int createRoiData(void);
		int createMolData(SMLMS::MoleculeList &);
		int createJudiData(SMLMS::JumpDistanceList &);
		int createStatisticData(void);
		int createWeightMatData(SMLMS::HMMSequence &);
		int createTransMatData(SMLMS::HMMSequence &);
		int createObsMatData(SMLMS::HMMSequence &);
		int createAlphabetData(SMLMS::HMMSequence &);
		int createModelWeightData(SMLMS::PhysicalModelBLD &);
		int createModelDiffData(SMLMS::PhysicalModelBLD &);
}; /* HDF5 */
} /* SMLMS */
#endif /* ERMINEHDF5_hpp */
