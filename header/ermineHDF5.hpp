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
		H5::DataSpace* _memSpace;
		H5::CompType* _compType;
		int _rank;
		hsize_t _dims[1];
		hsize_t _start[1];
		hsize_t _stride[1];
		hsize_t _count[1];
		hsize_t _block[1];
		hsize_t _entries;
		
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
		/* HDF functions */
		int archiveModel(const SMLMS::Microscope&,
					const SMLMS::MoleculeList&,
					const SMLMS::JumpDistanceList&,
					SMLMS::HMMSequence&,
					SMLMS::PhysicalModelBLD&);
		int extractModel(SMLMS::Microscope& mic,
					SMLMS::MoleculeList& mol,
					SMLMS::JumpDistanceList& judi,
					SMLMS::HMMSequence& hmm,
					SMLMS::PhysicalModelBLD& model);
	/* load hdf5 */
		int openDataSet(H5std_string groupName, H5std_string dataName);
		int tidy(void);
		/* write functions */
		int writeMicroscopeData(const SMLMS::Microscope &);
		int writeRoiData(const SMLMS::ROI &);
		int writeMolData(const SMLMS::MoleculeList &);
		int writeJudiData(const SMLMS::JumpDistanceList &);
		int writeStatisticData(SMLMS::HMMSequence &);
		int writeWeightMatData(SMLMS::HMMSequence &);
		int writeTransMatData(SMLMS::HMMSequence &);
		int writeObsMatData(SMLMS::HMMSequence &);
		int writeAlphabetData(SMLMS::HMMSequence &);
		int writeHmmData(SMLMS::HMMSequence &);
		int writePhysModData(SMLMS::PhysicalModelBLD &);
		/* read functions */
		int readMicroscopeData(SMLMS::Microscope &);
		int readRoiData(SMLMS::ROI &);
		int readMolData(SMLMS::MoleculeList &);
		int readJudiData(SMLMS::JumpDistanceList &);
		int readStatisticData(SMLMS::HMMSequence &);
		int readWeightMatData(SMLMS::HMMSequence &);
		int readTransMatData(SMLMS::HMMSequence &);
		int readObsMatData(SMLMS::HMMSequence &);
		int readAlphabetData(SMLMS::HMMSequence &);
		int readHmmData(SMLMS::HMMSequence &);
		int readPhysModData(SMLMS::PhysicalModelBLD &);
		/* create file function */
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
		int createMolData(const SMLMS::MoleculeList &);
		int createJudiData(const SMLMS::JumpDistanceList &);
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
