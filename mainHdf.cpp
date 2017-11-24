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


#include <iostream>
#include <string>
#include <vector>
#include "header/smlmsContainer.hpp"
#include "header/ermineFilenames.hpp"
#include "header/smlmsMicroscope.hpp"
#include "header/smlmsMolecules.hpp"
#include "header/smlmsJudi.hpp"
#include "header/smlmsHmmSequence.hpp"
#include "header/smlmsPhysModBrownLatDiff.hpp"
#include "header/smlmsExceptions.hpp"
#include "header/ermineHDF5.hpp"

int main(int argc, const char * argv[]){
	std::cout<<"Hello, World!\n"<<std::endl;
	SMLMS::Microscope mic;
	SMLMS::FileNames names;
	SMLMS::HDF5 file;
	std::string sourceName = "/home/malkusch/Dokumente/ermineTest/names.txt";
	names.setSourceFileName(sourceName);
	names.readNamesFromSourceFile();
	names.setFolderName("/home/malkusch/Dokumente/ermineTest/");	
	mic.loadMicroscope(names.microscopeName());
	SMLMS::ROI roi = {0.0, 100.0, 0.0, 100.0, 0.0, 1000.0, 0.0, 5000.0};
	SMLMS::MoleculeList mol;
	mol.readLocList(names.molListName());
	SMLMS::JumpDistanceList judi;
	judi.readJumpDistanceList(names.judiName());
	SMLMS::HMMSequence hmm;
	hmm.readHMM(names.hmmName());
	try{
	SMLMS::PhysicalModelBLD model;
	model.setMinValue(0.0);
	model.setMaxValue(500.0);
	model.setBinSize(1.0);
	model.setFolderName(names.folderName());
	model.readPhysMod(names.modelName());

	std::string fileName = "/home/malkusch/Dokumente/ermineTest/hdf.h5";
	
	file.setFileName(fileName);
	file.createFile();
	/* create groups */
	file.createDataGroup();
	file.createHmmGroup();
	file.createModelGroup();
	file.createSettingsGroup();
	/* create Data */
	file.createMicroscopeData();
	file.createRoiData();
	file.createMolData(mol);
	file.createJudiData(judi);
	file.createStatisticData();
	file.createWeightMatData(hmm);
	file.createTransMatData(hmm);
	file.createObsMatData(hmm);
	file.createAlphabetData(hmm);
	file.createModelWeightData(model);
	file.createModelDiffData(model);
	/* write data */
	file.writeMicroscopeData(mic);
	file.writeRoiData(roi);
	file.writeMolData(mol);
	file.writeJudiData(judi);
	file.writeStatisticData(hmm);
	file.writeWeightMatData(hmm);
	file.writeTransMatData(hmm);
	file.writeObsMatData(hmm);
	file.writeAlphabetData(hmm);
	}
	catch(SMLMS::SmlmsError &error){
		std::cout<<"error: "<<error.what()<<std::endl;
		return -1;
	}
	return 0;
}/* main */

