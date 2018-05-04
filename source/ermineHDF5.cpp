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
#include <fstream>
#include <boost/filesystem.hpp>
#include "header/smlmsContainer.hpp"
#include "header/smlmsMicroscope.hpp"
#include "header/smlmsMolecules.hpp"
#include "header/smlmsJudi.hpp"
#include "header/smlmsHmmSequence.hpp"
#include "header/smlmsMatrix.hpp"
#include "header/smlmsPhysModBrownLatDiff.hpp"
#include "header/ermineHDF5.hpp"

namespace SMLMS{

/* constructor */
HDF5::HDF5(){
	std::cout<<"HDF5 constructor called."<<std::endl;
	/* group names */
	_gnData = H5std_string("data");
	_gnHmm = H5std_string("hmm");
	_gnModel = H5std_string("physical model");
	_gnSettings = H5std_string("settings");
	/* data set names */
	_dsHmmStatistics = H5std_string("statistics");
	_dsJudi = H5std_string("judi");
	_dsMicroscope = H5std_string("microscope");
	_dsMol = H5std_string("molecules");
	_dsAlphabet = H5std_string("observation alphabet");
	_dsObsMat = H5std_string("observation matrix");
	_dsRoi = H5std_string("roi");
	_dsTransMat = H5std_string("transition matrix");
	_dsEquiMat = H5std_string("equilibrium matrix");
	_dsModelWeight = H5std_string("weight coefficient");
	_dsModelDiff = H5std_string("diffusion coefficient");
	/* HDF5 attributes */
	_dims[0] = 1;
	_start[0] = 1;
	_stride[0] = 1;
	_count[0] = 1;
	_block[0] = 1;
	_entries = 1;
}

/* destructor */
HDF5::~HDF5(){
	delete _file;
	delete _group;
	delete _data;
	delete _space;
	delete _memSpace;
	delete _compType;
	std::cout<<"HDF5 removed from Heap!"<<std::endl;
}

/* Copy Contructor */
HDF5::HDF5(const HDF5 &obj){
	std::cout<<"HDF5 copy constructor called."<<std::endl;
	_fileName = obj._fileName;
	_folderName = obj._folderName;
	/* group names */
	_gnData = obj._gnData;
	_gnHmm = obj._gnHmm;
	_gnModel = obj._gnModel;
	_gnSettings = obj._gnSettings;
	/* data set names */
	_dsHmmStatistics = obj._dsHmmStatistics;
	_dsJudi = obj._dsJudi;
	_dsMicroscope = obj._dsMicroscope;
	_dsMol = obj._dsMol;
	_dsRoi = obj._dsRoi;
	_dsTransMat = obj._dsTransMat;
	_dsEquiMat = obj._dsEquiMat;
	_dsObsMat = obj._dsObsMat;
	_dsAlphabet = obj._dsAlphabet;
	_dsModelWeight = obj._dsModelWeight;
	_dsModelDiff = obj._dsModelDiff;
	/* HDF5 types */
	_file = obj._file;
	_group = obj._group;
	_data = obj._data;
	_space = obj._space;
	_memSpace = obj._memSpace;
	_compType = obj._compType;
	_dims[0] = obj._dims[0];
	_start[0] = obj._start[0];
	_stride[0] = obj._stride[0];
	_count[0] = obj._count[0];
	_block[0] = obj._block[0];
	_entries = obj._entries;
}

/* elementary functions */
void HDF5::setFileName(std::string name){
	_fileName = name;
	std::size_t index = name.find_last_of("/");
	_folderName = name.substr(0,index);
}

std::string HDF5::fileName(void){
	return _fileName;
}

void HDF5::setFolderName(std::string name){
	_folderName = name;
	_fileName = name.append("/archive.h5");
}

std::string HDF5::folderName(void){
	return _folderName;
}
/* HDF functions */
int HDF5::archiveModel(const SMLMS::Microscope& mic,
			const SMLMS::MoleculeList& mol,
			const SMLMS::JumpDistanceList& judi,
			SMLMS::HMMSequence& hmm,
			SMLMS::PhysicalModelBLD& model){

	SMLMS::ROI roi = mol.roi();
	/* create file */
	if (createFile()<0){return -1;}
	/* create groups */
	if (createDataGroup()<0){return -1;}
	if (createHmmGroup()<0){return -1;}
	if (createModelGroup()<0){return -1;}
	if (createSettingsGroup()<0){return -1;}
	/* create Data */
 	if (createMicroscopeData()<0){return -1;}
	if (createRoiData()<0){return -1;}
	if (createMolData(mol)<0){return -1;}
	if (createJudiData(judi)<0){return -1;}
	if (createStatisticData()<0){return -1;}
	if (createEquiMatData(hmm)<0){return -1;}
	if (createTransMatData(hmm)<0){return -1;}
	if (createObsMatData(hmm)<0){return -1;}
	if (createAlphabetData(hmm)<0){return -1;}
	if (createModelWeightData(model)<0){return -1;}
	if (createModelDiffData(model)<0){return -1;}
	/* write data */
	if (writeMicroscopeData(mic)<0){return -1;}
	if (writeRoiData(roi)<0){return -1;}
	if (writeMolData(mol)<0){return -1;}
	if (writeJudiData(judi)<0){return -1;}
	if (writeHmmData(hmm)<0){return -1;}
	if (writePhysModData(model)<0){return -1;}
	/* read data */
	return 0;
}

int HDF5::extractModel(SMLMS::Microscope& mic,
			SMLMS::MoleculeList& mol,
			SMLMS::JumpDistanceList& judi,
			SMLMS::HMMSequence& hmm,
			SMLMS::PhysicalModelBLD& model){
	/* load hdf5 */
	if (readMicroscopeData(mic)<0){return -1;}
	SMLMS::ROI roi;
	if (readRoiData(roi)<0){return -1;}
	if (readMolData(mol)<0){return-1;}
	mol.setRoi(roi);
	if (readJudiData(judi)<0){return -1;}
	if (readHmmData(hmm)<0){return -1;}
	if (readPhysModData(model)<0){return -1;}
	/* create out file name base*/
	std::string outName;
	/* write data*/
	outName = _folderName;
	outName.append("/names.txt");
	/* open names file */
	std::ofstream outFile;
	outFile.open(outName.data());

	outName = _folderName;
	outName.append("/microscope.txt");
	outFile<<outName.data()<<std::endl;
	mic.saveMicroscope(_folderName);
	
	outName = _folderName;
	outName.append("/roi.txt");
	outFile<<outName.data()<<std::endl;
	mol.writeROI(_folderName);

	outName = _folderName;
	outName.append("/mol.txt");
	outFile<<outName.data()<<std::endl;
	mol.writeLocList(_folderName);

	outName = _folderName;
	outName.append("/judi.txt");
	outFile<<outName.data()<<std::endl;
	judi.writeJumpDistanceList(_folderName);
	
	outName = _folderName;
	outName.append("/hmm.txt");
	outFile<<outName.data()<<std::endl;
	hmm.writeHMM(_folderName);

	outName = _folderName;
	outName.append("/physMod.txt");
	outFile<<outName.data()<<std::endl;
	model.writePhysMod(_folderName);

	/* close name file */
	outFile.close();
	return 0;
}

int HDF5::openDataSet(H5std_string groupName, H5std_string dataName){
	_group = new H5::Group(_file->openGroup(groupName));
	_data =  new H5::DataSet(_group->openDataSet(dataName));
	_compType = new H5::CompType(_data->getCompType());
	_space = new H5::DataSpace(_data->getSpace());
	_rank = _space->getSimpleExtentNdims();
	_memSpace = new H5::DataSpace(_rank, _dims);
	_entries =  _space->getSimpleExtentNpoints(); 

	return 0;
}

int HDF5::tidy(void){
	_space->close();
	_memSpace->close();
	_compType->close();
	_data->close();
	_group->close();
	_file->close();
	return 0;
}
/* write functions */
int HDF5::writeMicroscopeData(const SMLMS::Microscope &microscope){
	try{
	H5::Exception::dontPrint();
	/* open dataset from file */
	_file = new H5::H5File (_fileName, H5F_ACC_RDWR);
	openDataSet(_gnSettings, _dsMicroscope);
	/* data transfer */
	SMLMS::Settings tempData={microscope.pxlSize(), microscope.intTime(), microscope.locPrec()};
	/* write data to file */
	_data->write(&tempData, *_compType, *_memSpace, *_space);
	/* tidy */
	tidy();
	}
	catch(H5::FileIException& error){
		std::cout<<"could not open file "<<_fileName<<std::endl;
		return -1;
	}
	catch(H5::GroupIException& error){
		std::cout<<"could not open Group "<<_gnSettings<<std::endl;
		return -1;
	}
	catch(H5::DataSetIException& error){
		std::cout<<"could not write data to dataset "<<_dsMicroscope<<std::endl;;
		return -1;
	}
	return 0;
}

int HDF5::writeRoiData(const SMLMS::ROI &roi){
	try{
	H5::Exception::dontPrint();
	/* open data set from file */
	_file = new H5::H5File (_fileName, H5F_ACC_RDWR);
	openDataSet(_gnSettings, _dsRoi);
	/* write data to file */
	_data->write(&roi, *_compType, *_memSpace, *_space);
	/* tidy */
	tidy();
	}
	catch(H5::FileIException& error){
		std::cout<<"could not open file "<<_fileName<<std::endl;
		return -1;
	}
	catch(H5::GroupIException& error){
		std::cout<<"could not open Group "<<_gnSettings<<std::endl;
		return -1;
	}
	catch(H5::DataSetIException& error){
		std::cout<<"could not write data to dataset "<<_dsRoi<<std::endl;
		return -1;
	}
	return 0;
}

int HDF5::writeMolData(const SMLMS::MoleculeList &mol){
	try{
	H5::Exception::dontPrint();
	/* open data set from file */
	_file = new H5::H5File (_fileName, H5F_ACC_RDWR);
	openDataSet(_gnData, _dsMol);
	SMLMS::Molecule tempData;
	for(int i=0; i<_entries; i++){
		/* data transfer*/
		tempData = mol.getMolecule(i);
		/* write data to file */
		_start[0]  = i;
		_space->selectHyperslab(H5S_SELECT_SET, _count, _start, _stride, _block);
		_data->write(&tempData, *_compType, *_memSpace, *_space);
	}
	/* tidy */
	tidy();
	}
	catch(H5::FileIException& error){
		std::cout<<"could not open file "<<_fileName<<std::endl;
		return -1;
	}
	catch(H5::GroupIException& error){
		std::cout<<"could not open Group "<<_gnData<<std::endl;
		return -1;
	}
	catch(H5::DataSetIException& error){
		std::cout<<"could not write data to dataset "<<_dsMol<<std::endl;
		return -1;
	}
	return 0;
}

int HDF5::writeJudiData(const SMLMS::JumpDistanceList &judi){
	try{
	H5::Exception::dontPrint();
	/* open data set from file */
	_file = new H5::H5File (_fileName, H5F_ACC_RDWR);
	openDataSet(_gnData, _dsJudi);
	SMLMS::Jump tempData;
	for(int i=0; i<_entries; i++){
		/* data transfer */
		tempData = judi.getJump(i);
		/* write data to file file */
		_start[0]  = i;
		_space->selectHyperslab( H5S_SELECT_SET, _count, _start, _stride, _block);
		_data->write(&tempData, *_compType, *_memSpace, *_space);
	}
	/* tidy */
	tidy();
	}
	catch(H5::FileIException& error){
		std::cout<<"could not open file "<<_fileName<<std::endl;
		return -1;
	}
	catch(H5::GroupIException& error){
		std::cout<<"could not open Group "<<_gnData<<std::endl;
		return -1;
	}
	catch(H5::DataSetIException& error){
		std::cout<<"could not write data to dataset "<<_dsJudi<<std::endl;
		return -1;
	}
	return 0;
}

int HDF5::writeStatisticData(SMLMS::HMMSequence& hmm){
	try{
	H5::Exception::dontPrint();
	/* open data set from file */
	_file = new H5::H5File (_fileName, H5F_ACC_RDWR);
	openDataSet(_gnHmm, _dsHmmStatistics);
	/* data transfer */
	SMLMS::HmmStatistics tempData={	hmm.stateNumber(),
					hmm.symbolNumber(),
					hmm.logLikelihood(),
					hmm.dof(),
					hmm.bic(),
					hmm.aic()};
	/* write data to file */
	_data->write(&tempData, *_compType, *_memSpace, *_space);
	/* tidy */
	tidy();
	}
	catch(H5::FileIException& error){
		std::cout<<"could not open file "<<_fileName<<std::endl;
		return -1;
	}
	catch(H5::GroupIException& error){
		std::cout<<"could not open Group "<<_gnHmm<<std::endl;
		return -1;
	}
	catch(H5::DataSetIException& error){
		std::cout<<"could not read data from dataset "<<_dsHmmStatistics<<std::endl;
		return -1;
	}
	return 0;
}

int HDF5::writeEquiMatData(SMLMS::HMMSequence& hmm){
	try{
	H5::Exception::dontPrint();
	/* open data set from file */
	_file = new H5::H5File (_fileName, H5F_ACC_RDWR);
	openDataSet(_gnHmm, _dsEquiMat);
	/* transfer data */
	SMLMS::Matrix matrix = hmm.equiPDF();
	/* write data to file */
	_data->write(matrix.data(), H5::PredType::NATIVE_DOUBLE);
	/* tidy */
	tidy();
	}
	catch(H5::FileIException& error){
		std::cout<<"could not open file "<<_fileName<<std::endl;
		return -1;
	}
	catch(H5::GroupIException& error){
		std::cout<<"could not open Group "<<_gnHmm<<std::endl;
		return -1;
	}
	catch(H5::DataSetIException& error){
		std::cout<<"could not write to dataset "<<_dsEquiMat<<std::endl;
		return -1;
	}
	return 0;
}

int HDF5::writeTransMatData(SMLMS::HMMSequence& hmm){
	try{
	H5::Exception::dontPrint();
	/* open data set from file */
	_file = new H5::H5File (_fileName, H5F_ACC_RDWR);
	openDataSet(_gnHmm, _dsTransMat);
	/* transfer data */
	SMLMS::Matrix matrix = hmm.transPDF();
	/* write data to file */
	_data->write(matrix.data(), H5::PredType::NATIVE_DOUBLE);
	/* tidy */
	tidy();
	}
	catch(H5::FileIException& error){
		std::cout<<"could not open file "<<_fileName<<std::endl;
		return -1;
	}
	catch(H5::GroupIException& error){
		std::cout<<"could not open Group "<<_gnHmm<<std::endl;
		return -1;
	}
	catch(H5::DataSetIException& error){
		std::cout<<"could not write to dataset "<<_dsTransMat<<std::endl;
		return -1;
	}
	return 0;
}

int HDF5::writeObsMatData(SMLMS::HMMSequence& hmm){
	try{
	H5::Exception::dontPrint();
	/* open data set from file */
	_file = new H5::H5File (_fileName, H5F_ACC_RDWR);
	openDataSet(_gnHmm, _dsObsMat);
	/* transfer data */
	SMLMS::Matrix matrix = hmm.obsPDF();
	/* write data to file */
	_data->write(matrix.data(), H5::PredType::NATIVE_DOUBLE);
	/* tidy */
	tidy();
	}
	catch(H5::FileIException& error){
		std::cout<<"could not open file "<<_fileName<<std::endl;
		return -1;
	}
	catch(H5::GroupIException& error){
		std::cout<<"could not open Group "<<_gnHmm<<std::endl;
		return -1;
	}
	catch(H5::DataSetIException& error){
		std::cout<<"could not write to dataset "<<_dsObsMat<<std::endl;
		return -1;
	}
	return 0;
}

int HDF5::writeAlphabetData(SMLMS::HMMSequence& hmm){
	try{
	H5::Exception::dontPrint();
	/* open data set from file */
	_file = new H5::H5File (_fileName, H5F_ACC_RDWR);
	openDataSet(_gnHmm, _dsAlphabet);
	/* transfer data */
	std::vector<double> alphabet = hmm.obsAlphabet();
	/* write data to file */
	_data->write(alphabet.data(), H5::PredType::NATIVE_DOUBLE);
	/* tidy */
	tidy();
	}
	catch(H5::FileIException& error){
		std::cout<<"could not open file "<<_fileName<<std::endl;
		return -1;
	}
	catch(H5::GroupIException& error){
		std::cout<<"could not open Group "<<_gnHmm<<std::endl;
		return -1;
	}
	catch(H5::DataSetIException& error){
		std::cout<<"could not write to dataset "<<_dsAlphabet<<std::endl;
		return -1;
	}
	return 0;
}

int HDF5::writeHmmData(SMLMS::HMMSequence& hmm){
	if(writeStatisticData(hmm)<0){
		return -1;
	}
	if(writeEquiMatData(hmm)<0){
		return -1;
	}
	if(writeTransMatData(hmm)<0){
		return -1;
	}
	if(writeObsMatData(hmm)<0){
		return -1;
	}
	if(writeAlphabetData(hmm)<0){
		return -1;
	}
	return 0;
}

int HDF5::writePhysModData(SMLMS::PhysicalModelBLD& pMod){
	try{
	H5::Exception::dontPrint();
	/* open data set from file */
	_file = new H5::H5File (_fileName, H5F_ACC_RDWR);
	openDataSet(_gnModel, _dsModelWeight);
	SMLMS::Matrix tempMat = pMod.paraMat();
	SMLMS::ModelState tempData;
	for(int i=0; i<_entries; i++){
		/* data transfer */
		tempData = 	{tempMat.at(i,0),
				unsigned(tempMat.at(i,1)),
				tempMat.at(i,2),
				tempMat.at(i,3),
				};
		/* write data to file */
		_start[0]  = i;
		_space->selectHyperslab( H5S_SELECT_SET, _count, _start, _stride, _block);
		_data->write(&tempData, *_compType, *_memSpace, *_space);
	}
	/* tidy */
	tidy();
	/* open data set from file */
	_file = new H5::H5File (_fileName, H5F_ACC_RDWR);
	openDataSet(_gnModel, _dsModelDiff);
	for(int i=0; i<_entries; i++){
		/* data transfer */
		tempData = 	{tempMat.at(i,4),
				unsigned(tempMat.at(i,5)),
				tempMat.at(i,6),
				tempMat.at(i,7),
				};
		/* write data to file */
		_start[0]  = i;
		_space->selectHyperslab( H5S_SELECT_SET, _count, _start, _stride, _block);
		_data->write(&tempData, *_compType, *_memSpace, *_space);
	}
	/* tidy */
	tidy();
	/* transfer data */
	}
	catch(H5::FileIException& error){
		std::cout<<"could not open file "<<_fileName<<std::endl;
		return -1;
	}
	catch(H5::GroupIException& error){
		std::cout<<"could not open Group "<<_gnModel<<std::endl;
		return -1;
	}
	catch(H5::DataSetIException& error){
		std::cout<<"could not write data to datasets "<<_dsModelWeight<<" and "<<_dsModelDiff<<std::endl;
		return -1;
	}
	return 0;
}

/* read functions */
int HDF5::readMicroscopeData(SMLMS::Microscope &microscope){
	try{
	H5::Exception::dontPrint();
	/* open dataset from file */
	_file = new H5::H5File (_fileName, H5F_ACC_RDONLY);
	openDataSet(_gnSettings, _dsMicroscope);
	/* read data from file */
	SMLMS::Settings tempData;
	_data->read(&tempData, *_compType, *_memSpace, *_space);
	/* transfer data */
	microscope.setPxlSize(tempData.pxl);
	microscope.setIntTime(tempData.dt);
	microscope.setLocPrec(tempData.prec);
	/* tidy */
	tidy();
	}
	catch(H5::FileIException& error){
		std::cout<<"could not open file "<<_fileName<<std::endl;
		return -1;
	}
	catch(H5::GroupIException& error){
		std::cout<<"could not open Group "<<_gnSettings<<std::endl;
		return -1;
	}
	catch(H5::DataSetIException& error){
		std::cout<<"could not read data from dataset "<<_dsMicroscope<<std::endl;;
		return -1;
	}
	return 0;
}

int HDF5::readRoiData(SMLMS::ROI &roi){
	try{
	H5::Exception::dontPrint();
	/* open data set from file */
	_file = new H5::H5File (_fileName, H5F_ACC_RDONLY);
	openDataSet(_gnSettings, _dsRoi);
	/* read data from file */
	_data->read(&roi, *_compType, *_memSpace, *_space);
	/* tidy */
	tidy();
	}
	catch(H5::FileIException& error){
		std::cout<<"could not open file "<<_fileName<<std::endl;
		return -1;
	}
	catch(H5::GroupIException& error){
		std::cout<<"could not open Group "<<_gnSettings<<std::endl;
		return -1;
	}
	catch(H5::DataSetIException& error){
		std::cout<<"could not read data from dataset "<<_dsRoi<<std::endl;
		return -1;
	}
	return 0;
}

int HDF5::readMolData(SMLMS::MoleculeList &mol){
	try{
	H5::Exception::dontPrint();
	/* open data set from file */
	_file = new H5::H5File (_fileName, H5F_ACC_RDONLY);
	openDataSet(_gnData, _dsMol);
	/* read data from file */
	SMLMS::Molecule tempData;
	mol.clearLocList();
	for(int i=0; i<_entries; i++){
		_start[0]  = i;
		_space->selectHyperslab(H5S_SELECT_SET, _count, _start, _stride, _block);
		_data->read(&tempData, *_compType, *_memSpace, *_space);
		/* transfer data */
		mol.addMoleculeToEnd(tempData);
	}
	/* tidy */
	tidy();
	}
	catch(H5::FileIException& error){
		std::cout<<"could not open file "<<_fileName<<std::endl;
		return -1;
	}
	catch(H5::GroupIException& error){
		std::cout<<"could not open Group "<<_gnData<<std::endl;
		return -1;
	}
	catch(H5::DataSetIException& error){
		std::cout<<"could not read  data from dataset "<<_dsMol<<std::endl;
		return -1;
	}
	return 0;
}

int HDF5::readJudiData(SMLMS::JumpDistanceList &judi){
	try{
	H5::Exception::dontPrint();
	/* open data set from file */
	_file = new H5::H5File (_fileName, H5F_ACC_RDONLY);
	openDataSet(_gnData, _dsJudi);
	/* read data from file */
	SMLMS::Jump tempData;
	judi.clearJumpDistanceList();
	for(int i=0; i<_entries; i++){
		_start[0]  = i;
		_space->selectHyperslab( H5S_SELECT_SET, _count, _start, _stride, _block);
		_data->read(&tempData, *_compType, *_memSpace, *_space);
		/* transfer data */
		judi.addJumpToEnd(tempData);
	}
	/* tidy */
	tidy();
	}
	catch(H5::FileIException& error){
		std::cout<<"could not open file "<<_fileName<<std::endl;
		return -1;
	}
	catch(H5::GroupIException& error){
		std::cout<<"could not open Group "<<_gnData<<std::endl;
		return -1;
	}
	catch(H5::DataSetIException& error){
		std::cout<<"could not read data from data Set "<<_dsJudi<<std::endl;
		return -1;
	}
	return 0;
}

int HDF5::readStatisticData(SMLMS::HMMSequence& hmm){
	try{
	H5::Exception::dontPrint();
	/* open data set from file */
	_file = new H5::H5File (_fileName, H5F_ACC_RDONLY);
	openDataSet(_gnHmm, _dsHmmStatistics);
	/* read data from file */
	SMLMS::HmmStatistics tempData;
	_data->read(&tempData, *_compType, *_memSpace, *_space);
	/* transfer data */
	hmm.setStateNumber(tempData.states);
	hmm.setSymbolNumber(tempData.symbols);
	hmm.setDof(tempData.dof);
	hmm.setBic(tempData.bic);
	hmm.setAic(tempData.aic);
	/* tidy */
	tidy();
	}
	catch(H5::FileIException& error){
		std::cout<<"could not open file "<<_fileName<<std::endl;
		return -1;
	}
	catch(H5::GroupIException& error){
		std::cout<<"could not open Group "<<_gnHmm<<std::endl;
		return -1;
	}
	catch(H5::DataSetIException& error){
		std::cout<<"could not read data from dataset "<<_dsHmmStatistics<<std::endl;
		return -1;
	}
	return 0;
}

int HDF5::readEquiMatData(SMLMS::HMMSequence& hmm){
	try{
	H5::Exception::dontPrint();
	/* open data set from file */
	_file = new H5::H5File (_fileName, H5F_ACC_RDONLY);
	openDataSet(_gnHmm, _dsEquiMat);
	/* read data from file */
	std::vector<double>  tempData(hmm.stateNumber(),0.0);
	_data->read(tempData.data(), H5::PredType::NATIVE_DOUBLE);//, *_memSpace, *_space);
	/* transfer data */
	SMLMS::Matrix tempMat(1,hmm.stateNumber());
	for (int i=0; i<hmm.stateNumber(); i++){
		tempMat.at(0,i,tempData.at(i));
	}
	hmm.setEquiPDF(tempMat);
	/* tidy */
	tidy();
	}
	catch(H5::FileIException& error){
		std::cout<<"could not open file "<<_fileName<<std::endl;
		return -1;
	}
	catch(H5::GroupIException& error){
		std::cout<<"could not open Group "<<_gnHmm<<std::endl;
		return -1;
	}
	catch(H5::DataSetIException& error){
		std::cout<<"could not read data from dataset "<<_dsEquiMat<<std::endl;
		return -1;
	}
	return 0;
}

int HDF5::readTransMatData(SMLMS::HMMSequence& hmm){
	try{
	H5::Exception::dontPrint();
	/* open data set from file */
	_file = new H5::H5File (_fileName, H5F_ACC_RDONLY);
	openDataSet(_gnHmm, _dsTransMat);
	/* read data from file */
	std::vector<double>  tempData(hmm.stateNumber()*hmm.stateNumber(),0.0);
	_data->read(tempData.data(), H5::PredType::NATIVE_DOUBLE);//, *_memSpace, *_space);
	/* transfer data */
	SMLMS::Matrix tempMat(hmm.stateNumber(),hmm.stateNumber());
	for (int i=0; i<hmm.stateNumber(); i++){
		for (int j=0; j<hmm.stateNumber(); j++){
			tempMat.at(i,j,tempData.at((i*hmm.stateNumber())+j));
		}
	}
	hmm.setTransPDF(tempMat);
	/* tidy */
	tidy();
	}
	catch(H5::FileIException& error){
		std::cout<<"could not open file "<<_fileName<<std::endl;
		return -1;
	}
	catch(H5::GroupIException& error){
		std::cout<<"could not open Group "<<_gnHmm<<std::endl;
		return -1;
	}
	catch(H5::DataSetIException& error){
		std::cout<<"could not read data from dataset "<<_dsTransMat<<std::endl;
		return -1;
	}
	return 0;
}

int HDF5::readObsMatData(SMLMS::HMMSequence& hmm){
	try{
	H5::Exception::dontPrint();
	/* open data set from file */
	_file = new H5::H5File (_fileName, H5F_ACC_RDONLY);
	openDataSet(_gnHmm, _dsObsMat);
	/* read data from file */
	std::vector<double>  tempData(hmm.symbolNumber()*hmm.stateNumber(),0.0);
	_data->read(tempData.data(), H5::PredType::NATIVE_DOUBLE);//, *_memSpace, *_space);
	/* transfer data */
	SMLMS::Matrix tempMat(hmm.stateNumber(),hmm.symbolNumber());
	for (int i=0; i<hmm.stateNumber(); i++){
		for (int j=0; j<hmm.symbolNumber(); j++){
			tempMat.at(i,j,tempData.at((i*hmm.symbolNumber())+j));
		}
	}
	hmm.setObsPDF(tempMat);
	/* tidy */
	tidy();
	}
	catch(H5::FileIException& error){
		std::cout<<"could not open file "<<_fileName<<std::endl;
		return -1;
	}
	catch(H5::GroupIException& error){
		std::cout<<"could not open Group "<<_gnHmm<<std::endl;
		return -1;
	}
	catch(H5::DataSetIException& error){
		std::cout<<"could not read data from dataset "<<_dsObsMat<<std::endl;
		return -1;
	}
	return 0;
}

int HDF5::readAlphabetData(SMLMS::HMMSequence& hmm){
	try{
	H5::Exception::dontPrint();
	/* open data set from file */
	_file = new H5::H5File (_fileName, H5F_ACC_RDONLY);
	openDataSet(_gnHmm, _dsAlphabet);
	/* read data from file */
	std::vector<double>  tempData(hmm.symbolNumber(),0.0);
	_data->read(tempData.data(), H5::PredType::NATIVE_DOUBLE);//, *_memSpace, *_space);
	/* transfer data */
	hmm.setObsAlphabet(tempData);
	/* tidy */
	tidy();
	}
	catch(H5::FileIException& error){
		std::cout<<"could not open file "<<_fileName<<std::endl;
		return -1;
	}
	catch(H5::GroupIException& error){
		std::cout<<"could not open Group "<<_gnHmm<<std::endl;
		return -1;
	}
	catch(H5::DataSetIException& error){
		std::cout<<"could not read data from dataset "<<_dsAlphabet<<std::endl;
		return -1;
	}
	return 0;
}

int HDF5::readHmmData(SMLMS::HMMSequence& hmm){
	hmm.clearHMM();
	if(readStatisticData(hmm)<0){
		return -1;
	}
	hmm.initEqui();
	if(readEquiMatData(hmm)<0){
		return -1;
	}
	hmm.initTrans();
	if(readTransMatData(hmm)<0){
		return -1;
	}
	hmm.initObs();
	if(readObsMatData(hmm)<0){
		return -1;
	}
	hmm.initObsAlphabet();
	if(readAlphabetData(hmm)<0){
		return -1;
	}
	return 0;
}

int HDF5::readPhysModData(SMLMS::PhysicalModelBLD& pMod){
	try{
	H5::Exception::dontPrint();
	/* open data set from file */
	_file = new H5::H5File (_fileName, H5F_ACC_RDONLY);
	openDataSet(_gnModel, _dsModelWeight);
	/* read data from file */
	SMLMS::Matrix tempMat(_entries, 8);
	pMod.setStateNumber(_entries);
	pMod.initParaMat();
	SMLMS::ModelState tempData;
	for(int i=0; i<_entries; i++){
		_start[0]  = i;
		_space->selectHyperslab( H5S_SELECT_SET, _count, _start, _stride, _block);
		_data->read(&tempData, *_compType, *_memSpace, *_space);
		/* transfer data */
		tempMat.at(i, 0, tempData.value);
		tempMat.at(i, 1, tempData.fix);
		tempMat.at(i, 2, tempData.min);
		tempMat.at(i, 3, tempData.max);
	}
	/* tidy */
	tidy();
	/* open data set from file */
	_file = new H5::H5File (_fileName, H5F_ACC_RDONLY);
	openDataSet(_gnModel, _dsModelDiff);
	/* read data from file */
	for(int i=0; i<_entries; i++){
		_start[0]  = i;
		_space->selectHyperslab( H5S_SELECT_SET, _count, _start, _stride, _block);
		_data->read(&tempData, *_compType, *_memSpace, *_space);
		/* transfer data */
		tempMat.at(i, 4, tempData.value);
		tempMat.at(i, 5, tempData.fix);
		tempMat.at(i, 6, tempData.min);
		tempMat.at(i, 7, tempData.max);
	}
	/* tidy */
	tidy();
	/* transfer data */
	pMod.setParaMat(tempMat);
	}
	catch(H5::FileIException& error){
		std::cout<<"could not open file "<<_fileName<<std::endl;
		return -1;
	}
	catch(H5::GroupIException& error){
		std::cout<<"could not open Group "<<_gnModel<<std::endl;
		return -1;
	}
	catch(H5::DataSetIException& error){
		std::cout<<"could not read data from data Sets "<<_dsModelWeight<<" and "<<_dsModelDiff<<std::endl;
		return -1;
	}
	return 0;
}

/* create functions */
int HDF5::createFile(void){
	if ( boost::filesystem::exists(_fileName)){
		std::cout<<_fileName<<" already exists and will be replaced."<<std::endl;

	}
	try{
		H5::Exception::dontPrint();
		std::cout<<"creaitng hdf5 file "<<_fileName<<std::endl;
		_file = new H5::H5File(_fileName, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
		/*
 		* tidy
 		*/
		_file->close();
	}
	catch(H5::FileIException& error){
		std::cout<<"could not create file "<<_fileName<<std::endl;
		return -1;
	}
	return 0;
}

/* create groups */
int HDF5::createGroup(std::string &groupName){
	try{
		H5::Exception::dontPrint();
		_file = new H5::H5File(_fileName, H5F_ACC_RDWR, H5P_DEFAULT, H5P_DEFAULT);
   		_group=new H5::Group(_file->createGroup(groupName.data()));
		/*
 		* tidy
 		*/
		_group->close();
		_file->close();
		std::cout<<"Created Group "<<groupName<<" in HDF5 file "<<_fileName<<std::endl;
	}
	catch(H5::GroupIException& error){
		std::cout<<"could not open file "<<_fileName<<std::endl;
		//error.printError();
		return -1;
	}
	catch(H5::FileIException& error){
		std::cout<<"Group "<<groupName<<"  already exists in "<<_fileName<<std::endl;
		//error.printError();
		return -1;
	}
	return 0;
}

int HDF5::createDataGroup(void){
	int success = createGroup(_gnData);
	return success;
}

int HDF5::createHmmGroup(void){
	int success = createGroup(_gnHmm);
	return success;
}

int HDF5::createModelGroup(void){
	int success = createGroup(_gnModel);
	return success;
}

int HDF5::createSettingsGroup(void){
	int success = createGroup(_gnSettings);
	return success;
}

/* create Data */
int HDF5::createStatisticData(void){
	/*
 	* define the data-type
 	*/
	try{
	H5::Exception::dontPrint();
	/*
 	* create CompType
 	*/
	_compType = new H5::CompType(sizeof(SMLMS::HmmStatistics));
	_compType->insertMember("states", HOFFSET(SMLMS::HmmStatistics, states), H5::PredType::NATIVE_INT);
	_compType->insertMember("symbols", HOFFSET(SMLMS::HmmStatistics, symbols), H5::PredType::NATIVE_INT);
	_compType->insertMember("logLikelihood", HOFFSET(SMLMS::HmmStatistics, logLikelihood), H5::PredType::NATIVE_DOUBLE);
	_compType->insertMember("dof", HOFFSET(SMLMS::HmmStatistics, dof), H5::PredType::NATIVE_INT);
	_compType->insertMember("bic", HOFFSET(SMLMS::HmmStatistics, bic), H5::PredType::NATIVE_DOUBLE);
	_compType->insertMember("aic", HOFFSET(SMLMS::HmmStatistics, aic), H5::PredType::NATIVE_DOUBLE);
	/*
 	* set dataset dimensions
 	*/
	hsize_t dim[1]={1};
	int rank = 1;
	/*
 	* create data set
 	*/
	_file = new H5::H5File(_fileName, H5F_ACC_RDWR);
	_group = new H5::Group(_file->openGroup(_gnHmm));
	_space = new H5::DataSpace(rank, dim);
	_data = new H5::DataSet(_group->createDataSet(_dsHmmStatistics, *_compType, *_space));
	/*
 	* tidy
 	*/
	_compType->close();
	_data->close();
	_space->close();
	_group->close();
	_file->close();
	std::cout<<"Created dataset "<<_dsHmmStatistics<<" in HDF5 file "<<_fileName<<std::endl;
	}
	catch(H5::GroupIException& error){
		//error.printError();
		std::cout<<"data set "<<_dsHmmStatistics<<" already exists in "<<_fileName<<std::endl;
		return -1;
	}
	
	return 0;
}

int HDF5::createJudiData(const SMLMS::JumpDistanceList &judi){
	hsize_t length = judi.getNumberOfJumps();
	/*
 	* define the data-type
 	*/
	try{
	H5::Exception::dontPrint();
	/*
 	* create CompType
 	*/
	_compType = new H5::CompType(sizeof(SMLMS::Jump));
	_compType->insertMember("trace", HOFFSET(SMLMS::Jump, trace), H5::PredType::NATIVE_INT);
	_compType->insertMember("distance [nm]", HOFFSET(SMLMS::Jump, jumpDistance), H5::PredType::NATIVE_DOUBLE);
	_compType->insertMember("state", HOFFSET(SMLMS::Jump, state), H5::PredType::NATIVE_INT);
	/*
 	* set dataset dimensions
 	*/
	hsize_t dim[1]={length};
	int rank = 1;
	/*
 	* create data set
 	*/
	_file = new H5::H5File(_fileName, H5F_ACC_RDWR);
	_group = new H5::Group(_file->openGroup(_gnData));
	_space = new H5::DataSpace(rank, dim);
	_data = new H5::DataSet(_group->createDataSet(_dsJudi, *_compType, *_space));
	/*
 	* tidy
 	*/
	_compType->close();
	_data->close();
	_space->close();
	_group->close();
	_file->close();
	std::cout<<"Created dataset "<<_dsJudi<<" in HDF5 file "<<_fileName<<std::endl;
	}
	catch(H5::GroupIException& error){
		//error.printError();
		std::cout<<"data set "<<_dsJudi<<" already exists in "<<_fileName<<std::endl;
		return -1;
	}
	
	return 0;
}

int HDF5::createMicroscopeData(void){
	/*
 	* define the data-type
 	*/
	try{
	H5::Exception::dontPrint();
	/*
 	* create CompType
 	*/
	_compType = new H5::CompType(sizeof(SMLMS::Settings));
	_compType->insertMember("pxl size [nm]", HOFFSET(SMLMS::Settings, pxl), H5::PredType::NATIVE_DOUBLE);
	_compType->insertMember("dt [s]", HOFFSET(SMLMS::Settings, dt), H5::PredType::NATIVE_DOUBLE);
	_compType->insertMember("precision [nm]", HOFFSET(SMLMS::Settings, prec), H5::PredType::NATIVE_DOUBLE);
	/*
 	* set dataset dimensions
 	*/
	hsize_t dim[1]={1};
	int rank = 1;
	/*
 	* create data set
 	*/
	_file = new H5::H5File(_fileName, H5F_ACC_RDWR);
	_group = new H5::Group(_file->openGroup(_gnSettings));
	_space = new H5::DataSpace(rank, dim);
	_data = new H5::DataSet(_group->createDataSet(_dsMicroscope, *_compType, *_space));
	/*
 	* tidy
 	*/
	_compType->close();
	_data->close();
	_space->close();
	_group->close();
	_file->close();
	std::cout<<"Created dataset "<<_dsMicroscope<<"  in HDF5 file "<<_fileName<<std::endl;
	}
	catch(H5::GroupIException& error){
		//error.printError();
		std::cout<<"dataset "<<_dsMicroscope<<" already exists in "<<_fileName<<std::endl;
		return -1;
	}
	
	return 0;
}

int HDF5::createMolData(const SMLMS::MoleculeList &mol){
	/*
 	* define the data-type
 	*/
	try{
	H5::Exception::dontPrint();
	/*
 	* create CompType
 	*/
	_compType = new H5::CompType(sizeof(SMLMS::Molecule));
	_compType->insertMember("trace", HOFFSET(SMLMS::Molecule, trace), H5::PredType::NATIVE_INT);
	_compType->insertMember("frame", HOFFSET(SMLMS::Molecule, frame), H5::PredType::NATIVE_INT);
	_compType->insertMember("x [nm]", HOFFSET(SMLMS::Molecule, x), H5::PredType::NATIVE_DOUBLE);
	_compType->insertMember("y [nm]", HOFFSET(SMLMS::Molecule, y), H5::PredType::NATIVE_DOUBLE);
	_compType->insertMember("state", HOFFSET(SMLMS::Molecule, state), H5::PredType::NATIVE_INT);
	_compType->insertMember("Intensity", HOFFSET(SMLMS::Molecule, intensity), H5::PredType::NATIVE_DOUBLE);
	_compType->insertMember("Precision [nm]", HOFFSET(SMLMS::Molecule, precision), H5::PredType::NATIVE_DOUBLE);
	/*
 	* set dataset dimensions
 	*/
	hsize_t length = mol.getNumberOfMolecules();
	hsize_t dim[1]={length};
	int rank = 1;
	/*
 	* create data set
 	*/
	_file = new H5::H5File(_fileName, H5F_ACC_RDWR);
	_group = new H5::Group(_file->openGroup(_gnData));
	_space = new H5::DataSpace(rank, dim);
	_data = new H5::DataSet(_group->createDataSet(_dsMol, *_compType, *_space));
	/*
 	* tidy
 	*/
	_compType->close();
	_data->close();
	_space->close();
	_group->close();
	_file->close();
	std::cout<<"Created dataset "<<_dsMol<<" in HDF5 file "<<_fileName<<std::endl;
	}
	catch(H5::GroupIException& error){
		//error.printError();
		std::cout<<"data set "<<_dsMol<<" already exists in "<<_fileName<<std::endl;
		return -1;
	}
	
	return 0;
}

int HDF5::createObsMatData(SMLMS::HMMSequence &hmm){
	/*
 	* define the data-type
 	*/
	try{
	H5::Exception::dontPrint();
	/*
 	* set dataset dimensions
 	*/
	hsize_t states = hmm.stateNumber();
	hsize_t obs = hmm.symbolNumber();
	hsize_t dim[2]={states, obs};
	int rank = 2;
	/*
 	* create data set
 	*/
	_file = new H5::H5File(_fileName, H5F_ACC_RDWR);
	_group = new H5::Group(_file->openGroup(_gnHmm));
	_space = new H5::DataSpace(rank, dim);
	_data = new H5::DataSet(_group->createDataSet(_dsObsMat, H5::PredType::NATIVE_DOUBLE, *_space));
	/*
 	* tidy
 	*/
	_data->close();
	_space->close();
	_group->close();
	_file->close();
	std::cout<<"Created dataset "<<_dsObsMat<<" in HDF5 file "<<_fileName<<std::endl;
	}
	catch(H5::GroupIException& error){
		//error.printError();
		std::cout<<"data set "<<_dsObsMat<<" already exists in "<<_fileName<<std::endl;
		return -1;
	}
	
	return 0;
}

int HDF5::createRoiData(void){
	/*
 	* define the data-type
 	*/
	try{
	H5::Exception::dontPrint();
	/*
 	* create CompType
 	*/
	_compType = new H5::CompType(sizeof(SMLMS::ROI));
	_compType->insertMember("minX [nm]", HOFFSET(SMLMS::ROI, minX), H5::PredType::NATIVE_DOUBLE);
	_compType->insertMember("maxX [nm]", HOFFSET(SMLMS::ROI, maxX), H5::PredType::NATIVE_DOUBLE);
	_compType->insertMember("minY [nm]", HOFFSET(SMLMS::ROI, minY), H5::PredType::NATIVE_DOUBLE);
	_compType->insertMember("maxY [nm]", HOFFSET(SMLMS::ROI, maxY), H5::PredType::NATIVE_DOUBLE);
	_compType->insertMember("minT [frames]", HOFFSET(SMLMS::ROI, minFrame), H5::PredType::NATIVE_DOUBLE);
	_compType->insertMember("maxT [frames]", HOFFSET(SMLMS::ROI, maxFrame), H5::PredType::NATIVE_DOUBLE);
	_compType->insertMember("minI [a.u.]", HOFFSET(SMLMS::ROI, minIntensity), H5::PredType::NATIVE_DOUBLE);
	_compType->insertMember("maxI [a.u.]", HOFFSET(SMLMS::ROI, maxIntensity), H5::PredType::NATIVE_DOUBLE);
	/*
 	* set dataset dimensions
 	*/
	hsize_t dim[1]={1};
	int rank = 1;
	/*
 	* create data set
 	*/
	_file = new H5::H5File(_fileName, H5F_ACC_RDWR);
	_group = new H5::Group(_file->openGroup(_gnSettings));
	_space = new H5::DataSpace(rank, dim);
	_data = new H5::DataSet(_group->createDataSet(_dsRoi, *_compType, *_space));
	/*
 	* tidy
 	*/
	_compType->close();
	_data->close();
	_space->close();
	_group->close();
	_file->close();
	std::cout<<"Created dataset "<<_dsRoi<<" in HDF5 file "<<_fileName<<std::endl;
	}
	catch(H5::GroupIException& error){
		//error.printError();
		std::cout<<"data set "<<_dsRoi<<" already exists in "<<_fileName<<std::endl;
		return -1;
	}
	
	return 0;
}

int HDF5::createTransMatData(SMLMS::HMMSequence &hmm){
	/*
 	* define the data-type
 	*/
	try{
	H5::Exception::dontPrint();
	/*
 	* set dataset dimensions
 	*/
	hsize_t states = hmm.stateNumber();
	hsize_t dim[2]={states, states};
	int rank = 2;
	/*
 	* create data set
 	*/
	_file = new H5::H5File(_fileName, H5F_ACC_RDWR);
	_group = new H5::Group(_file->openGroup(_gnHmm));
	_space = new H5::DataSpace(rank, dim);
	_data = new H5::DataSet(_group->createDataSet(_dsTransMat, H5::PredType::NATIVE_DOUBLE, *_space));
	/*
 	* tidy
 	*/
	_data->close();
	_space->close();
	_group->close();
	_file->close();
	std::cout<<"Created dataset "<<_dsTransMat<<" in HDF5 file "<<_fileName<<std::endl;
	}
	catch(H5::GroupIException& error){
		//error.printError();
		std::cout<<"data set "<<_dsTransMat<<" already exists in "<<_fileName<<std::endl;
		return -1;
	}
	
	return 0;
}

int HDF5::createEquiMatData(SMLMS::HMMSequence &hmm){
	/*
 	* define the data-type
 	*/
	try{
	H5::Exception::dontPrint();
	/*
 	* set dataset dimensions
 	*/
	hsize_t states = hmm.stateNumber();
	hsize_t dim[2]={1, states};
	int rank = 2;
	/*
 	* create data set
 	*/
	_file = new H5::H5File(_fileName, H5F_ACC_RDWR);
	_group = new H5::Group(_file->openGroup(_gnHmm));
	_space = new H5::DataSpace(rank, dim);
	_data = new H5::DataSet(_group->createDataSet(_dsEquiMat, H5::PredType::NATIVE_DOUBLE, *_space));
	/*
 	* tidy
 	*/
	_data->close();
	_space->close();
	_group->close();
	_file->close();
	std::cout<<"Created dataset "<<_dsEquiMat<<" in HDF5 file "<<_fileName<<std::endl;
	}
	catch(H5::GroupIException& error){
		//error.printError();
		std::cout<<"data set "<<_dsEquiMat<<" already exists in "<<_fileName<<std::endl;
		return -1;
	}
	
	return 0;
}

int HDF5::createAlphabetData(SMLMS::HMMSequence &hmm){
	/*
 	* define the data-type
 	*/
	try{
	H5::Exception::dontPrint();
	/*
 	* set dataset dimensions
 	*/
	hsize_t symbols = hmm.symbolNumber();
	hsize_t dim[2]={1, symbols};
	int rank = 2;
	/*
 	* create data set
 	*/
	_file = new H5::H5File(_fileName, H5F_ACC_RDWR);
	_group = new H5::Group(_file->openGroup(_gnHmm));
	_space = new H5::DataSpace(rank, dim);
	_data = new H5::DataSet(_group->createDataSet(_dsAlphabet, H5::PredType::NATIVE_DOUBLE, *_space));
	/*
 	* tidy
 	*/
	_data->close();
	_space->close();
	_group->close();
	_file->close();
	std::cout<<"Created dataset "<<_dsAlphabet<<" in HDF5 file "<<_fileName<<std::endl;
	}
	catch(H5::GroupIException& error){
		//error.printError();
		std::cout<<"data set "<<_dsAlphabet<<" already exists in "<<_fileName<<std::endl;
		return -1;
	}
	
	return 0;
}

int HDF5::createModelWeightData(SMLMS::PhysicalModelBLD &model){
	/*
 	* define the data-type
 	*/
	try{
	H5::Exception::dontPrint();
	/*
 	* create CompType
 	*/
	_compType = new H5::CompType(sizeof(SMLMS::ModelState));
	_compType->insertMember("weight", HOFFSET(SMLMS::ModelState, value), H5::PredType::NATIVE_DOUBLE);
	_compType->insertMember("fix weight", HOFFSET(SMLMS::ModelState, fix), H5::PredType::NATIVE_INT);
	_compType->insertMember("min weight", HOFFSET(SMLMS::ModelState, min), H5::PredType::NATIVE_DOUBLE);
	_compType->insertMember("max weight", HOFFSET(SMLMS::ModelState, max), H5::PredType::NATIVE_DOUBLE);
	/*
 	* set dataset dimensions
 	*/
	hsize_t length = model.stateNumber();
	hsize_t dim[1]={length};
	int rank = 1;
	/*
 	* create data set
 	*/
	_file = new H5::H5File(_fileName, H5F_ACC_RDWR);
	_group = new H5::Group(_file->openGroup(_gnModel));
	_space = new H5::DataSpace(rank, dim);
	_data = new H5::DataSet(_group->createDataSet(_dsModelWeight, *_compType, *_space));
	/*
 	* tidy
 	*/
	_compType->close();
	_data->close();
	_space->close();
	_group->close();
	_file->close();
	std::cout<<"Created dataset "<<_dsModelWeight<<" in HDF5 file "<<_fileName<<std::endl;
	}
	catch(H5::GroupIException& error){
		//error.printError();
		std::cout<<"data set "<<_dsModelWeight<<" already exists in "<<_fileName<<std::endl;
		return -1;
	}
	return 0;
}

int HDF5::createModelDiffData(SMLMS::PhysicalModelBLD &model){
	/*
 	* define the data-type
 	*/
	try{
	H5::Exception::dontPrint();
	/*
 	* create CompType
 	*/
	_compType = new H5::CompType(sizeof(SMLMS::ModelState));
	_compType->insertMember("D [nm^2 s^-1]", HOFFSET(SMLMS::ModelState, value), H5::PredType::NATIVE_DOUBLE);
	_compType->insertMember("fix", HOFFSET(SMLMS::ModelState, fix), H5::PredType::NATIVE_INT);
	_compType->insertMember("min D [nm^2 s^-1]", HOFFSET(SMLMS::ModelState, min), H5::PredType::NATIVE_DOUBLE);
	_compType->insertMember("max D [nm^2 s^-1]", HOFFSET(SMLMS::ModelState, max), H5::PredType::NATIVE_DOUBLE);
	/*
 	* set dataset dimensions
 	*/
	hsize_t length = model.stateNumber();
	hsize_t dim[1]={length};
	int rank = 1;
	/*
 	* create data set
 	*/
	_file = new H5::H5File(_fileName, H5F_ACC_RDWR);
	_group = new H5::Group(_file->openGroup(_gnModel));
	_space = new H5::DataSpace(rank, dim);
	_data = new H5::DataSet(_group->createDataSet(_dsModelDiff, *_compType, *_space));
	/*
 	* tidy
 	*/
	_compType->close();
	_data->close();
	_space->close();
	_group->close();
	_file->close();
	std::cout<<"Created dataset "<<_dsModelWeight<<" in HDF5 file "<<_fileName<<std::endl;
	}
	catch(H5::GroupIException& error){
		//error.printError();
		std::cout<<"data set "<<_dsModelWeight<<" already exists in "<<_fileName<<std::endl;
		return -1;
	}
	return 0;
}
}/* SMLMS */

