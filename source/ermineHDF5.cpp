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
	_dsWeightMat = H5std_string("weight matrix");
	_dsModelWeight = H5std_string("weigth coefficient");
	_dsModelDiff = H5std_string("diffusion coefficient");
}

/* destructor */
HDF5::~HDF5(){
	delete _file;
	delete _group;
	delete _data;
	delete _space;
	delete _compType;
	std::cout<<"HDF5 removed from Heap!"<<std::endl;
}

/* Copy Contructor */
HDF5::HDF5(const HDF5 &obj){
	std::cout<<"HDF5 copy constructor called."<<std::endl;
	_fileName = obj._fileName;
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
	_dsWeightMat = obj._dsWeightMat;
	_dsObsMat = obj._dsObsMat;
	_dsAlphabet = obj._dsAlphabet;
	_dsModelWeight = obj._dsModelWeight;
	_dsModelDiff = obj._dsModelDiff;
	/* HDF5 types */
	_file = obj._file;
	_group = obj._group;
	_data = obj._data;
	_space = obj._space;
	_compType = obj._compType;
}

/* elementary functions */
void HDF5::setFileName(std::string name){
	_fileName = name;
}

std::string HDF5::fileName(void){
	return _fileName;
}

/* read functions */

/* write functions */
int HDF5::writeMicroscopeData(SMLMS::Microscope &microscope){
	try{
	H5::Exception::dontPrint();
	_file = new H5::H5File (_fileName, H5F_ACC_RDWR);
	_group = new H5::Group(_file->openGroup(_gnSettings));
	_data =  new H5::DataSet(_group->openDataSet(_dsMicroscope));
	_compType = new H5::CompType(_data->getCompType());
	_space = new H5::DataSpace(_data->getSpace());
	int rank = _space->getSimpleExtentNdims();
	hsize_t dims[1]={1}; 
	SMLMS::Settings data={microscope.pxlSize(), microscope.intTime(), microscope.locPrec()};
	H5::DataSpace *dspace = new H5::DataSpace(rank, dims);

	hsize_t start[]  = {0};
	hsize_t stride[] = {1};
	hsize_t count[]  = {1};
	hsize_t block[]  = {1};

	_space->selectHyperslab( H5S_SELECT_SET, count, start, stride, block );
	_data->write(&data, *_compType, *dspace, *_space);

	_space->close();
	_compType->close();
	_data->close();
	_group->close();
	_file->close();
	}
	catch(H5::FileIException& error){
		std::cout<<"could not open file "<<_fileName<<std::endl;
		return -1;
	}
	catch(H5::GroupIException& error){
		std::cout<<"could not open Group settings"<<std::endl;
		return -1;
	}
	return 0;
}

int HDF5::writeRoiData(SMLMS::ROI &roi){
	try{
	H5::Exception::dontPrint();
	_file = new H5::H5File (_fileName, H5F_ACC_RDWR);
	_group = new H5::Group(_file->openGroup(_gnSettings));
	_data =  new H5::DataSet(_group->openDataSet(_dsRoi));
	_compType = new H5::CompType(_data->getCompType());
	_space = new H5::DataSpace(_data->getSpace());
	int rank = _space->getSimpleExtentNdims();
	hsize_t dims[1]={1}; 
	SMLMS::ROI data={roi.minX, roi.maxX, roi.minY, roi.maxY, roi.minFrame, roi.maxFrame, roi.minIntensity, roi.maxIntensity};
	H5::DataSpace *dspace = new H5::DataSpace(rank, dims);

	hsize_t start[]  = {0};
	hsize_t stride[] = {1};
	hsize_t count[]  = {1};
	hsize_t block[]  = {1};

	_space->selectHyperslab( H5S_SELECT_SET, count, start, stride, block );
	_data->write(&data, *_compType, *dspace, *_space);

	_space->close();
	_compType->close();
	_data->close();
	_group->close();
	_file->close();
	}
	catch(H5::FileIException& error){
		std::cout<<"could not open file "<<_fileName<<std::endl;
		return -1;
	}
	catch(H5::GroupIException& error){
		std::cout<<"could not open Group settings"<<std::endl;
		return -1;
	}
	return 0;
}

int HDF5::writeMolData(SMLMS::MoleculeList &mol){
	try{
	H5::Exception::dontPrint();
	_file = new H5::H5File (_fileName, H5F_ACC_RDWR);
	_group = new H5::Group(_file->openGroup(_gnData));
	_data =  new H5::DataSet(_group->openDataSet(_dsMol));
	_compType = new H5::CompType(_data->getCompType());
	_space = new H5::DataSpace(_data->getSpace());
	int rank = _space->getSimpleExtentNdims();
	hsize_t dims[1]={1}; 
	H5::DataSpace *dspace = new H5::DataSpace(rank, dims);
	SMLMS::Molecule data;
	
	for(hsize_t i=0; i<mol.getNumberOfMolecules(); i++){
		data = mol.getMolecule(i);

		hsize_t start[]  = {i};
		hsize_t stride[] = {1};
		hsize_t count[]  = {1};
		hsize_t block[]  = {1};

		_space->selectHyperslab( H5S_SELECT_SET, count, start, stride, block );
		_data->write(&data, *_compType, *dspace, *_space);
	}

	_space->close();
	_compType->close();
	_data->close();
	_group->close();
	_file->close();
	}
	catch(H5::FileIException& error){
		std::cout<<"could not open file "<<_fileName<<std::endl;
		return -1;
	}
	catch(H5::GroupIException& error){
		std::cout<<"could not open Group data"<<std::endl;
		return -1;
	}
	catch(H5::DataSetIException& error){
		std::cout<<"could not write all molecules.";
		std::cout<<" Make sure "<<_fileName.data()<<" can take up "<<mol.getNumberOfMolecules()<<" molecules."<<std::endl;
		return -1;
	}
	return 0;
}

int HDF5::writeJudiData(SMLMS::JumpDistanceList &judi){
	try{
	H5::Exception::dontPrint();
	_file = new H5::H5File (_fileName, H5F_ACC_RDWR);
	_group = new H5::Group(_file->openGroup(_gnData));
	_data =  new H5::DataSet(_group->openDataSet(_dsJudi));
	_compType = new H5::CompType(_data->getCompType());
	_space = new H5::DataSpace(_data->getSpace());
	int rank = _space->getSimpleExtentNdims();
	hsize_t dims[1]={1}; 
	H5::DataSpace *dspace = new H5::DataSpace(rank, dims);
	SMLMS::Jump data;
	for(hsize_t i=0; i<judi.getNumberOfJumps(); i++){
		data = judi.getJump(i);

		hsize_t start[]  = {i};
		hsize_t stride[] = {1};
		hsize_t count[]  = {1};
		hsize_t block[]  = {1};

		_space->selectHyperslab( H5S_SELECT_SET, count, start, stride, block );
		_data->write(&data, *_compType, *dspace, *_space);
	}

	_space->close();
	_compType->close();
	_data->close();
	_group->close();
	_file->close();
	}
	catch(H5::FileIException& error){
		std::cout<<"could not open file "<<_fileName<<std::endl;
		return -1;
	}
	catch(H5::GroupIException& error){
		std::cout<<"could not open Group data"<<std::endl;
		return -1;
	}
	catch(H5::DataSetIException& error){
		std::cout<<"could not write all jumps to dataset judi.";
		std::cout<<" Make sure "<<_fileName.data()<<" can take up "<<judi.getNumberOfJumps()<<" jumps."<<std::endl;
		return -1;
	}
	return 0;
}

int HDF5::writeStatisticData(SMLMS::HMMSequence &hmm){
	try{
	H5::Exception::dontPrint();
	_file = new H5::H5File (_fileName, H5F_ACC_RDWR);
	_group = new H5::Group(_file->openGroup(_gnHmm));
	_data =  new H5::DataSet(_group->openDataSet(_dsHmmStatistics));
	_compType = new H5::CompType(_data->getCompType());
	_space = new H5::DataSpace(_data->getSpace());
	int rank = _space->getSimpleExtentNdims();
	hsize_t dims[1]={1}; 
	SMLMS::HmmStatistics data={	hmm.stateNumber(),
					hmm.symbolNumber(),
					hmm.logLikelihood(),
					hmm.dof(),
					hmm.bic(),
					hmm.aic()};
	H5::DataSpace *dspace = new H5::DataSpace(rank, dims);

	hsize_t start[]  = {0};
	hsize_t stride[] = {1};
	hsize_t count[]  = {1};
	hsize_t block[]  = {1};

	_space->selectHyperslab( H5S_SELECT_SET, count, start, stride, block );
	_data->write(&data, *_compType, *dspace, *_space);

	_space->close();
	_compType->close();
	_data->close();
	_group->close();
	_file->close();
	}
	catch(H5::FileIException& error){
		std::cout<<"could not open file "<<_fileName<<std::endl;
		return -1;
	}
	catch(H5::GroupIException& error){
		std::cout<<"could not open Group "<<_gnHmm<<std::endl;
		return -1;
	}
	return 0;
}

int HDF5::writeWeightMatData(SMLMS::HMMSequence &hmm){
	try{
	H5::Exception::dontPrint();
	_file = new H5::H5File (_fileName, H5F_ACC_RDWR);
	_group = new H5::Group(_file->openGroup(_gnHmm));
	_data =  new H5::DataSet(_group->openDataSet(_dsWeightMat));
	SMLMS::Matrix matrix = hmm.equiPDF();
	_data->write(matrix.data(), H5::PredType::NATIVE_DOUBLE);
	/* tidy */
	_data->close();
	_group->close();
	_file->close();
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
		std::cout<<"could not write all states to dataset "<<_dsWeightMat;
		std::cout<<" Make sure "<<_fileName.data()<<" can take up "<<hmm.stateNumber()<<" states."<<std::endl;
		return -1;
	}
	return 0;
}

int HDF5::writeTransMatData(SMLMS::HMMSequence &hmm){
	try{
	H5::Exception::dontPrint();
	_file = new H5::H5File (_fileName, H5F_ACC_RDWR);
	_group = new H5::Group(_file->openGroup(_gnHmm));
	_data =  new H5::DataSet(_group->openDataSet(_dsTransMat));
	SMLMS::Matrix matrix = hmm.transPDF();
	_data->write(matrix.data(), H5::PredType::NATIVE_DOUBLE);
	/* tidy */
	_data->close();
	_group->close();
	_file->close();
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
		std::cout<<"could not write all states to dataset "<<_dsTransMat;
		std::cout<<" Make sure "<<_fileName.data()<<" can take up "<<hmm.stateNumber()<<" states."<<std::endl;
		return -1;
	}
	return 0;
}

int HDF5::writeObsMatData(SMLMS::HMMSequence &hmm){
	try{
	H5::Exception::dontPrint();
	_file = new H5::H5File (_fileName, H5F_ACC_RDWR);
	_group = new H5::Group(_file->openGroup(_gnHmm));
	_data =  new H5::DataSet(_group->openDataSet(_dsObsMat));
	SMLMS::Matrix matrix = hmm.obsPDF();
	_data->write(matrix.data(), H5::PredType::NATIVE_DOUBLE);
	/* tidy */
	_data->close();
	_group->close();
	_file->close();
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
		std::cout<<"could not write all states to dataset "<<_dsObsMat;
		std::cout<<" Make sure "<<_fileName.data()<<" can take up "<<hmm.stateNumber()<<" x "<<hmm.symbolNumber()<<" observations."<<std::endl;
		return -1;
	}
	return 0;
}

int HDF5::writeAlphabetData(SMLMS::HMMSequence &hmm){
	try{
	H5::Exception::dontPrint();
	_file = new H5::H5File (_fileName, H5F_ACC_RDWR);
	_group = new H5::Group(_file->openGroup(_gnHmm));
	_data =  new H5::DataSet(_group->openDataSet(_dsAlphabet));
	std::vector<double> alphabet = hmm.obsAlphabet();
	_data->write(alphabet.data(), H5::PredType::NATIVE_DOUBLE);
	/* tidy */
	_data->close();
	_group->close();
	_file->close();
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
		std::cout<<"could not write all states to dataset "<<_dsAlphabet;
		std::cout<<" Make sure "<<_fileName.data()<<" can take up "<<hmm.symbolNumber()<<" symbols."<<std::endl;
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

int HDF5::createJudiData(SMLMS::JumpDistanceList &judi){
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

int HDF5::createMolData(SMLMS::MoleculeList &mol){
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

int HDF5::createWeightMatData(SMLMS::HMMSequence &hmm){
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
	_data = new H5::DataSet(_group->createDataSet(_dsWeightMat, H5::PredType::NATIVE_DOUBLE, *_space));
	/*
 	* tidy
 	*/
	_data->close();
	_space->close();
	_group->close();
	_file->close();
	std::cout<<"Created dataset "<<_dsWeightMat<<" in HDF5 file "<<_fileName<<std::endl;
	}
	catch(H5::GroupIException& error){
		//error.printError();
		std::cout<<"data set "<<_dsWeightMat<<" already exists in "<<_fileName<<std::endl;
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
	_compType->insertMember("fix weight", HOFFSET(SMLMS::ModelState, fix), H5::PredType::NATIVE_INT);
	_compType->insertMember("weight", HOFFSET(SMLMS::ModelState, value), H5::PredType::NATIVE_DOUBLE);
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
	_compType->insertMember("fix", HOFFSET(SMLMS::ModelState, fix), H5::PredType::NATIVE_INT);
	_compType->insertMember("D [nm^2 s^-1]", HOFFSET(SMLMS::ModelState, value), H5::PredType::NATIVE_DOUBLE);
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

