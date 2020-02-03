/* ######################################################################
* File Name: ermineDelegate.cpp
* Project: ermine
* Version: 19.02
* Creation Date: 25.04.2018
* Created By Sebastian Malkusch
* <malkusch@chemie.uni-frankfurt.de>
* Goethe University of Frankfurt
* Physical and Theoretical Chemistry
* Single Molecule Biophysics
###################################################################### */


#include <iostream>
#include "header/ermineExceptions.hpp"
#include "header/smlmsExceptions.hpp"
#include "header/ermineStatement.hpp"
#include "header/ermineHDF5.hpp"
#include "header/ermineDelegate.hpp"

namespace po=boost::program_options;
namespace SMLMS{
/* Constructor */
Delegate::Delegate(void){
	//std::cout<<"Delegate constructor called."<<std::endl;
	_errors = 0;
}

Delegate::Delegate(po::variables_map &vm){
	//std::cout<<"Delegate constructor called."<<std::endl;
	_errors = 0;
	initDelegate(vm);
}

/* Destructor */
Delegate::~Delegate(){
	//std::cout<<"Delegate removed from heap."<<std::endl;
}

/* Copy Constructor */
Delegate::Delegate(const Delegate &obj){
	//std::cout<<"Delegate copy constructor called."<<std::endl;
	_errors = obj._errors;
	_eVar = obj._eVar;
	_fileNames = obj._fileNames;
	_microscope = obj._microscope;
	_molList = obj._molList;
	_judi = obj._judi;
	_physMod = obj._physMod;
	_hmm = obj._hmm;
	//_dwellTime = obj._dwellTime;
}

/* parse functions */
int Delegate::initDelegate(po::variables_map &vm){
	try{
		_eVar.parseArguments(vm);
		_eVar.printArguments();
	}
	catch(SMLMS::ErmineParserError& error){
		std::cout<<error.what()<<std::endl;
		_errors += 1;
		return SMLMS::SMLMS_FAILURE;
	}
	catch(SMLMS::SMLMSFolderError& error){
		std::cout<<error.what()<<std::endl;
		_errors += 1;
		return SMLMS::SMLMS_FAILURE;
	}
	catch (...){
		std::cout<<"Oops, the ermine discovered an unexpected error during argument parsing and is going to rest."<<std::endl;
		_errors += 1;
		return SMLMS::SMLMS_FAILURE;
	}
	return SMLMS::SMLMS_SUCCESS;
}

/* load functions */
int Delegate::loadFileNames(void){
	_fileNames.clearFileNames();
	_fileNames.setSourceFileName(_eVar.fileNameArgument());
	_fileNames.setFolderName(_eVar.folderNameArgument());	
	try{
		std::cout<<"reading file names."<<std::endl;
		_fileNames.readNamesFromSourceFile();
	}
	catch(SMLMS::ErmineFileNameError& error){
		std::cout<<error.what()<<std::endl;
		return SMLMS::SMLMS_FAILURE;
	}
	catch(...){
		std::cout<<"Oops, the ermine discovered an unexpected error during parsing of fileName file."<<std::endl;
		return SMLMS::SMLMS_FAILURE;
	}
	return SMLMS::SMLMS_SUCCESS;
}

int Delegate::loadMicroscope(void){
	try{
		std::cout<<"reading microscope."<<std::endl;
		_microscope.loadMicroscope(_fileNames.microscopeName());
	}
	catch(SMLMS::SmlmsError& error){
		std::cout<<error.what()<<std::endl;
		return SMLMS::SMLMS_FAILURE;
	}
	catch(...){
		std::cout<<"Oops, the ermine discovered an unexpected error during argument parsing and is going to rest."<<std::endl;
	return SMLMS::SMLMS_FAILURE;
	}
	std::cout<<"int time: "<<_microscope.intTime()<<std::endl;
	std::cout<<"pxl size: "<<_microscope.pxlSize()<<std::endl;
	std::cout<<"loc prec: "<<_microscope.locPrec()<<std::endl;
	return SMLMS::SMLMS_SUCCESS;
}

int Delegate::loadTrcList(void){
	SMLMS::MoleculeList tempList;
		try{
			std::cout<<"reading ROI."<<std::endl;
			_molList.readROI(_fileNames.roiName());
			std::cout<<"reading trc file 1."<<std::endl;
			_molList.readTrcList(_microscope, _fileNames.getTrcName(0));
		}
		catch(SMLMS::SmlmsError& error){
			std::cout<<error.what()<<std::endl;
			return SMLMS::SMLMS_FAILURE;
		}
		catch(std::out_of_range &error){
			std::cout<<error.what()<<std::endl;
			return SMLMS::SMLMS_FAILURE;
		}
		catch(...){
			std::cout<<"Oops, the ermine discovered an unexpected error while loading molecues."<<std::endl;
			return SMLMS::SMLMS_FAILURE;
		}
		/* load trc and append files */
		for (int i=1; i<_fileNames.trcNumber(); i++){
			try{
				std::cout<<"reading trc file "<<i+1<<"."<<std::endl;
				tempList.readTrcList(_microscope, _fileNames.getTrcName(i));
			}
			catch(SMLMS::SmlmsError& error){
				std::cout<<error.what()<<std::endl;
				return SMLMS::SMLMS_FAILURE;
			}
			catch(std::out_of_range &error){
				std::cout<<error.what()<<std::endl;
				return SMLMS::SMLMS_FAILURE;
			}
			catch(...){
				std::cout<<"Oops, the ermine discovered an unexpected error while loading molecules."<<std::endl;
				return SMLMS::SMLMS_FAILURE;
			}
			_molList.addMoleculeList(tempList);
		}
		/* filter trc files */
		std::cout<<"filtering localizations regarding to ROI."<<std::endl;
		_molList.filterMoleculeList();
		return SMLMS::SMLMS_SUCCESS;
}

int Delegate::loadMoleculeList(void){
	try{
		std::cout<<"reading molecules."<<std::endl;
		_molList.readMoleculeList(_fileNames.molListName(), _fileNames.roiName());
	}
	catch(SMLMS::SmlmsError& error){
		std::cout<<error.what()<<std::endl;
		return SMLMS::SMLMS_FAILURE;
	}
	catch(std::out_of_range &error){
		std::cout<<error.what()<<std::endl;
		return SMLMS::SMLMS_FAILURE;
	}
	catch(...){
		std::cout<<"Oops, the ermine discovered an unexpected error while loading molecules."<<std::endl;
		return SMLMS::SMLMS_FAILURE;
	}
	_molList.filterMoleculeList();
	return SMLMS::SMLMS_SUCCESS;
}

int Delegate::loadJumpDistanceList(void){
	try{
		std::cout<<"reading judi."<<std::endl;
		_judi.readJumpDistanceList(_fileNames.judiName());
	}
	catch(SMLMS::SmlmsError &error){
		std::cout<<error.what()<<std::endl;
		return SMLMS::SMLMS_FAILURE;
	}
	catch(std::out_of_range &error){
		std::cout<<error.what()<<std::endl;
		return SMLMS::SMLMS_FAILURE;
	}
	catch(...){
		std::cout<<"Oops, the ermine discovered an unexpected error while loading judi."<<std::endl;
		return SMLMS::SMLMS_FAILURE;
	}
	return SMLMS::SMLMS_SUCCESS;
}

int Delegate::loadPhysMod(void){
	_physMod.setMinValue(_eVar.minDistArgument());
	_physMod.setMaxValue(_eVar.maxDistArgument());
	_physMod.setBinSize(_eVar.jumpIntervalArgument());
	_physMod.setMicroscope(_microscope);
	try{
		/* loading Model */
		std::cout<<"reading probability density function parameters."<<std::endl;
		_physMod.readPhysMod(_fileNames.modelName());
	}
	catch(SMLMS::SmlmsError &error){
		std::cout<<error.what()<<std::endl;
		return SMLMS::SMLMS_FAILURE;
	}
	catch(...){
		std::cout<<"Oops, the ermine discovered an unexpected error while loading probability density function parameters."<<std::endl;
		return SMLMS::SMLMS_FAILURE;
	}
	return SMLMS::SMLMS_SUCCESS;
}

int Delegate::loadHmm(void){
	try{
		std::cout<<"loading hidden markov model."<<std::endl;
		_hmm.readHMM(_fileNames.hmmName());
	}
	catch(SMLMS::SmlmsError &error){
		std::cout<<error.what()<<std::endl;
		return SMLMS::SMLMS_FAILURE;
	}
	catch(std::out_of_range &error){
		std::cout<<error.what()<<std::endl;
		return SMLMS::SMLMS_FAILURE;
	}
	catch(...){
		std::cout<<"Oops, the ermine discovered an unexpected error while loading hidden markov model."<<std::endl;
		return SMLMS::SMLMS_FAILURE;
	}
	_hmm.setStopCrit(_eVar.stopCritArgument());
	_hmm.setMaxIt(_eVar.maxItArgument());
	_hmm.initLoadedHMM();
	_hmm.checkHMM();
	_hmm.printHMM();
	return SMLMS::SMLMS_SUCCESS;
}

/* create destination folder */
int Delegate::createFolder(void){
	try{
		_eVar.makeFolder();
	}
	catch(SMLMS::SMLMSFolderError& error){
		std::cout<<error.what()<<std::endl;
		return SMLMS::SMLMS_FAILURE;
	}
	catch(...){
		std::cout<<"Oops, the ermine discovered an unexpected error during destination folder creation."<<std::endl;
		return SMLMS::SMLMS_FAILURE;
	}
	return SMLMS::SMLMS_SUCCESS;
}

/* write functions */
int Delegate::writeParser(void){
	try{
		_eVar.writeErmineParser();
	}
	catch(...){
		std::cout<<"Oops, the ermine discovered an unexpected error while writing the parser arguments."<<std::endl;
		return SMLMS::SMLMS_FAILURE;
	}
	return SMLMS::SMLMS_SUCCESS;
}

int Delegate::writeMicroscope(void){
	try{
		_microscope.saveMicroscope(_fileNames.folderName());
	}
	catch(SMLMS::SmlmsError& error){
		std::cout<<error.what()<<std::endl;
		return SMLMS::SMLMS_FAILURE;
	}
	catch(...){
		std::cout<<"Oops, the ermine discovered an unexpected error while saving the microscope."<<std::endl;
		return SMLMS::SMLMS_FAILURE;
	}

	return SMLMS::SMLMS_SUCCESS;
}

int Delegate::writeMoleculeList(void){
	try{
		_molList.writeMoleculeList(_fileNames.folderName());
	}
	catch(SMLMS::SmlmsError& error){
		std::cout<<error.what()<<std::endl;
		return SMLMS::SMLMS_FAILURE;
	}
	catch(...){
		std::cout<<"Oops, the ermine discovered an unexpected error while saving the molecule list."<<std::endl;
		return SMLMS::SMLMS_FAILURE;
	}

	return SMLMS::SMLMS_SUCCESS;
}

int Delegate::writeJumpDistanceList(void){
	try{
		_judi.writeJumpDistanceList(_fileNames.folderName());
	}
	catch(SMLMS::SmlmsError& error){
		std::cout<<error.what()<<std::endl;
		return SMLMS::SMLMS_FAILURE;
	}
	catch(...){
		std::cout<<"Oops, the ermine discovered an unexpected error while saving the jump distance list."<<std::endl;
		return SMLMS::SMLMS_FAILURE;
	}

	return SMLMS::SMLMS_SUCCESS;
}

int Delegate::writePhysMod(void){
	try{
		_physMod.writePhysMod(_fileNames.folderName());
		_physMod.writePdfSuperPos(_fileNames.folderName());
		_physMod.plotPhysicalModel(_fileNames.folderName());
	}
	catch(SMLMS::SmlmsError& error){
		std::cout<<error.what()<<std::endl;
		return SMLMS::SMLMS_FAILURE;
	}
	catch(...){
		std::cout<<"Oops, the ermine discovered an unexpected error while saving the probability density function parameter."<<std::endl;
		return SMLMS::SMLMS_FAILURE;
	}

	return SMLMS::SMLMS_SUCCESS;
}

int Delegate::writeHmm(void){
	try{
		_hmm.writeHMM(_fileNames.folderName());
	}
	catch(SMLMS::SmlmsError& error){
		std::cout<<error.what()<<std::endl;
		return SMLMS::SMLMS_FAILURE;
	}
	catch(std::out_of_range &error){
		std::cout<<error.what()<<std::endl;
		return SMLMS::SMLMS_FAILURE;
	}
	catch(...){
		std::cout<<"Oops, the ermine discovered an unexpected error while saving the hidden markov model."<<std::endl;
		return SMLMS::SMLMS_FAILURE;
	}

	return SMLMS::SMLMS_SUCCESS;
}

/* run algorithm functions */
int Delegate::run(void){
	SMLMS::Statement statement;
	statement.printStart();
	loadFileNames();
	/* analysis */
	if (_eVar.algorithmArgument()=="batch"){
		statement.printBatch();
		createFolder();
		_eVar.writeErmineParser();
		runBatchAlgorithm();
	}
	else if(_eVar.algorithmArgument()=="mol2judi"){
		statement.printMol2Judi();
		createFolder();
		_eVar.writeErmineParser();
		runMol2JudiAlgorithm();
	}
	else if(_eVar.algorithmArgument()=="initPhysMod"){
		statement.printInitPhysMod();
		createFolder();
		_eVar.writeErmineParser();
		runInitPhysModAlgorithm();
	}
	else if(_eVar.algorithmArgument()=="fitPhysMod"){
		statement.printFitPhysMod();
		createFolder();
		_eVar.writeErmineParser();
		runFitPhysModAlgorithm();
	}
	else if(_eVar.algorithmArgument()=="initHMM"){
		statement.printInitHMM();
		createFolder();
		_eVar.writeErmineParser();
		if(_fileNames.proofModel()){runInitHmmWithPhysModAlgorithm();}
		else{runInitHmmAlgorithm();}
	}
	else if(_eVar.algorithmArgument()=="evaluate"){
		statement.printEvaluate();
		createFolder();
		_eVar.writeErmineParser();
		runEvaluateAlgorithm();
	}
	else if(_eVar.algorithmArgument()=="train"){
		statement.printTrain();
		createFolder();
		_eVar.writeErmineParser();
		if(_fileNames.proofModel()){runTrainWithPhysModAlgorithm();}
		else{runTrainAlgorithm();}
	}
	else if(_eVar.algorithmArgument()=="bestPath"){
		statement.printBestPath();
		createFolder();
		_eVar.writeErmineParser();
		runBestPathAlgorithm();
	}
	else if(_eVar.algorithmArgument()=="transferStates"){
		statement.printTransferStates();
		createFolder();
		_eVar.writeErmineParser();
		runTransferStatesAlgorithm();
	}
	else if(_eVar.algorithmArgument()=="archive"){
		statement.printArchive();
		createFolder();
		_eVar.writeErmineParser();
		runArchiveAlgorithm();
	}
	else if(_eVar.algorithmArgument()=="extract"){
		statement.printExtract();
		runExtractAlgorithm();
	}
	else if(_eVar.algorithmArgument()=="simulate"){
		statement.printSimulate();
		createFolder();
		_eVar.writeErmineParser();
		runSimulateAlgorithm();
	}
	else if(_eVar.algorithmArgument()=="wholeCell"){
		statement.printWholeCell();
		createFolder();
		_eVar.writeErmineParser();
		runWholeCellAnalysis();
	}
	else{
		std::cout<<_eVar.algorithmArgument();
		std::cout<<" is not a valid algorithm identifier. (see ermine -h for help)"<<std::endl;
		_errors += 1;
	}
	statement.printFinished();
	std::cout<<"Result:\nProceeding of "<<_eVar.algorithmArgument()<<" analysis completed with "<<_errors<<" errors."<<std::endl;
	statement.printTidy(); 
	return SMLMS::SMLMS_SUCCESS;
}

int Delegate::runBatchAlgorithm(void){
	/* load functions */
	_errors += loadMicroscope();
	_errors += loadTrcList();
	if (_errors){return SMLMS::SMLMS_FAILURE;}
	/* write functions */
	_errors += writeMicroscope();
	_errors += writeMoleculeList();
	if (_errors){return SMLMS::SMLMS_FAILURE;}
	return SMLMS::SMLMS_SUCCESS;
}

int Delegate::runMol2JudiAlgorithm(void){
	/* load functions */
	_errors += loadMicroscope();
	_errors += loadMoleculeList();
	if(_errors){return SMLMS::SMLMS_FAILURE;}
	/* analysis functions */
	_judi.calcJumpDistanceList(_molList);
	/* write functions */
	_errors += writeMicroscope();
	_errors += writeMoleculeList();
	_errors += writeJumpDistanceList();
	if(_errors){return SMLMS::SMLMS_FAILURE;}
	return SMLMS::SMLMS_SUCCESS;
}

int Delegate::runInitPhysModAlgorithm(void){
	unsigned int stateNumber = 0;
	std::cout<<"Number of states: ";
	std::cin>>stateNumber;
	std::cout<<std::endl;
	if(std::cin.fail()){
		std::cout << "User Error: State number needs to be of type positive integer.";
		std::cin.clear();
		return SMLMS::SMLMS_FAILURE;
	}
	/* load functions */
	_errors += loadMicroscope();
	if(_errors){return SMLMS::SMLMS_FAILURE;}
	/* analysis functions */
	try{
		_physMod.setStateNumber(stateNumber);
		_physMod.setMinValue(_eVar.minDistArgument());
		_physMod.setMaxValue(_eVar.maxDistArgument());
		_physMod.setBinSize(_eVar.jumpIntervalArgument());
		_physMod.setMicroscope(_microscope);
		_physMod.initModelByParameter();
		_physMod.initModelBLD();
	}
	catch(SMLMS::SmlmsError& error){
		std::cout<<error.what()<<std::endl;
		_errors += 1;
	}
	catch(...){
		std::cout<<"Oops, the ermine discovered an unexpected error during the initialization of the PDF parameters."<<std::endl;
		_errors += 1;
	}
	if(_errors){return SMLMS::SMLMS_FAILURE;}
	/* write functions */
	_errors += writeMicroscope();
	_errors += writePhysMod();
	if(_errors){return SMLMS::SMLMS_FAILURE;}
	return SMLMS::SMLMS_SUCCESS;
}

int Delegate::runFitPhysModAlgorithm(void){
	/* load functions */
	_errors += loadMicroscope();
	_errors += loadPhysMod();
	_errors += loadJumpDistanceList();
	if(_errors){return SMLMS::SMLMS_FAILURE;}
	/* analysis functions */
	_physMod.calcPdf(_judi);
	_physMod.normalizeModel();
	_physMod.fitPdfSuperPos();
	_physMod.calcFitSuperPosFromPara();
	_physMod.calcFitMatrixFromPara();
	_physMod.calcFitSuperPos();
	_physMod.calcResSuperPos();
	_physMod.calcChiSquareSuperPos();
	_physMod.printParaMat();
	_physMod.printChiSquareSuperPos();
	/* write functions */
	_errors += writeMicroscope();
	_errors += writePhysMod();
	_errors += writeJumpDistanceList();
	if(_errors){return SMLMS::SMLMS_FAILURE;}
	return SMLMS::SMLMS_SUCCESS;
}

int Delegate::runInitHmmWithPhysModAlgorithm(void){
	/* load functions */
	_errors += loadMicroscope();
	_errors += loadPhysMod();
	if(_errors){return SMLMS::SMLMS_FAILURE;}
	/* analysis functions */
	_physMod.calcFitMatrixFromPara();
	_physMod.calcFitSuperPosFromPara();
	_physMod.setPdfSuperPos(_physMod.fitSuperPos());
	_physMod.calcChiSquareSuperPos();
	_hmm.initFromPhysMod(_physMod);
	/* write functions */
	_errors += writeMicroscope();
	_errors += writePhysMod();
	_errors += writeHmm();
	if(_errors){return SMLMS::SMLMS_FAILURE;}
	return SMLMS::SMLMS_SUCCESS;
}


int Delegate::runInitHmmAlgorithm(void){
	/* analysis functions  */
	unsigned int stateNumber = 0;
	std::cout<<std::endl<<"Number of states: ";
	std::cin>>stateNumber;
	std::cout<<std::endl;
	if(std::cin.fail()){
		std::cout << "User Error: State number needs to be of type positive integer.";
		std::cin.clear();
		return SMLMS::SMLMS_FAILURE;
	}
	/* initHMM */
	_hmm.setStateNumber(stateNumber);
	_hmm.setMinValue(_eVar.minDistArgument());
	_hmm.setMaxValue(_eVar.maxDistArgument());
	_hmm.setSymbolInterval(_eVar.jumpIntervalArgument());
	try{
		_hmm.calcSymbolNumber();
		_hmm.initHMM();
		_hmm.calcObsAlphabetFromParas();
	}
	catch(SMLMS::SmlmsError &error){
		std::cout<<error.what()<<std::endl;
	 	_errors += 1;
	}
	catch(std::out_of_range &error){
		std::cout<<error.what()<<std::endl;
		_errors += 1;;
	}
	catch(...){
		std::cout<<"Oops, the ermine discovered an unexpected error during the initialization of the HMM."<<std::endl;
		_errors += 1;
	}
	if(_errors){return SMLMS::SMLMS_FAILURE;}
	/* write functions */
	_errors += writeHmm();
	if(_errors){return SMLMS::SMLMS_FAILURE;}
	return SMLMS::SMLMS_SUCCESS;
}

int Delegate::runEvaluateAlgorithm(void){
	/* load functions */
	_errors += loadJumpDistanceList();
	_errors += loadHmm();
	if(_errors){return SMLMS::SMLMS_FAILURE;}
	/* analysis functions */
	try{
		_hmm.setTraceNumber(_judi.traceNumber());
		_hmm.calcDof();
		_hmm.estimateSeqLikelihood(_judi);
	}
	catch(SMLMS::SmlmsError &error){
		std::cout<<error.what()<<std::endl;
		_errors += 1;
	}
	catch(std::out_of_range &error){
		std::cout<<error.what()<<std::endl;
		_errors += 1;
	}
	catch(...){
		std::cout<<"Oops, the ermine discovered an unexpected error during the evaluation of the HMM."<<std::endl;
		_errors += 1;
	}
	if(_errors){return SMLMS::SMLMS_FAILURE;}
	/* write functions */
	_errors += writeJumpDistanceList();
	_errors += writeHmm();
	if(_errors){return SMLMS::SMLMS_FAILURE;}
	return SMLMS::SMLMS_SUCCESS;
}

int Delegate::runTrainWithPhysModAlgorithm(void){
	/* load functions */
	_errors += loadMicroscope();
	_errors += loadJumpDistanceList();
	_errors += loadPhysMod();
	_errors += loadHmm();
	if(_errors){return SMLMS::SMLMS_FAILURE;}
	/* analysis functions */
	try{
		_physMod.calcPdf(_judi);
		_hmm.trainPhysModSequence(_judi, _physMod);
		_hmm.printHMM();
	}
	catch(SMLMS::SmlmsError &error){
		std::cout<<error.what()<<std::endl;
		_errors += 1;
	}
	catch(std::out_of_range &error){
		std::cout<<error.what()<<std::endl;
		_errors += 1;
	}
	catch(...){
		std::cout<<"Oops, the ermine discovered an unexpected error during the training of the HMM."<<std::endl;
		_errors +=  1;
	}
	if(_errors){return SMLMS::SMLMS_FAILURE;}
	/* write functions */
	_errors += writeMicroscope();
	_errors += writeJumpDistanceList();
	_errors += writePhysMod();
	_errors += writeHmm();
	if(_errors){return SMLMS::SMLMS_FAILURE;}
	return SMLMS::SMLMS_SUCCESS;
}

int Delegate::runTrainAlgorithm(void){
	/* load functions */
	_errors += loadJumpDistanceList();
	_errors += loadHmm();
	if(_errors){return SMLMS::SMLMS_FAILURE;}
	/* analysis functions */
	try {
		_hmm.trainSequence(_judi);
	}
	catch(SMLMS::SmlmsError& error){
		std::cout<<error.what()<<std::endl;
		_errors += 1;
	}
	if(_errors){return SMLMS::SMLMS_FAILURE;}
	/* write functions */
	_errors += writeJumpDistanceList();
	_errors += writeHmm();
	if(_errors){return SMLMS::SMLMS_FAILURE;}
	return SMLMS::SMLMS_SUCCESS;
}

int Delegate::runBestPathAlgorithm(void){
	/* load functions */
	_errors += loadJumpDistanceList();
	_errors += loadHmm();
	if(_errors){return SMLMS::SMLMS_FAILURE;}
	/* analysis functions */
	try{
		_hmm.estimateStateSequence(_judi);
	}
	catch(SMLMS::SmlmsError &error){
		std::cout<<error.what()<<std::endl;
		_errors += 1;
	}
	catch(...){
		std::cout<<"Oops, the ermine discovered an unexpected error during the optimization of the state sequence."<<std::endl;
		_errors +=  1;
	}
	if(_errors){return SMLMS::SMLMS_FAILURE;}
	/* write functions */
	_errors += writeJumpDistanceList();
	_errors += writeHmm();
	if(_errors){return SMLMS::SMLMS_FAILURE;}
	return SMLMS::SMLMS_SUCCESS;
}

int Delegate::runTransferStatesAlgorithm(void){
	/* load functions */
	_errors += loadMoleculeList();
	_errors += loadJumpDistanceList();
	if(_errors){return SMLMS::SMLMS_FAILURE;}
	/* analysis functions */
	try{
		_judi.transferStatesToMoleculeList(_molList);
	}
	catch(SMLMS::SmlmsError& error){
		std::cout<<error.what()<<std::endl;
		_errors += 1;
	}
	catch(...){
		std::cout<<"Oops, the ermine discovered an unexpected error during the transfer of the state sequence."<<std::endl;
		_errors += 1;
	}
	if(_errors){return SMLMS::SMLMS_FAILURE;}
	/* write functions */
	_errors += writeMoleculeList();
	_errors += writeJumpDistanceList();
	if(_errors){return SMLMS::SMLMS_FAILURE;}
	return SMLMS::SMLMS_SUCCESS;
}

int Delegate::runArchiveAlgorithm(void){
	/* load functions */
	_errors += loadMicroscope();
	_errors += loadMoleculeList();
	_errors += loadJumpDistanceList();
	_errors += loadPhysMod();
	_errors += loadHmm();
	if(_errors){return SMLMS::SMLMS_FAILURE;}
	/* analysis functions */
	try{
		SMLMS::HDF5 archive;
		archive.setFolderName(_fileNames.folderName());
		archive.archiveModel(_microscope, _molList, _judi, _hmm, _physMod);
	}
	catch(H5::FileIException& error){
		error.printErrorStack();
		_errors += 1;
	}
	catch(H5::GroupIException& error){
		error.printErrorStack();
		_errors += 1;
	}
	catch(H5::DataSetIException& error){
		error.printErrorStack();
		_errors += 1;;
	}
	if(_errors){return SMLMS::SMLMS_FAILURE;}
	/* write functions */
	_errors += writeMicroscope();
	_errors += writeMoleculeList();
	_errors += writeJumpDistanceList();
	_errors += writePhysMod();
	_errors += writeHmm();
	if(_errors){return SMLMS::SMLMS_FAILURE;}
	return SMLMS::SMLMS_SUCCESS;
}

int Delegate::runExtractAlgorithm(void){
	/* init functions */
	try{
		_physMod.setMinValue(_eVar.minDistArgument());
		_physMod.setMaxValue(_eVar.maxDistArgument());
		_physMod.setBinSize(_eVar.jumpIntervalArgument());
		_physMod.setMicroscope(_microscope);
	}
	catch(SMLMS::SmlmsError& error){
		std::cout<<error.what()<<std::endl;
		_errors += 1;
	}
	catch(...){
		std::cout<<"Oops, the ermine discovered an unexpected error during the initialization of the PDF parameters."<<std::endl;
		_errors += 1;
	}
	if(_errors){return SMLMS::SMLMS_FAILURE;}
	/* analysis functions */
	try{
		SMLMS::HDF5 archive;
		archive.setFileName(_fileNames.archiveName());
		archive.extractModel(_microscope, _molList, _judi, _hmm, _physMod);
	}
	catch(H5::FileIException& error){
		error.printErrorStack();
		_errors +=1;
	}
	catch(H5::GroupIException& error){
		error.printErrorStack();
		_errors += 1;
	}
	catch(H5::DataSetIException& error){
		error.printErrorStack();
		_errors += 1;
	}
	catch(SMLMS::SmlmsError& error){
		std::cout<<error.what()<<std::endl;
		_errors += 1;
	}
	catch(...){
		std::cout<<"Oops, the ermine discovered an unexpected error during archive extraction."<<std::endl;
		_errors += 1;
	}
	if(_errors){return SMLMS::SMLMS_FAILURE;}
	return SMLMS::SMLMS_SUCCESS;
}

int Delegate::runSimulateAlgorithm(void){
	/* load functions */
	_errors += loadHmm();
	if(_errors){return SMLMS::SMLMS_FAILURE;}
	/* analysis functions */
	_hmm.setTraceNumber(_eVar.particleArgument());
	try{
		_hmm.simulateSequence(_eVar.traceLengthArgument(), _judi);
	}
	catch(SMLMS::SmlmsError &error){
		std::cout<<error.what()<<std::endl;
		_errors += 1;
	}
	catch(std::out_of_range &error){
		std::cout<<error.what()<<std::endl;
		_errors += 1;
	}
	catch(...){
		std::cout<<"Oops, the ermine discovered an unknown error!"<<std::endl;
		_errors += 1;
	}
	if(_errors){return SMLMS::SMLMS_FAILURE;}
	/* write functions */
	_errors += writeJumpDistanceList();
	_errors += writeHmm();
	if(_errors){return SMLMS::SMLMS_FAILURE;}
	return SMLMS::SMLMS_SUCCESS;
}

int Delegate::runWholeCellAnalysis(void){
	/* load functions */
	_errors += loadMicroscope();
	_errors += loadMoleculeList();
	_errors += loadJumpDistanceList();
	_errors += loadPhysMod();
	_errors += loadHmm();
	if(_errors){return SMLMS::SMLMS_FAILURE;}
	/* analysis functions */
	/* train */
	try{
		_physMod.calcPdf(_judi);
		_hmm.trainPhysModSequence(_judi, _physMod);
	}
	catch(SMLMS::SmlmsError &error){
		std::cout<<error.what()<<std::endl;
		_errors += 1;
	}
	catch(std::out_of_range &error){
		std::cout<<error.what()<<std::endl;
		_errors += 1;
	}
	catch(...){
		std::cout<<"Oops, the ermine discovered an unexpected error during HMM training."<<std::endl;
		_errors +=  1;
	}
	if(_errors){return SMLMS::SMLMS_FAILURE;}
	/* fix D */
	_physMod.fixDiffusionCoefficients();
	/* fit model with fixed D */
	_physMod.calcPdf(_judi);
	_physMod.normalizeModel();
	_physMod.fitPdfSuperPos();
	_physMod.calcFitSuperPosFromPara();
	_physMod.calcFitMatrixFromPara();
	_physMod.calcFitSuperPos();
	_physMod.calcResSuperPos();
	_physMod.calcChiSquareSuperPos();
	/* loose D */
	_physMod.releaseDiffusionCoefficients();
	//evaluate
	try{
		//_hmm.setTraceNumber(_judi.traceNumber());
		_hmm.calcDof();
		_hmm.estimateSeqLikelihood(_judi);
	}
	catch(SMLMS::SmlmsError &error){
		std::cout<<error.what()<<std::endl;
		_errors += 1;
	}
	catch(std::out_of_range &error){
		std::cout<<error.what()<<std::endl;
		_errors += 1;
	}
	catch(...){
		std::cout<<"Oops, the ermine discovered an unexpected error during HMM evaluation."<<std::endl;
		_errors += 1;
	}
	if(_errors){return SMLMS::SMLMS_FAILURE;}
	/* best Path */
	try{
		_hmm.estimateStateSequence(_judi);
	}
	catch(SMLMS::SmlmsError &error){
		std::cout<<error.what()<<std::endl;
		_errors += 1;
	}
	catch(...){
		std::cout<<"Oops, the ermine discovered an unexpected error during state sequence optimization."<<std::endl;
		_errors +=  1;
	}
	if(_errors){return SMLMS::SMLMS_FAILURE;}

	/* transfer States */
	try{
		_judi.transferStatesToMoleculeList(_molList);
	}
	catch(SMLMS::SmlmsError& error){
		std::cout<<error.what()<<std::endl;
		_errors += 1;
	}
	catch(...){
		std::cout<<"Oops, the ermine discovered an unexpected error during state sequence transfer."<<std::endl;
		_errors += 1;
	}
	if(_errors){return SMLMS::SMLMS_FAILURE;}
	/* archive */
	try{
		SMLMS::HDF5 archive;
		archive.setFolderName(_fileNames.folderName());
		archive.archiveModel(_microscope, _molList, _judi, _hmm, _physMod);
	}
	catch(H5::FileIException& error){
		error.printErrorStack();
		_errors += 1;
	}
	catch(H5::GroupIException& error){
		error.printErrorStack();
		_errors += 1;
	}
	catch(H5::DataSetIException& error){
		error.printErrorStack();
		_errors += 1;;
	}
	if(_errors){return SMLMS::SMLMS_FAILURE;}
	/* write functions */
	_errors += writeMicroscope();
	_errors += writeMoleculeList();
	_errors += writeJumpDistanceList();
	_errors += writePhysMod();
	_errors += writeHmm();
	if(_errors){return SMLMS::SMLMS_FAILURE;}
	return SMLMS::SMLMS_SUCCESS;
}
} /* SMLMS */
