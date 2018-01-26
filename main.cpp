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
#include <iomanip>
#include <boost/program_options.hpp>
#include "header/ermineStatement.hpp"
#include "header/ermineParser.hpp"
#include "header/ermineFilenames.hpp"
#include "header/ermineExceptions.hpp"
#include "header/smlmsExceptions.hpp"
#include "header/smlmsMicroscope.hpp"
#include "header/smlmsMolecules.hpp"
#include "header/smlmsJudi.hpp"
#include "header/smlmsPhysModBrownLatDiff.hpp"
#include "header/smlmsHmmBase.hpp"
#include "header/smlmsHmmUnique.hpp"
#include "header/smlmsHmmSequence.hpp"
#include "header/smlmsDwellTime.hpp"
#include "header/ermineHDF5.hpp"

namespace po=boost::program_options;

int main(int argc, char *argv[]){
	/* initialize objects */
	SMLMS::Statement statement;
	statement.printStart();
	SMLMS::ErmineParser eVar;
	std::string ermineHeader( "\nEstimate Reaction-rates by Markov-based Investigation of Nanoscopy Experiments (ermine):\nsingle molecule biophysics\nIPTC\nGoethe University Frankfurt");


	/* parse boost programm options for ermine */
	po::options_description parserOpt(ermineHeader.data());
	parserOpt.add_options()
			("help,h", "show this help message.")
			("file,f", po::value<std::string>(), "input file")
			("algorithm,a", po::value<std::string>(), "analysis algorithm (see below for help)")
			("jumpInterval,j", po::value<double>()->default_value(10.0), "interval size of jump distances in pdf in [nm] (int)")
			("stopCrit,s", po::value<double>()->default_value(0.01), "stop criterion for model training (float)")
			("maxIt,i", po::value<int>()->default_value(100), "maximum number of iterations for model training (int)")
			("minDist", po::value<double>()->default_value(0.0), "minimal jump distance to analyze in [nm] (int)")
			("maxDist", po::value<double>()->default_value(100.), "maximal jump distance to analyze in [nm] (int)")
			("time,t", po::value<double>()->default_value(0.2), "time between jump measurements in [s] (float)")
			("duration,d", po::value< double>()->default_value(300), "duration of simulation in [s] (float)" )
			("particles,p", po::value<int>()->default_value(1000), "number of particles to simulate (int)")
		;
		po::variables_map vm;
	/* proof for required inputs */
	try{
		po::store(po::parse_command_line(argc, argv, parserOpt),vm);
		po::notify(vm);
	}
	catch(po::unknown_option& error){
		std::cout << error.what() << std::endl;
		return 1;
	}
	catch(...){
		std::cout<<"User Error: Unkown parameters given to ermine."<<std::endl;
		std::cout<<"Tyheader/ermineHDF5.hpppe --help (-h) to view help message."<<std::endl;
		return 1;
	}
	if(vm.count("help")){
		std::cout<<parserOpt<<std::endl;
		eVar.printAlgorithmHelp();
		return 0;
	}

	/* parse programm options to ErmineParser class */
	try{
		eVar.parseArguments(vm);
		eVar.printArguments();
	}
	catch(SMLMS::ErmineParserError& error){
		std::cout<<error.what()<<std::endl;
		return 1;
	}
	catch(SMLMS::SMLMSFolderError& error){
		std::cout<<error.what()<<std::endl;
		return 1;
	}
	catch (...){
		std::cout<<"oops, the ermine discovered an unexpected error during argument parsing and is going to rest"<<std::endl;
		return 1;
	}
	
	/* make directory */
	try{
		eVar.makeFolder();
	}
	catch(SMLMS::SMLMSFolderError& error){
		std::cout<<error.what()<<std::endl;
		return 1;
	}
	catch(...){
		std::cout<<"oops, the ermine discovered an unexpected error during folder creation and is going to rest"<<std::endl;
		return 1;
	}
	/* write ermine parser */
	try{
		eVar.writeErmineParser();
	}
	catch(...){
		std::cout<<"oops, the ermine discovered an unexpected error while writing the parser and is going to rest"<<std::endl;
		return 1;
	}
	
	/* read FileNames */
	SMLMS::FileNames fileNames;
	fileNames.clearFileNames();
	fileNames.setSourceFileName(eVar.fileNameArgument());
	fileNames.setFolderName(eVar.folderNameArgument());	
	try{
		fileNames.readNamesFromSourceFile();
	}
	catch(SMLMS::ErmineFileNameError& error){
		std::cout<<error.what()<<std::endl;
		return 1;
	}
	catch(...){
		std::cout<<"oops, the ermine discovered an unexpected error during argument parsing and is going to rest"<<std::endl;
		return 1;
	}

	/* choose algorithm */
	/* batch algorithm */
	if (eVar.algorithmArgument()=="batch"){
		statement.printBatch();
		/* load microscope */
		SMLMS::Microscope microscope;
		try{
			microscope.loadMicroscope(fileNames.microscopeName());
		}
		catch(SMLMS::SmlmsError& error){
			std::cout<<error.what()<<std::endl;
			return 1;
		}
		catch(...){
			std::cout<<"oops, the ermine discovered an unexpected error during argument parsing and is going to rest"<<std::endl;
			return 1;
		}
		std::cout<<"int time: "<<microscope.intTime()<<std::endl;
		std::cout<<"pxl size: "<<microscope.pxlSize()<<std::endl;
		std::cout<<"loc prec: "<<microscope.locPrec()<<std::endl;
		/* load file list and roi */
		SMLMS::MoleculeList tempList, molList;
		try{
			molList.readROI(fileNames.roiName());
			molList.readTrcList(microscope, fileNames.getTrcName(0));
		}
		catch(SMLMS::SmlmsError& error){
			std::cout<<error.what()<<std::endl;
			return 1;
		}
		catch(...){
			std::cout<<"oops, the ermine discovered an unexpected error during argument parsing and is going to rest"<<std::endl;
			return 1;
		}
		/* load trc and append files */
		for (int i=1; i<fileNames.trcNumber(); i++){
			try{
				tempList.readTrcList(microscope, fileNames.getTrcName(i));
			}
			catch(SMLMS::SmlmsError& error){
				std::cout<<error.what()<<std::endl;
				return 1;
			}
			catch(...){
				std::cout<<"oops, the ermine discovered an unexpected error during argument parsing and is going to rest"<<std::endl;
				return 1;
			}
			molList.addMoleculeList(tempList);
		}
		/* filter trc files */
		molList.filterMoleculeList();
		/* write results */
		try{
			microscope.saveMicroscope(fileNames.folderName().append("/microscope.txt"));
			molList.writeMoleculeList(fileNames.folderName().append("/mol.txt"), fileNames.folderName().append("/roi.txt"));
		}
		catch(SMLMS::SmlmsError& error){
			std::cout<<error.what()<<std::endl;
			return 1;
		}
		catch(...){
			std::cout<<"oops, the ermine discovered an unexpected error during argument parsing and is going to rest"<<std::endl;
			return 1;
		}
		/* start tidy */
		statement.printTidy();
	}
	/* calculate judi from mol: mol2judi algorithm */
	else if(eVar.algorithmArgument()=="mol2judi"){
		statement.printMol2Judi();
		SMLMS::Microscope microscope;
		try{
			microscope.loadMicroscope(fileNames.microscopeName());
		}
		catch(SMLMS::SmlmsError& error){
			std::cout<<error.what()<<std::endl;
			return 1;
		}
		catch(...){
			std::cout<<"oops, the ermine discovered an unexpected error during argument parsing and is going to rest"<<std::endl;
			return 1;
		}
		SMLMS::MoleculeList molList;
		try{
			molList.readMoleculeList(fileNames.molListName(), fileNames.roiName());
		}
		catch(SMLMS::SmlmsError& error){
			std::cout<<error.what()<<std::endl;
			return 1;
		}
		catch(...){
			std::cout<<"oops, the ermine discovered an unexpected error during argument parsing and is going to rest"<<std::endl;
			return 1;
		}
		/* filter molecule List */
		molList.filterMoleculeList();
		/* calculate judi from molecule list */
		SMLMS::JumpDistanceList judi;
		judi.calcJumpDistanceList(molList);
		/* write results */
		try{
			microscope.saveMicroscope(fileNames.folderName().append("/microscope.txt"));
			molList.writeMoleculeList(fileNames.folderName().append("/mol.txt"), fileNames.folderName().append("/roi.txt"));
			judi.writeJumpDistanceList(fileNames.folderName().append("/judi.txt"));
		}
		catch(SMLMS::SmlmsError& error){
			std::cout<<error.what()<<std::endl;
			return 1;
		}
		catch(...){
			std::cout<<"oops, the ermine discovered an unexpected error during argument parsing and is going to rest"<<std::endl;
			return 1;
		}
		/* start tidy */
		statement.printTidy();
	}
	/* Initialize Physical Model */
	else if(eVar.algorithmArgument()=="initPhysMod"){
		statement.printInitPhysMod();
		/* get state Number */
		std::cout<<std::endl<<"Number of states:"<<std::endl;
		int stateNumber = 0;
		std::cin>>stateNumber;
		if(std::cin.fail()){
			std::cout << "User Error: State number needs to be of type integer.";
			std::cin.clear();
			return 1;
		}
		/* Loading microscope */
		SMLMS::Microscope microscope;
		try{
			microscope.loadMicroscope(fileNames.microscopeName());
		}
		catch(SMLMS::SmlmsError& error){
			std::cout<<error.what()<<std::endl;
			return 1;
		}
		catch(...){
			std::cout<<"oops, the ermine discovered an unexpected error during argument parsing and is going to rest"<<std::endl;
			return 1;
		}
		/* init PhysMod */
		SMLMS::PhysicalModelBLD physMod;
		try{
			physMod.setStateNumber(stateNumber);
			physMod.setMinValue(eVar.minDistArgument());
			physMod.setMaxValue(eVar.maxDistArgument());
			physMod.setBinSize(eVar.jumpIntervalArgument());
			physMod.setFolderName(fileNames.folderName());
			physMod.setMicroscope(microscope);
			physMod.initModelByParameter();
			physMod.initModelBLD();
		}
		catch(SMLMS::SmlmsError& error){
			std::cout<<error.what()<<std::endl;
			return 1;
		}
		catch(...){
			std::cout<<"oops, the ermine discovered an unexpected error during argument parsing and is going to rest"<<std::endl;
			return 1;
		}
		/* save model */
		physMod.printMinValue();
		physMod.printMaxValue();
		physMod.printBinSize();
		physMod.printIncNumber();
		physMod.printAlphabet();
		try{
			physMod.writePhysMod();
			microscope.saveMicroscope(fileNames.folderName().append("/microscope.txt"));
		}
		catch(SMLMS::SmlmsError& error){
			std::cout<<error.what()<<std::endl;
			return 1;
		}
		catch(...){
			std::cout<<"oops, the ermine discovered an unexpected error during argument parsing and is going to rest"<<std::endl;
			return 1;
		}
		/* start tidy */
		statement.printTidy();
	}
	/* fit physMod to judi */
	else if(eVar.algorithmArgument()=="fitPhysMod"){
		statement.printFitPhysMod();
		/* init physical model */
		SMLMS::PhysicalModelBLD physMod;
		physMod.setMinValue(eVar.minDistArgument());
		physMod.setMaxValue(eVar.maxDistArgument());
		physMod.setBinSize(eVar.jumpIntervalArgument());
		physMod.setFolderName(fileNames.folderName());
		/* load microscope */
		SMLMS::Microscope microscope;
		try{
			microscope.loadMicroscope(fileNames.microscopeName());
		}
		catch(SMLMS::SmlmsError& error){
			std::cout<<error.what()<<std::endl;
			return 1;
		}
		catch(...){
			std::cout<<"oops, the ermine discovered an unexpected error during argument parsing and is going to rest"<<std::endl;
			return 1;
		}
		physMod.setMicroscope(microscope);
		try{
			/* loading Model */
			std::cout<<"reading physical model"<<std::endl;
			physMod.readPhysMod(fileNames.modelName());
		}
		catch(SMLMS::SmlmsError &error){
			std::cout<<error.what()<<std::endl;
			return 1;
		}
		catch(...){
			std::cout<<"oops, the ermine discovered an unexpected error while loading the physical model and is going to rest."<<std::endl;
			return 1;
		}
		/* load judi */
		SMLMS::JumpDistanceList judi;
		try{
			std::cout<<"reading judi."<<std::endl;
			judi.readJumpDistanceList(fileNames.judiName());
		}
		catch(SMLMS::SmlmsError &error){
			std::cout<<error.what()<<std::endl;
			return 1;
		}
		catch(...){
			std::cout<<"oops, the ermine discovered an unexpected error while loading judi data and is going to rest."<<std::endl;
			return 1;
		}
		physMod.printMinValue();
		physMod.printMaxValue();
		physMod.printBinSize();
		physMod.printIncNumber();
		physMod.printAlphabet();
		/* fit physMod to judi*/
		physMod.calcPdf(judi);
		physMod.normalizeModel();
		physMod.fitPdfSuperPos();
		physMod.calcFitSuperPosFromPara();
		physMod.calcFitMatrixFromPara();
		physMod.calcFitSuperPos();
		physMod.calcResSuperPos();
		physMod.calcChiSquareSuperPos();
		physMod.printParaMat();
		physMod.printChiSquareSuperPos();
		/* save model */
		try{
			physMod.writePhysMod();
			physMod.writePdfSuperPos();
			physMod.plotPhysicalModel();
		}
		catch(SMLMS::SmlmsError& error){
			std::cout<<error.what()<<std::endl;
			return 1;
		}
		catch(...){
			std::cout<<"oops, the ermine discovered an unexpected error while writing the physical model and is going to rest"<<std::endl;
			return 1;
		}
		// start tidy	
		statement.printTidy();
	}
	/* Initialize HMM */
	else if(eVar.algorithmArgument()=="initHMM"){
		statement.printInitHMM();
		SMLMS::HMMSequence hmm;	
		hmm.setFolderName(fileNames.folderName());
		/* test if model is provided */
		if(fileNames.proofModel()){
			/* init physical model */
			SMLMS::PhysicalModelBLD physMod;
			physMod.setMinValue(eVar.minDistArgument());
			physMod.setMaxValue(eVar.maxDistArgument());
			physMod.setBinSize(eVar.jumpIntervalArgument());
			physMod.setFolderName(fileNames.folderName());
			/* load microscope */
			SMLMS::Microscope microscope;
			try{
				microscope.loadMicroscope(fileNames.microscopeName());
			}
			catch(SMLMS::SmlmsError& error){
				std::cout<<error.what()<<std::endl;
				return 1;
			}
			catch(...){
				std::cout<<"oops, the ermine discovered an unexpected error during argument parsing and is going to rest"<<std::endl;
				return 1;
			}
			physMod.setMicroscope(microscope);

			/* loading Model */
			try{
				std::cout<<"reading physical model"<<std::endl;
				physMod.readPhysMod(fileNames.modelName());
			}
			catch(SMLMS::SmlmsError &error){
				std::cout<<error.what()<<std::endl;
				return 1;
			}
			catch(std::out_of_range &error){
				std::cout<<error.what()<<std::endl;
				return 1;
			}
			catch(...){
				std::cout<<"oops, the ermine discovered an unexpected error while loading the physical model and is going to rest."<<std::endl;
				return 1;
			}
			/* calc PDF */
			physMod.calcFitMatrixFromPara();
			physMod.calcFitSuperPosFromPara();
			physMod.setPdfSuperPos(physMod.fitSuperPos());
			physMod.calcChiSquareSuperPos();
			/* write model */
			try{
				physMod.writePhysMod();
				physMod.writePdfSuperPos();
				physMod.plotPhysicalModel();
			}
			catch(SMLMS::SmlmsError& error){
				std::cout<<error.what()<<std::endl;
				return 1;
			}
			catch(...){
				std::cout<<"oops, the ermine discovered an unexpected error while writing the physical model and is going to rest"<<std::endl;
				return 1;
			}
			/* init hmm */
			hmm.initFromPhysMod(physMod);
		}
		/* else if no model provided */
		else{
			/* get state Number */
			std::cout<<std::endl<<"Number of states:"<<std::endl;
			int stateNumber = 0;
			std::cin>>stateNumber;
			if(std::cin.fail()){
				std::cout << "User Error: State number needs to be of type integer.";
				std::cin.clear();
				return 1;
			}
			/* initHMM */
			hmm.setStateNumber(stateNumber);
			hmm.setMinValue(eVar.minDistArgument());
			hmm.setMaxValue(eVar.maxDistArgument());
			hmm.setSymbolInterval(eVar.jumpIntervalArgument());
			try{
				hmm.calcSymbolNumber();
				hmm.initHMM(); //calcSymbolNumber, calcObsAlphabet to init HMM
				hmm.calcObsAlphabetFromParas();
			}
			catch(SMLMS::SmlmsError &error){
				std::cout<<error.what()<<std::endl;
				return 1;
			}
			catch(std::out_of_range &error){
				std::cout<<error.what()<<std::endl;
				return 1;
			}
			catch(...){
				std::cout<<"oops, the ermine discovered an unexpected error while initializing the HMM and is going to rest."<<std::endl;
				return 1;
			}
		}
		/* save model */
		try{
			hmm.printHMM();
			hmm.writeHMM();
		}
		catch(SMLMS::SmlmsError &error){
			std::cout<<error.what()<<std::endl;
			return 1;
		}
		catch(std::out_of_range &error){
			std::cout<<error.what()<<std::endl;
			return 1;
		}
		catch(...){
			std::cout<<"oops, the ermine discovered an unexpected error while saving the HMM and is going to rest."<<std::endl;
			return 1;
		}
		/* start tidy */
		statement.printTidy();
	}
	/* Simulate */
	else if(eVar.algorithmArgument()=="simulate"){
		statement.printSimulate();
		/* create hmm instance*/
		SMLMS::HMMSequence hmm;	
		/* load hmm */
		try{
			hmm.readHMM(fileNames.hmmName());
			hmm.setTraceNumber(eVar.particleArgument());
			hmm.initLoadedHMM();
			hmm.checkHMM();
			hmm.printHMM();
		}
		catch(SMLMS::SmlmsError &error){
			std::cout<<error.what()<<std::endl;
			return 1;
		}
		catch(std::out_of_range &error){
			std::cout<<error.what()<<std::endl;
			return 1;
		}
		catch(...){
			std::cout<<"unkown error!"<<std::endl;
			return 1;
		}
		/* create judi instance */
		SMLMS::JumpDistanceList judi;
		/* simulate */
		try{
			hmm.simulateSequence(eVar.traceLengthArgument(), judi);
		}
		catch(SMLMS::SmlmsError &error){
			std::cout<<error.what()<<std::endl;
			return 1;
		}
		catch(std::out_of_range &error){
			std::cout<<error.what()<<std::endl;
			return 1;
		}
		catch(...){
			std::cout<<"unkown error!"<<std::endl;
			return 1;
		}
		try{
			judi.writeJumpDistanceList(fileNames.folderName().append("/judi.txt"));
		}
		catch(SMLMS::SmlmsError &error){
			std::cout<<error.what()<<std::endl;
			return 1;
		}
		catch(std::out_of_range &error){
			std::cout<<error.what()<<std::endl;
			return 1;
		}
		catch(...){
			std::cout<<"unkown error!"<<std::endl;
			return 1;
		}
		/* start tidy */
		statement.printTidy();
	}
	// evaluate
	else if(eVar.algorithmArgument()=="evaluate"){
		statement.printEvaluate();
		/* create hmm instance*/
		SMLMS::HMMSequence hmm;	
		/* load hmm */
		try{
			hmm.readHMM(fileNames.hmmName());
			hmm.setFolderName(eVar.folderNameArgument());
			hmm.initLoadedHMM();
			hmm.checkHMM();
			hmm.printHMM();
		}
		catch(SMLMS::SmlmsError &error){
			std::cout<<error.what()<<std::endl;
			return 1;
		}
		catch(std::out_of_range &error){
			std::cout<<error.what()<<std::endl;
			return 1;
		}
		catch(...){
			std::cout<<"unkown error!"<<std::endl;
			return 1;
		}
		/* create judi instance */
		SMLMS::JumpDistanceList judi;
		/* load judi */
		try{
			judi.readJumpDistanceList(fileNames.judiName());
		}
		catch(SMLMS::SmlmsError &error){
			std::cout<<error.what()<<std::endl;
			return 1;
		}
		catch(std::out_of_range &error){
			std::cout<<error.what()<<std::endl;
			return 1;
		}
		catch(...){
			std::cout<<"unkown error!"<<std::endl;
			return 1;
		}
		/* estimate hmm */
		try{
			std::cout<<"The ermine is evaluating the likelihood of model fitting the given sequence."<<std::endl;
			hmm.setTraceNumber(judi.traceNumber());
			hmm.calcDof();
			hmm.estimateSeqLikelihood(judi);
		}
		catch(SMLMS::SmlmsError &error){
			std::cout<<error.what()<<std::endl;
			return 1;
		}
		catch(std::out_of_range &error){
			std::cout<<error.what()<<std::endl;
			return 1;
		}
		catch(...){
			std::cout<<"unkown error!"<<std::endl;
			return 1;
		}
		/* write hmm */
		try{
			hmm.printHMM();
			hmm.writeHMM();
		}
		catch(SMLMS::SmlmsError &error){
			std::cout<<error.what()<<std::endl;
			return 1;
		}
		catch(std::out_of_range &error){
			std::cout<<error.what()<<std::endl;
			return 1;
		}
		catch(...){
			std::cout<<"oops, the ermine discovered an unexpected error while saving the HMM and is going to rest."<<std::endl;
			return 1;
		}
		/* start tidy */
		statement.printTidy();
	}
	// train
	else if(eVar.algorithmArgument()=="train"){
		statement.printTrain();
		/* create hmm instance*/
		SMLMS::HMMSequence hmm;	
		/* load hmm */
		try{
			hmm.readHMM(fileNames.hmmName());
			hmm.setFolderName(eVar.folderNameArgument());
			hmm.setStopCrit(eVar.stopCritArgument());
			hmm.setMaxIt(eVar.maxItArgument());
			hmm.initLoadedHMM();
			hmm.checkHMM();
			hmm.normalizeHMM();
			hmm.printHMM();
		}
		catch(SMLMS::SmlmsError &error){
			std::cout<<error.what()<<std::endl;
			return 1;
		}
		catch(std::out_of_range &error){
			std::cout<<error.what()<<std::endl;
			return 1;
		}
		catch(...){
			std::cout<<"unkown error!"<<std::endl;
			return 1;
		}
		/* create judi instance */
		SMLMS::JumpDistanceList judi;
		/* load judi */
		try{
			judi.readJumpDistanceList(fileNames.judiName());
		}
		catch(SMLMS::SmlmsError &error){
			std::cout<<error.what()<<std::endl;
			return 1;
		}
		catch(std::out_of_range &error){
			std::cout<<error.what()<<std::endl;
			return 1;
		}
		catch(...){
			std::cout<<"unkown error!"<<std::endl;
			return 1;
		}
		/* test if model is provided */
		if(fileNames.proofModel()){
			std::cout<<"moin"<<std::endl;
			/* init physical model */
			SMLMS::PhysicalModelBLD physMod;
			physMod.setMinValue(eVar.minDistArgument());
			physMod.setMaxValue(eVar.maxDistArgument());
			physMod.setBinSize(eVar.jumpIntervalArgument());
			physMod.setFolderName(fileNames.folderName());
			/* load microscope */
			SMLMS::Microscope microscope;
			try{
				microscope.loadMicroscope(fileNames.microscopeName());
			}
			catch(SMLMS::SmlmsError& error){
				std::cout<<error.what()<<std::endl;
				return 1;
			}
			catch(...){
				std::cout<<"oops, the ermine discovered an unexpected error during argument parsing and is going to rest"<<std::endl;
				return 1;
			}
			physMod.setMicroscope(microscope);

			/* loading Model */
			try{
				std::cout<<"reading physical model"<<std::endl;
				physMod.readPhysMod(fileNames.modelName());
			}
			catch(SMLMS::SmlmsError &error){
				std::cout<<error.what()<<std::endl;
				return 1;
			}
			catch(std::out_of_range &error){
				std::cout<<error.what()<<std::endl;
				return 1;
			}
			catch(...){
				std::cout<<"oops, the ermine discovered an unexpected error while loading the physical model and is going to rest."<<std::endl;
				return 1;
			}
			try{
				physMod.calcPdf(judi);
				hmm.trainPhysModSequence(judi, physMod);
			}
			catch(SMLMS::SmlmsError &error){
				std::cout<<error.what()<<std::endl;
				return 1;
			}
			catch(std::out_of_range &error){
				std::cout<<error.what()<<std::endl;
				return 1;
			}
			catch(...){
				std::cout<<"oops, the ermine discovered an unexpected error while loading the physical model and is going to rest."<<std::endl;
				return 1;
			}
			/* write model */
			try{
				physMod.writePhysMod();
				physMod.writePdfSuperPos();
				physMod.plotPhysicalModel();
			}
			catch(SMLMS::SmlmsError& error){
				std::cout<<error.what()<<std::endl;
				return 1;
			}
			catch(...){
				std::cout<<"oops, the ermine discovered an unexpected error while writing the physical model and is going to rest"<<std::endl;
				return 1;
			}
			
		}
		/* train hmm without physical model */
		else{
			try {
				hmm.trainSequence(judi);
			}
			catch(SMLMS::SmlmsError& error){
				std::cout<<error.what()<<std::endl;
				return 1;
			}
		}
		/* write hmm */
		try{
			hmm.printHMM();
			hmm.writeHMM();
		}
		catch(SMLMS::SmlmsError &error){
			std::cout<<error.what()<<std::endl;
			return 1;
		}
		catch(std::out_of_range &error){
			std::cout<<error.what()<<std::endl;
			return 1;
		}
		catch(...){
			std::cout<<"oops, the ermine discovered an unexpected error while saving the HMM and is going to rest."<<std::endl;
			return 1;
		}
		/* start tidy */
		statement.printTidy();
	}
	// calculate best path
	else if(eVar.algorithmArgument()=="bestPath"){
		statement.printBestPath();
		/* create hmm instance*/
		SMLMS::HMMSequence hmm;	
		/* load hmm */
		try{
			hmm.readHMM(fileNames.hmmName());
			hmm.setFolderName(eVar.folderNameArgument());
			hmm.setStopCrit(eVar.stopCritArgument());
			hmm.setMaxIt(eVar.maxItArgument());
			hmm.initLoadedHMM();
			hmm.checkHMM();
			hmm.normalizeHMM();
			hmm.printHMM();
		}
		catch(SMLMS::SmlmsError &error){
			std::cout<<error.what()<<std::endl;
			return 1;
		}
		catch(std::out_of_range &error){
			std::cout<<error.what()<<std::endl;
			return 1;
		}
		catch(...){
			std::cout<<"unkown error!"<<std::endl;
			return 1;
		}
		/* create judi instance */
		SMLMS::JumpDistanceList judi;
		/* load judi */
		try{
			judi.readJumpDistanceList(fileNames.judiName());
		}
		catch(SMLMS::SmlmsError &error){
			std::cout<<error.what()<<std::endl;
			return 1;
		}
		catch(std::out_of_range &error){
			std::cout<<error.what()<<std::endl;
			return 1;
		}
		catch(...){
			std::cout<<"unkown error!"<<std::endl;
			return 1;
		}
		/* Viterbi */
		hmm.estimateStateSequence(judi);
		/* write result */
		try{
			hmm.writeHMM();
			judi.writeJumpDistanceList(fileNames.folderName().append("/judi.txt"));
		}
		catch(SMLMS::SmlmsError &error){
			std::cout<<error.what()<<std::endl;
			return 1;
		}
		catch(std::out_of_range &error){
			std::cout<<error.what()<<std::endl;
			return 1;
		}
		catch(...){
			std::cout<<"oops, the ermine discovered an unexpected error while saving the HMM and is going to rest."<<std::endl;
			return 1;
		}
		/* start tidy */
		statement.printTidy();
	}
	else if(eVar.algorithmArgument()=="dwellTime"){
		statement.printDwellTime();
		/* create HMM instance*/
		SMLMS::HMMBase hmm;	
		/* load hmm */
		try{
			hmm.readHMM(fileNames.hmmName());
			hmm.setFolderName(eVar.folderNameArgument());
			hmm.setStopCrit(eVar.stopCritArgument());
			hmm.setMaxIt(eVar.maxItArgument());
			hmm.initLoadedHMM();
			hmm.checkHMM();
			hmm.normalizeHMM();
		}
		catch(SMLMS::SmlmsError &error){
			std::cout<<error.what()<<std::endl;
			return 1;
		}
		catch(std::out_of_range &error){
			std::cout<<error.what()<<std::endl;
			return 1;
		}
		catch(...){
			std::cout<<"unkown error!"<<std::endl;
			return 1;
		}
		/* create judi instance */
		SMLMS::JumpDistanceList judi;
		/* load judi */
		try{
			judi.readJumpDistanceList(fileNames.judiName());
		}
		catch(SMLMS::SmlmsError &error){
			std::cout<<error.what()<<std::endl;
			return 1;
		}
		catch(std::out_of_range &error){
			std::cout<<error.what()<<std::endl;
			return 1;
		}
		catch(...){
			std::cout<<"unkown error!"<<std::endl;
			return 1;
		}

		/* create dwell time instance */
		SMLMS::DwellTimeAnalysis dwellTime;
		dwellTime.setFolderName(eVar.folderNameArgument());
		dwellTime.setDt(eVar.jumpIntervalArgument());
		try{
			dwellTime.analyzeJudi(hmm, judi);
		}
		catch(SMLMS::SmlmsError &error) {std::cout<<error.what()<<std::endl;}
		catch(std::out_of_range &error) {std::cout<<error.what()<<std::endl;}
		catch(...) {std::cout<<"unkown error!"<<std::endl;}
		/* start tidy */
		statement.printTidy();
	}
	else if(eVar.algorithmArgument()=="transferStates"){
		statement.printTransferStates();
		/* create mol list */
		SMLMS::MoleculeList tempList, molList;
		/* read mol list */
		try{
			std::cout<<"reading molecule list."<<std::endl;
			molList.readLocList(fileNames.getTrcName(0));
		}
		catch(SMLMS::SmlmsError& error){
			std::cout<<error.what()<<std::endl;
			return 1;
		}
		catch(...){
			std::cout<<"oops, the ermine discovered an unexpected error during argument parsing and is going to rest"<<std::endl;
			return 1;
		}
		/* create judi */
		SMLMS::JumpDistanceList judi;
		/* read judi */
		try{
			std::cout<<"reading judi."<<std::endl;
			judi.readJumpDistanceList(fileNames.judiName());
		}
		catch(SMLMS::SmlmsError &error){
			std::cout<<error.what()<<std::endl;
			return 1;
		}
		catch(...){
			std::cout<<"oops, the ermine discovered an unexpected error while loading judi data and is going to rest."<<std::endl;
			return 1;
		}
		/* transfer states*/
		try{
			judi.transferStatesToMoleculeList(molList);
		}
		catch(SMLMS::SmlmsError& error){
			std::cout<<error.what()<<std::endl;
			return 1;
		}
		catch(...){
			std::cout<<"oops, the ermine discovered an unexpected error during argument parsing and is going to rest"<<std::endl;
			return 1;
		}
		/* write mol list */
		try{
			molList.writeLocList(fileNames.folderName().append("/mol.txt"));
		}
		catch(SMLMS::SmlmsError& error){
			std::cout<<error.what()<<std::endl;
			return 1;
		}
		catch(...){
			std::cout<<"oops, the ermine discovered an unexpected error during argument parsing and is going to rest"<<std::endl;
			return 1;
		}
		/* start tidy */
		statement.printTidy();
	}
	/* archive Model */
	else if(eVar.algorithmArgument()=="archive"){
		statement.printArchive();
		/* create instances */
		SMLMS::Microscope microscope;
		SMLMS::MoleculeList molList;
		SMLMS::JumpDistanceList judi;
		SMLMS::HMMSequence hmm;
		SMLMS::PhysicalModelBLD physMod;
		physMod.setMinValue(eVar.minDistArgument());
		physMod.setMaxValue(eVar.maxDistArgument());
		physMod.setBinSize(eVar.jumpIntervalArgument());
		physMod.setFolderName(fileNames.folderName());
		/* load data */
		try{
			microscope.loadMicroscope(fileNames.microscopeName());
			physMod.setMicroscope(microscope);
			molList.readMoleculeList(fileNames.molListName(), fileNames.roiName());
			judi.readJumpDistanceList(fileNames.judiName());
			hmm.readHMM(fileNames.hmmName());
			physMod.readPhysMod(fileNames.modelName());
		}
		catch(SMLMS::SmlmsError& error){
			std::cout<<error.what()<<std::endl;
			return 1;
		}
		catch(...){
			std::cout<<"oops, the ermine discovered an unexpected error during argument parsing and is going to rest"<<std::endl;
			return 1;
		}
		
		try{
			SMLMS::HDF5 archive;
			archive.setFileName(fileNames.archiveName());
			archive.archiveModel(microscope, molList, judi, hmm, physMod);
		}
		catch(H5::FileIException& error){
			error.printError();
			return -1;
		}
		catch(H5::GroupIException& error){
			error.printError();
			return -1;
		}
		catch(H5::DataSetIException& error){
			error.printError();
			return -1;
		}
		statement.printTidy();
	}
	/* archive Model */
	else if(eVar.algorithmArgument()=="extract"){
		statement.printExtract();
		/* create instances */
		SMLMS::Microscope microscope;
		SMLMS::MoleculeList molList;
		SMLMS::JumpDistanceList judi;
		SMLMS::HMMSequence hmm;
		SMLMS::PhysicalModelBLD physMod;
		try{
			physMod.setMinValue(eVar.minDistArgument());
			physMod.setMaxValue(eVar.maxDistArgument());
			physMod.setBinSize(eVar.jumpIntervalArgument());
			physMod.setFolderName(fileNames.folderName());
			physMod.setMicroscope(microscope);
		}
		catch(SMLMS::SmlmsError& error){
			std::cout<<error.what()<<std::endl;
			return 1;
		}
		catch(...){
			std::cout<<"oops, the ermine discovered an unexpected error during argument parsing and is going to rest"<<std::endl;
			return 1;
		}

		try{
			SMLMS::HDF5 archive;
			archive.setFileName(fileNames.archiveName());
			archive.extractModel(microscope, molList, judi, hmm, physMod);
		}
		catch(H5::FileIException& error){
			error.printError();
			return -1;
		}
		catch(H5::GroupIException& error){
			error.printError();
			return -1;
		}
		catch(H5::DataSetIException& error){
			error.printError();
			return -1;
		}
		catch(SMLMS::SmlmsError& error){
			std::cout<<error.what()<<std::endl;
			return 1;
		}
		catch(...){
			std::cout<<"oops, the ermine discovered an unexpected error during argument parsing and is going to rest"<<std::endl;
			return 1;
		}

		statement.printTidy();
	}
	/* no matching algorithm */
	else{
		std::cout<<std::endl<<eVar.algorithmArgument()<<" is not a valid argument for ermine"<<std::endl;
		std::cout<<"type -h for help"<<std::endl;
	}
	statement.printFinished();
	return 0;
}/* main */
