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
			("jumpInterval,j", po::value<int>()->default_value(10), "interval size of jump distances in pdf in [nm] (int)")
			("stopCrit,s", po::value<double>()->default_value(0.01), "stop criterion for model training (float)")
			("minDist", po::value<int>()->default_value(0), "minimal jump distance to analyze in [nm] (int)")
			("maxDist", po::value<int>()->default_value(100), "maximal jump distance to analyze in [nm] (int)")
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
		std::cout<<"Type --help (-h) to view help message."<<std::endl;
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
			microscope.saveMicroscope(fileNames.folderName().append("/microscope.mic"));
			molList.writeMoleculeList(fileNames.folderName().append("/molecule_list.mol"), fileNames.folderName().append("/region_of_interest.roi"));
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
		microscope.saveMicroscope(fileNames.folderName().append("/microscope.mic"));
		molList.writeMoleculeList(fileNames.folderName().append("/molecule_list.mol"), fileNames.folderName().append("/region_of_interest.roi"));
		judi.writeJumpDistanceList(fileNames.folderName().append("/jumpDistance_list.jud"));
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
		physMod.setStateNumber(stateNumber);
		physMod.setMinValue(eVar.minDistArgument());
		physMod.setMaxValue(eVar.maxDistArgument());
		physMod.setIncNumber(eVar.jumpIntervalArgument());
		physMod.setFolderName(fileNames.folderName());
		physMod.initModelByParameter();
		physMod.initModelBLD();
		physMod.updateFixModelParameter(microscope);
		/* save model */
		try{
			physMod.writePhysMod();
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
	/* Initialize HMM */
	else if(eVar.algorithmArgument()=="initHMM"){
		SMLMS::HMMSequence hmm;	
		hmm.setFolderName(fileNames.folderName());
		/* test if model is provided */
		if(fileNames.proofModel()){
			/* init physical model */
			SMLMS::PhysicalModelBLD physMod;
			physMod.setMinValue(eVar.minDistArgument());
			physMod.setMaxValue(eVar.maxDistArgument());
			physMod.setIncNumber(eVar.jumpIntervalArgument());
			physMod.setFolderName(fileNames.folderName());
		try{
			/* loading Model */
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
			std::cout<<"unkown error!"<<std::endl;
			return 1;
		}

			/* calc PDF */
			physMod.calcFitMatrixFromPara();
			physMod.writePhysMod();
			physMod.writePdfMatrix();
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
				std::cout<<"unkown error!"<<std::endl;
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
			std::cout<<"unkown error!"<<std::endl;
			return 1;
		}
		/* start tidy */
		statement.printTidy();
	}
	/* Simulate */
	else if(eVar.algorithmArgument()=="simulate"){
		/* create hmm instance*/
		SMLMS::HMMSequence hmm;	
		hmm.setFolderName(fileNames.folderName()); // muss weg
		/* load hmm */
		try{
			hmm.readHMM(fileNames.hmmName());
			hmm.setTraceNumber(eVar.particleArgument());
			hmm.initLoadedHMM();
			hmm.checkHMM(); // fehler bei der berechnung vom Jump interval
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
		std::cout<<"\nunder construction"<<std::endl;
		// start tidy	
		statement.printTidy();
	}
	// train
	else if(eVar.algorithmArgument()=="train"){
		std::cout<<"\nunder construction"<<std::endl;
		// start tidy	
		statement.printTidy();
	}
	// path
	else if(eVar.algorithmArgument()=="path"){
		std::cout<<"\nunder construction"<<std::endl;
		// start tidy	
		statement.printTidy();
	}
	// no matching algorithm
	else{
		std::cout<<std::endl<<eVar.algorithmArgument()<<" is still under construction"<<std::endl;
	}
	return 0;
}/* main */
