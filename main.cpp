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
#include "header/ermineExceptions.hpp"
#include "header/ermineParser.hpp"
#include "header/ermineFilenames.hpp"
#include "header/smlmsMicroscope.hpp"
#include "header/smlmsMolecules.hpp"
#include "header/ermineJudi.hpp"

namespace po=boost::program_options;

int main(int argc, char *argv[]){
	// initialize objects
	SMLMS::Statement statement;
	statement.printStart();
	SMLMS::ErmineParser eVar;
	std::string ermineHeader( "\nEstimate Reaction-rates by Markov-based Investigation of Nanoscopy Experiments (ermine):\nsingle molecule biophysics\nIPTC\nGoethe University Frankfurt");


	// parse boost programm options for ermine
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
	// proof for required inputs
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

	// parse programm options to ErmineParser class
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
	
	// make directory
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
	// write ermine parser
	try{
		eVar.writeErmineParser();
	}
	catch(...){
		std::cout<<"oops, the ermine discovered an unexpected error while writing the parser and is going to rest"<<std::endl;
		return 1;
	}
	
	// read FileNames
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

	// choose algorithm
	// batch algorithm
	if (eVar.algorithmArgument()=="batch"){
		statement.printBatch();
		//load microscope
		SMLMS::Microscope microscope;
		try{
			microscope.loadMicroscope(fileNames.microscopeName());
		}
		catch(SMLMS::SMLMSMicroscopeError& error){
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
		//load file list and roi
		SMLMS::MoleculeList tempList, molList;
		try{
			molList.readROI(fileNames.roiName());
			molList.readTrcList(microscope, fileNames.getTrcName(0));
		}
		catch(SMLMS::SMLMSMoleculesError& error){
			std::cout<<error.what()<<std::endl;
			return 1;
		}
		catch(...){
			std::cout<<"oops, the ermine discovered an unexpected error during argument parsing and is going to rest"<<std::endl;
			return 1;
		}
		// load trc and append files
		for (int i=1; i<fileNames.trcNumber(); i++){
			try{
				tempList.readTrcList(microscope, fileNames.getTrcName(i));
			}
			catch(SMLMS::SMLMSMoleculesError& error){
				std::cout<<error.what()<<std::endl;
				return 1;
			}
			catch(...){
				std::cout<<"oops, the ermine discovered an unexpected error during argument parsing and is going to rest"<<std::endl;
				return 1;
			}
			molList.addMoleculeList(tempList);
		}
		//filter trc files
		molList.filterMoleculeList();
		// write results
		microscope.saveMicroscope(fileNames.folderName().append("/microscope.mic"));
		molList.writeMoleculeList(fileNames.folderName().append("/molecule_list.mol"), fileNames.folderName().append("/region_of_interest.roi"));
		// start tidy	
		statement.printTidy();
	}
	// calculate judi from mol: mol2judi algorithm
	else if(eVar.algorithmArgument()=="mol2judi"){
		statement.printMol2Judi();
		SMLMS::Microscope microscope;
		try{
			microscope.loadMicroscope(fileNames.microscopeName());
		}
		catch(SMLMS::SMLMSMicroscopeError& error){
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
		catch(SMLMS::SMLMSMoleculesError& error){
			std::cout<<error.what()<<std::endl;
			return 1;
		}
		catch(...){
			std::cout<<"oops, the ermine discovered an unexpected error during argument parsing and is going to rest"<<std::endl;
			return 1;
		}
		//filter molecule List
		molList.filterMoleculeList();
		// calculate judi from molecule list
		SMLMS::JumpDistanceList judi;
		judi.calcJumpDistanceList(molList);
		// write results
		microscope.saveMicroscope(fileNames.folderName().append("/microscope.mic"));
		molList.writeMoleculeList(fileNames.folderName().append("/molecule_list.mol"), fileNames.folderName().append("/region_of_interest.roi"));
		judi.writeJumpDistanceList(fileNames.folderName().append("/jumpDistance_list.jud"));
		//start tidy
		statement.printTidy();
	}
	// Initialize HMM
	else if(eVar.algorithmArgument()=="initialize"){
		statement.printInitialize();
		// load microscope
		SMLMS::Microscope microscope;
		microscope.loadMicroscope(fileNames.microscopeName());
		// load judi
		SMLMS::JumpDistanceList judi;
		judi.readJumpDistanceList(fileNames.judiName());
		// read HMM
		/*
		SMLMS::HMM hmm;
		hmm.initObsAlphabet((double) parser.jumpIntervalArgument(), (double) parser.minDistArgument(), (double) parser.maxDistArgument());
		hmm.readHMM(fileNames.hmmName());
		// init hmm
		std::string outfilename= fileNames.sourceFolder();
		hmm.initialize(judi, outfilename);
		// write results
		hmm.writeHMM(fileNames.sourceFolder().append("HMM_result.hmm"));
		*/
		// start tidy
		statement.printTidy();
	}
	else{
		std::cout<<std::endl<<eVar.algorithmArgument()<<" is still under construction"<<std::endl;
	}
	return 0;
}/* main */
