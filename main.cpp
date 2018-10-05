/* ######################################################################
* File Name: main.cpp
* Project: ermine
* Version: 1809
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
#include "header/ermineDelegate.hpp"

namespace po=boost::program_options;

int main(int argc, char *argv[]){
	SMLMS::ErmineParser eVar;
	/* manage programm options */
	std::string ermineHeader( "\nEstimate Reaction-rates by Markov-based Investigation of Nanoscopy Experiments (ermine):\nsingle molecule biophysics\nIPTC\nGoethe University Frankfurt");
	/* parse boost programm options for ermine */
	po::options_description parserOpt(ermineHeader.data());
	parserOpt.add_options()
			("help,h", "show this help message.")
			("file,f", po::value<std::string>(), "input file")
			("algorithm,a", po::value<std::string>(), "analysis algorithm (see below for help)")
			("jumpInterval,j", po::value<double>()->default_value(1.0), "interval size of jump distances in pdf in [nm] (int)")
			("stopCrit,s", po::value<double>()->default_value(0.0001), "stop criterion for model training (float)")
			("maxIt,i", po::value<int>()->default_value(300), "maximum number of iterations for model training (int)")
			("minDist, n", po::value<double>()->default_value(0.0), "minimal jump distance to analyze in [nm] (int)")
			("maxDist, x", po::value<double>()->default_value(800.0), "maximal jump distance to analyze in [nm] (int)")
			("time,t", po::value<double>()->default_value(0.02), "time between jump measurements in [s] (float)")
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
	/* run Program */
	SMLMS::Delegate delegate(vm);
	delegate.run();
	return 0;
}/* main */
