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
#include <stdexcept>
#include <boost/program_options.hpp>
#include "header/ermineParser.hpp"

namespace po=boost::program_options;

int main(int argc, char *argv[]){

	std::string ermineHeader( "\nEstimate Reaction-rates by Markov-based Investigation of Nanoscopy Experiments (ermine):\nsingle molecule biophysics\nIPTC\nGoethe University Frankfurt");
	po::options_description parserOpt(ermineHeader.data());
	parserOpt.add_options()
			("help,h", "show this help message.")
			("file,f", po::value<std::string>(), "input file")
			("algorithm,a", po::value<std::string>(), "analysis algorithm (type -h for help)")
			("jumpInterval,j", po::value<int>(), "interval size of jump distances in pdf in [nm]")
			("stopCrit,s", po::value<double>(), "stop criterion for model training")
			("minDist", po::value<int>(), "minimal jump distance to analyze in [nm]")
			("maxDist", po::value<int>(), "maximal jump distance to analyze in [nm]")
			("time,t", po::value<double>(), "time between jump measurements in [s]")
			("duration,d", po::value<int>(), "duration of simulation in [s]")
			("particles,p", po::value<int>(), "number of particles to simulate")
		;
		po::variables_map vm;
	try{
		po::store(po::parse_command_line(argc, argv, parserOpt),vm);
		po::notify(vm);
	}
	catch(...){
		std::cout<<"User Error: Unkown."<<std::endl;
		std::cout<<"Type --help (-h) to view help message."<<std::endl;
		return 1;
	}
	if(vm.count("help")){
		std::cout<<parserOpt<<std::endl;
		return 0;
	}
	if(vm.count("algorithm")<1){
		std::cout<<"Parser Error: No algorithm selected."<<std::endl;
		std::cout<<"Type --help (-h) to view help message."<<std::endl;
		return 1;
	}
	if(vm.count("file")<1){
		std::cout<<"Parser Error: No filename given."<<std::endl;
		std::cout<<"Type --help (-h) to view help message."<<std::endl;
        	return 1;
	}
	else{
		std::cout<<vm["file"].as<std::string>()<<std::endl;
	}
	SMLMS::ErmineParser eVar;
	eVar.parseArguments(vm);
	std::cout<<eVar.fileNameArgument()<<std::endl;
	std::cout<<eVar.algorithmArgument()<<std::endl;
	return 0;
}

