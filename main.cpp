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
#include "header/ermineExceptions.hpp"
#include "header/ermineParser.hpp"


namespace po=boost::program_options;

int main(int argc, char *argv[]){
	SMLMS::ErmineParser eVar;
	std::string ermineHeader( "\nEstimate Reaction-rates by Markov-based Investigation of Nanoscopy Experiments (ermine):\nsingle molecule biophysics\nIPTC\nGoethe University Frankfurt");
	po::options_description parserOpt(ermineHeader.data());
	parserOpt.add_options()
			("help,h", "show this help message.")
			("file,f", po::value<std::string>(), "input file")
			("algorithm,a", po::value<std::string>(), "analysis algorithm (see below for help)")
			("jumpInterval,j", po::value<int>(), "interval size of jump distances in pdf in [nm] (int)")
			("stopCrit,s", po::value<double>(), "stop criterion for model training (float)")
			("minDist", po::value<int>(), "minimal jump distance to analyze in [nm] (int)")
			("maxDist", po::value<int>(), "maximal jump distance to analyze in [nm] (int)")
			("time,t", po::value<double>(), "time between jump measurements in [s] (float)")
			("duration,d", po::value<double>(), "duration of simulation in [s] (float)" )
			("particles,p", po::value<int>(), "number of particles to simulate (int)")
		;
		po::variables_map vm;
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

	try{
		eVar.parseArguments(vm);
	}
 	catch(SMLMS::NoFileName& error){
		std::cout<<error.returnError()<<std::endl;
		return 1;
	}
	catch (SMLMS::NoAlgorithm& error){
		std::cout<<error.returnError()<<std::endl;
		return 1;
	}
	catch (SMLMS::WrongAlgorithm& error){
		std::cout<<error.returnError()<<std::endl;
		eVar.printAlgorithmHelp();
		return 1;
	}
	catch (...){
		std::cout<<"oops, the ermine discovered an unexpected error and is going to rest"<<std::endl;
		return 1;
	}
return 0;
}

