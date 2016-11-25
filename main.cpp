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
			("file,f", po::value<std::string>(), "input file needs to be .txt")
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
	if(vm.count("file")<1){
		std::cout<<"User Error: No filename given."<<std::endl;
		std::cout<<"Type --help (-h) to view help message."<<std::endl;
        	return 1;
	}
	else{
		std::cout<<vm["file"].as<std::string>()<<std::endl;
	}
	SMLMS::ErmineParser eVar;
	eVar.printHelp();
	return 0;
}

