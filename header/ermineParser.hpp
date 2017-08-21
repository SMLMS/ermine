/* ######################################################################
* File Name: ermineParser.hpp
* Project: SMLMS
* Version: 1611
* Creation Date: 04.11.2016
* Created By Sebastian Malkusch
* <malkusch@chemie.uni-frankfurt.de>
* Goethe University of Frankfurt
* Physical and Theoretical Chemistry
* Single Molecule Biophysics
###################################################################### */


#ifndef Parser_hpp
#define Parser_hpp
#include <string>
#include <map>
#include <boost/program_options.hpp>
#include "header/smlmsFolder.hpp"

namespace po=boost::program_options;
namespace SMLMS{

class ErmineParser{
	private:
		// ErmineParser arguments
		std::map <std::string, int> _algorithmAlphabet;
		std::string _algorithmArgument; // a
		double _stopCritArgument; // c
		std::string _fileNameArgument; // f
		SMLMS::SMLMSFolder _folderArgument;
		std::string _folderNameArgument;
		double _jumpIntervalArgument; // i
		double _minDistArgument; // m
		double _maxDistArgument; // M
		double  _timeIntervalArgument; // t
		double _durationArgument; // d
		int _traceLengthArgument;
		int _particleArgument;// n
	public:
		// Constructor
		ErmineParser();
		// Destructor
		~ErmineParser();
		// Assessor functions of programs arguments
		void setAlgorithmAlphabet(std::map<std::string, int>);
		std::map<std::string, int> algorithmAlphabet();
		void setAlgorithmArgument(std::string);
		std::string algorithmArgument();
		void setStopCritArgument(double);
		double stopCritArgument();
		void setFileNameArgument(std::string);
		std::string fileNameArgument();
		void setFolderArgument(SMLMS::SMLMSFolder);
		SMLMS::SMLMSFolder folderArgument();
		void setFolderNameArgument(std::string);
		std::string folderNameArgument();
		void setJumpIntervalArgument(double);
		double jumpIntervalArgument();
		void setMinDistArgument(double);
		double minDistArgument();
		void setMaxDistArgument(double);
		double maxDistArgument();
		void setTimeIntervalArgument(double);
		double timeIntervalArgument();
		void setDurationArgument(double);
		double durationArgument();
		void setTraceLengthArgument(int);
		int traceLengthArgument();
		void setParticleArgument(int);
		int particleArgument();
		// Functions of class InputParameter
		void printAlgorithmHelp();
		void proofFilename(po::variables_map &);
		void proofAlgorithm(po::variables_map &);
		void proofAlgorithmArgument();
		void proofStopCrit(po::variables_map &);
		void proofJumpInterval(po::variables_map &);
		void proofMinDist(po::variables_map &);
		void proofMaxDist(po::variables_map &);
		void proofTime(po::variables_map &);
		void proofDuration(po::variables_map &);
		void proofTraceLengthRest(double);
		void proofParticles(po::variables_map &);
		void proofJumpIntervalValidity();
		void calcTraceLength();
		void parseArguments(po::variables_map &);
		void extractFolderName();
		void makeFolder();
		void writeErmineParser();
		void printArguments();
}; // ErmineParser
} // SMLMS
#endif /* Parser_hpp */
