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

namespace po=boost::program_options;
namespace SMLMS{

class ErmineParser{
	private:
		// ErmineParser arguments
		std::map <std::string, int> _algorithmAlphabet;
		std::string _algorithmArgument; // a
		double _stopCritArgument; // c
		std::string _fileNameArgument; // f
		std::string _folderNameArgument;
		int _jumpIntervalArgument; // i
		int _minDistArgument; // m
		int _maxDistArgument; // M
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
		void setFolderNameArgument(std::string);
		std::string folderNameArgument();
		void setJumpIntervalArgument(int);
		int jumpIntervalArgument();
		void setMinDistArgument(int);
		int minDistArgument();
		void setMaxDistArgument(int);
		int maxDistArgument();
		void setTimeIntervalArgument(double);
		double timeIntervalArgument();
		void setDurationArgument(double);
		double durationArgument();
		void setTraceLengthArgument(int);
		int traceLengthArgument();
		void setParticleArgument(int);
		int particleArgument();
		// Functions of class InputParameter
		void printHelp();
		void parseArguments(po::variables_map &);
		int proofArguments(std::string);
		int proofAlgorithmArgument();
		void calcTraceLength();
		void writeErmineParser();
}; // ErmineParser
} // SMLMS
#endif /* Parser_hpp */
