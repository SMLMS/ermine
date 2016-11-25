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
#include <boost/program_options.hpp>

namespace po=boost::program_options;
namespace SMLMS{

class ErmineParser{
	private:
		// ErmineParser indicators
		int _algorithmIndicator;
		int _stopCritIndicator;
		int _filenameIndicator;
		int _jumpIntervalIndicator;
		int _minDistIndicator;
		int _maxDistIndicator;
		int _lengthIndicator;
		int _particleIndicator;
		// ErmineParser arguments
		std::string _algorithmArgument; // a
		double _stopCritArgument; // c
		std::string _fileNameArgument; // f
		std::string _folderNameArgument;
		int _jumpIntervalArgument; // i
		int _minDistArgument; // m
		int _maxDistArgument; // M
		int _lengthArgument; // t
		int _particleArgument;// n
	public:
		// Constructor
		ErmineParser();
		// Destructor
		~ErmineParser();
		// Assessor functions of parser indicatiors
		void setAlgorithmIndicator(int);
		int algorithmIndicator();
		void setStopCritIndicator(int);
		int stopCritIndicator();
		void setFilenameIndicator(int);
		int filenameIndicator();
		void setJumpIntervalIndicator(int);
		int jumpIntervalIndicator();
		void setMinDistIndicator(int);
		int minDistIndicator();
		void setMaxDistIndicator(int);
		int maxDistIndicator();
		void setLengthIndicator(int);
		int lengthIndicator();
		void setParticleIndicator(int);
		int particleIndicator();
		// Assessor functions of programs arguments
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
		void setLengthArgument(int);
		int lengthArgument();
		void setParticleArgument(int);
		int particleArgument();
		// Functions of class InputParameter
		void printHelp();
		void parseArguments(po::variables_map &);
		void proofArguments();
		void writeErmineParser();
}; // ErmineParser
} // SMLMS
#endif /* Parser_hpp */
