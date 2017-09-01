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


#ifndef SMLMSDWELLTIME_hpp
#define SMLMSDWELLTIME_hpp

#include <vector>
#include <string>
#include "header/smlmsJudi.hpp"
#include "header/smlmsHmmBase.hpp"

namespace SMLMS{
class DwellTimeAnalysis{
	private:
		std::string _folderName;
		unsigned _stateNumber;
		double _dt;
		unsigned _originalState;
		unsigned _targetState;
		std::vector<double> _transTimes;
		std::vector<double> _time;
		std::vector<double> _transProbData;
		std::vector<double> _transProbFit;
		std::vector<double> _transProbRes;
		double _amplitude;
		double _transRate;
		double _transRateError;
		double _chiSquare;
	public:
		/* constructor */
		DwellTimeAnalysis();
		/* destructor */
		~DwellTimeAnalysis();
		/* copy-constructor */
		DwellTimeAnalysis(const DwellTimeAnalysis &);
		/* elementary functions */
		void setFolderName(std::string name);
		std::string folderName();
		void setStateNumber(unsigned);
		unsigned stateNumber();	
		void setDt(double);
		double dt();
		void setOriginalState(unsigned);
		unsigned originalState();
		void setTargetState(unsigned);
		unsigned targetState();
		void setTransTimes(std::vector<double>);
		std::vector<double> transTimes();
		void setTime(std::vector<double>);
		std::vector<double> time();
		void setTransProbData(std::vector<double>);
		std::vector<double> transProbData();
		void setTransProbFit(std::vector<double>);
		std::vector<double> transProbFit();
		void setTransProbRes(std::vector<double>);
		std::vector<double> transProbRes();
		void setAmplitude(double);
		double amplitude();
		void setTransRate(double);
		double transRate();
		void setTransRateError(double);
		double transRateError();
		void setChiSquare(double);
		double chiSquare();
		/* proof functions */
		void proofHistogram();
		/* print functions */
		void printFolderName();
		void printStateNumber();
		void printDt();
		void printOriginalState();
		void printTargetState();
		void printTransTimes();
		void printHist();
		void printAmplitude();
		void printTransRate();
		void printChiSquare();
		void printResult();
		/* write functions */
		void writeDwellTime();
		void plotDwellTime();
		/* help functions */
		bool startIndex(int &start, std::vector<int> &stateList);
		bool stopIndex(int &stop, std::vector<int> &stateList);
		void analyzeTraceTransitions(std::vector<int> &stateList);
		/* special functions */
		void estimateTransTimes(const JumpDistanceList &judi);
		void transTimesToHist();
		void fitHist();
		void estimateTransProbFit();
		void estimateResiduals();
		void estimateChiSquare();
		void analyzeJudi(HMMBase &hmm, const JumpDistanceList &Judi);
}; /* DwellTime */
} /* SMLMS */
#endif /* SMLMSDWELLTIME_hpp */

