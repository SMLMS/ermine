/* ######################################################################
* File Name: smlmsHmmBase.hpp
* Project: SMLMS
* Version: 16.03
* Creation Date: 22.03.2016
* Created By Sebastian Malkusch
* <malkusch@chemie.uni-frankfurt.de>
* Goethe University of Frankfurt
* Physical and Theoretical Chemistry
* Single Molecule Biophysics
###################################################################### */


#ifndef SmlmsHmmBase_hpp
#define SmlmsHmmBase_hpp

#include <vector>
#include <string>
#include "header/smlmsMatrix.hpp"
#include "header/smlmsPhysModBrownLatDiff.hpp"

namespace SMLMS{
class HMMBase{
	protected:
		unsigned _stateNumber;
		unsigned _symbolNumber;
		unsigned _minValue;
		unsigned _maxValue;
		double _symbolInterval;
		SMLMS::Matrix _equiPDF;
		SMLMS::Matrix _transPDF;
		SMLMS::Matrix _obsPDF;
		SMLMS::Matrix _equiCDF;
		SMLMS::Matrix _transCDF;
		SMLMS::Matrix _obsCDF;
		std::vector<double> _obsAlphabet;
		double _logLikelihood;
		double _bic;
		double _aic;
		double _stopCrit;
		int _maxIt;

	public:
		/* Constructor */
		HMMBase();
		HMMBase(unsigned states, unsigned symbols);	
		/* Destructor */
		virtual ~HMMBase();
		/* copy Constructor */
		HMMBase(const HMMBase &);
		/* elementary functions */
		void setStateNumber(unsigned);
		unsigned stateNumber();
		void setSymbolNumber(unsigned);
		unsigned symbolNumber();
		void setMinValue(unsigned);
		unsigned minValue();
		void setMaxValue(unsigned);
		unsigned maxValue();
		void setSymbolInterval(double);
		double symbolInterval();
		void setEquiPDF(SMLMS::Matrix);
		SMLMS::Matrix equiPDF();
		void setTransPDF(SMLMS::Matrix);
		SMLMS::Matrix transPDF();
		void setObsPDF(SMLMS::Matrix);
		SMLMS::Matrix obsPDF(void);
		void setObsAlphabet(std::vector<double>);
		std::vector<double> obsAlphabet();
		double logLikelihood();
		double bic();
		double aic();
		void setStopCrit(double);
		double stopCrit();
		void setMaxIt(int);
		int maxIt();
        	/* load/save funCtions */
		void readHMM(std::string const &);
		void writeHMM(std::string const &);
		/* clear functions */
		void clearHMM();
		/* proof functions */
		void checkStateNumber();
		void checkSymbolNumber();
		void checkEqui();
		void checkTrans();
		void checkObs();
		void checkAlphabet();
		void checkProbDen(SMLMS::Matrix &);
		void checkHMM();
		void checkPhysicalModelCompatibility(SMLMS::PhysicalModelBLD &);
		/* init functions */
		void initObsAlphabet();
		void initEqui();
		void initEquiPDF();
		void initEquiCDF();
		void initTrans();
		void initTransPDF();
		void initTransCDF();
		void initObs();
		void initObsPDF();
		void initObsCDF();
		void initHMM();
		/* init from file */
		void initLoadedHMM();	
		void initFromPhysMod(SMLMS::PhysicalModelBLD &);
		/* normalization */
		void normalizePDF(SMLMS::Matrix &);
		void normalizeHMM();
		/* calc functions */
		void calcSymbolInterval();
		void calcMinValue();
		void calcMaxValue();
		void calcObsAlphabetFromParas();
		void calcParasFromObsAlphabet();
		void calcCDF(const SMLMS::Matrix& pdf, SMLMS::Matrix& cdf);
		/* print functions */
		void printStateNumber();
		void printSymbolNumber();
		void printEquiPDF();
		void printTransPDF();
		void printObsPDF();
		void printEquiCDF();
		void printTransCDF();
		void printObsCDF();
		void printObsAlphabet();
		void printLogLikelihood();
		void printBic();
		void printAic();
		void printHMM();
};/* HMMBase*/
}/* SMLMS */
#endif /* smlmsHmmBase_hpp */

