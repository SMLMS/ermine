/* ######################################################################
* File Name: smlmsHmmBase.hpp
* Project: SMLMS
* Version: 18.09
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
#include <boost/mpi.hpp>
#include <boost/mpi/collectives.hpp>
#include "header/smlmsMatrix.hpp"
#include "header/smlmsPhysModBrownLatDiff.hpp"

namespace SMLMS{
class HMMBase{
	protected:
		friend class boost::serialization::access;
		template<class Archive>
		void serialize(Archive & ar, const unsigned int version){
			ar & _stateNumber;
			ar & _symbolNumber;
			ar & _minValue;
			ar & _maxValue;
			ar & _symbolInterval;
			ar & _equiPDF;
			ar & _transPDF;
			ar & _obsPDF;
			ar & _equiCDF;
			ar & _transCDF;
			ar & _obsCDF;
			ar & _obsAlphabet;
			ar & _logLikelihood;
			ar & _dof;
			ar & _bic;
			ar & _aic;
			ar & _stopCrit;
			ar & _maxIt;
		}
		unsigned _stateNumber;
		unsigned _symbolNumber;
		double _minValue;
		double _maxValue;
		double _symbolInterval;
		SMLMS::Matrix _equiPDF;
		SMLMS::Matrix _transPDF;
		SMLMS::Matrix _obsPDF;
		SMLMS::Matrix _equiCDF;
		SMLMS::Matrix _transCDF;
		SMLMS::Matrix _obsCDF;
		std::vector<double> _obsAlphabet;
		double _logLikelihood;
		unsigned _dof;
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
		void setMinValue(double);
		double minValue();
		void setMaxValue(double);
		double maxValue();
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
		void setLogLikelihood(double);
		double logLikelihood();
		void setDof(unsigned);
		unsigned dof();
		void setBic(double);
		double bic();
		void setAic(double);
		double aic();
		void setStopCrit(double);
		double stopCrit();
		void setMaxIt(int);
		int maxIt();
        	/* load/save funCtions */
		void readHMM(const std::string &folderName);
		void writeHMM(const std::string &folderName);
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
		void calcSymbolNumber();
		void calcSymbolInterval();
		void extractMinValue();
		void extractMaxValue();
		void extractSymbolInterval();
		void calcObsAlphabetFromParas();
		void extractParasFromObsAlphabet();
		void calcCDF(const SMLMS::Matrix& pdf, SMLMS::Matrix& cdf);
		void calcDof(void);
		void calcDof(SMLMS::PhysicalModelBLD &);
		void calcBic(unsigned);
		void calcAic(unsigned);
		void calcModelSelection(unsigned);
		/* print functions */
		void printStateNumber();
		void printMinValue();
		void printMaxValue();
		void printSymbolNumber();
		void printSymbolInterval();
		void printEquiPDF();
		void printTransPDF();
		void printObsPDF();
		void printEquiCDF();
		void printTransCDF();
		void printObsCDF();
		void printObsAlphabet();
		void printLogLikelihood();
		void printDof();
		void printBic();
		void printAic();
		void printHMM();
};/* HMMBase*/
}/* SMLMS */
#endif /* smlmsHmmBase_hpp */

