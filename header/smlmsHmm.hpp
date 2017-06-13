/* ######################################################################
* File Name: hmm.hpp
* Project: SMLMS
* Version: 16.03
* Creation Date: 22.03.2016
* Created By Sebastian Malkusch
* <malkusch@chemie.uni-frankfurt.de>
* Goethe University of Frankfurt
* Physical and Theoretical Chemistry
* Single Molecule Biophysics
###################################################################### */


#ifndef Hmm_hpp
#define Hmm_hpp

#include <vector>
#include <string>
#include "header/smlmsMatrix.hpp"
#include "header/smlmsRandom.hpp"

namespace SMLMS{
class HMM{
	protected:
		SMLMS::SMLMSRandom _randGen;
		unsigned _stateNumber;
		unsigned _symbolNumber;
		SMLMS::Matrix _equiPDF;
		SMLMS::Matrix _transPDF;
		SMLMS::Matrix _obsPDF;
		SMLMS::Matrix _equiCDF;
		SMLMS::Matrix _transCDF;
		SMLMS::Matrix _obsCDF;
		std::vector<double> _obsAlphabet;
		SMLMS::Matrix _alpha;
		std::vector<double> _norm;
		SMLMS::Matrix _normAlpha;
		std::vector<int> _scaleAlpha;
		SMLMS::Matrix _beta;
		SMLMS::Matrix _normBeta;
		std::vector<int> _scaleBeta;
		SMLMS::Matrix _obsProb;
		std::vector<double> _numer;
		std::vector<double> _denom;
		SMLMS::Matrix _xi;
		std::vector<SMLMS::Matrix> _Xi;
		SMLMS::Matrix _gamma;
		double _minl;
		double _invMinl;
		double _likelihood;
		double _logLikelihood;
		bool _fbDone;
	public:
		/* Constructor */
		HMM();
		HMM(unsigned states, unsigned symbols);	
		/* Destructor */
		~HMM();
		/* copy Constructor */
		HMM(const HMM &);
		/* elementary functions */
		void setStateNumber(unsigned);
		unsigned stateNumber();
		void setSymbolNumber(unsigned);
		unsigned symbolNumber();
		void setEquiPDF(SMLMS::Matrix);
		SMLMS::Matrix equiPDF();
		void setTransPDF(SMLMS::Matrix);
		SMLMS::Matrix transPDF();
		void setObsPDF(SMLMS::Matrix);
		SMLMS::Matrix obsPDF(void);
		void setObsAlphabet(std::vector<double>);
		std::vector<double> obsAlphabet();
		void setObservations(std::vector<double>);
		SMLMS::Matrix obsProb();
		double likelihood();
		double logLikelihood();
        	/* load/save funCtions */
		void readHMM(std::string const &);
		void writeHMM(std::string const &);
		/* init functions */
		void clearHMM();
		void checkStateNumber();
		void checkSymbolNumber();
		void checkEqui();
		void checkTrans();
		void checkObs();
		void checkProbDen(SMLMS::Matrix &);
		void checkHMM();
		void initEqui();
		void initTrans();
		void initObs();
		void initHMM();
		void initAlpha(int);
		void initNorm(int);
		void initNormAlpha(int);
		void initScaleAlpha(int);
		void initBeta(int);
		void initNormBeta(int);
		void initScaleBeta(int);
		void initAnalysis(int);
		void initObsProb(int);
		void initXi();
		void initXiVector(int);
		void initGamma(int);
		/* normalization */
		void normalizePDF(SMLMS::Matrix&);
		void normalizeHMM();
		/* calc functions */
		void calcObsAlphabet(double);
		void calcCDF(SMLMS::Matrix& pdf, SMLMS::Matrix& cdf);
		/* print functions */
		void printAlpha();
		void printBeta();
		void printNorm();
		void printXi();
		void printGamma();
		void printHMM();
		/* help functionc */
		int findMatch(SMLMS::Matrix& cdf,int state, double event);
		int obsPDFMatch(double event);
		bool obsSymbolMatch(int symbol, double event);
		void estimateLikelihood(std::vector<double>& observations);
		void estimateLikelihood2(std::vector<double>& observations);
		void calcXi();
		void calcXi2(std::vector<double>& observations);
		void calcGamma();
		void calcGamma2(std::vector<double>& observations);
		void reestimateEquiPDF();
		void reestimateEquiPDF2();
		void reestimateTransPDF();
		void reestimateTransPDF2();
		void reestimateObsPDF2(std::vector<double>& observations);
		/* core functions: simulation */
		void simulate(std::vector<double>& observations, std::vector<int>& states);
		void simulateStates(std::vector<int>&);
		void simulateObservations(std::vector<double>& observations, std::vector<int>& states);
		/* core functions likelihood estimation */
		void calcObsProb(std::vector<double>& observations);
		void forwardAlgorithm(std::vector<double>& observations);
		void forwardAlgorithm2(std::vector<double>& observations);
		void backwardAlgorithm(std::vector<double>& observations);
		void backwardAlgorithm2(std::vector<double>& observations);
		void forwardBackward(std::vector<double>& observations);
		void forwardBackward2(std::vector<double>& observations);
		void baumWelch(std::vector<double>& observations);
		void baumWelch2(std::vector<double>& observations);
		//void simulate(unsigned stepNumber, unsigned traceNumber, SMLMS::JumpDistanceList &);
		//void estimate(JumpDistanceList &);/* judi übergeben */ 
		//void train(); /* Abbruch parameter und judi übergeben */
		//void bestPath(); /* judi übergeben und zurück */

};/* HMM */
}/* SMLMS */
#endif /* hmm_hpp */

