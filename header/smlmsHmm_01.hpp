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
	private:
		SMLMS::SMLMSRandom _randGen;
		//std::vector<double> _obsSeq;
		//std::vector<int> _stateSeq;
		//SMLMS::Matrix _obsProb;
		unsigned _obsNumber;
		unsigned _stateNumber;
		unsigned _symbolNumber;
		double _symbolInterval;
		SMLMS::Matrix _equiPDF;
		SMLMS::Matrix _transPDF;
		SMLMS::Matrix _obsPDF;
		SMLMS::Matrix _equiCDF;
		SMLMS::Matrix _transCDF;
		SMLMS::Matrix _obsCDF;
		std::vector<double> _obsAlphabet;
		std::vector<double> _normLH;
		SMLMS::Matrix _alpha;
		SMLMS::Matrix _beta;
		std::vector<SMLMS::Matrix> _xi;
		SMLMS::Matrix _gamma;
		double _logLikelihood;
		bool _fbDone;
	public:
		/* Constructor */
		HMM();
		HMM(unsigned states, unsigned symbols, unsigned obs);	
		/* Destructor */
		~HMM();
		/* copy Constructor */
		HMM(const HMM &);
		/* elementary functions */
		//void setObsSeq(std::vector<double>);
		//std::vector<double> obsSeq();
		//void setStateSeq(std::vector<int>);
		//std::vector<int> stateSeq();
		void setObsNumber(unsigned);
		unsigned obsNumber();
		void setStateNumber(unsigned);
		unsigned stateNumber();
		void setSymbolNumber(unsigned);
		unsigned symbolNumber();
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
		std::vector<double> normLH();
		double logLikelihood();
        	/* load/save funCtions */
		void readHMM(std::string const &);
		void writeHMM(std::string const &);
		/* clear functions */
		void clearHMM();
		/* proof functions */
		void checkObsNumber();
		void checkObsLength(unsigned seqLength);
		void checkStateNumber();
		void checkSymbolNumber();
		void checkEqui();
		void checkTrans();
		void checkObs();
		void checkProbDen(SMLMS::Matrix &);
		void checkHMM();
		/* init functions */
		// hier gehts weiter
		// void initObsSeq();
		// void initStateSeq();	
		void initEqui();
		void initTrans();
		void initObs();
		//void initObsProb();
		//void resetObsProb();
		void initHMM();
		void initAlpha();
		void resetAlpha();
		void initNormLH();
		void resetNormLH();
		void initBeta();
		void resetBeta();
		void initXi();
		void resetXi();
		void initGamma();
		void resetGamma();
		void initAnalysis();
		void resetAnalysis();	
		/* normalization */
		void normalizePDF(SMLMS::Matrix&);
		void normalizeHMM();
		/* calc functions */
		void calcObsAlphabet(double);
		void calcCDF(const SMLMS::Matrix& pdf, SMLMS::Matrix& cdf);
		/* print functions */
		void printObsNumber();
		void printStateNumber();
		void printSymbolNumber();
		//void printObsSeq();
		//void printObsProb();
		//void printStateSeq();
		void printEquiPDF();
		void printTransPDF();
		void printObsPDF();
		void printEquiCDF();
		void printTransCDF();
		void printObsCDF();
		void printObsAlphabet();
		void printAlpha();
		void printBeta();
		void printNormLH();
		void printXi();
		void printGamma();
		void printLogLikelihood();
		void printHMM();
		/* help functionc */
		int findMatch(const SMLMS::Matrix& cdf,int state, double event);
		int obsPDFMatch(double event);
		bool obsSymbolMatch(int symbol, double event);
		void calcXi(const std::vector<double> &obsSeq);
		void calcGamma();
		void reestimateEquiPDF();
		void reestimateTransPDF();
		void reestimateObsPDF(const std::vector<double> &obsSeq);
		/* core functions: simulation */
		void simulate(std::vector<double> &obs, std::vector<int> &states);
		void simulateStates(std::vector<int> &stateSeq);
		void simulateObservations(std::vector<double> &obsSeq, const std::vector<int> &stateSeq);
		/* core functions likelihood estimation */
		void forwardAlgorithm(const std::vector<double> &obsSeq);
		void backwardAlgorithm(const std::vector<double> &obsSeq);
		void forwardBackward(const std::vector<double> &obsSeq);
		void estimateLikelihood();
		void baumWelch(const std::vector<double> &obsSeq);
		void train(const std::vector<double> &obsSeq, double stopCrit, int maxIt);
		//void bestPath(); /* judi übergeben und zurück */

};/* HMM */
}/* SMLMS */
#endif /* hmm_hpp */

