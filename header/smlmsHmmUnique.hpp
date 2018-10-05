/* ######################################################################
* File Name: smlmsHmmUnique.hpp
* Project: SMLMS
* Version: 18.09
* Creation Date: 22.03.2016
* Created By Sebastian Malkusch
* <malkusch@chemie.uni-frankfurt.de>
* Goethe University of Frankfurt
* Physical and Theoretical Chemistry
* Single Molecule Biophysics
###################################################################### */


#ifndef SmlmsHmmUnique_hpp
#define SmlmsHmmUnique_hpp

#include <vector>
#include <string>
#include <iostream>
#include "header/smlmsHmmBase.hpp"
#include "header/smlmsMatrix.hpp"
#include "header/smlmsRandom.hpp"

namespace SMLMS{
class HMMUnique: public HMMBase{
	private:
		SMLMS::SMLMSRandom _randGen;
		unsigned _obsNumber;
		std::vector<double> _normLH;
		SMLMS::Matrix _alpha;
		SMLMS::Matrix _beta;
		std::vector<SMLMS::Matrix> _xi;
		SMLMS::Matrix _gamma;
		SMLMS::Matrix _transPDFNumer; //new
		SMLMS::Matrix _obsPDFNumer; //new
		SMLMS::Matrix _pdfDenominator; //new
		bool _fbDone;
	public:
		/* Constructor */
		HMMUnique();
		HMMUnique(unsigned states, unsigned symbols);	
		HMMUnique(unsigned states, unsigned symbols, unsigned obs);	
		/* Destructor */
		//~HMMUnique();{std::cout<<"HMMUnique removed from Heap!"<<std::endl;}
		/* copy Constructor */
		HMMUnique(const HMMUnique &);
		/* elementary functions */
		void setObsNumber(unsigned);
		unsigned obsNumber();
		std::vector<double> normLH();
		SMLMS::Matrix transPDFNumer();
		SMLMS::Matrix obsPDFNumer();
		SMLMS::Matrix pdfDenominator();
		/* proof functions */
		void checkObsNumber();
		void checkObsLength(unsigned seqLength);
		/* init functions */
		void initNormLH();
		void resetNormLH();
		void initAlpha();
		void resetAlpha();
		void initBeta();
		void resetBeta();
		void initXi();
		void resetXi();
		void initGamma();
		void resetGamma();
		void initTransPDFNumer();
		void resetTransPDFNumer();
		void initObsPDFNumer();
		void resetObsPDFNumer();
		void initPDFDenominator();
		void resetPDFDenominator();
		void initAnalysis();
		void resetAnalysis();	
		/* print functions */
		void printObsNumber();
		void printNormLH();
		void printAlpha();
		void printBeta();
		void printXi();
		void printGamma();
		void printTransPDFNumer();
		void printObsPDFNumer();
		void printPDFDenominator();
		/* help functionc */
		int findMatch(const SMLMS::Matrix& cdf,int state, double event);
		int obsPDFMatch(double event);
		bool obsSymbolMatch(int symbol, double event);
		void calcXi(const std::vector<double> &obsSeq);
		void calcGamma();
		void reestimateEquiPDF();
		void estimateTransPDFNumer(); //new
		void reestimateTransPDF();
		//void reestimateTransPDF02(); //new
		void estimateObsPDFNumer(const std::vector<double> &obsSeq); //new
		void estimatePDFDenominator(); //new
		void reestimateObsPDF(const std::vector<double> &obsSeq);
		//void reestimateObsPDF02(const std::vector<double> &obsSeq); //new
		/* core functions: simulation */
		void simulate(std::vector<double> &obs, std::vector<int> &states);
		void simulateStates(std::vector<int> &stateSeq);
		void simulateObservations(std::vector<double> &obsSeq, const std::vector<int> &stateSeq);
		/* core functions likelihood estimation */
		void forwardAlgorithm(const std::vector<double> &obsSeq);
		void backwardAlgorithm(const std::vector<double> &obsSeq);
		void forwardBackward(const std::vector<double> &obsSeq);
		void estimateLikelihood();
		//void estimateBic();
		//void estimateAic();
		void baumWelch(const std::vector<double> &obsSeq);
		void train(const std::vector<double> &obsSeq);
		void viterbi(std::vector<int> &stateSeq, const std::vector<double> &obsSeq);
};/* HMMUnique */
}/* SMLMS */
#endif /* hmm_hpp */

