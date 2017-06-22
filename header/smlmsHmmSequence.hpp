/* ######################################################################
* File Name: smlmsHmmSequence.hpp
* Project: SMLMS
* Version: 17.02
* Creation Date: 23.02.2017
* Created By Sebastian Malkusch
* <malkusch@chemie.uni-frankfurt.de>
* Goethe University of Frankfurt
* Physical and Theoretical Chemistry
* Single Molecule Biophysics
###################################################################### */


#ifndef SmlmsHmmSequence_hpp
#define SmlmsHmmSequence_hpp

#include <vector>
#include <string>
#include <iostream>
#include "header/smlmsHmmBase.hpp"
#include "header/smlmsHmmUnique.hpp"
#include "header/smlmsMatrix.hpp"
#include "header/smlmsJudi.hpp"
//#include "header/smlmsPhysModBrownLatDiff.hpp"

namespace SMLMS{
class HMMSequence: public HMMBase{
	private:
		unsigned _traceNumber;
		bool _fbSeqDone;
		bool _modelAdjustInd;
		double _seqLogLikelihood;
		SMLMS::Matrix _equiPDFNumer;
		SMLMS::Matrix _transPDFNumer;
		SMLMS::Matrix _obsPDFNumer;
		SMLMS::Matrix _pdfDenominator;
	public:
		/* Constructor */
		HMMSequence();
		HMMSequence(unsigned states, unsigned symbols);	
		HMMSequence(unsigned states, unsigned symbols, unsigned trace);	
		/* Destructor */
		~HMMSequence(){std::cout<<"HMM removed from Heap!"<<std::endl;}
		/* copy Constructor */
		HMMSequence(const HMMSequence &);
		/* elementary functions */
		void setTraceNumber(unsigned);
		unsigned traceNumber();
		void setModelAdjustInd(bool);
		bool modelAdjustInd();
		double seqLogLikelihood();
		/* proof functions */
		void checkTraceNumber();
		void checkSimulationDimension(const SMLMS::JumpDistanceList &judi, unsigned obsVal);
		/* init functions */
		void initSeqLogLikelihood();
		void initEquiPDFNumer();
		void initTransPDFNumer();
		void initObsPDFNumer();
		void initPDFDenominator();
		void initTrainingSequences();
		/* clear functions */
		void resetEquiPDFNumer();
		void resetTransPDFNumer();
		void resetObsPDFNumer();
		void resetPDFDenominator();
		void resetTrainingSequences();
		/* print Functions */
		void printTraceNumber();
		void printSeqLogLikelihood();
		void printEquiPDFNumer();
		void printTransPDFNumer();
		void printObsPDFNumer();
		void printPDFDenominator();
		/* help functions */
		void initHelpUniqueHMM(SMLMS::HMMUnique &, unsigned obsNumber);
		void increaseSimulationSequence(SMLMS::JumpDistanceList &, unsigned trace, const std::vector<double> &obs, const std::vector<int> &states);
		// model adjustment
		void reestimateEquiPDF();
		void reestimateTransPDF();
		void reestimateObsPDF();
		void reestimateHMM();
		bool benchmarkTrainingsResult(double llResult, int itStep);
		/* core functions */
		void simulateSequence(unsigned obsVal, SMLMS::JumpDistanceList &);
		void estimateSeqLikelihood(const SMLMS::JumpDistanceList&);
		void estimateSeqBic(unsigned obsNumber);
		void estimateSeqAic();
		void baumWelchSequence(const SMLMS::JumpDistanceList &judi);
		void trainSequence(const SMLMS::JumpDistanceList& judi);
		void trainPhysModSequence(const SMLMS::JumpDistanceList &judi, SMLMS::PhysicalModelBLD &physMod);
		void estimateStateSequence(SMLMS::JumpDistanceList&);
}; /* HMMSequence */
} /* SMLMS */
#endif /* SmlmsHmmSequence_hpp */

