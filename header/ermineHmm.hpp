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
#include "header/ermineJudi.hpp"
//#include "header/probabilityCalculus.hpp"

namespace SMLMS{
class HMM{
	private:
		unsigned _stateNumber;
		double _minValue;
		double _maxValue;
		double _binSize;
		std::vector<double> _stateVector;
		std::vector<double> _stateVectorCdf;
		SMLMS::Matrix _transMatrix;
		SMLMS::Matrix _transMatrixCdf;
		std::vector<double> _obsAlphabet;
		SMLMS::Matrix _obsMatrix;
		SMLMS::Matrix _obsMatrixCdf;
        	SMLMS::Matrix _obsParaMatrix;
	public:
		/* Constructor */
		HMM();	
		//HMM(unsigned, double size, double min, double max);
		/* Destructor */
		~HMM();
		/* copy Constructor */
		HMM(const HMM &);
		/* elementary functions */
		void setStateNumber(unsigned);
		unsigned stateNumber();
		void setMinValue(double);
		double minValue();
		void setMaxValue(double);
		double maxValue();
		void setBinSize(double);
		double binSize();
		void setStateVector(std::vector<double> const &);
		std::vector<double> stateVector();
		void setStateVectorCdf(std::vector<double> const&);
		std::vector<double> stateVectorCdf(void);
		void setTransMatrix(SMLMS::Matrix const &);
		Matrix transMatrix();
		void setTransMatrixCdf(SMLMS::Matrix const &);
		SMLMS::Matrix transMatrixCdf(void);
		void setObsAlphabet(std::vector<double> const &);
		std::vector<double> obsAlphabet(void);
		void setObsMatrix(SMLMS::Matrix const &);
		SMLMS::Matrix obsMatrix();
		void setObsMatrixCdf(SMLMS::Matrix const &);
		SMLMS::Matrix obsMatrixCdf(void);
        	void setObsParaMatrix(SMLMS::Matrix const&);
        	SMLMS::Matrix obsParaMatrix();
        	/* load/save funCtions */
		void readHMM(std::string const &);
		void writeHMM(std::string const &);
		/* special functions */
		void clearHMM();
		void checkHMM();
		void checkAlphabet();
		void checkStateNumber();
		void initStateVector();
		void initTransMatrix();
		void initObsAlphabet(double size, double min, double max);
		void initObsParaMatrix();
		void initObsMatrix();
		void initStateVectorCdf();
		void initTransMatrixCdf();
		void initObsMatrixCdf();
		void initHMM();
		void calcObsMatrix(void);	
		void calcObsMatrixCdf(void);
		void calcStateVectorCdf(void);
		void calcTransMatrixCdf(void);
		void calcCdf(void);
		/* core functions */
		void initialize(SMLMS::JumpDistanceList &, std::string &);
		void simulate(unsigned stepNumber, unsigned traceNumber, SMLMS::JumpDistanceList &);
		//void estimate(JumpDistanceList &);/* judi übergeben */ 
		//void train(); /* Abbruch parameter und judi übergeben */
		//void bestPath(); /* judi übergeben und zurück */
		
		/* help functions */
		double judiPdf(double *, double *);
		double judiSuperPosPdf(double *, double *);
		/*normalize functions */
		void normalizeTransMatrix(void);
		void normalizeObsMatrix(void);
		void normalizeStateVector(void);
		void normalizeHMM(void);
		double integratePdf(std::vector<double> &);
		void normalizePdf(double area, std::vector<double> &pdf);
		/* Match probability functions */
		int findInitStateMatch(double event);
		int findTransStateMatch(double event);
		int findObseMatch(double event);

		std::vector<double> getSingleStateObsPara(int state);
		void calcPdfFromPara(std::vector<double> &pdf, std::vector<double> &para);
		void simStateSequence(unsigned stepNumber, unsigned traceNumber, SMLMS::JumpDistanceList &);
};/* HMM */
}/* SMLMS */
#endif /* hmm_hpp */

