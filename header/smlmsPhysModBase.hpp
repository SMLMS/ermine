/* ######################################################################
* File Name: smlmsPhysModBase.hpp
* Project: SMLMS
* Version: 18.09
* Creation Date: 22.03.2016
* Created By Sebastian Malkusch
* <malkusch@chemie.uni-frankfurt.de>
* Goethe University of Frankfurt
* Physical and Theoretical Chemistry
* Single Molecule Biophysics
###################################################################### */


#ifndef SMLMSPHYSMODBASE_HPP
#define SMLMSPHYSMODBASE_HPP

#include <vector>
#include <string>
#include "header/smlmsMatrix.hpp"
#include "header/smlmsJudi.hpp"

namespace SMLMS{
class PhysicalModelBase{
	protected:
		double _minValue;
		double _maxValue;
		double _binSize;
		unsigned _incNumber;
		unsigned _stateNumber;
		std::vector<double> _alphabet;
		double _areaPdfSuperPos;
		double _areaFitSuperPos;
		std::vector<double> _pdfSuperPos;
		std::vector<double> _fitSuperPos;
		std::vector<double> _resSuperPos;
		double _chiSquareSuperPos;
		std::vector<double> _areaPdf;
		std::vector<double> _areaFit;
		SMLMS::Matrix _pdfMatrix;
		SMLMS::Matrix _fitMatrix;
		SMLMS::Matrix _resMatrix;
		std::vector<double> _chiSquare;
		std::vector<double> _pdfWeight;
	public:
		/* Constructor */
		PhysicalModelBase();
		PhysicalModelBase(const std::vector<double> &xVal, int stateVal);
		PhysicalModelBase(double minVal, double maxVal, double incVal, int stateVal);
		/* Destructor */
		virtual ~PhysicalModelBase();
		/* Copy Constructor */
		PhysicalModelBase(const PhysicalModelBase&);
		/* elementary functions */
		void setMinValue(double);
		double minValue();
		void setMaxValue(double);
		double maxValue();
		void setBinSize(double);
		double binSize();
		void setIncNumber(unsigned);
		unsigned incNumber();
		void setStateNumber(unsigned);
		unsigned stateNumber();
		void setAlphabet(std::vector<double>);
		std::vector<double> alphabet();
		double areaPdfSuperPos();
		double areaFitSuperPos();
		void setPdfSuperPos(std::vector<double>);
		std::vector<double> pdfSuperPos();
		void setFitSuperPos(std::vector<double>);
		std::vector<double> fitSuperPos();
		void setResSuperPos(std::vector<double>);
		std::vector<double> resSuperPos();
		double chiSquareSuperPos();
		std::vector<double> areaPdf();
		std::vector<double> areaFit();
		void setPdfMatrix(SMLMS::Matrix);
		SMLMS::Matrix pdfMatrix();
		void setFitMatrix(SMLMS::Matrix);
		SMLMS::Matrix fitMatrix();
		void setResMatrix(SMLMS::Matrix);
		SMLMS::Matrix resMatrix();
		std::vector<double> chiSquare();
		void setPdfWeight(std::vector<double>&);
		std::vector<double> pdfWeight();
		/* description */
		void description();
		/* print functions */
		void printMinValue();
		void printMaxValue();
		void printBinSize();
		void printIncNumber();
		void printStateNumber();
		void printAlphabet();
		void printAreaPdfSuperPos();
		void printAreaFitSuperPos();
		void printPdfSuperPos();
		void printFitSuperPos();
		void printResSuperPos();
		void printChiSquareSuperPos();
		void printAreaPdf();
		void printAreaFit();
		void printPdfMatrix();
		void printFitMatrix();
		void printResMatrix();
		void printChiSquare();
		void printPdfWeight();
		/* check functions */
		void checkMinValue();
		void checkMaxValue();
		void checkBinSize();
		void checkIncNumber();
		void checkStateNumber();
		void checkAlphabetSize();
		void checkAlphabetInc();
		void checkAreaPdfSuperPos();
		void checkAreaFitSuperPos();
		void checkPdfSuperPos();
		void checkFitSuperPos();
		void checkResSuperPos();
		void checkChiSquareSuperPos();
		void checkAreaPdfSize();
		void checkAreaPdf();
		void checkAreaFitSize();
		void checkAreaFit();
		void checkPdfMatrix();
		void checkFitMatrix();
		void checkResMatrix();
		void checkChiSquare();
		void checkPdfWeight();
		void checkPhysicalModelBase();
		/* init functions */
		void initModel();
		void initModelByParameter();
		void initModelByAlphabet();
		/* write functions */
		void writePdfSuperPos(const std::string &folderName);
		void writePdfMatrix(const std::string &folderName);
		/* plot functions */
		void plotPhysicalModel(const std::string &folderName);
		/* norm functions */
		void intPdfSuperPos(double &area, const std::vector<double> &pdf);
		void normPdfSuperPos(double &area, std::vector<double> &pdf);
		void intPdfMatrix(std::vector<double> &area, const SMLMS::Matrix &pdf);
		void normPdfMatrix(std::vector<double> &area, SMLMS::Matrix &pdf);
		void normalizeModel();
		/* calc funcions */
		void calcPdf(const SMLMS::JumpDistanceList&);
		void calcPdfSuperPos(const SMLMS::JumpDistanceList&);
		void calcFitSuperPos();
		void calcResSuperPos();
		void calcPdfMatrix(const SMLMS::JumpDistanceList&);
		void calcPdfSingle(std::vector<double>&, const SMLMS::JumpDistanceList&);
		void calcResMatrix();
		double calcMin(const std::vector<double> &);
		double calcMax(const std::vector<double> &);
		void calcChiSquareSuperPos();
		void calcChiSquareMatrix();
}; /* PhysicalModelBase */
} /* SMLMS */
#endif /* SmlmsPhysModBase_hpp */

