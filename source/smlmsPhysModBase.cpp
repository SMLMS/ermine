/* ######################################################################
* File Name: smlmsPhysModBase.cpp
* Project: ermine
* Version: 19.02
* Creation Date: 28.03.2017
* Created By Sebastian Malkusch
* <malkusch@chemie.uni-frankfurt.de>
* Goethe University of Frankfurt
* Physical and Theoretical Chemistry
* Single Molecule Biophysics
###################################################################### */

#include <iostream>
#include <string>
#include <vector>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <limits>
#include "TH1.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TGraph.h"
#include "header/smlmsPhysModBase.hpp"
#include "header/smlmsMatrix.hpp"
#include "header/smlmsContainer.hpp"
#include "header/smlmsExceptions.hpp"
#include "header/smlmsHmmUnique.hpp"
#include "header/smlmsHmmSequence.hpp"

namespace SMLMS{
/* constructor */
PhysicalModelBase::PhysicalModelBase(){
	//std::cout<<"Physical Model Base constructor called"<<std::endl;
}

PhysicalModelBase::PhysicalModelBase(const std::vector<double> &xVal, int stateVal){
	//std::cout<<"Physical Model Base constructor called"<<std::endl;
	_alphabet = xVal;
	_stateNumber = stateVal;
	initModelByAlphabet();
}


PhysicalModelBase::PhysicalModelBase(double minVal, double maxVal, double incVal, int stateVal){
	//std::cout<<"Physical Model Base constructor called"<<std::endl;
	_minValue = minVal;
	_maxValue = maxVal;
	_incNumber = incVal;
	_stateNumber = stateVal;
	initModelByParameter();
}

/* Destructor */
PhysicalModelBase::~PhysicalModelBase(){
	//std::cout<<"Physical Model Base destructor called"<<std::endl;
}

/* Copy Constructor */
PhysicalModelBase::PhysicalModelBase(const PhysicalModelBase &obj){
	//std::cout<<"Physical Model Base copy constructor called"<<std::endl;
	_minValue = obj._minValue;
	_maxValue = obj._maxValue;
	_binSize = obj._binSize;
	_incNumber = obj._incNumber;
	_stateNumber = obj._stateNumber;
	_alphabet = obj._alphabet;
	_areaPdfSuperPos = obj._areaPdfSuperPos;
	_areaFitSuperPos = obj._areaFitSuperPos;
	_pdfSuperPos = obj._pdfSuperPos;
	_fitSuperPos = obj._fitSuperPos;
	_resSuperPos = obj._resSuperPos;
	_chiSquareSuperPos = obj._chiSquareSuperPos;
	_areaPdf = obj._areaPdf;
	_areaFit = obj._areaFit;
	_pdfMatrix = obj._pdfMatrix;
	_fitMatrix = obj._fitMatrix;
	_resMatrix = obj._resMatrix;
	_chiSquare = obj._chiSquare;
	_pdfWeight = obj._pdfWeight;
}

/* elementary functions */
void PhysicalModelBase::setMinValue(double val){
	_minValue = val;
}

double PhysicalModelBase::minValue(){
	return _minValue;
}

void PhysicalModelBase::setMaxValue(double val){
	_maxValue = val;
}

double PhysicalModelBase::maxValue(){
	return _maxValue;
}

void PhysicalModelBase::setBinSize(double val){
	_binSize = val;
}

double PhysicalModelBase::binSize(){
	return _binSize;
}

void PhysicalModelBase::setIncNumber(unsigned val){
	_incNumber = val;
}

unsigned PhysicalModelBase::incNumber(){
	return _incNumber;
}

void PhysicalModelBase::setStateNumber(unsigned val){
	_stateNumber = val;
}

unsigned PhysicalModelBase::stateNumber(){
	return _stateNumber;
}

void PhysicalModelBase::setAlphabet(std::vector<double> val){
	_alphabet = val;
}

std::vector<double> PhysicalModelBase::alphabet(){
	return _alphabet;
}

double PhysicalModelBase::areaPdfSuperPos(){
	return _areaPdfSuperPos;
}

double PhysicalModelBase::areaFitSuperPos(){
	return _areaFitSuperPos;
}

void PhysicalModelBase::setPdfSuperPos(std::vector<double> val){
	_pdfSuperPos = val;
}

std::vector<double> PhysicalModelBase::pdfSuperPos(){
	return _pdfSuperPos;
}

void PhysicalModelBase::setFitSuperPos(std::vector<double> val){
	_fitSuperPos = val;
}

std::vector<double> PhysicalModelBase::fitSuperPos(){
	return _fitSuperPos;
}

void PhysicalModelBase::setResSuperPos(std::vector<double> val){
	_resSuperPos = val;
}

std::vector<double> PhysicalModelBase::resSuperPos(){
	return _resSuperPos;
}

double PhysicalModelBase::chiSquareSuperPos(){
	return _chiSquareSuperPos;
}

std::vector<double> PhysicalModelBase::areaPdf(){
	return _areaPdf;
}

std::vector<double> PhysicalModelBase::areaFit(){
	return _areaFit;
}

void PhysicalModelBase::setPdfMatrix(SMLMS::Matrix val){
	_pdfMatrix = val;
}

SMLMS::Matrix PhysicalModelBase::pdfMatrix(){
	return _pdfMatrix;
}

void PhysicalModelBase::setFitMatrix(SMLMS::Matrix val){
	_fitMatrix = val;
}

SMLMS::Matrix PhysicalModelBase::fitMatrix(){
	return _fitMatrix;
}

void PhysicalModelBase::setResMatrix(SMLMS::Matrix val){
	_resMatrix = val;
}

SMLMS::Matrix PhysicalModelBase::resMatrix(){
	return _resMatrix;
}

std::vector<double> PhysicalModelBase::chiSquare(){
	return _chiSquare;
}

void PhysicalModelBase::setPdfWeight(std::vector<double> &val){
	_pdfWeight = val;
}

std::vector<double> PhysicalModelBase::pdfWeight(){
	return _pdfWeight;
}

/* description */
void PhysicalModelBase::description(){
	std::stringstream descMessage;
	descMessage<<"\nSMLMS Function:"<<std::endl;
	descMessage<<"Name: PhysicalModelBase:"<<std::endl;
	descMessage<<"Inheritance: None declared"<<std::endl;
	descMessage<<"Purpose: Translates numerical physical models into observation alphabets and vice versa."<<std::endl;
	std::cout<<descMessage.str()<<std::endl;
}

/* print functions */
void PhysicalModelBase::printMinValue(){
	std::cout<<"min alphabet value: "<<_minValue<<std::endl;
}

void PhysicalModelBase::printMaxValue(){
	std::cout<<"max alphabet value: "<<_maxValue<<std::endl;
}

void PhysicalModelBase::printBinSize(){
	std::cout<<"alphabet bin size: "<<_binSize<<std::endl;
}

void PhysicalModelBase::printIncNumber(){
	std::cout<<"number of increments within alphabet: "<<_incNumber<<std::endl;
}

void PhysicalModelBase::printStateNumber(){
	std::cout<<"number of states: "<<_stateNumber<<std::endl;
}

void PhysicalModelBase::printAlphabet(){
	unsigned i;
	std::cout<<"observation alphabet: "<<std::endl;
	for (i=0; i<_incNumber; i++)std::cout<<_alphabet.at(i)<<"\t";
	std::cout<<std::endl;
}

void PhysicalModelBase::printAreaPdfSuperPos(){
	std::cout<<"area of pdf super position: "<<_areaPdfSuperPos<<std::endl;
}

void PhysicalModelBase::printAreaFitSuperPos(){
	std::cout<<"area of fit super position: "<<_areaFitSuperPos<<std::endl;
}

void PhysicalModelBase::printPdfSuperPos(){
	unsigned i;
	std::cout<<"pdf super position: "<<std::endl;
	for(i=0; i<_incNumber; i++) std::cout<<_pdfSuperPos.at(i)<<"\t";
	std::cout<<std::endl;
}

void PhysicalModelBase::printFitSuperPos(){
	unsigned i;
	std::cout<<"fit super position: "<<std::endl;
	for(i=0; i<_incNumber; i++) std::cout<<_fitSuperPos.at(i)<<"\t";
	std::cout<<std::endl;
}

void PhysicalModelBase::printResSuperPos(){
	unsigned i;
	std::cout<<"residues super position: "<<std::endl;
	for(i=0; i<_incNumber; i++) std::cout<<_resSuperPos.at(i)<<"\t";
	std::cout<<std::endl;
}

void PhysicalModelBase::printChiSquareSuperPos(){
	std::cout<<"chi squared of super position: "<<_chiSquareSuperPos<<std::endl;
}

void PhysicalModelBase::printAreaPdf(){
	unsigned i;
	std::cout<<"area of single state pdf: "<<std::endl;
	for(i=0; i<_stateNumber; i++) std::cout<<_areaPdf.at(i)<<"\t";
	std::cout<<std::endl;
}

void PhysicalModelBase::printAreaFit(){
	unsigned i;
	std::cout<<"area of single state pdf fit: "<<std::endl;
	for(i=0; i<_stateNumber; i++) std::cout<<_areaFit.at(i)<<"\t";
	std::cout<<std::endl;
}

void PhysicalModelBase::printPdfMatrix(){
	unsigned i,j;
	std::cout<<"pdf matrix: "<<std::endl;
	for (i=0; i<_incNumber; i++){
		for (j=0; j<_stateNumber; j++)std::cout<<_pdfMatrix.at(j,i)<<"\t";
		std::cout<<std::endl;
	}
}

void PhysicalModelBase::printFitMatrix(){
	unsigned i,j;
	std::cout<<"pdf fit matrix: "<<std::endl;
	for (i=0; i<_incNumber; i++){
		for (j=0; j<_stateNumber; j++)std::cout<<_fitMatrix.at(j,i)<<"\t";
		std::cout<<std::endl;
	}
}

void PhysicalModelBase::printResMatrix(){
	unsigned i,j;
	std::cout<<"residues matrix: "<<std::endl;
	for (i=0; i<_incNumber; i++){
		for (j=0; j<_stateNumber; j++)std::cout<<_resMatrix.at(j,i)<<"\t";
		std::cout<<std::endl;
	}
}

void PhysicalModelBase::printChiSquare(){
	unsigned i;
	std::cout<<"single state chi suqare: "<<std::endl;
	for(i=0; i<_stateNumber; i++)std::cout<<_chiSquare.at(i)<<"\t";
	std::cout<<std::endl;
}

void PhysicalModelBase::printPdfWeight(){
	unsigned i;
	std::cout<<"pdf matrix weight factors: "<<std::endl;
	for(i=0; i<_stateNumber; i++)std::cout<<_pdfWeight.at(i)<<"\t";
	std::cout<<std::endl;
}

/* check functions */
void PhysicalModelBase::checkMinValue(){
	if (_minValue>=_maxValue){
		std::stringstream errorMessage;
		errorMessage<<"Physical model base instance:"<<std::endl;
		errorMessage<<"minValue needs to be smaller than maxValue."<<std::endl;
		SMLMS::SmlmsError error(errorMessage.str());
		throw error;
	}
}

void PhysicalModelBase::checkMaxValue(){
	if (_maxValue<=_minValue){
		std::stringstream errorMessage;
		errorMessage<<"Physical model base instance:"<<std::endl;
		errorMessage<<"minValue needs to be smaller than maxValue."<<std::endl;
		SMLMS::SmlmsError error(errorMessage.str());
		throw error;
	}
}

void PhysicalModelBase::checkBinSize(){
	if (_binSize<=0.0){
		std::stringstream errorMessage;
		errorMessage<<"Physical model base instance:"<<std::endl;
		errorMessage<<"BinSize needs to be a positive float."<<std::endl;
		SMLMS::SmlmsError error(errorMessage.str());
		throw error;
	}
	double doubleIncNumber = (_maxValue - _minValue) / _binSize;
	double tempIncNumber;
	double rest = std::modf(doubleIncNumber, &tempIncNumber); 
	if(rest >  0.0){
		std::stringstream errorMessage;
		errorMessage<<"The container size does not correspond to the observation interval. Thus, no valid observation number can be calculated.!"<<std::endl;
		SMLMS::SmlmsError error(errorMessage.str());
		throw error;
	}
}

void PhysicalModelBase::checkIncNumber(){
	if (_incNumber<=0){
		std::stringstream errorMessage;
		errorMessage<<"Physical model base instance:"<<std::endl;
		errorMessage<<"incNumber needs to be a positive integer."<<std::endl;
		SMLMS::SmlmsError error(errorMessage.str());
		throw error;
	}
}

void PhysicalModelBase::checkStateNumber(){
	if (_stateNumber<=0){
		std::stringstream errorMessage;
		errorMessage<<"Physical model base instance:"<<std::endl;
		errorMessage<<"stateNumber needs to be a positive integer."<<std::endl;
		SMLMS::SmlmsError error(errorMessage.str());
		throw error;
	}
}

void PhysicalModelBase::checkAlphabetSize(){
	if (_alphabet.size()!=_incNumber){
		std::stringstream errorMessage;
		errorMessage<<"Physical model base instance:"<<std::endl;
		errorMessage<<"incNumber differs from alphabet size."<<std::endl;
		SMLMS::SmlmsError error(errorMessage.str());
		throw error;
	}
}

void PhysicalModelBase::checkAlphabetInc(){
	double inc;
	unsigned size=_alphabet.size();
	if (size>1){
		inc = std::abs(_alphabet.at(size-1)-_alphabet.at(size-2));
		if (inc!=_binSize){
			std::stringstream errorMessage;
			errorMessage<<"Physical model base instance:"<<std::endl;
			errorMessage<<"alphabet bin size differs from binSize."<<std::endl;
			SMLMS::SmlmsError error(errorMessage.str());
			throw error;
		}
	}
}

void PhysicalModelBase::checkAreaPdfSuperPos(){
	int prec;
	prec = -1 * std::numeric_limits<double>::max_digits10;
	if ((_areaPdfSuperPos < 1.0-std::pow(1,prec)) or (_areaPdfSuperPos > 1.0 + std::pow(1,prec))){
		std::stringstream errorMessage;
		errorMessage<<"Physical model base instance:"<<std::endl;
		errorMessage<<"areaPdfSuperPos needs to be normalized to 1."<<std::endl;
		SMLMS::SmlmsError error(errorMessage.str());
		throw error;
	}
}

void PhysicalModelBase::checkAreaFitSuperPos(){
	int prec;
	prec = -1 * std::numeric_limits<double>::max_digits10;
	if ((_areaFitSuperPos < 1.0-std::pow(1,prec)) or (_areaFitSuperPos > 1.0 + pow(1,prec))){
		std::stringstream errorMessage;
		errorMessage<<"Physical model base instance:"<<std::endl;
		errorMessage<<"areaFitSuperPos needs to be normalized to 1."<<std::endl;
		SMLMS::SmlmsError error(errorMessage.str());
		throw error;
	}
}

void PhysicalModelBase::checkPdfSuperPos(){
	if(_pdfSuperPos.size() != _incNumber){
		std::stringstream errorMessage;
		errorMessage<<"Physical model base instance:"<<std::endl;
		errorMessage<<"size of pdfSuperPos differs from incNumber."<<std::endl;
		SMLMS::SmlmsError error(errorMessage.str());
		throw error;
	}
}

void PhysicalModelBase::checkFitSuperPos(){
	if(_fitSuperPos.size() != _incNumber){
		std::stringstream errorMessage;
		errorMessage<<"Physical model base instance:"<<std::endl;
		errorMessage<<"size of fitSuperPos differs from incNumber."<<std::endl;
		SMLMS::SmlmsError error(errorMessage.str());
		throw error;
	}
}

void PhysicalModelBase::checkResSuperPos(){
	if(_resSuperPos.size() != _incNumber){
		std::stringstream errorMessage;
		errorMessage<<"Physical model base instance:"<<std::endl;
		errorMessage<<"size of resSuperPos differs from incNumber."<<std::endl;
		SMLMS::SmlmsError error(errorMessage.str());
		throw error;
	}
}

void PhysicalModelBase::checkChiSquareSuperPos(){
	if(_chiSquareSuperPos < 0){
		std::stringstream errorMessage;
		errorMessage<<"Physical model base instance:"<<std::endl;
		errorMessage<<"chi square of super position pdf needs to be an unsigned integer."<<std::endl;
		SMLMS::SmlmsError error(errorMessage.str());
		throw error;
	}
}

void PhysicalModelBase::checkAreaPdfSize(){
	if(_areaPdf.size() != _stateNumber){
		std::stringstream errorMessage;
		errorMessage<<"Physical model base instance:"<<std::endl;
		errorMessage<<"size of areaPdf differs from stateNumber."<<std::endl;
		SMLMS::SmlmsError error(errorMessage.str());
		throw error;
	}
}

void PhysicalModelBase::checkAreaPdf(){
	unsigned i, prec;
	prec = -1 * std::numeric_limits<double>::max_digits10;
	for (i=0; i<_stateNumber; i++){
		if((_areaPdf.at(i) < 1.0 - std::pow(1,prec)) or (_areaPdf.at(i) > 1.0 + std::pow(1,prec))){
			std::stringstream errorMessage;
			errorMessage<<"Physical model base insatnce:"<<std::endl;
			errorMessage<<"at least one member of areaPdf indicates that pdfMatrix is not normalized to 1."<<std::endl;
			SMLMS::SmlmsError error(errorMessage.str());
			throw error;
		}
	}
}

void PhysicalModelBase::checkAreaFitSize(){
	if(_areaFit.size() != _stateNumber){
		std::stringstream errorMessage;
		errorMessage<<"Physical model base instance:"<<std::endl;
		errorMessage<<"size of areaFit differs from stateNumber."<<std::endl;
		SMLMS::SmlmsError error(errorMessage.str());
		throw error;
	}
}

void PhysicalModelBase::checkAreaFit(){
	unsigned i, prec;
	prec = -1 * std::numeric_limits<double>::max_digits10;
	for (i=0; i<_stateNumber; i++){
		if((_areaFit.at(i) < 1.0 - std::pow(1,prec)) or (_areaFit.at(i) > 1.0 + std::pow(1,prec))){
			std::stringstream errorMessage;
			errorMessage<<"Physical model base insatnce:"<<std::endl;
			errorMessage<<"at least one member of areaFit indicates that fitMatrix is not normalized to 1."<<std::endl;
			SMLMS::SmlmsError error(errorMessage.str());
			throw error;
		}
	}
}

void PhysicalModelBase::checkPdfMatrix(){
	if(_pdfMatrix.numberOfRows() != _stateNumber){
		std::stringstream errorMessage;
		errorMessage<<"Physical model base instance:"<<std::endl;
		errorMessage<<"row number of pdf matrix  differs from stateNumber."<<std::endl;
		SMLMS::SmlmsError error(errorMessage.str());
		throw error;
	}
	if(_pdfMatrix.numberOfColumns() != _incNumber){
		std::stringstream errorMessage;
		errorMessage<<"Physical model base instance:"<<std::endl;
		errorMessage<<"column number of pdf matrix differs from incNumber."<<std::endl;
		SMLMS::SmlmsError error(errorMessage.str());
		throw error;
	}
}

void PhysicalModelBase::checkFitMatrix(){
	if(_fitMatrix.numberOfRows() != _stateNumber){
		std::stringstream errorMessage;
		errorMessage<<"Physical model base instance:"<<std::endl;
		errorMessage<<"row number of fit matrix  differs from stateNumber."<<std::endl;
		SMLMS::SmlmsError error(errorMessage.str());
		throw error;
	}
	if(_fitMatrix.numberOfColumns() != _incNumber){
		std::stringstream errorMessage;
		errorMessage<<"Physical model base instance:"<<std::endl;
		errorMessage<<"column number of fit matrix differs from incNumber."<<std::endl;
		SMLMS::SmlmsError error(errorMessage.str());
		throw error;
	}
}

void PhysicalModelBase::checkResMatrix(){
	if(_resMatrix.numberOfRows() != _stateNumber){
		std::stringstream errorMessage;
		errorMessage<<"Physical model base instance:"<<std::endl;
		errorMessage<<"row number of res matrix  differs from stateNumber."<<std::endl;
		SMLMS::SmlmsError error(errorMessage.str());
		throw error;
	}
	if(_resMatrix.numberOfColumns() != _incNumber){
		std::stringstream errorMessage;
		errorMessage<<"Physical model base instance:"<<std::endl;
		errorMessage<<"column number of res matrix differs from incNumber."<<std::endl;
		SMLMS::SmlmsError error(errorMessage.str());
		throw error;
	}
}

void PhysicalModelBase::checkChiSquare(){
	if(_chiSquare.size() != _stateNumber){
		std::stringstream errorMessage;
		errorMessage<<"Physical model base instance:"<<std::endl;
		errorMessage<<"number of chi square values differs from stateNumber."<<std::endl;
		SMLMS::SmlmsError error(errorMessage.str());
		throw error;
	}
}

void PhysicalModelBase::checkPdfWeight(){
	if(_pdfWeight.size() != _stateNumber){
		std::stringstream errorMessage;
		errorMessage<<"Physical model base instance:"<<std::endl;
		errorMessage<<"number of pdf weight factors differs from stateNumber."<<std::endl;
		SMLMS::SmlmsError error(errorMessage.str());
		throw error;
	}
}

void PhysicalModelBase::checkPhysicalModelBase(){
	checkMinValue();
	checkMaxValue();
	checkBinSize();
	checkIncNumber();
	checkStateNumber();
	checkAlphabetSize();
	checkAlphabetInc();
	checkAreaPdfSuperPos();
	checkAreaFitSuperPos();
	checkPdfSuperPos();
	checkFitSuperPos();
	checkResSuperPos();
	checkChiSquareSuperPos();
	checkAreaPdfSize();
	checkAreaPdf();
	checkAreaFitSize();
	checkAreaFit();
	checkPdfMatrix();
	checkFitMatrix();
	checkResMatrix();
	checkChiSquare();
	checkPdfWeight();
}

/* init functions */
void PhysicalModelBase::initModel(){
	unsigned i, j;
	/* init superPos */
	_areaPdfSuperPos = 1.0;
	_areaFitSuperPos = 1.0;
	std::vector<double> tempVec1(_incNumber, _areaPdfSuperPos/double(_incNumber));
	_pdfSuperPos = tempVec1;
	_fitSuperPos = tempVec1;
	std::vector<double> tempVec2(_incNumber, 0.0);
	_resSuperPos = tempVec2;
	_chiSquareSuperPos = 0.0;
	/* init matrices */
	std::vector<double> tempArea(_stateNumber, 1.0);
	SMLMS::Matrix tempMatrix(_stateNumber, _incNumber);
	_resMatrix = tempMatrix;
	_areaPdf = tempArea;
	_areaFit = tempArea;
	for (j=0; j<_stateNumber; j++){
		for (i=0; i<_incNumber; i++){
			tempMatrix(j,i) = 1.0/(double(_incNumber)*_binSize);
		}
	}
	_pdfMatrix = tempMatrix;
	_fitMatrix = tempMatrix;
	std::vector<double> tempChi(_stateNumber, 0.0);
	_chiSquare = tempChi;
	std::vector<double> tempWeight(_stateNumber, 1.0/double(_stateNumber));
	_pdfWeight = tempWeight;
	/* check stuff */
	checkPhysicalModelBase();
}

void PhysicalModelBase::initModelByParameter(){
	unsigned i;
	/* check stuff */
	checkBinSize();
	checkMinValue();
	checkMaxValue();
	/* init alphabet */
	_incNumber =int((_maxValue - _minValue)/double(_binSize));
	std::vector<double> tempVec(_incNumber, 0.0);
	for (i=0; i<_incNumber; i++)tempVec.at(i) = _minValue+((i+1)*_binSize);
	_alphabet = tempVec;
	/* init model */
	initModel();
}

void PhysicalModelBase::initModelByAlphabet(){
	/* init parameter */
	_incNumber = _alphabet.size();
	_maxValue = _alphabet.at(_incNumber-1);
	_binSize = 0.0;
	if (_incNumber>1)_binSize = _alphabet.at(1)-_alphabet.at(0);
	_minValue = _alphabet.at(0)-_binSize;
	/* init model */
	initModel();
}

/* write functions */
void PhysicalModelBase::writePdfSuperPos(const std::string &folderName){
	unsigned i,j;
	std::string name;
	name = folderName;
	name.append("/pdfSuperPos.txt");
	std::ofstream outFile(name.data());
	if (outFile.is_open()){
		outFile<<std::scientific;
		outFile<<std::setprecision(6);
		/* header line */
		outFile<<"# PDF Super Position "<<std::endl;
		/* State Number */
		outFile<<"# Number of states: "<<_stateNumber<<std::endl;
		/* chi Square */
		outFile<<"# chi square: "<<_chiSquareSuperPos<<std::endl;
		/* data */
		outFile<<"r[nm]\tsPDF\tsFIT\tsRES\t";
		for (j=0; j<_stateNumber; j++) outFile<<j<<"PDF\t";
		outFile<<std::endl;
		for (i=0; i<_incNumber; i++){
			outFile<<_alphabet.at(i)<<"\t"<<_pdfSuperPos.at(i)<<"\t"<<_fitSuperPos.at(i)<<"\t"<<_resSuperPos.at(i)<<"\t";
			for(j=0; j<_stateNumber; j++) outFile<<_fitMatrix.at(j,i)*_pdfWeight.at(j)<<"\t";
			outFile<<std::endl;
		}
	}
	outFile.close();
}

void PhysicalModelBase::writePdfMatrix(const std::string &folderName){
	unsigned i,j;
	std::string name;
	name = folderName;
	name.append("/pdfMatrix.txt");
	std::ofstream outFile(name.data());
	if (outFile.is_open()){
		outFile<<std::scientific;
		outFile<<std::setprecision(6);
		/* header line */
		outFile<<"# PDF Matrix"<<std::endl;
		/* State Number */
		outFile<<"# Number of states: "<<_stateNumber<<std::endl;
		/* chi Square */
		outFile<<"# chi square: "<<std::endl;
		for(j=0; j<_stateNumber;j++)outFile<<_chiSquare.at(j)<<"\t";
		outFile<<std::endl;
		/* data */
		outFile<<"r[nm]\t";
		for (j=0; j<_stateNumber; j++) outFile<<j<<"PDF\t"<<j<<"Fit\t"<<j<<"Res\t";
		outFile<<std::endl;
		for (i=0; i<_incNumber; i++){
			outFile<<_alphabet.at(i)<<"\t";
			for(j=0; j<_stateNumber; j++) outFile<<_pdfMatrix.at(j,i)<<"\t"<<_fitMatrix.at(j,i)<<"\t"<<_resMatrix.at(j,i)<<"\t";
			outFile<<std::endl;
		}
	}
	outFile.close();
}

/* plot functions */
void PhysicalModelBase::plotPhysicalModel(const std::string &folderName){
	unsigned i,j;
	/* create Filename */
	std::string name;
	name = folderName;
	name.append("/ProbabilityDensityDistribution.pdf");
	/* define and set style */
	TStyle eStyle("Plain","ermine style");
	eStyle.SetOptStat(0);
	eStyle.SetTitleFontSize(0.05);
	eStyle.SetTitleOffset(1.0, "x");
	eStyle.SetTitleOffset(1.0, "y");
   	eStyle.SetCanvasBorderMode(0);
   	eStyle.SetPadBorderMode(0);
   	eStyle.SetPadColor(0);
   	eStyle.SetCanvasColor(0);
   	eStyle.SetTitleColor(0);
   	eStyle.SetStatColor(0);
	eStyle.cd();
	/* create canvas */
	TCanvas c1("c1", "PDF Matrix", 700, 500);
	c1.Divide(1,2);
	/* change TPad */
	c1.cd(1);
	/* draw pdf to Histogram */
	TH1F pdfHist("PDF", "PDF", _incNumber, _minValue, _maxValue);
	for (i=0; i<_incNumber; i++) pdfHist.SetBinContent(i, _pdfSuperPos.at(i));
	pdfHist.SetLineColor(17);
	pdfHist.SetFillColor(17);
	pdfHist.SetTitle("pdf fit");
	pdfHist.GetYaxis()->SetTitle("P(r) [a.u.]");
	pdfHist.GetXaxis()->SetTitle("r [nm]");
	pdfHist.GetXaxis()->SetTitleSize(0.04);
	pdfHist.GetXaxis()->SetTickLength(0.04);
	pdfHist.GetXaxis()->SetLabelSize(0.04);
	pdfHist.GetYaxis()->SetTitleSize(0.04);
	pdfHist.GetYaxis()->SetTickLength(0.01);
	pdfHist.GetYaxis()->SetLabelSize(0.04);
	pdfHist.GetXaxis()->SetLimits(_minValue, _maxValue);
	pdfHist.SetBarOffset(+1);
	pdfHist.Draw("b");
	pdfHist.SetBit(TObject::kCanDelete); //Delegate Ownership to Canvas
	/* set legend */
	TLegend leg1(0.6,0.7,0.9,0.9);
   	leg1.SetHeader("jump distance dirtsibutions");
	leg1.AddEntry(&pdfHist, "judi distribution pdf", "L");
	std::stringstream legendState;
	/* draw single state judi */
	std::vector<double> ySingleFit(_incNumber);
	std::vector<TGraph> grSingleFit(_stateNumber);
	for (j=0; j<_stateNumber;j++){
		for (i=0; i<_incNumber; i++) ySingleFit.at(i) = _fitMatrix.at(j,i)*_pdfWeight.at(j);
		grSingleFit.at(j) = TGraph(_incNumber, _alphabet.data(), ySingleFit.data());
		grSingleFit.at(j).SetLineColor(j+2);
		grSingleFit.at(j).SetLineWidth(2);
		grSingleFit.at(j).Draw("Csame");
		grSingleFit.at(j).SetBit(TObject::kCanDelete); //Delegate Ownership to Canvas
		legendState.str("");
		legendState<<"model pdf state "<<j+1;
   		leg1.AddEntry(&grSingleFit.at(j),legendState.str().data(), "L");
	}
   	leg1.Draw();
	leg1.SetBit(TObject::kCanDelete); //Delegate Ownership to Canvas
	/* draw pdf super pos fit */
	TGraph grSuperFit(_incNumber, _alphabet.data(), _fitSuperPos.data());
	grSuperFit.SetLineColor(1);
	grSuperFit.SetLineWidth(2);
	grSuperFit.SetLineStyle(2);
	grSuperFit.Draw("Csame");
	grSuperFit.SetBit(TObject::kCanDelete); //Delegate Ownership to Canvas
   	leg1.AddEntry(&grSuperFit,"model pdf super position", "L");
	/* change TPad */
	c1.cd(2);
	/* draw Res Expextation value */
	std::vector<double> expVal(_incNumber, 0.0); 
        TGraph grExpVal(_incNumber,_alphabet.data(),expVal.data()); 
	grExpVal.SetTitle("residuals between model and judi pdf");
	grExpVal.GetYaxis()->SetTitle("P(r) [a.u.]");
	grExpVal.GetXaxis()->SetTitle("r [nm]");
	grExpVal.GetXaxis()->SetTitleSize(0.04);
	grExpVal.GetXaxis()->SetTickLength(0.04);
	grExpVal.GetXaxis()->SetLabelSize(0.04);
	grExpVal.GetYaxis()->SetTitleSize(0.04);
	grExpVal.GetYaxis()->SetTickLength(0.01);
	grExpVal.GetYaxis()->SetLabelSize(0.04);
        grExpVal.SetLineColor(1);   
        grExpVal.SetLineWidth(2);
	grExpVal.GetXaxis()->SetLimits(_minValue, _maxValue);
	double min = calcMin(_resSuperPos);
	double minFac = 0.9;
	if (min<0.0) minFac = 1.1;
	grExpVal.SetMinimum(minFac*min);
	double max = calcMax(_resSuperPos);
	double maxFac = 1.1;
	if (min<0.0) minFac = 0.9;
	grExpVal.SetMaximum(maxFac*max);
	grExpVal.Draw();
	grExpVal.SetBit(TObject::kCanDelete); //Delegate Ownership to Canvas
	/* draw residuals */
	TGraph grResVal(_incNumber, _alphabet.data(), _resSuperPos.data());
        grResVal.SetLineColor(1);   
        grResVal.SetLineWidth(1);
        grResVal.SetMarkerColor(1);
        grResVal.SetMarkerStyle(3);
	grResVal.SetMarkerSize(0.1);
	grResVal.Draw("Psame");
	grResVal.SetBit(TObject::kCanDelete); //Delegate Ownership to Canvas
	/* set legend */
	std::stringstream chiStream;
	chiStream<<"Chi Square: "<<_chiSquareSuperPos<<std::endl;
	//TLegend *leg2 = new TLegend(0.6,0.7,0.9,0.9);
	TLegend leg2(0.6,0.7,0.9,0.9);
   	leg2.SetHeader("residual between model and judi pdf");
	leg2.AddEntry(&grExpVal, "expected residuals", "L");
   	leg2.AddEntry(&grSuperFit,"residuals", "P");
	leg2.AddEntry((TObject*)0, chiStream.str().data(), "");
	leg2.Draw();
	leg2.SetBit(TObject::kCanDelete); //Delegate Ownership to Canvas
	/* save */
	c1.cd();
	c1.Update();
	c1.Print(name.data());
	/* delete root class instances */
	c1.Close();
}

/* norm Functions */
void PhysicalModelBase::intPdfSuperPos(double &area, const std::vector<double> &pdf){
	unsigned i;
	area = 0.0;
	for (i=0; i<_incNumber; i++) area += pdf.at(i);
	area *= _binSize; 
}

void PhysicalModelBase::normPdfSuperPos(double &area, std::vector<double> &pdf){
	unsigned i;
	intPdfSuperPos(area, pdf);
	for (i=0; i<_incNumber; i++) pdf.at(i) /= area;
	intPdfSuperPos(area, pdf);
}

void PhysicalModelBase::intPdfMatrix(std::vector<double> &area, const SMLMS::Matrix &pdf){
	unsigned i,j;
	for (j=0; j<_stateNumber; j++){
		area.at(j)=0.0;
		for (i=0; i<_incNumber; i++) area.at(j) += pdf.at(j,i);
		area.at(j) *= _binSize;
	}
}

void PhysicalModelBase::normPdfMatrix(std::vector<double> &area, SMLMS::Matrix &pdf){
	unsigned i,j;
	intPdfMatrix(area, pdf);
	for (j=0; j<_stateNumber; j++){
		for (i=0; i<_incNumber; i++) pdf(j,i) /= area.at(j);
	}
	intPdfMatrix(area, pdf);
}

void PhysicalModelBase::normalizeModel(){
	/* check stuff */
	checkPdfSuperPos();
	checkFitSuperPos();
	checkPdfMatrix();
	checkFitMatrix();
	/* norm stuff */
	normPdfSuperPos(_areaPdfSuperPos, _pdfSuperPos);
	normPdfSuperPos(_areaFitSuperPos, _fitSuperPos);
	normPdfMatrix(_areaPdf, _pdfMatrix);
	normPdfMatrix(_areaFit, _fitMatrix);
	/* check */
	checkPhysicalModelBase();
}

/* calc functions */
void PhysicalModelBase::calcPdf(const SMLMS::JumpDistanceList &judi){
	calcPdfSuperPos(judi);
	calcPdfMatrix(judi);
}

void PhysicalModelBase::calcPdfSuperPos(const SMLMS::JumpDistanceList &judi){
	std::vector<double> tempPdf(_incNumber, 0.0);
	calcPdfSingle(tempPdf, judi);
	_pdfSuperPos = tempPdf;
}

void PhysicalModelBase::calcFitSuperPos(){
	unsigned i,j;
	/* check stuff */
	checkFitSuperPos();
	checkFitMatrix();
	checkPdfWeight();
	/* calc Fit super pos */
	for(i=0; i<_incNumber; i++){
		_fitSuperPos.at(i)=0.0;
		for(j=0; j<_stateNumber; j++)_fitSuperPos.at(i) += _fitMatrix.at(j,i) * _pdfWeight.at(j);
	}
	normPdfSuperPos(_areaFitSuperPos, _fitSuperPos);
}

void PhysicalModelBase::calcResSuperPos(){
	unsigned i;
	/* check stuff */
	checkPdfSuperPos();
	checkFitSuperPos();
	checkResSuperPos();
	/* calc res super pos */
	for (i=0; i<_incNumber; i++) _resSuperPos.at(i) = _pdfSuperPos.at(i) - _fitSuperPos.at(i);
}

void PhysicalModelBase::calcPdfMatrix(const SMLMS::JumpDistanceList &judi){
	unsigned i,j, jumpNumber;
	jumpNumber = judi.getNumberOfJumps();
	std::vector<double> tempPdf(_incNumber, 0.0);
	SMLMS::JumpDistanceList tempJudi;
	SMLMS::Jump tempJump;
	for (j=0; j<_stateNumber; j++){
		for (i=0; i<jumpNumber; i++){
			tempJump = judi.getJump(i);
			if (unsigned(tempJump.state)==j) tempJudi.addJumpToEnd(tempJump);	
		}
		calcPdfSingle(tempPdf, tempJudi);
		for (i=0; i<_incNumber; i++) _pdfMatrix(j,i) = tempPdf.at(i);
		tempJudi.clearJumpDistanceList();
	}
}

void PhysicalModelBase::calcPdfSingle(std::vector<double> &pdf, const SMLMS::JumpDistanceList &judi){
	unsigned i,j,jumpNumber;
	int prec;
	double tempDist;
	SMLMS::Jump tempJump;
	prec = -1 * std::numeric_limits<float>::max_digits10;
	// make hist
	jumpNumber = judi.getNumberOfJumps();
	for (i=0; i<pdf.size(); i++) pdf.at(i)=0;
	for (i=0; i<jumpNumber; i++){
		tempJump = judi.getJump(i);
		tempDist = tempJump.jumpDistance;
		if ((tempDist<_maxValue)&&(tempDist>=_minValue)){
			j = unsigned (std::floor((tempDist-(std::pow(10.0, prec)))/_binSize));
			j -= unsigned(std::floor(_minValue/_binSize));
			pdf.at(j) +=1;
		}
		
	}
}

void PhysicalModelBase::calcResMatrix(){
	unsigned i,j;
	/* check stuff */
	checkPdfMatrix();
	checkFitMatrix();
	checkResMatrix();
	/* clac res matrix */
	for (j=0; j<_stateNumber; j++){
		for (i=0; i<_incNumber; i++) _resMatrix(j,i) = _pdfMatrix.at(j,i) - _fitMatrix.at(j,i);
	}
}

double PhysicalModelBase::calcMin(const std::vector<double> &pdf){
	unsigned i;
	double min=1;
	for (i=0; i<pdf.size(); i++){
		if(pdf.at(i)<min) min = pdf.at(i); 
	}
	return min;
}

double PhysicalModelBase::calcMax(const std::vector<double> &pdf){
	unsigned i;
	double max=0;
	for (i=0; i<pdf.size(); i++){
		if(pdf.at(i)>max) max = pdf.at(i); 
	}
	return max;
}

void PhysicalModelBase::calcChiSquareSuperPos(){
	unsigned i;
	/* check stuff */
	checkFitSuperPos();
	checkResSuperPos();
	/* calc chi */
	_chiSquareSuperPos = 0.0;
	for (i=0; i<_incNumber; i++) _chiSquareSuperPos += std::pow(_resSuperPos.at(i),2)/ _fitSuperPos.at(i);
}

void PhysicalModelBase::calcChiSquareMatrix(){
	unsigned i,j;
	/* check stuff */
	checkFitMatrix();
	checkResMatrix();
	for (j=0; j<_stateNumber; j++){
		_chiSquare.at(j) = 0.0;
		for (i=0; i<_incNumber; i++) _chiSquare.at(j) += std::pow(_resMatrix.at(j,i),2)/_fitMatrix.at(j,i);
	}
}
}/* SMLMS */
