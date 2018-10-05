/* ######################################################################
* File Name: smlmsPhysModBrownLatDiff.cpp* Project: SMLMS
* Version: 18.09
* Creation Date: 30.03.2017
* Created By Sebastian Malkusch
* <malkusch@chemie.uni-frankfurt.de>
* Goethe University of Frankfurt
* Physical and Theoretical Chemistry
* Single Molecule Biophysics
###################################################################### */

#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <cmath>
#include "TH1.h"
#include "TF1.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "header/smlmsPhysModBrownLatDiff.hpp"
#include "header/smlmsPhysModBase.hpp"
#include "header/smlmsMatrix.hpp"
#include "header/smlmsExceptions.hpp"
#include "header/smlmsMicroscope.hpp"
#include "header/smlmsPdfFunctions.hpp"
#include "header/smlmsHmmUnique.hpp"
#include "header/smlmsHmmSequence.hpp"

namespace SMLMS{
/* Constructor */
PhysicalModelBLD::PhysicalModelBLD(): SMLMS::PhysicalModelBase(){
	std::cout<<"Physical Model (brownian lateral diffusion) constructor called."<<std::endl;
}

PhysicalModelBLD::PhysicalModelBLD(const std::vector<double> &xVal, int stateVal, SMLMS::Microscope microscope): SMLMS::PhysicalModelBase(xVal, stateVal){
	std::cout<<"Physical Model (brownian lateral diffusion) constructor called."<<std::endl;
	_microscope = microscope;
}


PhysicalModelBLD::PhysicalModelBLD(double minVal, double maxVal, int incVal, int stateVal, SMLMS::Microscope microscope): SMLMS::PhysicalModelBase(minVal, maxVal, incVal, stateVal){
	std::cout<<"Physical Model (brownian lateral diffusion) constructor called."<<std::endl;
	_microscope = microscope;
}

/* Destructor */
PhysicalModelBLD::~PhysicalModelBLD(){
	std::cout<<"Physical Model (brownian lateral diffusion) removed from Heap!"<<std::endl;
}

/* Copy Constructor */
PhysicalModelBLD::PhysicalModelBLD(const PhysicalModelBLD &obj){
	std::cout<<"Physical Model (brownian lateral diffusion) copy constructor called"<<std::endl;
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
	_microscope = obj._microscope;
	_paraMat = obj._paraMat;
	_paraVect = obj._paraVect;
	_contAreaSuperPos = obj._contAreaSuperPos;
	_contArea = obj._contArea;
}

/* elementary functions */
void PhysicalModelBLD::setMicroscope(SMLMS::Microscope microscope){
	_microscope = microscope;
}

SMLMS::Microscope PhysicalModelBLD::microscope(){
	return _microscope;
}

void PhysicalModelBLD::setParaMat(SMLMS::Matrix val){
	_paraMat = val;
}

SMLMS::Matrix PhysicalModelBLD::paraMat(){
	return _paraMat;
}

std::vector<double> PhysicalModelBLD::paraVect(){
	return _paraVect;
}

double PhysicalModelBLD::contAreaSuperPos(){
	return _contAreaSuperPos;
}

std::vector<double> PhysicalModelBLD::contArea(){
	return _contArea;
}

/* init functions */
void PhysicalModelBLD::initModelBLD(){
	initParaMat();
	initParaVect();
	initContArea();
}

void PhysicalModelBLD::initParaMat(){
	/* check stuff */
	checkStateNumber();
	SMLMS::Matrix tempMat(_stateNumber, 8);
	_paraMat = tempMat;
}

void PhysicalModelBLD::initParaVect(){
	/* check stuff */
	checkStateNumber();
	std::vector<double> tempVect((_stateNumber*2)+1,0.0);
	_paraVect = tempVect;
}

void PhysicalModelBLD::initContArea(){
	// check stuff
	checkStateNumber();
	unsigned i;
	_contAreaSuperPos = 0.0;
	_areaPdfSuperPos = 0.0;
	std::vector<double> tempVec(_stateNumber,0.0);
	for(i=0; i<_stateNumber; i++)tempVec.at(i)=0.0;
	_contArea = tempVec;
	//_areaPdf = tempVec;
}


/* check functions */
void PhysicalModelBLD::checkMicroscope(){
	if (_microscope.pxlSize() <= 0){
		std::stringstream errorMessage;
		errorMessage<<"Physical Model Error: Pixel size in microscope class needs to be > 0."<<std::endl;
		SMLMS::SmlmsError error(errorMessage.str());
		throw error;
	}
	if (_microscope.intTime() <= 0){
		std::stringstream errorMessage;
		errorMessage<<"Physical Model Error: Integreation time in microscope class needs to be > 0."<<std::endl;
		SMLMS::SmlmsError error(errorMessage.str());
		throw error;
	}
	if (_microscope.locPrec() <= 0){
		std::stringstream errorMessage;
		errorMessage<<"Physical Model Error: Localization precision in microscope class needs to be > 0."<<std::endl;
		SMLMS::SmlmsError error(errorMessage.str());
		throw error;
	}
}

void PhysicalModelBLD::checkParaMat(){
	if(_paraMat.numberOfRows() != _stateNumber){
		std::stringstream errorMessage;
		errorMessage<<"Physical Model Error: Number of states in brownian lateral diffusion parameter Matrix is not consistent with HMM state Number."<<std::endl;
		SMLMS::SmlmsError error(errorMessage.str());
		throw error;
	}
	if(_paraMat.numberOfColumns() != 8){
		std::stringstream errorMessage;
		errorMessage<<"Physical Model Error: Brownial lateral diffusion parameter matrix needs 8 entries for each state."<<std::endl;
		SMLMS::SmlmsError error(errorMessage.str());
		throw error;
	}
}

void PhysicalModelBLD::checkParaVect(){
	if(_paraVect.size() != (2*_stateNumber)+1){
		std::stringstream errorMessage;
		errorMessage<<"Physical Model Error: Parmeter vector is of insufficient size."<<std::endl;
		SMLMS::SmlmsError error(errorMessage.str());
		throw error;
	}
}

void PhysicalModelBLD::checkContArea(){
	if(_contArea.size() != (_stateNumber)){
		std::stringstream errorMessage;
		errorMessage<<"Physical Model Error: Cont area vector is of insufficient size."<<std::endl;
		SMLMS::SmlmsError error(errorMessage.str());
		throw error;
	}
}


/* print functions */
void PhysicalModelBLD::printParaMat(){
	unsigned i,j;
	std::cout<<"\nBrownian lateral diffusion parameter matrix:"<<std::endl;
	std::cout<<"Number of states: "<<_stateNumber<<std::endl<<_stateNumber<<std::endl;;
	std::cout<<"Weight\tfix\tmin\tmax\tD[nm^2/s]\tfix\tmin\tmax"<<std::endl;
	for (j=0; j<_stateNumber; j++){
		for (i=0; i<8; i++)std::cout<<_paraMat.at(j,i)<<"\t";
		std::cout<<std::endl;
	}
	std::cout<<std::endl;
}

void PhysicalModelBLD::printParaVect(){
	unsigned i, vSize;
	vSize = _paraVect.size();
	std::cout<<"\nBrownian lateral diffusion parameter vector:"<<std::endl;
	for (i=0; i<vSize; i++) std::cout<<_paraVect.at(i)<<"\t";
	std::cout<<std::endl;
}

void PhysicalModelBLD::printContAreaSuperPos(){
	std::cout<<"\nIntegral of super pos function: "<<_contAreaSuperPos<<std::endl;
}

void PhysicalModelBLD::printContArea(){
	unsigned i;
	std::cout<<"\nIntegrals of pdf Matrix functions:"<<std::endl;
	for (i=0; i<_stateNumber; i++) std::cout<<_contArea.at(i)<<"\t";
	std::cout<<std::endl;
}


/* read and write functions */
void PhysicalModelBLD::writePhysMod(const std::string &folderName){
	unsigned i,j;
	std::string name;
	name = folderName;
	name.append("/physMod.txt");
	std::ofstream outFile(name.data());
	if (outFile.is_open()){
		outFile<<std::scientific;
		outFile<<std::setprecision(6);
		/* header line */
		outFile<<"# Physical Model for Probaility Density Function:"<<std::endl;
		outFile<<"# Lateral Brownian Diffusion"<<std::endl;
		/* state Number */
		outFile<<"# Number of states:"<<std::endl;
		outFile<<_stateNumber<<std::endl;
		/* Model Parameter */
		outFile<<"# Model parameter:\n# pi\tfix\tmin\tmax\tD[nm^2/s]\tfix\tmin\tmax"<<std::endl;
		for(j=0; j<_stateNumber; j++){
			for(i=0; i<8; i++)outFile<<_paraMat.at(j,i)<<"\t";
			outFile<<std::endl;
		}
		outFile.close();
	}
	else{
		std::stringstream errorMessage;
		errorMessage<<"Physical Model Error: The ermine can not write physMod to file."<<std::endl;
		SMLMS::SmlmsError error(errorMessage.str());
		throw error;
	}
}

void PhysicalModelBLD::readPhysMod(const std::string &name){
	unsigned i,j,n;
	std::string line;
	std::fstream inFile(name.data());
	n=0;
	if (inFile.is_open()){
		while(std::getline(inFile, line)){
			/* skip comments */
			if(!line.find("#")) continue;
			/* read line */
			std::stringstream lineContent(line);
			/* read state Number */
			if (n==0){
				lineContent>>_stateNumber;
				initModelBLD();
				initModelByParameter();
			}
			/* read model parameter to para mat */
			if(n>0 && n<(1+_stateNumber)){
				j=n-1;
				for (i=0; i<8; i++)lineContent>>_paraMat(j, i);
			}
			n+=1;
		}
		inFile.close();
		/* calc para Vector */
		paraMat2paraVect();
		updatePdfWeight();
	}
	else{
		std::stringstream errorMessage;
		errorMessage<<"Physical Model Error: The ermine could not read physMod from file."<<std::endl;
		SMLMS::SmlmsError error(errorMessage.str());
		throw error;
	}
}

/* contineous normalization function */

void PhysicalModelBLD::contIntSuperPos(double &area, const std::vector<double> &pdf){
	unsigned i;
	area=0.0;
 	// integrate TrapZ way
 	/*
	for (i=0; i<_incNumber-1; i++){
		 area+= (pdf.at(i)*_binSize) - (0.5*_binSize*(pdf.at(i)-pdf.at(i+1)));
	}
	*/
	for(i=1; i<_incNumber; i++){
		area += 0.5*_binSize*(pdf.at(i-1)+pdf.at(i));
	}
}

void PhysicalModelBLD::contNormSuperPos(double &area, std::vector<double> &pdf){
	unsigned i;
	contIntSuperPos(area, pdf);
	for (i=0; i<_incNumber; i++) pdf.at(i) /= area;
	contIntSuperPos(area, pdf);
}

void PhysicalModelBLD::contIntPdfMat(std::vector<double> &area, const SMLMS::Matrix &pdf){
	unsigned i,j;
	for (j=0; j<_stateNumber; j++){
		area.at(j)=0.0;
 		// integrate TrapZ way
 		/*
		for (i=0; i<_incNumber-1; i++){
		 	area.at(j)+= (pdf.at(j,i)*_binSize) - (0.5*_binSize*(pdf.at(j,i)-pdf.at(j,i+1)));
		}
		*/
		for(i=1; i<_incNumber; i++){
			area.at(j) += 0.5*_binSize*(pdf.at(j,i-1)+pdf.at(j,i));
		}
	}
}

void PhysicalModelBLD::contNormPdfMat(std::vector<double> &area, SMLMS::Matrix &pdf){
	unsigned i,j;
	contIntPdfMat(area, pdf);
	for (j=0; j<_stateNumber; j++){
		for (i=0; i<_incNumber; i++) pdf(j,i) /= area.at(j);
	}
	contIntPdfMat(area, pdf);
}

void PhysicalModelBLD::contNormModel(){
	// fit Matrix 
	contNormPdfMat(_contArea, _fitMatrix);
	contNormSuperPos(_contAreaSuperPos, _fitSuperPos);
	// pdf Matrix
	contNormPdfMat(_contArea, _pdfMatrix);
	contNormSuperPos(_contAreaSuperPos, _pdfSuperPos);
}

/* model parameter functions */
void PhysicalModelBLD::paraMat2paraVect(){
	/* check stuff */
	checkIncNumber();
	checkBinSize();
	checkStateNumber();
	checkParaMat();
	checkMicroscope();
	/* calculate expDist and transfer to paraVect */
	unsigned j;
	double tempPi, tempDist;
	_paraVect.clear();
	_paraVect.push_back(_stateNumber);
	for (j=0; j<_paraMat.numberOfRows(); j++){
		tempPi = _paraMat(j,0);
		_paraVect.push_back(tempPi);
		tempDist = SMLMS::expectedDistance(_paraMat(j,4), _microscope.intTime(), _microscope.locPrec()); 
		_paraVect.push_back(tempDist);
	}
}

void PhysicalModelBLD::paraVect2paraMat(){
	/* check stuff */
	checkIncNumber();
	checkBinSize();
	checkStateNumber();
	checkParaVect();
	checkMicroscope();
	/* calculate expDiff and transfer to paraMat */
	unsigned j;
	double tempPi, tempDiff;
	for (j=0; j<_stateNumber; j++){
		tempPi=_paraVect.at(1+(j*2));
		_paraMat(j,0) = tempPi;
		tempDiff = SMLMS::expectedDiffCoeff(_paraVect.at(1+(j*2)+1), _microscope.intTime(), _microscope.locPrec());
		_paraMat(j,4) = tempDiff;
	}
}

void PhysicalModelBLD::calcFitSuperPosFromPara(){
	/* check stuff */
	checkParaVect();
	checkFitSuperPos();
	checkAlphabetSize();
	/* pdf by para */
	for (unsigned i=0; i<_incNumber; i++){
		_fitSuperPos.at(i)=judiSuperPosPdf(&_alphabet.at(i), _paraVect.data());
	}
	/* normalize */
	normPdfSuperPos(_areaPdfSuperPos, _fitSuperPos);
	//contNormSuperPos(_contAreaSuperPos, _fitSuperPos);
}

void PhysicalModelBLD::calcFitMatrixFromPara(){
	/* check stuff */
	checkParaVect();
	checkFitMatrix();
	checkAlphabetSize();
	/* pdf by para */
	unsigned i,j;
	for (j=0; j<_stateNumber; j++){
		for ( i=0; i<_incNumber; i++){
			_fitMatrix(j,i)=judiPdf(&_alphabet.at(i), &_paraVect[1+(j*2)]);
		}
	}
	/* normalize */
	//normPdfMatrix(_areaPdf, _fitMatrix);
	contNormPdfMat(_contArea, _fitMatrix);
}

void PhysicalModelBLD::updatePdfWeight(void){
	checkParaMat();
	for (unsigned i=0; i<_stateNumber; i++) _pdfWeight.at(i) = _paraMat.at(i,0);
}

void PhysicalModelBLD::updateWeight(const SMLMS::Matrix &pi){
	unsigned j;
	/* check stuff */
	checkParaMat();
	checkParaVect();
	for (j=0; j<_stateNumber; j++) _paraMat(j,0) = pi.at(0,j);
	paraMat2paraVect();
}

void PhysicalModelBLD::updatePi(SMLMS::Matrix &pi){
	unsigned j;
	for (j=0; j<_stateNumber; j++) pi(0,j)=_paraMat.at(j,0);
}

void PhysicalModelBLD::fixDiffusionCoefficients(void){
	unsigned i;
	for (i=0; i<_stateNumber; i++) _paraMat(i,5)=1;
	paraMat2paraVect();
}

void PhysicalModelBLD::releaseDiffusionCoefficients(void){
	unsigned i;
	for (i=0; i<_stateNumber; i++) _paraMat(i,5)=0;
	paraMat2paraVect();
}

/* fit functions */
void PhysicalModelBLD::fitPdfSuperPos(){
	unsigned i;
	/* check stuff */
	checkParaMat();
	checkParaVect();
	checkPdfSuperPos();
	checkFitSuperPos();
	checkPdfMatrix();
	checkFitMatrix();
	checkIncNumber();
	checkMinValue();
	checkMaxValue();
	checkBinSize();
	/* normalize pdf SuperPos */
	normPdfSuperPos(_areaPdfSuperPos, _pdfSuperPos);
	paraMat2paraVect();
	/* create canvas */
	TCanvas c1 ("c1", "PDF Matrix", 700, 500);
	/* create Histogramm */
	TH1F pdfHist("PDF", "PDF", _incNumber,_minValue+(_binSize),_maxValue+(_binSize));
	for (i=0; i<_incNumber; i++){
		pdfHist.SetBinContent(i, _pdfSuperPos.at(i));
	}
	/* normalize */
	double histIntegral;
	histIntegral = pdfHist.Integral();
	double scale = 1/histIntegral;
	pdfHist.Scale(scale);
	pdfHist.Draw();
	pdfHist.SetBit(TObject::kCanDelete); //Delegate Ownership to Canvas
	/* create fit function */
	TF1 fitFcn("fitFcn", judiSuperPosPdf ,_minValue+_binSize,_maxValue+_binSize,_paraVect.size());
	fitFcn.SetParName(0, "states");
	fitFcn.FixParameter(0,double(_stateNumber));
	std::string weightName, weightBaseName("weight");
	std::string DistName, DistBaseName("distance");
	printParaMat();
	printParaVect();	
	double minDist, maxDist;
	for (i=0; i<_stateNumber; i++){
		weightName=weightBaseName;
		weightName.append(std::to_string(i));
		DistName=DistBaseName;
		DistName.append(std::to_string(i));
		fitFcn.SetParName(1+(i*2), weightName.data());
 		fitFcn.SetParName(2+(i*2), DistName.data());
          	if (_paraMat.at(i,1)>0){
        		fitFcn.FixParameter(1+(i*2),_paraVect.at(1+(i*2)));
		}
		else{
	        	fitFcn.SetParameter(1+(i*2),_paraVect.at(1+(i*2)));
			fitFcn.SetParLimits(1+(i*2), _paraMat.at(i,2), _paraMat.at(i,3));
		}
		if (_paraMat.at(i,5)>0){
			fitFcn.FixParameter(2+(i*2),_paraVect.at(2+(i*2)));
		}
		else{
			fitFcn.SetParameter(2+(i*2),_paraVect.at(2+(i*2)));
			minDist = SMLMS::expectedDistance(_paraMat(i,6), _microscope.intTime(), _microscope.locPrec()); 
			maxDist = SMLMS::expectedDistance(_paraMat(i,7), _microscope.intTime(), _microscope.locPrec()); 
			fitFcn.SetParLimits(2+(i*2), minDist, maxDist);
		}
	}
	pdfHist.Fit("fitFcn", "I,WL,M,Q");
	/* parse fit parameter to new rf1 */
	TF1 fit = *pdfHist.GetFunction("fitFcn");
	for (i=0; i<_paraVect.size(); i++){
		_paraVect.at(i)=std::abs(fit.GetParameter(i));
	}
	paraVect2paraMat();
	double normFac = 0.0;
	for (i=0; i<_stateNumber; i++) normFac += _paraMat.at(i,0);
	for (i=0; i<_stateNumber; i++) _paraMat(i,0) /= normFac;
	updatePdfWeight();
	calcFitSuperPosFromPara();
	calcFitMatrixFromPara();
	/* delete root class instances */
	c1.Close();
}

void PhysicalModelBLD::fitPdfMatState(int j, SMLMS::Matrix &pdf){
	unsigned i;
	/* check stuff */
	checkParaVect();
	checkPdfMatrix();
	checkFitMatrix();
	checkContArea();
	checkAreaPdf();
	checkIncNumber();
	checkMinValue();
	checkMaxValue();
	checkBinSize();
	/* create canvas */	
	TCanvas c1("c1", "PDF Matrix", 700, 500);
	/* create histogram */
	TH1F pdfHist("PDF", "PDF", _incNumber,(_minValue+_binSize),(_maxValue+_binSize));
	for (i=0; i<_incNumber; i++){
		pdfHist.SetBinContent(i, pdf.at(j,i));
	}
	/* normalize */
	double histIntegral;
	histIntegral = pdfHist.Integral();
	double scale = 1/histIntegral;
	pdfHist.Scale(scale);
	pdfHist.Draw();
	pdfHist.SetBit(TObject::kCanDelete); //Delegate Ownership to Canvas
	/* create fit function */
	TF1 fitFcn("fitFcn", judiPdf ,_minValue,_maxValue,2);
	std::string weightName("weight");
	std::string distName("dist");
	fitFcn.SetParName(0, weightName.data());
	fitFcn.SetParName(1, distName.data());
	fitFcn.SetParameter(0,1.0);
	if (_paraMat.at(j, 5) > 0){
		fitFcn.FixParameter(1,_paraVect.at(2+(j*2)));
	}
	else{
		double minDist, maxDist;
		fitFcn.SetParameter(1,_paraVect.at(2+(j*2)));
		minDist = SMLMS::expectedDistance(_paraMat(j,6), _microscope.intTime(), _microscope.locPrec()); 
		maxDist = SMLMS::expectedDistance(_paraMat(j,7), _microscope.intTime(), _microscope.locPrec()); 
		fitFcn.SetParLimits(2+(j*2), minDist, maxDist);
	}
	/* fit */
	pdfHist.Fit("fitFcn", "I,WL,Q");
	/* parse fit parameter to new rf1 */
	TF1 fit = *pdfHist.GetFunction("fitFcn");
	_paraVect.at(2+(j*2)) = fit.GetParameter(1);
	/* delete root class instances */
	c1.Close();
}

void PhysicalModelBLD::fitPdf(SMLMS::Matrix &pdf){
	/* normalize pdf SuperPos */
	contNormPdfMat(_contArea, pdf);
	unsigned j;
	for (j=0; j<_stateNumber; j++) fitPdfMatState(j, pdf);
	paraVect2paraMat();
}

void PhysicalModelBLD::initFit(SMLMS::Matrix &pdf){
	fitPdfSuperPos();
	pdf = _fitMatrix;
	calcResSuperPos();
	}

void PhysicalModelBLD::baumWelchFit(const SMLMS::Matrix &pi, SMLMS::Matrix &pdf){
	updateWeight(pi);
	paraMat2paraVect();
	_pdfMatrix = pdf;
	contNormPdfMat(_contArea, _pdfMatrix);
	fitPdf(_pdfMatrix);
	updateWeight(pi);
	updatePdfWeight();
	calcFitMatrixFromPara();
	calcFitSuperPosFromPara();
	/* normalize histogram */
	normalizeModel();
	/* calc statistics */
	calcResMatrix();
	calcChiSquareMatrix();
	calcResSuperPos();
	calcChiSquareSuperPos();
	/* return pdf fit */
	pdf = _fitMatrix;
}
} /* SMLMS */
