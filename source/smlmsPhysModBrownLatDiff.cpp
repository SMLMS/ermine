/* ######################################################################
* File Name: smlmsPhysModBrownLatDiff.cpp* Project: SMLMS
* Version: 16.03
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

PhysicalModelBLD::PhysicalModelBLD(const std::vector<double> &xVal, int stateVal, std::string &name): SMLMS::PhysicalModelBase(xVal, stateVal, name){
	std::cout<<"Physical Model (brownian lateral diffusion) constructor called."<<std::endl;
}


PhysicalModelBLD::PhysicalModelBLD(double minVal, double maxVal, int incVal, int stateVal, std::string &name): SMLMS::PhysicalModelBase(minVal, maxVal, incVal, stateVal, name){
	std::cout<<"Physical Model (brownian lateral diffusion) constructor called."<<std::endl;
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
	_fileNameBase = obj._fileNameBase;
	_paraMat = obj._paraMat;
	_paraVect = obj._paraVect;
	_contAreaSuperPos = obj._contAreaSuperPos;
	_contArea = obj._contArea;
}

/* elementary functions */
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
	SMLMS::Matrix tempMat(_stateNumber, 4);
	_paraMat = tempMat;
}

void PhysicalModelBLD::initParaVect(){
	/* check stuff */
	checkStateNumber();
	std::vector<double> tempVect((_stateNumber*4)+1,0.0);
	_paraVect = tempVect;
}

void PhysicalModelBLD::initContArea(){
	/* check stuff */
	checkStateNumber();
	int i;
	_contAreaSuperPos = 0.0;
	std::vector<double> tempVec(_stateNumber,0.0);
	for(i=0; i<_stateNumber; i++)tempVec.at(i)=0.0;
	_contArea = tempVec;
}

/* check functions */
void PhysicalModelBLD::checkParaMat(){
	if(_paraMat.numberOfRows() != _stateNumber){
		std::stringstream errorMessage;
		errorMessage<<"Physical Model Error: Number of states in brownian lateral diffusion parameter Matrix is not consistent with HMM state Number."<<std::endl;
		SMLMS::SmlmsError error(errorMessage.str());
		throw error;
	}
	if(_paraMat.numberOfColumns() != 4){
		std::stringstream errorMessage;
		errorMessage<<"Physical Model Error: Brownial lateral diffusion parameter matrix needs 4 entries for each state."<<std::endl;
		SMLMS::SmlmsError error(errorMessage.str());
		throw error;
	}
}

void PhysicalModelBLD::checkParaVect(){
	if(_paraVect.size() != (4*_stateNumber)+1){
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
	int i,j;
	std::cout<<"\nBrownian lateral diffusion parameter matrix:"<<std::endl;
	std::cout<<"Number of states: "<<_stateNumber<<std::endl<<_stateNumber<<std::endl;;
	std::cout<<"Weight\tD[nm^2/s]\tdt[s]\tsigma[nm]"<<std::endl;
	for (j=0; j<_stateNumber; j++){
		for (i=0; i<4; i++)std::cout<<_paraMat.at(j,i)<<"\t";
		std::cout<<std::endl;
	}
	std::cout<<std::endl;
}

void PhysicalModelBLD::printParaVect(){
	int i, vSize;
	vSize = _paraVect.size();
	std::cout<<"\nBrownian lateral diffusion parameter vector:"<<std::endl;
	for (i=0; i<vSize; i++) std::cout<<_paraVect.at(i)<<"\t";
	std::cout<<std::endl;
}

void PhysicalModelBLD::printContAreaSuperPos(){
	std::cout<<"\nIntegral of super pos function: "<<_contAreaSuperPos<<std::endl;
}

void PhysicalModelBLD::printContArea(){
	int i;
	std::cout<<"\nIntegrals of pdf Matrix functions:"<<std::endl;
	for (i=0; i<_stateNumber; i++) std::cout<<_contArea.at(i)<<"\t";
	std::cout<<std::endl;
}

/* read and write functions */
void PhysicalModelBLD::writePhysMod(){
	int i,j;
	std::string name;
	name = _fileNameBase;
	name.append("physicalModelBLD.txt");
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
		outFile<<"# Model parameter:\n# pi\tD[nm^2/s]\tdt[s]\tsigma"<<std::endl;
		for(j=0; j<_stateNumber; j++){
			for(i=0; i<4; i++)outFile<<_paraMat.at(j,i)<<"\t";
			outFile<<std::endl;
		}
		outFile.close();
	}
}

void PhysicalModelBLD::readPhysMod(const std::string &name){
	int i,j,n;
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
			}
			/* read model parameter to para mat */
			if(n>0 && n<(1+_stateNumber)){
				j=n-1;
				for (i=0; i<4; i++)lineContent>>_paraMat(j, i);
			}
			n+=1;
		}
		inFile.close();
		/* calc para Vector */
		paraMat2paraVect();
	}
}

/* contineous normalization function */
void PhysicalModelBLD::contIntSuperPos(double &area, const std::vector<double> &pdf){
	int i;
	area=0.0;
 	/* integrate TrapZ way */
	for (i=0; i<_incNumber-1; i++){
		 area+= (pdf.at(i)*_binSize) - (0.5*_binSize*(pdf.at(i)-pdf.at(i+1)));
	}
}

void PhysicalModelBLD::contNormSuperPos(double &area, std::vector<double> &pdf){
	int i;
	contIntSuperPos(area, pdf);
	for (i=0; i<_incNumber; i++) pdf.at(i) /= area;
	contIntSuperPos(area, pdf);
}

void PhysicalModelBLD::contIntPdfMat(std::vector<double> &area, const SMLMS::Matrix &pdf){
	int i,j;
	for (j=0; j<_stateNumber; j++){
		area.at(j)=0.0;
 		/* integrate TrapZ way */
		for (i=0; i<_incNumber-1; i++){
		 	area.at(j)+= (pdf.at(j,i)*_binSize) - (0.5*_binSize*(pdf.at(j,i)-pdf.at(j,i+1)));
		}
	}
}

void PhysicalModelBLD::contNormPdfMat(std::vector<double> &area, SMLMS::Matrix &pdf){
	int i,j;
	contIntPdfMat(area, pdf);
	for (j=0; j<_stateNumber; j++){
		for (i=0; i<_incNumber; i++) pdf(j,i) /= area.at(j);
	}
	contIntPdfMat(area, pdf);
}

void PhysicalModelBLD::contNormModel(){
	/* fit Matrix */
	contNormPdfMat(_contArea, _fitMatrix);
	contNormSuperPos(_contAreaSuperPos, _fitSuperPos);
	/* pdf Matrix */
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
	int i,j;
	_paraVect.clear();
	_paraVect.push_back(_stateNumber);
	SMLMS::Matrix tempMat = _paraMat;
	for (j=0; j<_stateNumber; j++) tempMat(j,0) /= (_incNumber*_binSize);
	for (j=0; j<tempMat.numberOfRows(); j++){
		for (i=0; i<tempMat.numberOfColumns(); i++){
			_paraVect.push_back(tempMat.at(j,i));
		}
	}
}

void PhysicalModelBLD::paraVect2paraMat(){
	/* check stuff */
	checkIncNumber();
	checkBinSize();
	checkStateNumber();
	checkParaVect();
	int i,j;
	_paraMat.clearMatrix();
	SMLMS::Matrix tempMat(_stateNumber, 4);
	for (j=0; j<_stateNumber; j++){
		for(i=0; i<4; i++) tempMat(j,i)=_paraVect.at(1+(j*4)+i);
	}
	for (j=0; j<_stateNumber; j++) tempMat(j,0) *= (_incNumber*_binSize);
	_paraMat = tempMat;
}

void PhysicalModelBLD::calcFitSuperPosFromPara(){
	/* check stuff */
	checkParaVect();
	checkFitSuperPos();
	checkAlphabetSize();
	/* pdf by para */
	for (int i=0; i<_incNumber; i++){
		_fitSuperPos.at(i)=judiSuperPosPdf(&_alphabet.at(i), _paraVect.data());
	}
	/* normalize */
	contNormSuperPos(_contAreaSuperPos, _fitSuperPos);
}

void PhysicalModelBLD::calcFitMatrixFromPara(){
	/* check stuff */
	checkParaVect();
	checkFitMatrix();
	checkAlphabetSize();
	/* pdf by para */
	int i,j;
	for (j=0; j<_stateNumber; j++){
		for (int i=0; i<_incNumber; i++){
			_fitMatrix(j,i)=judiPdf(&_alphabet.at(i), &_paraVect[1+(j*4)]);
		}
	}
	/* normalize */
	contNormPdfMat(_contArea, _fitMatrix);
}

void PhysicalModelBLD::updateWeight(const SMLMS::Matrix &pi){
	int j;
	/* check stuff */
	checkParaMat();
	checkParaVect();
	for (j=0; j<_stateNumber; j++) _paraMat(j,0) = pi.at(0,j);
	paraMat2paraVect();
}

void PhysicalModelBLD::updatePi(SMLMS::Matrix &pi){
	int j;
	for (j=0; j<_stateNumber; j++) pi(0,j)=_paraMat.at(j,0);
}

/* fit functions */
void PhysicalModelBLD::fitPdfSuperPos(){
	int i;
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
	contNormSuperPos(_contAreaSuperPos, _pdfSuperPos);
	paraMat2paraVect();
	/* create canvas */
	TCanvas *c1 = new TCanvas("c1", "PDF Matrix", 700, 500);
	/* create Histogramm */
	TH1F* pdfHist = new TH1F("PDF", "PDF", _incNumber,_minValue+(_binSize),_maxValue+(_binSize));
	for (i=0; i<_incNumber; i++){
		pdfHist->SetBinContent(i, _pdfSuperPos.at(i));
	}
	/* normalize */
	double histIntegral;
	histIntegral = pdfHist->GetEntries()*pdfHist->GetBinWidth(0);
	double scale = 1/histIntegral;
	pdfHist->Scale(scale);
	pdfHist->Draw();
	/* create fit function */
	TF1 *fitFcn = new TF1("fitFcn", judiSuperPosPdf ,_minValue,_maxValue,_paraVect.size());
	fitFcn->SetParName(0, "states");
	fitFcn->FixParameter(0,(double) _stateNumber);
	std::string weightName, weightBaseName("weight");
	std::string DName, DBaseName("D");
	std::string dtName, dtBaseName("dt");
	std::string sigmaName, sigmaBaseName("sigma");
	for (i=0; i<_stateNumber; i++){
		weightName=weightBaseName;
		weightName.append(std::to_string(i));
		DName=DBaseName;
		DName.append(std::to_string(i));
		dtName=dtBaseName;
		dtName.append(std::to_string(i));
		sigmaName=sigmaBaseName;
		sigmaName.append(std::to_string(i));
		fitFcn->SetParName(1+(i*4), weightName.data());
 		fitFcn->SetParName(2+(i*4), DName.data());
          	fitFcn->SetParName(3+(i*4), dtName.data());
          	fitFcn->SetParName(4+(i*4), sigmaName.data());
	        fitFcn->SetParameter(1+(i*4),_paraVect.at(1+(i*4)));
        	fitFcn->SetParameter(2+(i*4),_paraVect.at(2+(i*4)));
        	fitFcn->FixParameter(3+(i*4),_paraVect.at(3+(i*4)));
        	fitFcn->FixParameter(4+(i*4),_paraVect.at(4+(i*4)));
	}
	pdfHist->Fit("fitFcn", "I,L,M");
	/* parse fit parameter to new rf1 */
	TF1 *fit = pdfHist->GetFunction("fitFcn");
	for (i=0; i<_paraVect.size(); i++){
		_paraVect.at(i)=fit->GetParameter(i);
	}
	paraVect2paraMat();
	double normFac = 0.0;
	for (i=0; i<_stateNumber; i++) normFac += _paraMat.at(i,0);
	for (i=0; i<_stateNumber; i++) _paraMat(i,0) /= normFac;
	calcFitSuperPosFromPara();
	calcFitMatrixFromPara();
	/* delete root class instances */
	delete fit;
	delete fitFcn;
	delete pdfHist;
	delete c1;
}

void PhysicalModelBLD::fitPdfMatState(int j, SMLMS::Matrix &pdf){
	int i;
	/* check stuff */
	checkParaVect();
	checkPdfMatrix();
	checkFitMatrix();
	checkContArea();
	checkIncNumber();
	checkMinValue();
	checkMaxValue();
	checkBinSize();
	/* create canvas */	
	TCanvas *c1 = new TCanvas("c1", "PDF Matrix", 700, 500);
	/* create histogram */
	TH1F* pdfHist = new TH1F("PDF", "PDF", _incNumber,_minValue+(_binSize),_maxValue+(_binSize));
	for (i=0; i<_incNumber; i++){
		pdfHist->SetBinContent(i, pdf.at(j,i));
	}
	/* normalize */
	double histIntegral;
	histIntegral = pdfHist->GetEntries()*pdfHist->GetBinWidth(0);
	double scale = 1/histIntegral;
	pdfHist->Scale(scale);
	pdfHist->Draw();
	/* create fit function */
	TF1 *fitFcn = new TF1("fitFcn", judiPdf ,_minValue,_maxValue,4);
	std::string weightName("weight");
	std::string DName("D");
	std::string dtName("dt");
	std::string sigmaName("sigma");
	fitFcn->SetParName(0, weightName.data());
	fitFcn->SetParName(1, DName.data());
	fitFcn->SetParName(2, dtName.data());
	fitFcn->SetParName(3, sigmaName.data());
	fitFcn->SetParameter(0,1.0/(_incNumber*_binSize));
	fitFcn->SetParameter(1,_paraVect.at(2+(j*4)));
	fitFcn->FixParameter(2,_paraVect.at(3+(j*4)));
	fitFcn->FixParameter(3,_paraVect.at(4+(j*4)));
	/* fit */
	pdfHist->Fit("fitFcn", "I,L,Q");
	/* parse fit parameter to new rf1 */
	TF1 *fit = pdfHist->GetFunction("fitFcn");
	for (i=1; i<4; i++){
		_paraVect.at(1+i+(j*4))=fit->GetParameter(i);
	}
	/* delete root class instances */
	delete fit;
	delete pdfHist;
	delete c1;
}

void PhysicalModelBLD::fitPdf(SMLMS::Matrix &pdf){
	/* normalize pdf SuperPos */
	contNormPdfMat(_contArea, pdf);
	int j;
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

void PhysicalModelBLD::viterbiFit(){
	//calc PdfMat from judi -> pdfMatrix
	//calc pdfSuperPos from Judi -> pdfSuperPos
	//calc pdfFitMat from para
	//calc fitSuperPos by sum pdfFitMat
	//calc resSuperPos
	//calc resPdfMat
	//calc chi Square
	//plot
	//save
}

} /* SMLMS */
