/* ######################################################################
* File Name:
* Project: 
* Version:
* Creation Date:
* Created By Sebastian Malkusch
* <malkusch@chemie.uni-frankfurt.de>
* Goethe University of Frankfurt
* Physical and Theoretical Chemistry
* Single Molecule Biophysics
###################################################################### */


#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <cmath>
#include <algorithm>
#include "TH1.h"
#include "TF1.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "header/smlmsDwellTime.hpp"
#include "header/smlmsExceptions.hpp"
#include "header/smlmsContainer.hpp"
#include "header/smlmsJudi.hpp"
#include "header/smlmsPdfFunctions.hpp"
#include "header/smlmsHmmBase.hpp"
#include "header/smlmsMatrix.hpp"

namespace SMLMS{
/* constructor */
DwellTimeAnalysis::DwellTimeAnalysis(){
	std::cout<<"DwellTimeAnalysis constructor called."<<std::endl;
}

/* destructor */
DwellTimeAnalysis::~DwellTimeAnalysis(){
	std::cout<<"DwellTimeAnalysis is removed from heap."<<std::endl;
}

/* copy-constructor */
DwellTimeAnalysis::DwellTimeAnalysis(const DwellTimeAnalysis &obj){
	std::cout<<"DwellTimeAnalysis copy-constructor called."<<std::endl;
	_folderName = obj._folderName;
	_stateNumber = obj._stateNumber;
	_dt = obj._dt;
	_originalState = obj._originalState;
	_targetState = obj._targetState;
	_transTimes = obj._transTimes;
	_time = obj._time;
	_transProbData = obj._transProbData;
	_transProbFit = obj._transProbFit;
	_transProbRes = obj._transProbRes;
	_amplitude = obj._amplitude;
	_transRate = obj._transRate;
	_transRateError = obj._transRateError;
	_chiSquare = obj._chiSquare;
}

/* elementary functions */
void DwellTimeAnalysis::setFolderName(std::string name){
	_folderName = name;
}

std::string DwellTimeAnalysis::folderName(){
	return _folderName;
}

void DwellTimeAnalysis::setStateNumber(unsigned number){
	_stateNumber = number;
}

unsigned DwellTimeAnalysis::stateNumber(){
	return _stateNumber;
}

void DwellTimeAnalysis::setDt(double t){
	_dt = t;
}

double DwellTimeAnalysis::dt(){
	return _dt;
}

void DwellTimeAnalysis::setOriginalState(unsigned state){
	_originalState = state;
}

unsigned DwellTimeAnalysis::originalState(){
	return _originalState;
}

void DwellTimeAnalysis::setTargetState(unsigned state){
	_targetState = state;
}

unsigned DwellTimeAnalysis::targetState(){
	return _targetState;
}

void DwellTimeAnalysis::setTransTimes(std::vector<double> time){
	_transTimes = time;
}

std::vector<double> DwellTimeAnalysis::transTimes(){
	return _transTimes;
}

void DwellTimeAnalysis::setTime(std::vector<double> time){
	_time = time;
}

std::vector<double> DwellTimeAnalysis::time(){
	return _time;
}

void DwellTimeAnalysis::setTransProbData(std::vector<double> data){
	_transProbData = data;
}

std::vector<double> DwellTimeAnalysis::transProbData(){
	return _transProbData;
}

void DwellTimeAnalysis::setTransProbFit(std::vector<double> data){
	_transProbFit = data;
}

void DwellTimeAnalysis::setAmplitude(double amp){
	_amplitude = amp;
}

double DwellTimeAnalysis::amplitude(){
	return _amplitude;
}

void DwellTimeAnalysis::setTransRate(double rate){
	_transRate = rate;
}

double DwellTimeAnalysis::transRate(){
	return _transRate;
}

void DwellTimeAnalysis::setTransRateError(double error){
	_transRateError = error;
}

double DwellTimeAnalysis::transRateError(){
	return _transRateError;
}

void DwellTimeAnalysis::setChiSquare(double chi){
	_chiSquare = chi;
}

double DwellTimeAnalysis::chiSquare(){
	return _chiSquare;
}

/* proof functions */
void DwellTimeAnalysis::proofHistogram(){
	if(_time.empty()){
		std::stringstream errorMessage;
		errorMessage<<"Histogram vectors in DwellTimeAnalysis are empty!"<<std::endl;
		SMLMS::SmlmsError error(errorMessage.str());
		throw error;
	}
	if((_time.size() != _transProbData.size() ) or (_time.size() != _transProbFit.size()) or (_time.size() != _transProbRes.size())){
		std::stringstream errorMessage;
		errorMessage<<"Histogram vectors in DwellTimeAnalysis must have the same size!"<<std::endl;
		SMLMS::SmlmsError error(errorMessage.str());
		throw error;
	}
}

/* print functions */
void DwellTimeAnalysis::printFolderName(){
	std::cout<<"folder name: "<<_folderName<<std::endl;
}

void DwellTimeAnalysis::printStateNumber(){
	std::cout<<"state number: "<<_stateNumber<<std::endl;
}

void DwellTimeAnalysis::printDt(){
	std::cout<<"time interval: "<<_dt<<std::endl;
}

void DwellTimeAnalysis::printOriginalState(){
	std::cout<<"original state: "<<_originalState<<std::endl;
}

void DwellTimeAnalysis::printTargetState(){
	std::cout<<"target state: "<<_targetState<<std::endl;
}

void DwellTimeAnalysis::printTransTimes(){
	int i;
	std::cout<<"transition times:"<<std::endl;
	for (i=0; i<_transTimes.size(); i++){
		std::cout<<_transTimes.at(i)<<std::endl;
	}
}

void DwellTimeAnalysis::printHist(){
	int i=0;
	std::cout<<"dwell time analysis histogram:"<<std::endl;
	std::cout<<"time[s]\tdata[a.u.]\tfit[a.u.]\tres[a.u.]"<<std::endl;
	for (i=0; i<_time.size(); i++){
		std::cout<<_time.at(i)<<"\t"<<_transProbData.at(i)<<"\t"<<_transProbFit.at(i)<<"\t"<<_transProbRes.at(i)<<std::endl;
	}
}

void DwellTimeAnalysis::printAmplitude(){
	std::cout<<"amplitude is: "<<_amplitude<<std::endl;
}

void DwellTimeAnalysis::printTransRate(){
	std::cout<<"transition rate: "<<_transRate<<std::endl;
	std::cout<<"transition rate error: "<<_transRateError<<std::endl;
}

void DwellTimeAnalysis::printChiSquare(){
	std::cout<<"chi square: "<<_chiSquare<<std::endl;
}

void DwellTimeAnalysis::printResult(){
	std::cout<<std::endl<<"dwell time analysis result for transition "<<_originalState<<" -> "<<_targetState<<std::endl;
	printTransRate();
	printChiSquare();
}

/* write functions */
void DwellTimeAnalysis::writeDwellTime(){
	int i;
	proofHistogram();
	std::string name;
	name = _folderName;
	std::stringstream fileName;
	fileName<<"/DwellTimeAnalysis_"<<_originalState<<"to"<<_targetState<<".txt";
	name.append(fileName.str());
	std::ofstream outFile(name.data());
	if (outFile.is_open()){
		outFile<<std::scientific;
		outFile<<std::setprecision(6);
		/* header line */
		outFile<<"# Dwell Time Analysis"<<std::endl;
		outFile<<"# ermine format"<<std::endl;
		/* number of states */
		outFile<<"# HMM state number: "<<_stateNumber<<std::endl;
		/* transition */
		outFile<<"# analyzed transition: "<<_originalState<<" -> "<<_targetState<<std::endl;
		/* amplitide */
		outFile<<"# amplitude [a.u.]: "<<_amplitude<<std::endl;
		/* transition rate */
		outFile<<"# transition rate [s^-1]: "<<_transRate<<" + "<<_transRateError<<std::endl;
		/* chi square */
		outFile<<"# chi square: "<<_chiSquare<<std::endl;
		/* data info */
		outFile<<"# time[s]\tdata[a.u.]\tfit[a.u.]\tres[a.u.]"<<std::endl;
		for(i=0; i<_time.size(); i++){
			outFile<<_time.at(i)<<"\t"<<_transProbData.at(i)<<"\t"<<_transProbFit.at(i)<<"\t"<<_transProbRes.at(i)<<std::endl;
		}
		outFile.close();
	}
}

void DwellTimeAnalysis::plotDwellTime(){
	int i;
	proofHistogram();
	/* create filename */
	std::string name;
	name =_folderName;
	std::stringstream fileName;
	fileName<<"/DwellTimeAnalysis_"<<_originalState<<"to"<<_targetState<<".pdf";
	name.append(fileName.str());
	/* define and set style */
	gStyle->SetOptStat(0);
	gStyle->SetTitleFontSize(0.05);
	gStyle->SetTitleOffset(1.0, "x");
	gStyle->SetTitleOffset(1.0, "y");
	/* create canvas */
	TCanvas *c1 = new TCanvas("c1", "Dwell time analysis", 700, 500);
	c1->Divide(1,2);
	/* change TPad */
	c1->cd(1);
	/* draw dwell time data histogram */
	TH1F* dtHist = new TH1F("Dwell Time Analysis", "Dwell Time Analysis", _time.size(), _time.at(0), _time.back());
	for(i=0; i<_time.size(); i++) dtHist->SetBinContent(i, _transProbData.at(i));
	dtHist->SetLineColor(17);
	dtHist->SetFillColor(17);
	dtHist->SetTitle("dwell time");
	dtHist->GetYaxis()->SetTitle("fraction of molecules [a.u]");
	dtHist->GetXaxis()->SetTitle("t [s]");
	dtHist->GetXaxis()->SetTitleSize(0.04);
	dtHist->GetXaxis()->SetTickLength(0.04);
	dtHist->GetXaxis()->SetLabelSize(0.04);
	dtHist->GetYaxis()->SetTitleSize(0.04);
	dtHist->GetYaxis()->SetTickLength(0.01);
	dtHist->GetYaxis()->SetLabelSize(0.04);
	dtHist->SetBarOffset(0.5);
	dtHist->Draw("b");
	/* draw dwell time fit */
	TGraph *dtFit = new TGraph(_time.size(), _time.data(), _transProbFit.data());
	dtFit->SetLineColor(2);
	dtFit->SetLineWidth(2);
	dtFit->SetMinimum(0.0);
	dtFit->SetMaximum(1.0);
	dtFit->Draw("Csame");
	/* set legend */
	TLegend *leg1 = new TLegend(0.6, 0.7, 0.9, 0.9);
	std::stringstream header;
	header<<"dwell time analysis for state transition "<<_originalState<<" to "<<_targetState;
	leg1->SetHeader(header.str().data());
	leg1->AddEntry(dtHist, "dwell time data", "L");
	leg1->AddEntry(dtFit, "dwell time Fit", "L");
	std::stringstream result;
	result<<"decay rate: "<<_transRate<<"+ "<<_transRateError<<" [s^-1]";
	leg1->AddEntry((TObject*)0, result.str().data(), "");
	leg1->Draw();
	/* change TPad */
	c1->cd(2);
	/* draw Res Expextation value */
	std::vector<double> expVal(_time.size(), 0.0); 
        TGraph *grExpVal = new TGraph(_time.size(),_time.data(),expVal.data()); 
	grExpVal->SetTitle("residuals between dwell time data and fit");
	grExpVal->GetYaxis()->SetTitle("fraction of molecules [a.u.]");
	grExpVal->GetXaxis()->SetTitle("t [s]");
	grExpVal->GetXaxis()->SetTitleSize(0.04);
	grExpVal->GetXaxis()->SetTickLength(0.04);
	grExpVal->GetXaxis()->SetLabelSize(0.04);
	grExpVal->GetYaxis()->SetTitleSize(0.04);
	grExpVal->GetYaxis()->SetTickLength(0.01);
	grExpVal->GetYaxis()->SetLabelSize(0.04);
        grExpVal->SetLineColor(1);   
        grExpVal->SetLineWidth(2);
	grExpVal->GetXaxis()->SetLimits(_time.at(0), _time.back());
	double min = *std::min_element(_transProbRes.begin(), _transProbRes.end());
	double minFac = 0.9;
	if (min<0.0) minFac = 1.1;
	grExpVal->SetMinimum(minFac*min);
	double max = *std::max_element(_transProbRes.begin(), _transProbRes.end());
	double maxFac = 1.1;
	if (min<0.0) minFac = 0.9;
	grExpVal->SetMaximum(maxFac*max);
	grExpVal->Draw();
	/* draw residuals */
	TGraph *grResVal = new TGraph(_time.size(), _time.data(), _transProbRes.data());
        grResVal->SetLineColor(1);   
        grResVal->SetLineWidth(1);
        grResVal->SetMarkerColor(1);
        grResVal->SetMarkerStyle(3);
	grResVal->SetMarkerSize(0.1);
	grResVal->Draw("Psame");
	/* set legend */
	std::stringstream chiStream;
	chiStream<<"Chi Square: "<<_chiSquare<<std::endl;
	TLegend *leg2 = new TLegend(0.6,0.7,0.9,0.9);
   	leg2->SetHeader("residual between dwell time data and fit");
	leg2->AddEntry(grExpVal, "expected residuals", "L");
   	leg2->AddEntry(grResVal,"residuals", "P");
	leg2->AddEntry((TObject*)0, chiStream.str().data(), "");
	leg2->Draw();
	/* save */
	c1->cd();
	c1->Update();
	c1->Print(name.data());
	/* delete root class instances */
	delete leg1;
	delete leg2;
	delete dtHist;
	delete grResVal;
	delete grExpVal;
	delete c1;
}

/* help functions */
bool DwellTimeAnalysis::startIndex(int &start, std::vector<int> &stateList){
	/* retuns index of first timepoint in right state to be measured */
	int i;
	bool success=false;
	for (i = start+1; i<stateList.size(); i++){
		if ((stateList.at(i-1)!=_originalState) && (stateList.at(i)==_originalState)){
			success = true;
			start = i;
			break;
		}
	}
	return success;
}

bool DwellTimeAnalysis::stopIndex(int &stop, std::vector<int> &stateList){
	/* retuns index of last time point in right state to be measured */
	int i;
	bool success = false;
	for (i = stop+1; i<stateList.size(); i++){
		if((stateList.at(i-1)==_originalState) && (stateList.at(i)!=_originalState)){
			stop = i-1;
			if (stateList.at(i) == _targetState) success = true;
			break;
		}
	}
	return success;
}

void DwellTimeAnalysis::analyzeTraceTransitions(std::vector<int> &stateList){
	int i; 
	int start=0;
	int stop;
	double timeInterval = 0.0;
	while (startIndex(start, stateList)){
		stop = start;
		if(stopIndex(stop, stateList)){
			timeInterval = (stop-start+1.0)*_dt;
			_transTimes.push_back(timeInterval);
			start = stop;
		}
	}
}

/* special functions */
void DwellTimeAnalysis::estimateTransTimes(const JumpDistanceList &judi){
	int i;
	std::vector<int> stateList;
	_transTimes.clear();
	for (i=0; i<judi.traceNumber(); i++){
		stateList.clear();
		stateList = judi.getTraceStates(i+1);
		analyzeTraceTransitions(stateList);
	}
}

void DwellTimeAnalysis::transTimesToHist(){
	int i, j;
	double tempTime;
	double tMin = _dt;
	double tMax;
	/* find tMax */
	tMax = *std::max_element(_transTimes.begin(), _transTimes.end());
	/* init histogram */
	std::vector<double> tempVec(((tMax/_dt)+10), 0.0);
	_time = tempVec;
	_transProbData = tempVec;
	_transProbFit = tempVec;
	_transProbRes = tempVec;
	/* calc time */
	for (i=0; i<_time.size(); i++){
		_time.at(i) = (i)*_dt;
	}
	/* estimate Data Histogram */
	for (i=0; i<_transTimes.size(); i++){
		tempTime = _transTimes.at(i);
		if ((tempTime<=tMax)&&(tempTime>=tMin)){
			j = unsigned (floor(tempTime/_dt));
			tempVec.at(j) +=1;
		}
	}
	/* transform to normalized decay Histogram */
	_amplitude = 0;
	for (i=0; i<tempVec.size(); i++){
		_amplitude += tempVec.at(i);
	}
	_transProbData.at(0) = _amplitude - tempVec.at(0);
	for (i=1; i<_transProbData.size(); i++){
		_transProbData.at(i) = _transProbData.at(i-1) - tempVec.at(i);
	}
	for (i=0; i<_transProbData.size(); i++){
		_transProbData.at(i) /= _amplitude;
	}
	_amplitude = _transProbData.at(0);
}


void DwellTimeAnalysis::fitHist(){
	int i, incNumber;
	/* check stuff */
	proofHistogram();
	/* get incNumber */
	incNumber = _time.size();
	/* create canvas */
	TCanvas *c1 = new TCanvas("c1", "dwell time analysis", 700, 500);
	/* create histogram */
	TH1F* decayHist = new TH1F("decayData", "decayData", incNumber,_time.at(0)+(0.5*_dt),_time.back()+(0.5*_dt));
	for (i=0; i<incNumber; i++){
		decayHist->SetBinContent(i, _transProbData.at(i));
	}
	decayHist->Draw();
	/* create fit function */
	TF1 *decayFit = new TF1("decayFit", monoExpDecFunc,_time.at(0),_time.back(),2);
	std::string AName("amplitude");
	std::string kName("rate");
	decayFit->SetParName(0, AName.data());
	decayFit->SetParName(1, kName.data());
	decayFit->FixParameter(0,_amplitude);
	decayFit->SetParameter(1,_transRate);
	/* fit */
	decayHist->Fit("decayFit", "I,L,Q");
	/* parse fit parameter to new rf1 */
	TF1 *fit = decayHist->GetFunction("decayFit");
	_amplitude = fit->GetParameter(0);
	_transRate = fit->GetParameter(1);
	_transRateError = fit->GetParError(1);
	/* delete root class instances */
	delete fit;
	delete decayHist;
	delete c1;
}

void DwellTimeAnalysis::estimateTransProbFit(){
	int i;
	/* proof hist */
	proofHistogram();
	/* estimate fit */
	std::vector<double> paraVect(2, 0.0);
	paraVect.at(0) = _amplitude;
	paraVect.at(1) = _transRate;
	for(i=0; i<_time.size(); i++){
		_transProbFit.at(i) = monoExpDecFunc(&_time.at(i), &paraVect[0]);
	}
}

void DwellTimeAnalysis::estimateResiduals(){
	int i;
	/* check stuff */
	proofHistogram();
	for (i=0; i<_time.size(); i++) _transProbRes.at(i) = _transProbData.at(i) - _transProbFit.at(i);
}

void DwellTimeAnalysis::estimateChiSquare(){
	int i;
	/* check stuff */
	proofHistogram();
	/* calc chi */
	_chiSquare = 0.0;
	for (i=0; i<_time.size(); i++) _chiSquare += std::pow((_transProbData.at(i)-_transProbFit.at(i)),2)/ _transProbFit.at(i);
}

void DwellTimeAnalysis::analyzeJudi(HMMBase &hmm, const JumpDistanceList &judi){
	int i,j;
	_stateNumber = hmm.stateNumber();
	SMLMS::Matrix transMat = hmm.transPDF();
	for (i=0; i<_stateNumber; i++){
		for(j=0; j<_stateNumber; j++){
			if (i!=j){
				_originalState = i;
				_targetState = j;
				_transRate = transMat.at(i,j)/_dt;
				estimateTransTimes(judi);
				transTimesToHist();
				fitHist();
				estimateTransProbFit();
				estimateResiduals();
				estimateChiSquare();
				writeDwellTime();
				plotDwellTime();
				printResult();
			}
		}
	}
}
}/* SMLMS */
