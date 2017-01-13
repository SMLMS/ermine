/* ######################################################################
* File Name: matrix.cpp
* Project: SMLMS
* Version: 16.03
* Creation Date: 21.03.2016
* Created By Sebastian Malkusch
* <malkusch@chemie.uni-frankfurt.de>
* Goethe University of Frankfurt
* Physical and Theoretical Chemistry
* Single Molecule Biophysics
###################################################################### */

#include <iostream>
#include <vector>
#include <stdexcept>
#include "header/smlmsMatrix.hpp"

namespace SMLMS{
/* constructor */
Matrix::Matrix (){
	std::cout<<"Matrix constructor called"<<std::endl;
}
Matrix::Matrix(unsigned rowNumber, unsigned columnNumber){
	std::cout<<"Matrix constructor called"<<std::endl;
	setNumberOfRows(rowNumber);
	setNumberOfColumns(columnNumber);
	calcNumberOfElements();
	initMatrix();
}

/* destructor */
Matrix::~Matrix(){
	std::cout<<"Matrix removed from Heap!"<<std::endl;
}

/* copy constructor */
Matrix::Matrix(const Matrix &obj){
	std::cout<<"Matrix copy constructor called"<<std::endl;
	_numberOfRows = obj._numberOfRows;
	_numberOfColumns = obj._numberOfColumns;
	calcNumberOfElements();
	_entries = obj._entries;
}

/* operator overload */
double& Matrix::operator()(unsigned rowNumber, unsigned columnNumber){
	if (rowNumber>=_numberOfRows || columnNumber >= _numberOfColumns){
		throw std::out_of_range("Matrix index out of bounds.");
	}
	return _entries.at(rowNumber*numberOfColumns()+columnNumber);
}

double Matrix::operator()(unsigned rowNumber, unsigned columnNumber)const{
        if (rowNumber>=_numberOfRows || columnNumber >= _numberOfColumns){
                throw std::out_of_range("Matrix index out of bounds.");
        }
        return _entries.at(rowNumber*_numberOfColumns+columnNumber);
}

Matrix& Matrix::operator=(const Matrix &obj){
	if (this != &obj){
		_numberOfRows = obj._numberOfRows;
		_numberOfColumns = obj._numberOfColumns;
		calcNumberOfElements();
		_entries = obj._entries;
	}
	return *this;
}
/* elementary functions */
void Matrix::setNumberOfRows(unsigned rowNumber){
	_numberOfRows = rowNumber;
}

unsigned Matrix::numberOfRows(){
	return _numberOfRows;
}

void Matrix::setNumberOfColumns(unsigned columnNumber){
	_numberOfColumns = columnNumber;
}

unsigned Matrix::numberOfColumns(){
	return _numberOfColumns;
}

unsigned Matrix::numberOfElements(){
	return _numberOfElements;
}

/* special functions */
void Matrix::calcNumberOfElements(){
	_numberOfElements = numberOfRows()*numberOfColumns();
}

void Matrix::initMatrix(){
	_entries.clear();
	for (unsigned i=0; i<numberOfElements(); i++){
		_entries.push_back(0.0);
	}	
}

void Matrix::clearMatrix(){
	setNumberOfRows(0);
	setNumberOfColumns(0);
	calcNumberOfElements();
	initMatrix();
}

void Matrix::at(unsigned rowNumber, unsigned columnNumber, double entry){
	if (rowNumber>=numberOfRows() || columnNumber>=numberOfColumns()){
		throw std::out_of_range(" Matrix entry out of bounds!");
	}
	_entries.at(rowNumber*numberOfColumns()+columnNumber)=entry; 
}

double Matrix::at(unsigned rowNumber, unsigned columnNumber){
	if (rowNumber>=numberOfRows() || columnNumber>=numberOfColumns()){
		throw std::out_of_range(" Matrix entry out of bounds!"); 
        }
	return  _entries.at(rowNumber*numberOfColumns()+columnNumber);
}
}/* SMLMS */

