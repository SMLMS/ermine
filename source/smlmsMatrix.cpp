/* ######################################################################
* File Name: smlmMatrix.cpp
* Project: ermine
* Version: 19.02
* Creation Date: 21.03.2016
* Created By Sebastian Malkusch
* <malkusch@chemie.uni-frankfurt.de>
* Goethe University of Frankfurt
* Physical and Theoretical Chemistry
* Single Molecule Biophysics
###################################################################### */

#include <iostream>
#include <string>
#include <sstream>
#include <vector>
#include <stdexcept>
#include "header/smlmsMatrix.hpp"

namespace SMLMS{
/* constructor */
Matrix::Matrix (){
	//std::cout<<"Matrix constructor called"<<std::endl;
}
Matrix::Matrix(unsigned rowNumber, unsigned columnNumber){
	//std::cout<<"Matrix constructor called"<<std::endl;
	setNumberOfRows(rowNumber);
	setNumberOfColumns(columnNumber);
	calcNumberOfElements();
	initMatrix();
}

/* destructor */
Matrix::~Matrix(){
	//std::cout<<"Matrix removed from Heap!"<<std::endl;
}

/* copy constructor */
Matrix::Matrix(const Matrix &obj){
	//std::cout<<"Matrix copy constructor called"<<std::endl;
	_numberOfRows = obj._numberOfRows;
	_numberOfColumns = obj._numberOfColumns;
	calcNumberOfElements();
	_entries = obj._entries;
}

/* operator overload */
double& Matrix::operator()(unsigned rowNumber, unsigned columnNumber){
	if (rowNumber>=_numberOfRows || columnNumber >= _numberOfColumns){
		std::stringstream errorMessage;
		errorMessage<<"Matrix out of range error: Matrix index out of bounds!"<<std::endl;
		std::out_of_range matrixError(errorMessage.str());
		throw matrixError;
	}
	return _entries.at(rowNumber*numberOfColumns()+columnNumber);
}

double Matrix::operator()(unsigned rowNumber, unsigned columnNumber)const{
        if (rowNumber>=_numberOfRows || columnNumber >= _numberOfColumns){
		std::stringstream errorMessage;
		errorMessage<<"Matrix out of range error: Matrix index out of bounds!"<<std::endl;
		std::out_of_range matrixError(errorMessage.str());		
                throw matrixError;
        }
        return _entries.at(rowNumber*_numberOfColumns+columnNumber);
}

Matrix Matrix::operator+(const Matrix &obj){
        if (obj.numberOfRows()!=_numberOfRows || obj.numberOfColumns() != _numberOfColumns){
		std::stringstream errorMessage;
		errorMessage<<"Matrix out of range error: Matrix dimensions need to be the same!"<<std::endl;
		std::out_of_range matrixError(errorMessage.str());		
                throw matrixError;
        }
	SMLMS::Matrix obj2(_numberOfRows, _numberOfColumns);
	for (unsigned i=0; i<_numberOfRows; i++){
		for (unsigned j=0; j<_numberOfColumns; j++) obj2.at(i,j,at(i,j)+obj.at(i,j));
	}
	return obj2;
}

Matrix Matrix::operator+(const Matrix &obj) const{
        if (obj.numberOfRows()!=_numberOfRows || obj.numberOfColumns() != _numberOfColumns){
		std::stringstream errorMessage;
		errorMessage<<"Matrix out of range error: Matrix dimensions need to be the same!"<<std::endl;
		std::out_of_range matrixError(errorMessage.str());		
                throw matrixError;
        }
	SMLMS::Matrix obj2(_numberOfRows, _numberOfColumns);
	for (unsigned i=0; i<_numberOfRows; i++){
		for (unsigned j=0; j<_numberOfColumns; j++) obj2.at(i,j,at(i,j)+obj.at(i,j));
	}
	return obj2;
}

Matrix& Matrix::operator+=(const Matrix &obj){
        if (obj.numberOfRows()!=_numberOfRows || obj.numberOfColumns() != _numberOfColumns){
		std::stringstream errorMessage;
		errorMessage<<"Matrix out of range error: Matrix dimensions need to be the same!"<<std::endl;
		std::out_of_range matrixError(errorMessage.str());		
                throw matrixError;
        }
	for (unsigned i=0; i<_numberOfRows; i++){
		for (unsigned j=0; j<_numberOfColumns; j++) at(i,j,(at(i,j))+obj.at(i,j));
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

unsigned Matrix::numberOfRows()const{
	return _numberOfRows;
}

void Matrix::setNumberOfColumns(unsigned columnNumber){
	_numberOfColumns = columnNumber;
}

unsigned Matrix::numberOfColumns(){
	return _numberOfColumns;
}

unsigned Matrix::numberOfColumns()const{
	return _numberOfColumns;
}

unsigned Matrix::numberOfElements(){
	return _numberOfElements;
}

unsigned Matrix::numberOfElements()const{
	return _numberOfElements;
}

double* Matrix::data(void){
	return _entries.data(); 
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
		std::stringstream errorMessage;
		errorMessage<<"Matrix out of range error: Matrix index outof bounds!"<<std::endl;
		std::out_of_range matrixError(errorMessage.str());
		throw matrixError;
	}
	_entries.at(rowNumber*numberOfColumns()+columnNumber)=entry; 
}

double Matrix::at(unsigned rowNumber, unsigned columnNumber){
	if (rowNumber>=numberOfRows() || columnNumber>=numberOfColumns()){
		std::stringstream errorMessage;
		errorMessage<<"Matrix out of range error: Matrix index outof bounds!"<<std::endl;
		std::out_of_range matrixError(errorMessage.str());
		throw matrixError; 
        }
	return  _entries.at(rowNumber*numberOfColumns()+columnNumber);
}

double Matrix::at(unsigned rowNumber, unsigned columnNumber)const{
	if (rowNumber>=numberOfRows() || columnNumber>=numberOfColumns()){
		std::stringstream errorMessage;
		errorMessage<<"Matrix out of range error: Matrix index outof bounds!"<<std::endl;
		std::out_of_range matrixError(errorMessage.str());
		throw matrixError; 
        }
	return  _entries.at(rowNumber*numberOfColumns()+columnNumber);
}

}/* SMLMS */

