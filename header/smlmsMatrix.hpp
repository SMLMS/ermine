/* ######################################################################
* File Name: smlmsMatrix.hpp
* Project: ermine
* Version: 19.02
* Creation Date: 21.03.2016
* Created By Sebastian Malkusch
* <malkusch@chemie.uni-frankfurt.de>
* Goethe University of Frankfurt
* Physical and Theoretical Chemistry
* Single Molecule Biophysics
###################################################################### */


#ifndef SMLMSMATRIX_HPP
#define SMLMSMATRIX_HPP
#include <vector>

namespace SMLMS{
class Matrix{
	private:		
		// private members
		unsigned _numberOfRows;
		unsigned _numberOfColumns;
		unsigned _numberOfElements;
		std::vector<double> _entries;
	public:
		/* constructor */
		Matrix();
		Matrix(unsigned rowNumber, unsigned columnNumber);
		/* destructor */
		~Matrix();
		/* copy constructor */
		Matrix(const Matrix &);
		/* operator overload */
		double& operator()(unsigned rowNumber, unsigned columnNumber);
		double operator()(unsigned rowNumber, unsigned columnNumber) const;
		Matrix operator+(const Matrix &);
		Matrix operator+(const Matrix &) const;
		Matrix& operator+=(const Matrix &);
		/* elementary functions */
		void setNumberOfRows(unsigned);
		unsigned numberOfRows();
		unsigned numberOfRows() const;
		void setNumberOfColumns(unsigned);
		unsigned numberOfColumns();
		unsigned numberOfColumns() const;
		unsigned numberOfElements();
		unsigned numberOfElements() const;
		double* data(void);
		/* special functions */
       	 	void calcNumberOfElements();
		void initMatrix();
        	void clearMatrix();
        	void at (unsigned rowNumber, unsigned columnNumber, double entry);
		double at (unsigned rowNumber, unsigned columnNumber);
		double at (unsigned rowNumber, unsigned columnNumber) const;
};/* Matrix */
}/* SMLMS */
#endif /* Matrix_hpp */

