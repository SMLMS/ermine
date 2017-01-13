/* ######################################################################
* File Name: matrix.hpp
* Project: SMLMS
* Version: 16.03
* Creation Date: 21.03.2016
* Created By Sebastian Malkusch
* <malkusch@chemie.uni-frankfurt.de>
* Goethe University of Frankfurt
* Physical and Theoretical Chemistry
* Single Molecule Biophysics
###################################################################### */


#ifndef Matrix_hpp
#define Matrix_hpp
#include <vector>

namespace SMLMS{
class Matrix{
	private:
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
		Matrix& operator=(const Matrix &);
		/* elementary functions */
		void setNumberOfRows(unsigned);
		unsigned numberOfRows();
		void setNumberOfColumns(unsigned);
		unsigned numberOfColumns();
		unsigned numberOfElements();
		/* special functions */
       	 	void calcNumberOfElements();
		void initMatrix();
        	void clearMatrix();
        	void at (unsigned rowNumber, unsigned columnNumber, double entry);
		double at (unsigned rowNumber, unsigned columnNumber);
};/* Matrix */
}/* SMLMS */
#endif /* Matrix_hpp */

