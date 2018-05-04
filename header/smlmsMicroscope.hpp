/* ######################################################################
* File Name: microscope.hpp
* Project: SMLMS
* Version: 16.02.
* Creation Date: 29.02.2016
* Created By Sebastian Malkusch
* <malkusch@chemie.uni-frankfurt.de>
* Goethe University of Frankfurt
* Physical and Theoretical Chemistry
* Single Molecule Biophysics
###################################################################### */


#ifndef Microscope_hpp
#define Microscope_hpp
#include <string>

namespace SMLMS{
class Microscope{
	private:
		double _pxlSize;
		double _intTime;
		double _locPrec;
	public:
		/* Constructor */
		Microscope();
		Microscope(double pxlSize, double intTime, double locPrec);
		/* Destructor */
		~Microscope();
		/* Copy Constructor*/
		Microscope(const Microscope &);
		/* Elementary functions */
		void setPxlSize(double);
		double pxlSize();
		double pxlSize()const;
		void setIntTime(double);
		double intTime();
		double intTime()const;
		void setLocPrec(double);
		double locPrec();
		double locPrec()const;
		/* Sprecial Functions */
		void clearMicroscope();
		void loadMicroscope(const std::string &name);
		void saveMicroscope(const std::string &folderName);
		
};/* Microscope */
}/* SMLMS */
#endif /* Name_hpp */

