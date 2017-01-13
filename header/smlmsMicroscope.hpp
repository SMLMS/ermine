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
		void setIntTime(double);
		double intTime();
		void setLocPrec(double);
		double locPrec();
		/* Sprecial Functions */
		void clearMicroscope();
		void loadMicroscope(std::string name);
		void saveMicroscope(std::string name);
		
};/* Microscope */
}/* SMLMS */
#endif /* Name_hpp */

