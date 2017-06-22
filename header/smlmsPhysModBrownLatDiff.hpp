/* ######################################################################
* File Name: smlmsPhysModBrownLat.hpp
* Project: SMLMS
* Version: 16.03
* Creation Date: 30.03.2017
* Created By Sebastian Malkusch
* <malkusch@chemie.uni-frankfurt.de>
* Goethe University of Frankfurt
* Physical and Theoretical Chemistry
* Single Molecule Biophysics
###################################################################### */


#ifndef smlmsPhysModBrownLatDiff_hpp
#define smlmsPhysModBrownLatDiff_hpp

#include <string>
#include <vector>
#include "header/smlmsMatrix.hpp"
#include "header/smlmsPhysModBase.hpp"
#include "header/smlmsMicroscope.hpp"

namespace SMLMS{
class PhysicalModelBLD: public PhysicalModelBase{
	private:
		SMLMS::Matrix _paraMat;
		std::vector<double> _paraVect;
		double _contAreaSuperPos;
		std::vector<double> _contArea;
	public:
		/* Constructor */
		PhysicalModelBLD();
		PhysicalModelBLD(const std::vector<double> &xVal, int stateVal, std::string&);
		PhysicalModelBLD(double minVal, double maxVal, int incVal, int stateVal, std::string&);
		/* Destructor */
		~PhysicalModelBLD();
		/* Copy Constructor */
		PhysicalModelBLD(const PhysicalModelBLD&);
		/* elementary functions */
		void setParaMat(SMLMS::Matrix);
		SMLMS::Matrix paraMat();
		std::vector<double> paraVect();
		double contAreaSuperPos();
		std::vector<double> contArea();
		/* init functions */
		void initModelBLD();
		void initParaMat();
		void initParaVect();
		void initContArea();
		/* check functions */
		void checkParaMat();
		void checkParaVect();
		void checkContArea();
		void checkContAreaSuperPos();
		/* print functions */
		void printParaMat();
		void printParaVect();
		void printContAreaSuperPos();
		void printContArea();
		/* read and write functions */
		void writePhysMod();
		void readPhysMod(const std::string&);
		/* contineous normalization  Functions */
		void contIntSuperPos(double &area, const std::vector<double> &pdf);
		void contNormSuperPos(double &area, std::vector<double> &pdf);
		void contIntPdfMat(std::vector<double> &area, const SMLMS::Matrix &pdf);
		void contNormPdfMat(std::vector<double> &area, SMLMS::Matrix &pdf);
		void contNormModel();
		/* model parameter functions */
		void paraMat2paraVect();
		void paraVect2paraMat();
		void calcFitSuperPosFromPara();
		void calcFitMatrixFromPara();
		void updateWeight(const SMLMS::Matrix &pi);
		void updatePi(SMLMS::Matrix &pi);
		void updateFixModelParameter(SMLMS::Microscope &);
		/* fit functions */
		void fitPdfSuperPos();
		void fitPdfMatState(int state, SMLMS::Matrix &pdf);
		void fitPdf(SMLMS::Matrix&);
		void initFit(SMLMS::Matrix&);
		void baumWelchFit(const SMLMS::Matrix &pi, SMLMS::Matrix &pdf);
		void viterbiFit();
}; /* PhysicalModelBLD*/
} /* SMLMS */
#endif /* smlmsPhysModBrownLatDiff_hpp */