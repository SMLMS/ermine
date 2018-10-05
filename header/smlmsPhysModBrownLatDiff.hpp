/* ######################################################################
* File Name: smlmsPhysModBrownLat.hpp
* Project: SMLMS
* Version: 18.09
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
		SMLMS::Microscope _microscope;
		SMLMS::Matrix _paraMat;
		std::vector<double> _paraVect;
		double _contAreaSuperPos;
		std::vector<double> _contArea;
	public:
		/* Constructor */
		PhysicalModelBLD();
		PhysicalModelBLD(const std::vector<double> &xVal, int stateVal, SMLMS::Microscope microscope);
		PhysicalModelBLD(double minVal, double maxVal, int incVal, int stateVal, SMLMS::Microscope microscope);
		/* Destructor */
		~PhysicalModelBLD();
		/* Copy Constructor */
		PhysicalModelBLD(const PhysicalModelBLD&);
		/* elementary functions */
		void setMicroscope(SMLMS::Microscope);
		SMLMS::Microscope microscope();
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
		void checkMicroscope();
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
		void writePhysMod(const std::string &folderName);
		void readPhysMod(const std::string &folderName);
		/* numerical normalization  Functions */
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
		void updatePdfWeight(void);
		void updateWeight(const SMLMS::Matrix &pi);
		void updatePi(SMLMS::Matrix &pi);
		void fixDiffusionCoefficients(void);
		void releaseDiffusionCoefficients(void);
		/* fit functions */
		void fitPdfSuperPos();
		void fitPdfMatState(int state, SMLMS::Matrix &pdf);
		void fitPdf(SMLMS::Matrix&);
		void initFit(SMLMS::Matrix&);
		void baumWelchFit(const SMLMS::Matrix &pi, SMLMS::Matrix &pdf);
}; /* PhysicalModelBLD*/
} /* SMLMS */
#endif /* smlmsPhysModBrownLatDiff_hpp */
