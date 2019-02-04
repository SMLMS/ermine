/* ######################################################################
* File Name: ermineDelegate.hpp
* Project: ermine
* Version: 19.02
* Creation Date: 2016
* Created By Sebastian Malkusch
* <malkusch@chemie.uni-frankfurt.de>
* Goethe University of Frankfurt
* Physical and Theoretical Chemistry
* Single Molecule Biophysics
###################################################################### */


#ifndef ERMINEDELEGATE_HPP
#define ERMINEDELEGATE_HPP

#include <boost/program_options.hpp>
#include "header/ermineParser.hpp"
#include "header/ermineFilenames.hpp"
#include "header/smlmsMicroscope.hpp"
#include "header/smlmsMolecules.hpp"
#include "header/smlmsJudi.hpp"
#include "header/smlmsPhysModBrownLatDiff.hpp"
#include "header/smlmsHmmSequence.hpp"


namespace po=boost::program_options;
namespace SMLMS{
class Delegate{
	private:
		unsigned int _errors;
		SMLMS::ErmineParser _eVar;
		SMLMS::FileNames _fileNames;
		SMLMS::Microscope _microscope;
		SMLMS::MoleculeList _molList;
		SMLMS::JumpDistanceList _judi;
		SMLMS::PhysicalModelBLD _physMod;
		SMLMS::HMMSequence _hmm;
	public:
		/* Constructor */
		Delegate();
		Delegate(po::variables_map &);
		/* Destructor */
		~Delegate();
		/* copy Constructor */
		Delegate(const SMLMS::Delegate &);
		/* init functions */
		int initDelegate(po::variables_map &);
		/* load functions */
		int loadFileNames(void);
		int loadMicroscope(void);
		int loadTrcList(void);
		int loadMoleculeList(void);
		int loadJumpDistanceList(void);
		int loadPhysMod(void);
		int loadHmm(void);
		/* create destination folder */
		int createFolder(void);
		/* write functions */
		int writeParser(void);
		int writeMicroscope(void);
		int writeMoleculeList(void);
		int writeJumpDistanceList(void);
		int writePhysMod(void);
		int writeHmm(void);
		/* proof functions */
		/* run algorithm functions */
		int run(void);
		int runBatchAlgorithm(void);
		int runMol2JudiAlgorithm(void);
		int runInitPhysModAlgorithm(void);
		int runFitPhysModAlgorithm(void);
		int runInitHmmAlgorithm(void);
		int runInitHmmWithPhysModAlgorithm(void);
		int runEvaluateAlgorithm(void);
		int runTrainAlgorithm(void);
		int runTrainWithPhysModAlgorithm(void);
		int runBestPathAlgorithm(void);
		int runTransferStatesAlgorithm(void);
		int runArchiveAlgorithm(void);
		int runExtractAlgorithm(void);
		int runWholeCellAnalysis(void);
		int runSimulateAlgorithm(void);
}; /* Delegate*/
} /* SMLMS */

#endif /* ERMINEDELEGATE_HPP */

