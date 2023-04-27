#ifndef GENERATELOW_H
#define GENERATELOW_H

#include "Pythia8/Pythia.h"
#include <vector>
#include <iostream>
#include <math.h>

#include <chrono>

using namespace Pythia8;

class GeneratorLow
{
	public:
		GeneratorLow();
		~GeneratorLow();

		// Set all the parameters to your liking
		void setNEtaBins(int n){nEtaBins = n; std::cout << "Number of eta bins set" << std::endl;}
		void setNPhiBins(int n){nPhiBins = n; std::cout << "Number of phi bins set" << std::endl;}
		void setEtaRange(double beg, double end){
			etaBegin = beg;
			etaEnd = end;
		}
		void setPhiRange(double beg, double end){
			phiBegin = beg;
			phiEnd = end;
		}
		
		void setNEBins(int n){nEBins = n;}
		void setERange(double beg, double end){
			Ebeg = beg;
			Eend = end;
		}

		void setProbePhotE(double probE){probePhotE = probE;}
		void setPhotStatus(int stat){photStatus = stat;}
		
		void setVerbose(bool verb){verbose = verb;}
		void setSeed(int a){seed = a;}

		// Then call Init
		void Init();

		// Then call generate photons, which adds the photons to your stack
		void generatePhotons(Pythia* pythia);
	
	private:
		void setPythiaStackPtr(Pythia* pythia){mPythia = pythia;}
		void getBinCenter1D(int EBin, double& EBinCenter);
		void getBinCenters2D(int etaBin, int phiBin, double& etaBinCenter, double& phiBinCenter);
		void getBinRange1D(int EBin, double& EBinBegin, double& EBinEnd);
		void getBinRanges2D(int etaBin, int phiBin, 
							double& etaBinBegin, double& etaBinEnd, 
							double& phiBinBegin, double& phiBinEnd);
		void getRandom(double* array, int& EBin);
		void getRandom2(double** array, int& etaBin, int& phiBin);

		void rejectSample(int targetFunc, double scale, double* range,
						double* sampledVars, int freeParams,
						int iterations);
		int poisson(double mean);
		// Define the friend functions
		// evalLow(GeneratorLow* GLow, double* x);
		// evalEDist(GeneratorLow* GLow, double* x);
		double evalLow(double* x);
		double evalEDist(double* x);

		double mNormalizationFactor;
		Rndm* rat; // Random number generator
		double Ebeg; // GeV
		double Eend; // GeV
		double EBW; // Energy Bin Width
		int nEBins;
		double probePhotE; // GeV
		double* energyArray;
		bool verbose;

		// Photon properties - coming from Low's theorem
		int photID; // PDG code of the photon
		int photStatus; // Statuses from 201 and above are free and can be user defined 
		int photMother1; // The "mother" is the event as a whole
		int photMother2; // Idem
		int photDaughter1; // No daughters (so far)
		int photDaughter2; // Idem
		int photCol; // 
		int photACol;
		double photM;
		double photScale;
		double photPol;

		// Container and paramters for the grid
		double** etaPhiGamma2DArray;
		int nEtaBins, nPhiBins;
		double etaBegin, etaEnd, etaBW;
		double phiBegin, phiEnd, phiBW;

		Pythia* mPythia;

		bool increasePhotons; // Artificially increase the number of photons emitted
		int seed; // Seed of the random number generator

		std::vector<int> chargedPartList; // List of all the charged particles

};

#endif