#include "low_Gen.h"

GeneratorLow::GeneratorLow(){ 
	// Default values. After changing them with the set functions, call Init();
	
	// Photon properties for PYTHIA
	// These are the photons coming from Low's theorem
	photID = 22;
	photStatus = 201;
	photMother1 = 0;
	photMother2 = 0;
	photDaughter1 = 0;
	photDaughter2 = 0;
	photCol = 0;
	photACol = 0;
	photM = 0.;
	photScale = 0.;
	photPol = 9.;

	// Energy range for the emitted photons
	Ebeg = 0.001; // GeV
	Eend = 1; // GeV
	probePhotE = 0.01; // GeV
	nEBins = 1000;
	
	// Set parameters for etaPHiGamma2DArray
	nEtaBins = 100;
	nPhiBins = 100;
	etaBegin = -5.;
	etaEnd = -3.;
	phiBegin = 0.;
	phiEnd = 2. * M_PI;

	verbose = false;

	seed = 0;

	std::cout << "Low Photon Generator Constructed" << std::endl;

}

GeneratorLow::~GeneratorLow(){
	for(int i = 0; i < nEtaBins; i++){
		delete [] etaPhiGamma2DArray[i];
	}
	delete [] etaPhiGamma2DArray;
	delete [] energyArray;
}


void GeneratorLow::Init(){
	// Calculating bin widths
	EBW = (Eend - Ebeg)/ nEBins;
	etaBW = (etaEnd - etaBegin) / (double)nEtaBins;
	phiBW = (phiEnd - phiBegin) / (double)nPhiBins;

	mNormalizationFactor = 1./137. / (4.*M_PI*M_PI);

	// Make array for etaPhiGamma
	etaPhiGamma2DArray = new double*[nEtaBins];
	for(int i = 0; i < nEtaBins; i++){
		etaPhiGamma2DArray[i] = new double[nPhiBins];
	}

	// Fill EnergyDistribution 1/E_gamma to draw photons from
	energyArray = new double[nEBins];
	for(int i = 0; i < nEBins; i++){
		double EBinCenter;
		getBinCenter1D(i, EBinCenter);
		energyArray[i] = 1./EBinCenter;
	}

	// Set Seed of the random number generator
	rat = new Rndm(0);

	std::cout << "Low Generator Initialized" << std::endl;

	return;
}

void GeneratorLow::getBinCenter1D(int EBin, double& EBinCenter){
	EBinCenter = Ebeg + EBW/2. + EBW * EBin;
}

void GeneratorLow::getBinCenters2D(int etaBin, int phiBin, double& etaBinCenter, double& phiBinCenter){
	etaBinCenter = etaBegin + etaBW/2. + etaBW * etaBin;
	phiBinCenter = phiBegin + phiBW/2. + phiBW * phiBin;
	return;
}

void GeneratorLow::getBinRange1D(int EBin, double& EBinBegin, double& EBinEnd){
	EBinBegin = Ebeg + EBW * EBin;
	EBinEnd = Ebeg + EBW * (EBin + 1);
	return;
}

void GeneratorLow::getBinRanges2D(int etaBin, int phiBin, double& etaBinBegin, double& etaBinEnd, double& phiBinBegin, double& phiBinEnd){
	etaBinBegin = etaBegin + etaBW * etaBin;
	etaBinEnd = etaBegin + etaBW * (etaBin + 1);
	phiBinBegin = phiBegin + phiBW * phiBin;
	phiBinEnd = phiBegin + phiBW * (phiBin + 1);
	return;
}

// Picks 1 random value from a 1D grid. Principle can be found here
// https://root.cern.ch/doc/master/TH2_8cxx_source.html#l01153
void GeneratorLow::getRandom(double* array, int& EBin){
	int nBins = nEBins;

	// Compute integral
	double* fIntegral = new double[nBins + 1];
	int iBin = 0; fIntegral[iBin] = 0.;
	for(int ix = 0; ix < nBins; ix++){
		++iBin;
		double y = array[ix];
		fIntegral[iBin] = fIntegral[iBin - 1] + y;
	}
	// Normalize integral to 1
	for(int i = 1; i <= nBins; ++i){ fIntegral[i] /= fIntegral[nBins]; }
	
	// Perform a binary search over the integral to match r1
	double r1 = rat->flat();
	auto pind = std::lower_bound(fIntegral, fIntegral + nBins, r1);
	int position;
	if ( (pind != fIntegral + nBins) && (*pind == r1) ){
		position = pind - fIntegral;
	} else {
		position = pind - fIntegral - 1;
	}

	// Calculate the positions of the bins (1D is trivial)
	EBin = position;
		
	delete [] fIntegral;
	return;
}

// Picks 2 random values from a 2D grid. Principle can be found here
// https://root.cern.ch/doc/master/TH2_8cxx_source.html#l01153
// It unwraps the 2D array into a 1D array, calculates the integral, normalizes it to 1
// and compares this to a random value taken from a uniform distribution
// to calculate the bin positions
void GeneratorLow::getRandom2(double** array, int& etaBin, int& phiBin){
	int nBins = nEtaBins * nPhiBins;
	
	// Compute integral
	double* fIntegral = new double[nBins + 1];
	int iBin = 0; fIntegral[iBin] = 0.;
	for(int ix = 0; ix < nEtaBins; ix++){
		for(int iy = 0; iy < nPhiBins; iy++){
			++iBin;
			double y = array[ix][iy];
			fIntegral[iBin] = fIntegral[iBin - 1] + y;
		}
	}
	// Normalize integral to 1
	for(int i = 1; i <= nBins; ++i){ fIntegral[i] /= fIntegral[nBins]; }

	// Perform a binary search over the integral to match r1
	double r1 = rat->flat();
	auto pind = std::lower_bound(fIntegral, fIntegral + nBins, r1);
	int position;
	if ( (pind != fIntegral + nBins) && (*pind == r1) ){
		position = pind - fIntegral;
	} else {
		position = pind - fIntegral - 1;
	}

	// Calculate the positions of the bins
	phiBin = position/nEtaBins;
	etaBin = position - nEtaBins*phiBin;
		
	delete [] fIntegral;
	return;
}

void GeneratorLow::rejectSample(int targetFunc, double scale, double* range,
				double* sampledVars, int freeParams,
				int iterations = 10000){
	
	// This used to be a template, but I could not get it to work, so now there's this
	// int targetFunc: 0: evalLow, 1: evalEDist
	
	// It is up to the user to choose the right scaling. The scale should be equal to or 
	// higher than the highest point of the target function
	
	// The parameters of range dictiate the area where you would like to sample.
	// They should be like so
	// x0 = range[0], x1 = range[1]
	// y0 = range[2], y1 = range[3]
	// z0 = range[4], z1 = range[5]
	// ...

	// freeParams is the number of parameters each the target function uses
	// Even though a function could be constant in one dimension
	// still take this into account when constructing the function.

	// sampledVars should be an empty array with length equivalent 
	// to the number of free parameters

	// With iterations the user can dictate how many times should be tried for a 
	// random numbers. If the limit is exceeded, consider using a different scaling
	// or a different envelop function
	
	double *testVar = new double[freeParams];

	bool found = false;
	int counter = 0;

	// Will loop until it has found a random number or breached the counter
	while(!found && counter < iterations){
		for(int i = 0; i < freeParams; i++){
			testVar[i] = range[2*i] + (range[2*i + 1] - range[2*i]) * rat->flat();
		}

		double targetVal = 0;
		// The original intention was that rejectSample was a template
		// which would take the eval functions as template parameter (targetfunc)
		// but this did work out. Feel free to make a pull request to fix it like that
		if(targetFunc == 0){
			targetVal = evalLow(testVar);
		} else if(targetFunc == 1){
			targetVal = evalEDist(testVar);
		}

		double u = rat->flat();
		if(u < (targetVal / scale)){
			found = true;
		}
		counter++;
	}

	for(int i = 0; i < freeParams; i++){
		sampledVars[i] = testVar[i];
	}

	delete [] testVar;

	return;

}

double GeneratorLow::evalLow(double* x){
	// x[0] : eta
 	// x[1] : phi
	double PhotPx, PhotPy, PhotPz, PhotPtot;

	Pythia8::Vec4 MainP;
	double theta = 2. * atan(exp(-x[0]));
	PhotPx = probePhotE * cos(x[1]) * sin(theta);
	PhotPy = probePhotE * sin(x[1]) * sin(theta);
	PhotPz = probePhotE * cos(theta);
	PhotPtot = sqrt(PhotPx*PhotPx + PhotPy*PhotPy + PhotPz*PhotPz);
	Pythia8::Vec4 pgamma(PhotPx, PhotPy, PhotPz, PhotPtot);
	MainP.reset();
	std::vector<int>::iterator iter;

	for(iter = chargedPartList.begin(); iter != chargedPartList.end(); iter++){
		double prefactorfast = mPythia->event[*iter].charge();
		if(*iter <= 2){ prefactorfast *= -1.; }
		MainP = MainP + prefactorfast * mPythia->event[*iter].p() / (mPythia->event[*iter].p()*pgamma);
	}

	double prod = MainP * MainP;
	return -prod * PhotPtot * mNormalizationFactor * 2. * cos(theta/2.) * sin(theta/2.) * sin(theta);
}

double GeneratorLow::evalEDist(double* x){
	return 1./x[0];
}

int GeneratorLow::poisson(double mean){
	// Prob(N) = exp(-mean)*mean^N/Factorial(N)
	int n;
	if (mean <= 0) return 0;
	double expmean = exp(-mean);
	double pir = 1;
	n = -1;
	while(1) {
		n++;
		pir *= rat->flat();
		if (pir <= expmean) break;
	}
	return n;
}

void GeneratorLow::generatePhotons(Pythia* pythia){
	cout.precision(10);

	setPythiaStackPtr(pythia);

	int ndp = mPythia->event.size(); // Number of particles
	if(verbose){
		std::cout << "Number of particles: " << ndp << std::endl;
	}

	// Create the particle list of particles to be considered. Doing this once speeds up the process
	for(int ip = 1; ip < ndp; ip++){
		if((abs(mPythia->event[ip].charge())>0.001 && abs(mPythia->event[ip].status())/10==8) || ip<=2){ // 
			chargedPartList.push_back(ip);
		}
	}

	double dNPhot = 0.; 	// d number of photons
	int NPhot = 0; 		// Final number of photons
	double photEta; 		// Eta of the emitted photon
	double photTheta;		// Theta of the emitted photon
	double photPhi; 		// Phi of the emitted photon
	double thisPhotE; 	// Energy of the emitted photon

	// Via Low's theorem, calculate how many photons are emitted
	double integralEtaPhiGamma = 0;
	for(int etaBin = 0; etaBin < nEtaBins; etaBin++){
		for(int phiBin = 0; phiBin < nPhiBins; phiBin++){
			double etaBinCenter, phiBinCenter;
			getBinCenters2D(etaBin, phiBin, etaBinCenter, phiBinCenter);
			double toEval[2] = {etaBinCenter, phiBinCenter};
			double prod = evalLow(toEval);
			etaPhiGamma2DArray[etaBin][phiBin] = prod;
			integralEtaPhiGamma += etaBW * phiBW * prod;
		}
	}

	// Number of photons of this (Px, Py, Pz) photon energy times delta phi delta eta
	dNPhot = integralEtaPhiGamma * probePhotE * (log(Eend) - log(Ebeg));
	NPhot = poisson(dNPhot);
	if(verbose){
		std::cout << "Expected number of photons dNPhot: " << dNPhot << std::endl;
		std::cout << "Number of photons after Poisson " << NPhot << std::endl;
	}

	// Add photons
	for(int i = 0; i < NPhot; i++){
		// Pick a random value from the Energy distribution for the energy of the photon
		int EBin = -1;

		getRandom(energyArray, EBin); // Get a random bin from the distribution

		double EBinBegin, EBinEnd;
		getBinRange1D(EBin, EBinBegin, EBinEnd); // Get the range of that bin
		double rangeE[2] = {EBinBegin, EBinEnd};
		double EtoEval[1] = {EBinBegin};
		double EMaxVal = evalEDist(EtoEval); // Calculate the max value in that bin
		double sampledE[1] = {0.};

		rejectSample(1, EMaxVal, rangeE, sampledE, 1); // Perform rejection sampling over this range

		// Pick a random value for eta and phi based on the etaPhiGamma2DArray
		int etaBin, phiBin;

		getRandom2(etaPhiGamma2DArray, etaBin, phiBin); // Get a random phi and eta bin based on the etaPhiGammaArray

		double etaBinCenter, phiBinCenter;
		getBinCenters2D(etaBin, phiBin, etaBinCenter, phiBinCenter); // Get the values for the bin center
		double etaBinBegin, etaBinEnd, phiBinBegin, phiBinEnd;
		getBinRanges2D(etaBin, phiBin, etaBinBegin, etaBinEnd, phiBinBegin, phiBinEnd); // Get the ranges this bin spans
		double rangeEtaPhi[4] = {etaBinBegin, etaBinEnd, phiBinBegin, phiBinEnd};
		double etaPhiToEval[2] = {etaBinCenter, phiBinCenter};
		double etaPhiMaxVal = evalLow(etaPhiToEval); // Evaluate Low at the Bin centers 
		double sampledEtaPhi[2] = {0., 0.};

		rejectSample(0, 3. * etaPhiMaxVal, rangeEtaPhi, sampledEtaPhi, 2); // Perform rejection sampling over this bin

		photTheta = 2 * atan(exp(-sampledEtaPhi[0]));
		photEta = sampledEtaPhi[0];
		photPhi = sampledEtaPhi[1];

		// Calculate px, py, pz of this photon
		thisPhotE = sampledE[0];
		double thisPhotPx = thisPhotE * cos(photPhi) * sin(photTheta);
		double thisPhotPy = thisPhotE * sin(photPhi) * sin(photTheta);
		double thisPhotPz = thisPhotE * cos(photTheta);
		
		// Add particle to the event (or stack)
		mPythia->event.append(photID, photStatus, photMother1, photMother2, photDaughter1, photDaughter2, 
							photCol, photACol, thisPhotPx, thisPhotPy, thisPhotPz, thisPhotE, photM,
							photScale, photPol);
		
		if(verbose){
			std::cout << "For photon " << i << std::endl;
			std::cout << "Energy: " << thisPhotE << " GeV" << std::endl;
			std::cout << "Px: " << 1000 * thisPhotPx << " MeV/c" << std::endl;
			std::cout << "Py: " << 1000 * thisPhotPy << " MeV/c" << std::endl;
			std::cout << "Pz: " << 1000 * thisPhotPz << " MeV/c" << std::endl;
			std::cout << "Phi: " << photPhi << " rad" << std::endl;
			std::cout << "Eta: " << photEta << " pseurap" << std::endl;
			std::cout << "pT: " << 1000. * thisPhotE / cosh(photEta) << " MeV/c" << std::endl;
			std::cout << " " << std::endl;
		}
	}

	// Prepare all the containers for the next event
	for(int i = 0; i < nEtaBins; i++){
		for(int j = 0; j < nPhiBins; j++){
			etaPhiGamma2DArray[i][j] = 0.;
		}
	}
	chargedPartList.clear();

	return;
}
