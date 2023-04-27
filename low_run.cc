// Example program how to run the Low Photon Generator.
// It takes the event and calculates the amount of soft photons generated via Low's theorem
// You initialize pythia, initialize the low photon generator, run an event, pass the event stack to the low_Gen class
// and the class calculates how many Low photons are generated
// Since these are soft photons, they do not influence the other particles of the event
// and are therefore added directly to the stack without altering the particles already present


#include "Pythia8/Pythia.h"
#include "low_Gen.h"
using namespace Pythia8;

int main(){
	Pythia pythia;
	pythia.readFile("low.cmnd");

	int nEvent = pythia.mode("Main:numberOfEvents");
	int nAbort = pythia.mode("Main:timesAllowErrors");

	pythia.init();

	int iAbort = 0;
	GeneratorLow lowGen;
	// Set all your parameters for the Low Generator
	lowGen.setNEtaBins(100);
	lowGen.setNPhiBins(100);
	lowGen.setVerbose(true);
	// Initialize the Low Photon Generator
	lowGen.Init();

	for (int iEvent = 0; iEvent < nEvent; ++iEvent) {
		// Generate events. Quit if too many failures.
		if (!pythia.next()) {
			if (++iAbort < nAbort) continue;
			cout << " Event generation aborted prematurely, owing to error!\n";
			break;
		}
		std::cout << "-------------------------------- " << std::endl;
		std::cout << "Generating event " << iEvent << std::endl;
		
		int ndp = pythia.event.size();

		// Now the pointer to the event stack will get passed to the low_Gen class and photons will get added
		lowGen.generatePhotons(&pythia);
		
		// Check out the photons
		int ndp_after_Low = pythia.event.size();
		int photon_counter = 0;
		for(int i = ndp; i < ndp_after_Low; i++){
			std::cout << "For photon " << photon_counter << std::endl;
			std::cout << "ParticlePDG : " << pythia.event[i].id() << std::endl;
			std::cout << "Px : " << 1000 * pythia.event[i].px() << " MeV/c" << std::endl;
			std::cout << "Py : " << 1000 * pythia.event[i].py() << " MeV/c" << std::endl;
			std::cout << "Pz : " << 1000 * pythia.event[i].pz() << " MeV/c" << std::endl;
			std::cout << "E : " << pythia.event[i].e() << " GeV" << std::endl;
			if(pythia.event[i].phi() < 0){
				std::cout << "Phi : " << 2 * M_PI + pythia.event[i].phi() << " rad" << std::endl;
			} else {
				std::cout << "Phi : " << pythia.event[i].phi() << " rad" << std::endl;
			}
			std::cout << "Eta : " << pythia.event[i].eta() << " pseurap" << std::endl;
			std::cout << "Statuscode: " << pythia.event[i].status() << std::endl;
			std::cout << "" << std::endl;
			photon_counter++;
		}

	}

	std::cout << "Generation finished" << std::endl;
	return 0;
}
