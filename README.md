# Low-Photon-Generator
An addition to PYTHIA 8.3 that will generate the soft photons from the inner bremsstrahlung predicted via [Low's theorem](https://doi.org/10.1103/PhysRev.110.974) (F. E. Low, Phys. Rev. 110, 974 – Published 15 May 1958).

## How to use

This program can be run by adding it to the exampels folder of pythia and adding the following lines to the makefile
```make
# Low Photon Gen
low_run: $(PYTHIA) low_run.cc low_Gen.cc
	$(CXX) $@.cc low_Gen.cc -o $@ $(CXX_COMMON)
```
