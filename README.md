# Low-Photon-Generator
An addition to PYTHIA 8.3 that will generate the soft photons from the inner bremsstrahlung predicted via [Low's theorem](https://doi.org/10.1103/PhysRev.110.974) (F. E. Low, Phys. Rev. 110, 974 – Published 15 May 1958).

## How to run

This program can be build by adding it to the exampels folder of pythia and adding the following lines to the makefile
```make
# Low Photon Gen
low_run: $(PYTHIA) low_run.cc low_Gen.cc
	$(CXX) $@.cc low_Gen.cc -o $@ $(CXX_COMMON)
```

and then doing
```sh
make low_run
```

To run, do
```sh
./low_run > low_output.txt
```
And have a look at the output

## How to use with Pythia
The example shows you how to use the Low Photon Generator with Pythia. The concept boils down to this
1) Initialize Pythia
2) Construct the Low Photon Generator
3) Set the parameters of the Low Photon Generator
4) Initialize the Low Photon Generator
5) After an event is generated with pythia, pass the stack to the Low PHoton Generator. It will then add the photons to the stack
