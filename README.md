# RICH Detector for HELIX

HELIX (High Energy Light Isotope eXperiment) is a balloon-borne experiment with the goal of measuring chemical and isotopic abundance of cosmic ray nuclei. In particular, the experiment aims to measure the 10Be/9Be ratio, a measurement needed to constrain cosmic-ray propagation models. The detector consists of a mass spectrometer, and time-of-flight counters and a ring-imaging Cherenkov detector (RICH) to measure particle velocities. The proximity-focused RICH consists of a radiator made of aerogel tiles (refractive index approximately 1.15) and a detector plane. 

## Required packages:
- cmake
- geant4.10.04

## Geant4 - an Object-oriented Toolkit for Simulating High Energy Physics

Geant4 is a toolkit for the simulation of the passage of particles through matter. Its areas of application include high energy, nuclear and accelerator physics, as well as studies in medical and space science.

### Geant4 Installation 

[Geant-4 Download](https://geant4.web.cern.ch/support/download)

[Geant-4 Installation Guide ](https://indico.cern.ch/event/679723/contributions/2792554/attachments/1559217/2453759/Geant4InstallationGuide.pdf)

## Usage
To run the simulation, change the current directory to the build directory. 
```bash
cd HR-build
```
To build the project and link the required data files, run the lines below in the build directory. 
```bash
cmake -DGeant4_DIR=../../../geant4-build ../HR
make
source ../../../geant4-install/bin/geant4.sh
```
To execute the simulation, run the line below in the build directory. 
```bash
./rich
```
Simulation parameters can changed in the [Simulation Constants file](../HR/src/SimulationConstants.cpp). 

After executing, the TSV output file (table of detected particles) can be found in the build directory. Each run can be parsed and loaded into a Pandas Dataframe for visual analysis. 

## License
[MIT](https://choosealicense.com/licenses/mit/)