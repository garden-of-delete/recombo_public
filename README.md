# RECOMBO 

## Synopsis
RECOMBO is the working title of the Topological Molecular Biology Recombination Tool Suite. The suite is made up of a number of different applications that share a common code base are designed to work well together. As a basic overview, these applications cover everything from generating SAPs of fixed topology (MMC), performing reconnection (recomboFromFile), computing homfly polynomials from extended gauss codes (homfly), identifying homfly polynomials with a topology (idknot). RECOMBO also includes an integrated unit testing suite, and a thoroughly experimental z-value finder (zAnalyzer). Each application in this suite is described in detail below.

##Installation
This software suite is intended for use on Linux platforms and may be compatible with most OSX platforms. Install by following these steps:

1. Download and extract this repository.

2. In a shell command line, navigate inside the repository's containing folder. 
3. execute `make all`.
4. Assuming compilation is successful and all unit tests pass, navigate to the /src/bin subdirectory to find the executables. 
5. If compilation is not successful, examine the error messages and ensure your compiler and environment are treating the code appropriately. 

## MMC
MMC is the largest and most complex piece of software in the RECOMBO suite. It uses an extension of the BFACF algorithm (CMC-BFACF) [x] to generate and sample self-avoiding polygons. Because of its complexity, there are a number of options and run-modes avaliable to the user to best meet your simulation needs. 

### Usage
Usage information goes here

## recomboFromFile
recomboFromFile takes a file of CUBE binary formatted polygons, and attempts to perform reconnection on each one according to user-specified reconnection criteria. It is capable of performing reconnection in direct ~~and invertet~~ repeat. 

## homfly
homfly

## idknot

## zAnalyzer

## Other Software

## Contributors (ordered chronologically)
Rob Schaerine, Reuben Brasher, Robert Stolz

## License
All rights reserved. 2017. 

This repository is currently private, and this source code is for internal use during development. A stable branch will be made public at an appropriate time.  
