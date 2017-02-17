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
MMC is the largest and most complex piece of software in the RECOMBO suite. It uses an extension of the BFACF algorithm (CMC-BFACF) [Stolz et. al., 2017] to generate and sample self-avoiding polygons. Because of its complexity, there are a number of options and run-modes avaliable to the user to best meet your simulation needs. 

#### Usage
Usage information goes here

## recomboFromFile
recomboFromFile takes an input file or stream of CUBE binary formatted polygons, and attempts to perform reconnection on each one according to user-specified reconnection criteria. It is capable of performing reconnection in direct ~~and invertet~~ repeat. 

## homfly
Homfly takes as an input a file of extended gauss codes (.egc), and outputs the HOMFLY-PT polynomial [x] for each one. 

#### Usage
`homfly` for interactive usage.

`homfly input_file output_file` 

`homfly -- |` or ` homfly -- > output_file` for use in a stream.

## idknot
Idknot takes a input file or stream of homfly polynomials and associates each one with a topology, or 'unknown' if unrecognized. 

#### Usage
Running `idknot` with no other arguments will print out detailed usage information. Some common usage scenarios:

`idknot input_file output_file -nz -u -P` to specify input output files.

`idknot input_file output_file -k` to print conformation identifications line-by-line. 

`idknot -- output_file -nz -u -P` or `idknot -- -nz -u -P > output_file` substitutes the input file for a data stream. 

## zAnalyzer
N/A

## Other Software
There are other useful applications that are not included as part of RECOMBO: 

Xinger attempts to remove reidermeister crossings at the extended gauss code level, improving the performance of Homfly. It is avaliable in limited circumstances upon request by emailing a request to: mariel@math.ucdavis.edu. 

Seqconvert is a generally useful application that efficiently converts files filled with self-avoiding polygons between different formats. For example, it can be used to convert CUBE binary formatted polygons to ascii formatted files. It can also perform and output some useful computations for all polygons in a file, such as polygon length and writhe. It is also avaliable in limited circumstances upon request by emailing a request to: mariel@math.ucdavis.edu.

Knotplot is a wonderful knot computation and visualization software developed and maintained by Rob Scharein. In the context of reconnection studies, it is useful for computing the extended gauss codes for each Polygon. Knotplot can be purchased for a small fee at http://knotplot.com/. 

## Contributors (ordered chronologically)
Rob Schaerine, Reuben Brasher, Robert Stolz

## License
All rights reserved. 2017. 

This repository is currently private, and this source code is for internal use during development. A stable branch will be made public at an appropriate time.  
