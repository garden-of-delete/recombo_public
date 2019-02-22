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
MMC is the largest and most complex piece of software in the RECOMBO suite. It uses an extension of the BFACF algorithm (CMC-BFACF) [Stolz et. al., 2017,Szafron,Orlandini] to generate and sample self-avoiding polygons. Because of its complexity, there are a number of options and run-modes avaliable to the user to best meet your simulation needs. 

#### Usage
For a more detailed explanation of the theory behind Composite Markov Chain Monte Carlo, I refer the reader to [Stolz PhD Thesis]. Here is described a number of the important arguments with examples. Determining which arguments are specifically best for their simulation is left to the reader, though suggestions appear below. (wording)

#### Argument List

`-zmin` the smallest z value set to one of the initial markov chains

`-zmax` the largest z value set to one of the initial markov chains

`-sr` the user-defined minimum acceptaple swap ratio that must be satisfied between adjacent markov chains

`-s` the swap interval, or the number of BFACF moves between swap attempts in each chain

`-c` the sampling interval, or the number of BFACF moves between samples (or attempted samples)

`-m` initial number of chains initialized on the interval defined by the `-zmin` and `-zmax` arguments. 

The sampling mode dictates how BFACF samples conformations. Basic options:

`-mode a` Analyze only mode. Perform warmup and calibration, measuring and outputting statistics for each chain. 

`-mode s` Sample mode. Sample conformations in at an interval specified by the `-c` option into a binary format file.

`-mode b` Combination of `a` and `s` mode.

Advanced operating modes work in conjunction with recombo to perform reconnection in mmc:

`-mode r` Reconnection mode. Samples like `s` mode but saves a post reconnection version of each sampled conformation in a separate file. Requires reconnection criteria to also be specified.

`-mode f` Filter sample mode. Samples like `r` mode but saves only conformations that have sites satisfying the defined reconnection criteria.

`-minarc n` Mandatory when in one of the two advanced operating modes. Defines the minimum arclength to identify a reconnection site as the integer `n`.

`-maxarc n` Mandatory when in one of the two advanced operating modes. Defines the maximum arclength to identify a reconnection site as the integer `n`.

`-targetlength n` Optional when in one of the two advanced operating modes. Def


`-recomboParas < >` is used in place of recomboFromFile.

`-recomboParas ds` Optional when in f mode. Do standard (on-lattice) recombination in direct repeat

`-recomboParas da` Optional when in f mode. Do writhe_based virtual (off-lattice) recombination in direct repeat

`-recomboParas dp` Optional when in f mode. Do positive virtual (off-lattice) recombination in direct repeat

`-recomboParas dn` Optional when in f mode. Do negative virtual (off-lattice) recombination in direct repeat 

`-recomboParas dsa` Optional when in f mode. Do a mix of standard and writhe-based recombination in direct repeat

`-recomboParas dsp` Optional when in f mode. Do a mix of standard and positive virtual recombination in direct repeat

`-recomboParas dsn` Optional when in f mode. Do a mix of standard and negative virtual recombination in direct repeat

`-recomboParas is` Optional when in f mode. Do standard recombination in inverted repeat

`-recomboParas ip` Optional when in f mode. Do positive virtual recombination in inverted repeat  (Not implemented)

`-recomboParas in` Optioanl when in f mode. Do positive virtual recombination in inverted repeat  (Not implemented)

`-recomboParas ia` Optional when in f mode. Do writhe-based virtual recombination in inverted repeat  (Not implemented)

`-recomboParas isp` Optioanl when in f mode. Do a mix of standard and positive virtual recombination in inverted repeat  (Not implemented)

`-recomboParas isn` Optional when in f mode. Do a mix of standard and negative virtual recombination in inverted repeat  (Not implemented)

`-recomboParas isa` Optional when in f mode. Do a mix of standard and writhe-based virtual recombination in inverted repeat  (Not implemented)

##### Optional Arguments
`-q` (Optional) the q value used for all markov chains in the mmc process. If unspecified, defaults to q=1.

`-seed` (Optional) Manually specify a seed for the random number generator. If unspecified, uses the system clock and records the seed to stdout (console output).

`+s` (Optional) Supresses progress output. Should be used when running long term jobs on a server that are not being directly monitored during the run.

`-info`(Optional) Outputs information needed to perform batch-mean analysis as described in George Fishman's "Monte Carlo: Concepts, Algorithms, and Applications" into a new file with the .info extension. 

##### Halting Conditions
At least one of these must be specified for the program to terminate. If more than one is specified, mmc will sample conformations until any halting condition is met. 

`-n` number of conformations to sample

`-t` integer number of hours to run before safe termination. 

#### Example Usage
`mmc initial/3_1 3_1 -zmin 0.20466 -zmax 0.20966 -q 1 -sr .8 -s 5 -n 800000 -c 20000 -m 5 -w 1000000 -mode m -minarc 57 -maxarc 63 -targetlength 120 -seed 42 +s > 3_1_log.txt &`

## recomboFromFile (Not working now. It is integrated to MMC by -recomboParas < >)
recomboFromFile takes an input file or stream of CUBE binary formatted polygons, and attempts to perform reconnection on each one according to user-specified reconnection criteria. It is capable of performing reconnection in direct ~~and inverted~~ repeat. 

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
