# RECOMBO 

### Synopsis
RECOMBO is the working title of the Topological Molecular Biology Recombination Tool Suite. The suite is made up of a number of different applications that share a mostly common code base and work well together. A basic overview: these applications cover generating and performing topological reconnection on SAP conformations of fixed topology (MMC), removing Reidemeister crossings in knot projections (xinger), computing homfly polynomials from extended gauss codes (homfly), and identifying homfly polynomials with a topology (idknot). RECOMBO also includes an integrated unit testing suite with ok-ish test coverage, and a thoroughly experimental fugacity parameter finder (zAnalyzer). A task-oriented walkthrough is provided by the 'Guide' section below. Below that, the 'Manual' section provides an application-oriented look at the optional arguments of each program. 

# Guide
Updated Dec 10, 2020

### Introduction
In this guide, we will walk through setup and all the steps in the RECOMBO workflow. These steps are: generating randomized self-avoiding polygons, modifying the geometry of existing conformations, identifying knot and link-type, and analysis of the resulting transition data.

For our example, we will start where every project starts, a question: "What linktypes are possible outcomes of topological reconnection on the 5_1 knot, and which  is most likely to occur?" If some of those words or concepts are unfamiliar to you, I would recommend pausing to read [this] and [this](https://www.nature.com/articles/s41598-017-12172-2).

### Installation (OSX)
0. [Install and configure git](https://git-scm.com/book/en/v2/Getting-Started-Installing-Git), and a C++ compiler (G++,GCC,etc). In OSX, both of these can be installed from Xcode (avaliable for free on the app store) [like so](https://www.embarcadero.com/starthere/xe5/mobdevsetup/ios/en/installing_the_commandline_tools.html).  

1. Set up the directory structure and pull the repository using Git.
Open the terminal. Start by navigating to the directory where you want to download the repository using the `cd` command. You can create a new directory using the `mkdir` command. 

    For example, if we wanted to create a new folder to host the repository called `recombo_public` in our home directory, we might type:
```bash
cd ~
mkdir recombo_public
cd recombo_public
git init
```

2. Download the repository using the `git pull` command, as described in the Git documentation. Continuing the above example, we might use:
```bash
git pull https://github.com/garden-of-delete/recombo_public master
```

3. Building requires `autotools`. From the project root directory type the following
```bash
./configure
make
```

4. Assuming compilation is successful, project root directory will contain these four executables: `mmc`, `xinger`, `homfly`, and `idknot`.

### Generating randomized self-avoiding polygons
The Multiple Markov Chain BFACF executable (`mmc`) is the most direct path to generating randomized self avoiding polygons of known knot or link-type from scratch. `mmc` is often the starting point for the RECOMBO workflow. `mmc` requires a coordinate file specifying an initial SAP conformation to randomize. Initial conformations for many knot and link-types can be found in the `/initial` folder. The initial conformations are named according to the following scheme, derived from [Robert Scharein's Knotplot](https://knotplot.com/):
- 2_2_1a : initial conformation in knotplot (type `load 3.1` in the KnotPlot console)  
- 2_2_1b : initial conformation in knotplot with one component reversed  
- 2_2_1c : initial cofnormation in knotplot reflected with one component reversed  
- 2_2_1d : initial conformation in knotplot reflected 
- 3_1 : initial conformtion in knotplot  
- 3_1s : initial conformation in knotplot reflected 

There is also a [new biologically motivated notation](https://pdfs.semanticscholar.org/aac2/2f9e6d967d83991c5797f90fae26525ccd53.pdf) which resolves ambiguity around the 'default' chirality and orientation. A great project for a future undergraduate student would be to update or re-generate initial conformation directory with a new naming scheme based on this notation.

`mmc` has many parameters that must be chosen, and the interested user should dig into the [RECOMBO appendix] to learn more about the theoretical significance of each one. Most importantly, we must choose:  
- the fugacity parameter (`-zmax`: larger values increase average length)  
- the number of conformations to be sampled (`-n`)
- the number of BFACF moves between samples (`-c`: larger values decrease the statistical dependence of sequential samples)  
- the number of BFACF moves to take as a warmup to randomize the initial conformation (`-w`: larger values decrease the statistical dependence between the first sample and the initial conformation)  
- the sampling mode (`-mode`: controls what is saved / output from `mmc`).  

The parameter that controls the number of Markov chains (`-m`) can be set to 1, which makes `mmc` run as an instance of serial BFACF. I recommend the beginner run `mmc` in this mode, since running more than one Markov chain requires a slew of other parameters to be specified (and understood). For the advanced user, these parameters are the lower fugacity parameter (`-zmin`), the swap ratio (`-sr`), the number of steps between attempted swaps (`-s`). I have provided reasonable choices of these parameters below as a guide

To start our example, we would like to generate a collection of random 5_1 knots. We will return to the project directory, create a directory for our results, and execute mmc like so:
```bash
cd ~/recombo_public
mkdir results
src/bin/mmc initial/5_1 results/5_1 -zmin 0.1800 -zmax 0.1800 -q 1 -sr .8 -s 5 -n 5000 -c 20000 -m 1 -w 100000 -mode s -bfs 5000 +s > 5_1_log.txt
```
The log for the run can be found (with non-verbose output)can be found in the file `5_1_log.txt`. The `+s` option should be used to supress progress output when sending stdout to a file. The `-bfs 5000` option tells `mmc` to save up to 5000 samples per output file, and since we are asking for 5000 samples with the `-n 5000` option, all the conformations we sample will go to a single output file.  

When `mmc` finishes, 5000 5_1 knot conformations should have been deposited in a file called `5_1n1.b` in our `results/` sub-directory. 

### CMC Mode (Optional)
The `mmc` executable also has the ability to generate conformations and sample as a composite Markov chain. One might be motivated to use this functionality if sampling a wide range of lengths, or if randomization of a particular knot or link-type is too slow in serial BFACF mode. The software orders the BFACF Markov chains by their fugacity parameter, which provides a concept of 'adjacency' between two chains. To run `mmc` in composite Markov chain mode, the `-m` parameter, which specifies the initial number of Markov chains, should be set to an integer value larger than 1. The `-zmin` parameter must also be specified as a value lower than `-zmax`. Any Markov chains beyond 2 will have their fugacity parameter specified automatically by equidistantly chosen points on a linear curve between `-zmin` and `-zmax`. During the warmup, `mmc` will check if the swap-ratio `-sr` (the proportion of attempted conformation exchanges between adjacent chains) is at or above the specified frequency. Failure to swap conformations frequently enough is a consequence of adjacent Markov chains with fugacity parameters too far apart, so new Markov chains will be added as needed with fugacity parameters between them. The user must also choose a number of BFACF moves between attempted conformation swaps with the `-s` option. Small integer values work best, as the goal of MMC is to swap information as frequently as possible within performance limitations.

This is <b>not</b> part of our example, but we might use a commmand like this to generate a similar dataset to the above in composite markov chain mode.
```bash
src/bin/mmc initial/5_1 results/5_1 -zmin 0.1750 -zmax 0.1850 -q 1 -sr .8 -s 5 -n 5000 -c 20000 -m 5 -w 100000 -mode s -bfs 5000 +s > 5_1_cmc_log.txt
```
NOTE: It is possible that up to `m-1` additional conformations will be saved above the `-n` conformations requested when running in CMC mode. Everything is ok. 

### Performing Topological Reconnection
It may be desirable to modify the geometry of the generated conformations in a controlled way. The RECOMBO suite in its current form was assembled to carry out one such modification, called topological reconnection (see [this](https://www.nature.com/articles/s41598-017-12172-2) for more details). The `mmc` executable supports topological reconnection at the time of polygon generation through the `-mode f` (filter mode) and `-mode r` (reconnection mode). Filter mode saves only conformations that have sites meeting the arclength specified through both `-minarc` and `-maxarc`, and/or `-targetlength`. Reconnection mode saves conformations on the specified interval like `-mode s`, but also saves a post-reconnection version of any conformations that has a site meeting the specified criteria.

To continue our example, suppose we wanted to generate a larger collection of 5_1 knots, but also perform reconnection in inverted repeat when able. For our arclength criteria, we will require that at least 6 edges fall on one side of the reconnection site, and we will require the conformation to be exactly 120 edges long.
```bash
src/bin/mmc initial/5_1 5_1 -zmin 0.1800 -zmax 0.1800 -q 1 -sr .8 -s 5 -n 5000 -c 20000 -m 1 -w 100000 -mode f -minarc 4 -maxarc 240 -seed 42 +s > 5_1_log.txt
```
Note: Conformations with valid reconnection sites range from relatively rare to supremely rare. The above command can be expted to take multiple days to finish. A pre-computed output file from a similar run can be found at `demo/5_1_after.b`. Copy it to the results directory with the following command:
```bash
cp demo/5_1_after.b results/
```

### Identifying Knot and Link-type
The conformation identification pathway relies on an existing installation of Knotplot with BFACF functionality enabled. Knotplot is commercial software that can be purchased [here](https://knotplot.com/download/). Its creator, Rob scharein, must be contacted via email to provide instructions and a special license which enables the BFACF engine. Once knotplot is installed, the we must locate the stand-alone knotplot executable included in the installation (default location is `/Applications/KnotPlot/utilities/knotplot`), and copy/paste it into the `/knotplot` subdirectory of our local repository.

```bash
cp /Applications/KnotPlot/utilities/knotplot knotplot/
```

The first step in the identification pathway uses knotplot to carry out a sequence of critical tasks. First, the conformation is loaded into knotplot's BFACF engine. Then, a number of BFACF moves are performed allowing only moves that maintain or reduce the size of the conformation. This has the effect of reducing the conformation's length . The conformation is taken off the lattice and embedded in R3, and an energy minimization algorithm straightens the conformation out to disambiguate the crossings in the projection. Finally, a projection of the conformation is taken onto the plane and the extended Gauss code is computed. We will carry out all these tasks in sequence using a knotplot script. It is neccisary to know the number of components in our conformations of interest, as knot and two-component link conformations require different treatment. 

To continue our example, the following knotplot script would be suitable for identifying the post-reconnection conformations we generated in the previous section:
```
alias idknot "bfacf load; bfacf prob 1 1 0; bfacf step 20000; ago 200; centre; align axes; gauss"
alias idlink "save |tmp; sequence +; load combine |tmp"
gauss noblank yes

gauss open 5_1.egc
sequence open 5_1_after.b
until seqend "idlink; idknot; sequence +"
sequence close
gauss close
```
To execute this script, we should make sure we are in the main directory of our local repository, then invoke knotplot like so:
```bash
cd ~/recombo_public
knotplot/knotplot -nog < knotplot/idlink_example.kps
```

The above knotplot script and another example for knotted conformations have been included in the `/knotplot` subdirectory of this repository. More information about knotplot's command line interface can be found [here](https://knotplot.com/KPman/commands.html).

After the extended Gauss codes are computed, we can remove Reidemeister crossings (if applicable), compute the HOMFLY-PT polynomial, and identify the results. Using the piping feature in the terminal, we can do this with one commaqnd. :
```bash
cat results/5_1_after.egc | src/bin/homfly -- | src/bin/idknot -- -nz -u -P > results/5_1_after.txt
```  

Now open `results/5_1_after.txt` and behold, a distribution of resulting link-types. 

### Analysis of statistical dependence in sampled conformations
Until now, we have not been particularly concerned about wether or not the SAP conformations we generate are correlated statistically. However, when the goal is to make an unbiased estimation that reflects the population of 5_1 conformations, we must either:
1. ensure that our samples are likely identically and independently distributed
2. compensate for the correlated samples through a batch-mean approach

Option 1 is best accessed through the `-mode a` or `-mode b` option in `mmc`. Both modes run an autocorrelation and length analysis on each chain after warmup. An example:

```bash
cd ~/recombo_public
src/bin/mmc initial/5_1 results/5_1 -zmin 0.1800 -zmax 0.1800 -q 1 -sr .8 -s 5 -n 5000 -c 20000 -m 1 -w 1000000 -mode a -bfs 5000 -seed 42 +s
```
Output:
```bash
s 5000 -seed 42 +s
filename= results/5_1
start_time= 2020-10-19-08-26
zmin= 0.18
zmax= 0.18
q= 1
target_swap_ratio= 0.8
swap_interval= 5
n_samples= 5000
steps_between_samples= 20000
initial_n_chains= 1
warmup= 1000000
sample_mode= a
seed= 42
block_file_size= 5000
time_limit= 0

Warming up 1 chains: (1000000 steps)

Starting calibration...
Testing chains... 

Starting sampling with 1 chains...

Estimating average lengths with swapping... (crude estimate)

0.18 0.5 0.00316228 58.8855 0.0517601

0.18 58.8855 133
```
First we see that the run parameters have been output. Here, `0.18 0.519492 0.00734673 58.9078 0.0532237` is the result. From left to right, the values displayed are: fugacity parameter, estimated integrated autocorrelation time, standard error in estimated integrated autocorrelation time, mean length, and standard error in the mean length. An estimated integrated autocorrelation time close to .5 indicates statistical independence. On the next line, `0.18 58.8855 133` is an the fugacity parameter, independent estimate of the average length made without using the autocorrelation module, and a display of the raw variance in length. The second line provides a sanity check on rather opaque code behind the first, and one should expect the average length values to agree exactly.  

Option 1 is most suitable for exploratory work, when the relationship between fugacity parameter and length is unclear, or when limited on disk space for generating a large dataset. 

Option 2 has the advantage of allowing the use of correlated samples with the drawback of requiring substantial additional information to be saved during the sampling process. It is therefore most suitable for large studies with a well explored BFACF parameter set and lots of storage. The `mmc` executable only saves this data when running in `-mode m`. See `scripts/batch_mean_analysis.py` to get a sense of what to do with the output from this mode. 

# Manual
### Installation
This software suite is intended for use on Linux and OSX. Install by following these steps:
1. Download and extract this repository.
2. In the terminal, navigate inside the repository's `src/` subdirectory. (`cd /path/to/extracted/directory')
3. Execute `make all`. Warnings may be generated from the use of depricated language features and syntax.
4. Assuming compilation is successful and all unit tests pass, navigate to the `/src/bin` subdirectory to find the executables.
5. If compilation is not successful, examine the error messages and ensure your compiler and environment are treating the code appropriately.

### MMC
MMC is the largest and most complex piece of software in the RECOMBO suite. It implements an extension of the BFACF algorithm, CMC-BFACF, [Stolz et. al., 2017](https://www.nature.com/articles/s41598-017-12172-2),Szafron,Orlandini to generate and sample self-avoiding polygons. Because of its complexity, there are a number of options and run-modes avaliable to the user to best meet your simulation needs. 

#### Usage
For a more detailed explanation of the theory behind Composite Markov Chain Monte Carlo, I refer the reader to [Stolz PhD Thesis]. Here is described a number of the important arguments with examples. Determining which arguments are specifically best for their simulation is left to the reader, though suggestions appear below. (wording)

#### Argument List

`-zmin` the smallest z value set to one of the initial markov chains

`-zmax` the largest z value set to one of the initial markov chains

`-sr` the user-defined minimum acceptaple swap ratio that must be satisfied between adjacent markov chains. If the swap-ratio is not satisfied, a new Markov chain will be added with a z value half-way between them on a linear curve. 

`-s` the swap interval, or the number of BFACF moves between swap attempts in each chain

`-c` the sampling interval, or the number of BFACF moves between samples (or attempted samples)

`-m` initial number of chains initialized on the interval defined by the `-zmin` and `-zmax` arguments. 

The sampling mode dictates how BFACF samples conformations. Basic options:

`-mode a` Analyze only mode. Perform warmup and calibration, measuring and outputting statistics for each chain. 

`-mode s` Sample mode. Sample conformations in at an interval specified by the `-c` option into a binary format file.

`-bfs n` Block-file size. Specifies the number of conformations to save per file, `n`. Useful for splitting up the dataset. With the `-bfs` option unused or if `n` is 1, each conformation will be output in its own CUBE binary formatted `.b` file. NOTE: The binary representation of the CUBE format should be thought of as a linked list. Retrieving a conformation at a particular index requires iterating through all the preceding conformations in the file. Therefore, it is advisable to think about how you will want to organize and analyze the output conformations when specifing the `-bfs` option.

`-mode b` Combination of `a` and `s` mode.

Advanced operating modes work in conjunction with recombo to perform reconnection in mmc:

`-mode r` Reconnection mode. Samples like `s` mode but saves a post reconnection version of each sampled conformation in a separate file. Requires reconnection criteria to also be specified.

`-mode f` Filter sample mode. Samples like `r` mode but saves only conformations that have sites satisfying the defined reconnection criteria.

`-minarc n` Mandatory when in one of the two advanced operating modes. Defines the minimum arclength to identify a reconnection site as the integer `n`.

`-maxarc n` Mandatory when in one of the two advanced operating modes. Defines the maximum arclength to identify a reconnection site as the integer `n`.

`-targetlength n` Optional when in one of the two advanced operating modes. For knots, the conformation must be of length `n` within a defined arclength tolerance. For two-component links, the components must add to `n` within a defined arclength tolerance.

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

#### Example Usage (revise)
`mmc initial/3_1 3_1 -zmin 0.20466 -zmax 0.20966 -q 1 -sr .8 -s 5 -n 800000 -c 20000 -m 5 -w 1000000 -mode m -minarc 57 -maxarc 63 -targetlength 120 -seed 42 +s > 3_1_log.txt &`

## homfly
Homfly takes as an input a file of extended gauss codes (.egc), and outputs the HOMFLY-PT polynomial [x] for each one. 

#### Usage
`homfly` for interactive usage. If you are a mathematician looking to 

`homfly input_file output_file` 

`homfly -- |` or ` homfly -- > output_file` for use in a stream.

## xinger
Usage: xinger [options/operators] input-xing-sequence output-xing-sequence
       or
       xinger [options/operators] input-xing-sequence > output-xing-sequence
       or
       some-xing-producing-program | xinger -- [options/operators] > output-xing-sequence
       Program reads a sequence of knots with extended gauss or dowker codes (dowker codes
       plus orientation information),
       then performs a (possibly null) series of operations on the crossings in the knot,
       then writes out 
       
       Can read / write:
                   --- old style triple file format (wastes disc space)
                   --- compact binary format
                   --- gauss code format (default format for writing
       
       xing operators:
       +r          remove crossings using R1 and R2 moves until minimized
       +R          remove crossings using R1 and R2 moves until minimized, ignoring LIMIT
       +1          remove crossings using only R1
       +x          pick a crossing uniformly at random and flip it
       +c          pick a clasp uniformly at random and flip it
       +p          pick a positive crossing uniformly at random and flip it
       +n          pick a negative crossing uniformly at random and flip it
       +P          pick a positive super crossing uniformly at random and flip it
       +N          pick a negative super crossing uniformly at random and flip it
       +i          pick the first isolated clasp (if it exists) and flip it
       +o          reverse the orientation of the knot
       +z          reflect the knot in the z-direction
       +Z          flip all crossings randomly
       +b          pick a crossing uniformly at random and recombo it
       +S          pick a super-crossing uniformly at random and recombo it
       +l          find standard labellings for crossings
       +k          write number of crossings (useful together with the -wn option)
       +K          write number of isolated clasps (useful together with the -wn option)
       +D          write debug / diagnostic information about each diagram
       -prog PROG  run program PROG (not implemented yet)
       
       
       xing options:
       -p PROB     set probability value to PROB (default is 1)
       -L LIMIT    minimum number of crossings during removal (default is 0, ignored using -R)
        
       
       read / write options:
       -r3         read old style from three files
       -rg         read gauss format (only valid with the `--' option, ignored otherwise)
       -rb         read binary format only valid with the `--' option, ignored otherwise)
       -rx         read KnotPlot X-file format (only valid with the `--' option, ignored otherwise)
       -wg         write a gauss code (the default)
       -wag        write all gauss codes consistent with knot shadow, ignoring chirality
       -wag2       write all gauss codes consistent with knot shadow, considering chirality
       -wNg N      write N gauss codes, choosing signs of crossing randomly
       -w3         write old style to three files, wasting disc space
       -wb         write binary
       -wj         write Jenkins HOMFLY format (default style)
       -wjo        write Jenkins HOMFLY format (old style)
       -wjn        write Jenkins HOMFLY format (new style)
       -wkg        write knot group representation
       -wem        write Ewing / Millett format
       -wn         write nothing (no knots), useful for checking for errors
       -wx         write number of crossings
       -anc        write average number of isolated clasps
       -nx MIN MAX only write knots that have number crossing >= MIN and <= MAX
       -v          produce verbose output
       -q          run quietly (the default)
       -s SEED     set random number seed to SEED (use this only for testing purposes)
       -f FORMAT   set 3-file format base name to FORMAT
       -zok        OK to output knots with zero crossings
       -znb        don't write blank lines for knots with zero crossings
       --          read data from standard input instead of opening a file
      
## idknot
Idknot takes a input file or stream of homfly polynomials and associates each one with a topology, or 'unknown' if unrecognized. It expects an input file of homfly polynomials, one on each line, and can output the identifications in a variety of ways as described below.

#### Usage
Running `idknot` with no other arguments will print out detailed usage information. Some common usage scenarios:

`idknot input_file output_file -nz -u -P` to specify input output files.

`idknot input_file output_file -k` to print conformation identifications line-by-line. 

`idknot -- output_file -nz -u -P` or `idknot -- -nz -u -P > output_file` substitutes the input file for a data stream. 

## Other Software
There are other useful applications that are not included as part of RECOMBO: 

Xinger tries to remove reidermeister crossings at the extended Gauss code level, which improves the performance of Homfly. It is not required to compute the Homfly polynomial on an extended Gauss code. It is avaliable in limited circumstances upon request by emailing a request to: mariel@math.ucdavis.edu. 

Seqconvert is a generally useful application that efficiently converts files filled with self-avoiding polygons between different formats. For example, it can be used to convert CUBE binary formatted polygons to ascii formatted files. It can also perform and output some useful computations for all polygons in a file, such as polygon length and writhe. It is also avaliable in limited circumstances upon request by emailing a request to: mariel@math.ucdavis.edu.

Knotplot is a wonderful knot computation and visualization software developed and maintained by Rob Scharein. In the context of reconnection studies, it is useful for computing the extended gauss codes for each Polygon. Knotplot can be purchased for a small fee at http://knotplot.com/ .

---

## Contributors (ordered chronologically)
Rob Scharein, Reuben Brasher, Robert Stolz, Michelle Flanner, Zihao Zhu, Diwen Lu

The code from `src/autoCorr.cpp` was extracted from a original C program by Buks van Rensburg and Enzo Orlandini in Stu Whittington's research group.

Repository currently maintained by Robert Stolz. 

## Citations

### Pathways of DNA unlinking: A story of stepwise simplification

```
@article{stolz2017pathways,
  title={Pathways of DNA unlinking: A story of stepwise simplification},
  author={Stolz, Robert and Yoshida, Masaaki and Brasher, Reuben and Flanner, Michelle and Ishihara, Kai and Sherratt, David J and Shimokawa, Koya and Vazquez, Mariel},
  journal={Scientific reports},
  volume={7},
  number={1},
  pages={1--11},
  year={2017},
  publisher={Nature Publishing Group}
}
```

### New biologically motivated knot table
Initial conformations in directory initial/ the convention from this paper.

```
@misc{brasher2013new,
  title={New biologically motivated knot table},
  author={Brasher, Reuben and Scharein, Rob G and Vazquez, Mariel},
  year={2013},
  publisher={Portland Press Ltd.}
}
```

### KnotPlot
Much of the work related to this project would not have been feasible without use of Rob Scharein's KnotPlot tool. Works using KnotPlot should cite the following two:

```
@PhdThesis{SchareinPhD,
  author={Robert G. Scharein},
  title={Interactive Topological Drawing},
  school={Department of Computer Science,
The University of British Columbia},
  year=1998}
```

```
@book{borwein2002multimedia,
  title={Multimedia Tools for Communicating Mathematics:[presentations at an International Workshop MTCM2000, Organized at the Centro de Matem{\'a}tica E Aplica{\c{c}}oes Fundamentais at the University of Lisbon, in November 2000]},
  author={Borwein, Jonathan and Morales, Maria H and Polthier, Konrad and Rodrigues, Jose F},
  year={2002},
  publisher={Springer Science \& Business Media}
}
```


## Funding
This research was supported by the following: Japan Society for the Promotion of Science KAKENHI grant numbers 25400080, 26310206, 16H03928, 16K13751, 17H06463(to K.S.), 26800081 (to K.I.); National Science Foundation DMS1716987 (MF, MV) and CAREER Grant DMS1057284 (MV, RS, MF, RB) and NIH-R01GM109457 (MV); Welcome Trust SIA 099204/Z/12Z and 200782/Z/16/Z (DJS). We are grateful to Rob Scharein for providing assistance with Knotplot and for his work on the first version of the reconnection software; C. Soteros, M. Szafron and M. Schmirler for contributing their statistical expertise; J. Arsuaga, D.W. Sumners and S. Witte for helpful discussions.

## License
This software is under a GPL 3.0 License. 
