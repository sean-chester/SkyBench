## SkyBench

Version 1.1

© 2015-2016 Darius Šidlauskas, Sean Chester, and Kenneth S. Bøgh

-------------------------------------------
### Table of Contents 

  * [Introduction](#introduction)
  * [Algorithms](#algorithms)
  * [Datasets](#datasets)
  * [Requirements](#requirements)
  * [Usage](#usage)
  * [License](#license)
  * [Contact](#contact)
  * [References](#references)
  

------------------------------------
### Introduction

The *SkyBench* software suite contains software for efficient main-memory 
computation of skylines. The state-of-the-art sequential (i.e., single-threaded) and 
multi-core (i.e., multi-threaded) algorithms are included. 

[The skyline operator](https://en.wikipedia.org/wiki/Skyline_operator) [1] identifies 
so-called pareto-optimal points in a multi-dimensional dataset. In two dimensions, the 
problem is often presented as 
[finding the silhouette of Manhattan](http://stackoverflow.com/q/1066234/2769271):  
if one has knows the position of the corner points of every building, what parts of 
which buildings are visible from across the river? 
The two-dimensional case is trivial to solve and not the focus of *SkyBench*.

In higher dimensions, the problem is formalised with the concept of _dominance_: a point 
_p_ is _dominated by_ another point _q_ if _q_ has better or equal values for every 
attribute and the points are distinct. All points that are not dominated are part of 
the skyline. For example, if the points correspond to hotels, then any hotel that is 
more expensive, farther from anything of interest, and lower-rated than another choice 
would _not_ be in the skyline.  In the table below, _Marge's Hotel_ is dominated by 
_Happy Hostel_, because it is more expensive, farther from Central Station, and lower 
rated, so it is not in the skyline. On the other hand, _The Grand_ has the best rating 
and _Happy Hostel_ has the best price. _Lovely Lodge_ does not have the best value for 
any one attribute, but neither _The Grand_ nor _Happy Hostel_ outperform it on every 
attribute, so it too is in the skyline and represents a good _balance_ of the attributes. 


|Name         |Price per Night|Rating|Distance to Central Station|In skyline?|
|:------------|--------------:|:----:|:-------------------------:|:---------:|
|The Grand    |           $325| ⋆⋆⋆⋆⋆|                      1.2km|          ✓|
|Marge's Motel|            $55|    ⋆⋆|                      3.6km|           |
|Happy Hostel |            $25|   ⋆⋆⋆|                      0.4km|          ✓|
|Lovely Lodge |           $100|  ⋆⋆⋆⋆|                      8.2km|          ✓|


As the number of dimensions/attributes increases, so too does the size of and difficulty 
in producing the skyline. Parallel algorithms, such as those implemented here, quickly 
become necessary. 

*SkyBench* is released in conjunction with our recent ICDE paper [2]. All of the 
code and scripts necessary to repeat experiments from that paper are available in 
this software suite. To the best of our knowledge, this is also the first publicly 
released C++ skyline software, which will hopefully be a useful resource for the 
academic and industry research communities.


------------------------------------
### Algorithms

The following algorithms have been implemented in SkyBench:

 * **Hybrid** [2]: Located in [src/hybrid](src/hybrid). 
 It is the state-of-the-art multi-core algorithm, based on two-level 
 quad-tree partitioning of the data and memoisation of point-to-point 
 relationships.
 
 * **Q-Flow** [2]: Located in [src/qflow](src/qflow). 
 It is a simplification of Hybrid to demonstrate control flow.
 
 * **PSkyline** [3]: Located in [src/pskyline](src/pskyline).
 It was the previous state-of-the-art multi-core algorithm, based 
 on a divide-and-conquer paradigm.
 
 * **BSkyTree** [4]: Located in [src/bskytree](src/bskytree). 
 It is the state-of-the-art sequential algorithm, based on a 
 quad-tree partitioning of the data and memoisation of point-to-point 
 relationships.
  
All four algorithms are implementations of the common interface defined in 
[common/skyline_i.h](common/skyline_i.h) and use common dominance tests from  
[common/common.h](common/common.h) and [common/dt_avx.h](common/dt_avx.h) 
(the latter when vectorisation is enabled).

------------------------------------
### Datasets

For reproducibility of the experiments in [2], we include three datasets.
The [WEATHER](workloads/elv_weather-U-15-566268.csv) dataset was originally obtained from 
[The University of East Anglia Climatic Research Unit](http://www.cru.uea.ac.uk/cru/data/hrg/tmc) 
and preprocessed for skyline computation.
We also include two classic skyline datasets, exactly as used in [2]: 
[NBA](workloads/nba-U-8-17264.csv) and 
[HOUSE](workloads/house-U-6-127931.csv).

The synthetic workloads can be generated with the standard benchmark skyline 
data generator [1] hosted on 
[pgfoundry](http://pgfoundry.org/projects/randdataset).
  

------------------------------------
### Requirements

*SkyBench* depends on the following applications:

 * A C++ compiler that supports C++11 and OpenMP (e.g., the newest 
 [GNU compiler](https://gcc.gnu.org/)) 
 
 * The GNU `make` program

 * AVX or AVX2 if vectorised dominance tests are to be used


------------------------------------
### Usage

To run, the code needs to be compiled with the given number of dimensions.^
For example, to compute the skyline of the 8-dimensional NBA data set located
in `workloads/nba-U-8-17264.csv`, do:

> make all DIMS=8
>
> ./bin/SkyBench -f workloads/nba-U-8-17264.csv

By default, it will compute the skyline with all algorithms. Running `./bin/SkyBench`
without parameters will provide more details about the supported options.

You can make use of the provided shell script (`/script/runExp.sh`) that does all of
the above automatically. For details, execute:
> ./script/runExp.sh

To reproduce the experiment with real datasets (Table II in [2]), do (assuming
a 16-core machine):
> ./scripts/realTest.sh 16 T "bskytree pbskytree pskyline qflow hybrid"

^For performance reasons, skyline implementations that we obtained from other 
authors compile their code for a specific number of dimensions. For a fair
comparison, we adopted the same approach.


------------------------------------
### License

This software is subject to the terms of 
[The MIT License](http://opensource.org/licenses/MIT), 
which [has been included in this repository](LICENSE.md).


------------------------------------
### Contact

This software suite will be expanded soon with new algorithms; so, you are 
encouraged to ensure that this is still the latest version. Please do not 
hesitate to contact the authors if you have comments, questions, or bugs to report.
>[SkyBench on GitHub](https://github.com/sean-chester/SkyBench) 


------------------------------------
### References

 1. 
S. Börzsönyi, D. Kossmann, and K. Stocker. 
(2001) 
"The Skyline Operator."
In _Proceedings of the 17th International Conference on Data Engineering (ICDE 2001)_, 
421--432.
http://infolab.usc.edu/csci599/Fall2007/papers/e-1.pdf

 2. 
S. Chester, D. Šidlauskas, I Assent, and K. S. Bøgh. 
(2015) 
"Scalable parallelization of skyline computation for multi-core processors."
In _Proceedings of the 31st IEEE International Conference on Data Engineering (ICDE 2015)_, 
1083--1094.
http://cs.au.dk/~schester/publications/chester_icde2015_mcsky.pdf

 3. 
H. Im, J. Park, and S. Park. 
(2011) 
"Parallel skyline computation on multicore architectures."
_Information Systems_ 36(4): 
808--823.
http://dx.doi.org/10.1016/j.is.2010.10.005

 4. 
J. Lee and S. Hwang. 
(2014) 
"Scalable skyline computation using a balanced pivot selection technique."
_Information Systems_ 39: 
1--21.
http://dx.doi.org/10.1016/j.is.2013.05.005

------------------------------------
