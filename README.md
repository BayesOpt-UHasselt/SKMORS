# SKORS
Multiobjective Ranking and Selection Using Stochastic Kriging
Source code for the SK-MORS algorithm

Copyright 2022 Sebastian Rojas Gonzalez

srojas33@gmail.com

sebastian.rojasgonzalez@ugent.be

Source code for the paper:

%--------------------------------------------------------------------------

Rojas Gonzalez, S., Branke, J., & van Nieuwenhuyse, I. (2022). Multiobjective Ranking and Selection Using Stochastic Kriging. arXiv e-prints, arXiv-2209.

%--------------------------------------------------------------------------

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the above copyright notice is retained.

To run any of the algorithms, execute this command in the terminal:

main('-algorithm', @SKMORS, '-problem', @WFG4, '-M', 2, '-D', 5, '-evaluation', 1000, '-mode', 1);

where

M: Number of objectives

D: Number of decision variables

evaluation: Number of algorithm iterations

mode:

1: Display data

2: Save data

-> This SK-MORS algorithm code uses partial implementations of PlatEMO to facilitate its usage:

%--------------------------------------------------------------------------

Ye Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB Platform for Evolutionary Multi-Objective Optimization [Educational Forum], IEEE Computational Intelligence Magazine, 2017, 12(4): 73-87".

%--------------------------------------------------------------------------

-> The GLOBAL class encodes the global properties of the experimental setting.

-> The INDIVIDUAL class encodes the individual properties of each design point.

System requirements:

Code has been tested using Matlab 2018b, 2019b, 2020a, 2020b, 2021a, 2021b and 2022a running on 64-bit Linux Debian/Ubuntu, with an Intel i7 VPro CPU with 12 cores in 2.60 GHz and 32 Gb of RAM.

The MOCBA source code is coded in C, based on this reference:

%--------------------------------------------------------------------------

Chen, C. H., & Lee, L. H. (2011). "Stochastic simulation optimization: an optimal computing budget allocation (Vol. 1)". World Scientific.

%--------------------------------------------------------------------------

-> C-Matlab bindings are provided to run the code in Matlab via mocba.mexa64

--> execute mocba(Mat_Obj',Mat_Var'), where

    --> Mat_Obj is the matrix of objective values transposed and
    
    --> Mat_Var is the matrix of variance values transposed 
-> These bindings are compiled for 64-bit Linux and haven't been tested on Matlab running on Windows or Mac.
