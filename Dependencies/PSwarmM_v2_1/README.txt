Copyright (C) 2009 A. Ismael F. Vaz and L. N. Vicente.

http://www.norg.uminho.pt/aivaz
http://www.mat.uc.pt/~lnv

This library is free software; you can redistribute it and/or modify
it  under  the  terms  of  the  GNU Lesser General Public License as
published by the Free Software Foundation; either version 2.1 of the
License, or (at your option) any later version.

This library is distributed in the hope that it will be useful,  but
WITHOUT  ANY  WARRANTY;  without  even  the  implied   warranty   of
MERCHANTABILITY  or  FITNESS  FOR  A PARTICULAR PURPOSE. See the GNU
Lesser General Public License for more details.

You should have received a copy of the  GNU  Lesser  General  Public
License  along with this library; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307,
USA.



                     =========================
                     PSwarm MATLAB Version 2.1
                     =========================

                          ---------------
                            Read me file
                          ---------------


Welcome to Pattern Swarm MATLAB version 2.1 (PSwarmM v2.1)

This code provides an implementation of the algorithm described in:

A.I.F.  Vaz  and L.N.Vicente, A particle swarm pattern search method
for   bound  constrained  global  optimization,  Journal  of  Global
Optimization, 39 (2007) 197-219.

A.I.F.   Vaz   and   L.N.  Vicente.  PSwarm:  A  hybrid  solver  for
linearly    constrained    global    derivative-free   optimization.
Optimization   Methods   and   Software,  Optimization  Methods  and
Software, 24 (4-5), 669-685.

Le Thi Hoai An, A.I.F. Vaz and L.N. Vicente. Optimizing Radial Basis
Functions  by  D.C.  Programming  and  its  use in Direct Search for
Global Derivative-Free Optimization, to be submitted.



This file describes the following topics:


1. System requirements
2. Installation
3. Contents of the pswarmm_v2_1 directory
4. The command line
5. Credits


1. System requirements
----------------------------------------

Your computer should have a MATLAB numerical compiler installed.


2. Installation
----------------------------------------

The code requires approximately 8Mb of hard disk space. The majority
of  disk  requirement  is  due to the AMPL interface and the problem
models. Only the .m files are required for a standalone run  of  the
solver.

Unzip  the  PSwarmM_v2_1.zip file in any location using your favorite
unzip  software.  A  directory PSwarmM_v2_1 will be created, having a
subdirectory  'models',  which  refers to the AMPL problems model and
the corresponding .nl files.

This directory must be the MATLAB working directory or  included  in
the MATLAB search path.


3. Contents of the PSwarmM_v2_1 directory
----------------------------------------

In the directory PSwarmM_v2_1 there are the following files:

------ README -------------
README.txt - This readme file

lgpl.txt - GNU LESSER GENERAL PUBLIC LICENSE

------ AMPL Related -------
matampl.dll  -  An AMPL mex file. This is system dependent of MATLAB
               v7 (R14) and Microsoft windows XP.
matampl.mexw32 - AMPL mex file.
ampl_obj.m   - Objective function computation using AMPL.

------ Examples -----------
RunPswarm.m  -  Runs  the  PSwarmM  algorithm  for  a  test problem.

Solve_all.m - A batch script to run a set of problems.

hs024.m - Run PSwarmM on the hs024 problem.

hs024_obj.m - The hs024 objective function.

draw_ellipse.m - Draws a two dimensional ellipse.

------ Solver files -------
InitPatternSearch.m  - Initializes the search directions used in the
                      poll   step.   One  can  edit  it  to  provide
                      additional search directions.

PenaltyEval.m   -  Computes  the  objective  function  value.  Since
                   projection  into  the  feasible  domain  is  made
                   before  calling  this routine, checking for bound
                   feasibility is not requested.

PollStep.m - Preforms a poll step at the population (swarm) leader.

Projection.m - Projects  a  given  point  into  the  bound  feasible
               domain.

PSwarm.m - The main algorithm. Applies  a  particle  swarm  strategy
           whenever progress in the leader  is  being  obtained.  If
           progress is no longer occurring a poll step is made.

mve_presolve.m - Computes the maximum volume ellipsoid center. See
           Credits.

mve_solve.m - Computes the maximum volume ellipsoid. See Credits.

outputfcn.m - An output function called each IPrint iterations.

CheckStoppingCriteria.m - Checks for the stopping criteria.

InitCache.m - Initializes the cache.

InitPopulation.m  -  Initializes  que population (it may be only one
                     particle).

ParticleSwarm.m - A particle swarm iteration (search step).

UpdateInfo.m  -  To  be  called after a search and poll step (in the
                 particle   swarm  case  it  updates  the  particles
                 velocities).

RBF.m - Search step with the RBF model.




------ Directories --------
models - AMPL models for global optimization problems. This  problem
         collection  was used to obtain the numerical results in the
         previous  cited  paper. The .nl intermediate AMPL files are
         also  provided. They are requested for the MATLAB interface
         with AMPL and were obtained from the problems in the models
         directory.
Wild  -  MATLAB  files  from Stefen Wild with test problems for DFO,
         obtained at http://www.mcs.anl.gov/~more/dfo/



4. The command line
----------------------------------------

The only requested files to run the algorithm are the "Solver files"
(see  previous  section) and a user provided objective function. The
file RunPswarm.m provides an example of a run  with  the  test_obj.m
objective function.

To obtain the available options just  type  in  the  MATLAB  command
window:

PSwarm('defaults')

To obtain help on using the PSwarm library then just type

help PSwarm



5. Credits
----------------------------------------

The MVE MATLAB interior point code is due to Zhang and Gao:

Y.~Zhang  and  L.~Gao,  "On numerical solution of the maximum volume
ellipsoid  problem",  SIAM  Journal  on  Optimization, 14(1):53--76,
2003. http://www.caam.rice.edu/~zhang/mve/index.html
