# pyADAPT

Analysis of Dynamic Adaptations in Parameter Trajectories

[http://bmi.bmt.tue.nl/sysbio/software/adapt.html](http://bmi.bmt.tue.nl/sysbio/software/adapt.html)

Implementation in Python

## Dev Plan

1. Reproduce a similar procedure in python
    - The toy model in MATLAB
    - trehalose model with van Heerden's glucose pulse response
    - design Lotka model, generate data
2. Read Dissertation chapter 2 analysis and write an analysis routine
3. read least squares optimization on wikipedia and differential equations on Julia Diffeq's documentation (improve writing and documentation)
4. Parameter sensitivity analysis
5. development on the `pyADAPT` package
    - the `k3` parameter of toy model (see the matlab results)
    - the SBML format support
        - this is to support trehalose model
    - TODOs in the source code
    - plotting and profiling
6. finish CLI functions. in functional programming fashion.
    - convert, analysis
    - maybe just leave the templating work to cookicutter.

## Let's split: What should be OOP and what should be functional?

This is a really simple project so please don't bother writing everything in OOP fashion which is unnecessary and tedious.

## Writing plans

The development has been halted due to a lack of motivation. Now I have a view of reproducing every case in the 2013 case. However, the problem is that only reproducing and showing result is not good enough as a master thesis.
