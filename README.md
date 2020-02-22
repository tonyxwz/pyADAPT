# pyADAPT

Analysis of Dynamic Adaptations in Parameter Trajectories

[http://bmi.bmt.tue.nl/sysbio/software/adapt.html](http://bmi.bmt.tue.nl/sysbio/software/adapt.html)

Implementation in Python

## Dev Plan

1. constant parameters simulation
2. toy model simulation without considering input
3. clamp model with input
4. David yeast models


## Tips for Speeding Up the Process

My progress on my thesis is really slow in the past few month. Of course there
is the problem with the supervisors, but I cannot be excused for not working
hard and being lazy.

I need to focus on the following things from now on:

1. Read David's MATLAB code
2. Reproduce a similar procedure in python
3. Read Dissertation chapter 2 analysis
4. learn more about least squares optimization
5. improvement on the `pyADAPT` package

  - why it gives wrong results
  - write convert function
  - write documentation along the way

6. finish CLI functions. in functional programming fashion.

## Let's split: What should be OOP and what should be functional?

This is a really simple project so please don't bother writing everything in OOP
fashion which is unnecessary and tedious.

### OOP
  - Model, SBML model, bio parts
  - app
  - data set
  - convert, a convert procedure

### Functional, as library modules (procedures)
  - all the examples
  - analysis, an analysis procedure
  - objective functions
  - io, an IO procedure. BTW should this cover convert?
  - optimize, 


## Project plans

The development has been halted due to a lack of motivation. Now I have a view of
reproducing every case in the 2013 case. However, the problem is that only reproducing
and showing result is not good enough as a master thesis.

### What should be included?

So far, the project has been very "programming" focus, such as writing python
package for `pyADAPT`. The package doesn't give the correct results yet. 
