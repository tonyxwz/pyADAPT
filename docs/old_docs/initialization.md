# Initialization method of pyADAPT

Over the course of years of development, there're mainly three variants of initialization in ADAPT algorithm.

1. Natal van Riel (Interface Focus 2013): this method repeatedly generate random initial parameters from a `range`, then use least squares method to simulate the model for `ss_time` so that the model reaches to a fake steady states. Then check the cost returned by the lsq method. Accept the random initial parameters only if the cost is less than `sseThres`.
2. Roderick Snel (https://github.com/nvanriel/ADAPT): not worth to mention because it is really poor performance.
3. Pascal van Beek (master thesis): similar to the first method. But instead of simulating a long time step, it simulate the entire course of the time range that the model lives in.

pyADAPT only support Pascal's initialization method as it is the result of the latest research on this topic (model initialization).

## Why special for 1st time step?

ADAPT is an iterative optimization procedure, which means it takes the information from the previous time step as the starting point of the next time step. Of course, for the first time step, there's no previous information to start the optimization routine (`fit_timestep`). That's why we need to provide the first time step with some additionally information.

Let's take a look the first iteration in function `fit_iteration`, where the first time step is handled. When `i_ts = 1`, the `initial_guess` and `begin_stats` provided for `fit_timestep` are the first row of the trajectories.

For a time step in the middle of a iteration, the information it needs are

1. `initial_guess`: starting point of lsq on the parameters
2. `begin_state`: the states (s1, s2)

We need to manually calculate the above information at time step 0 in order to feed time step 1.

1. Generate random values for those states who are not provided in the data set.
2. Generate random initial parameters.
3. use these parameters to simulate the model on the entire time range.


## Pascal van Beek's Model Initialization

First estimate the parameters using classical approach.


## Solution (discussion with David, May 6th, 2020)

Trust only States from the data. Don't compare the fluxes, assign 0 to the weight of fluxes in the objective function.