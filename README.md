# mvb-detector

Repository with the code used for the paper "The Multivariate Bernoulli
detector: Change point estimation in discrete survival analysis" by Willem van
den Boom, Maria De Iorio, Fang Qian and Alessandra Guglielmi (in preparation)


## Description of files

* [`mvb_detector.R`](mvb_detector.R) implements the local-global Markov chain
Monte Carlo for the multivariate Bernoulli detector.

* [`simulation_with_censoring.R`](simulation_with_censoring.R) produces the
results for the simulation study with censoring in the main text of the paper.
[`simulation_with_censoring.R`](simulation_with_censoring.R) loads
[`mvb_detector.R`](mvb_detector.R).

* [`simulation_no_change_points.R`](simulation_no_change_points.R) produces the
results for the simulation study without change points in the appendices of the
paper. [`simulation_no_change_points.R`](simulation_no_change_points.R) loads
[`mvb_detector.R`](mvb_detector.R).
