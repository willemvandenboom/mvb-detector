# mvb-detector

Repository with the code used for the paper "The Multivariate Bernoulli
detector: Change point estimation in discrete survival analysis" by Willem van
den Boom, Maria De Iorio, Fang Qian and Alessandra Guglielmi
([arXiv:2308.10583](https://arxiv.org/abs/2308.10583))


## Description of files

The repository provides an R package named `mvb.detector` that implements the
multivariate Bernoulli detector and which can be installed using
`remotes::install_github("willemvandenboom/mvb-detector")`.

Additionally, the folder [paper](paper/) contains code to reproduce the results
in the paper:

* [`paper/setup.R`] installs the required R packages including `mvb.detector`.

* The folder [paper/icu](paper/icu/) contains the code for the application to
ICU length of stay.

* The folder [paper/simulation](paper/simulation/) contains the code for the
simulation studies.
