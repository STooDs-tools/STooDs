# STooDs

STooDs is a framework to build and estimate probabilistic models for data varying in Space, Time or other Dimensions. A typical STooDs case study can be described as follows:

1. The target dataset contains the values taken by one of several predictand variables varying along one or more dimensions.
2. Each variable in the dataset is assumed to follow a distribution whose parameters may vary along the same dimensions as the dataset.
3. This variability is specified by regression formulas that may use parameters, covariates and processes.

The current code is a computational engine and mostly aims at generating MCMC samples from the posterior distribution associated with the model. Further analyses such as exploring MCMC samples, making predictions or plotting results can be performed with other tools such as the R package [RSTooDs](https://github.com/STooDs-tools/RSTooDs). 

### Usage
STooDs is implemented as an executable file controlled by a set of configuration text files. The easiest usage is to call STooDs at the command line by typing `./STooDs -wk "path/to/workspace/"` (Linux) or `STooDs -wk "path/to/workspace/"` (Windows), where the workspace is the folder containing all configuration files.

### Getting STooDs executable
A recent STooDs executable for your system (Windows or Linux) can be downloaded [here](https://github.com/STooDs-tools/RSTooDs/tree/main/inst/bin).

Alternatively, the executable can be recompiled from sources using the provided makefile or Code::Blocks project. Files from the following projects are needed for this purpose:

1. [BMSL](https://github.com/benRenard/BMSL)
2. [BaM](https://github.com/BaM-tools/BaM)
3. miniDMSL (available soon!)

### Managing STooDs configuration files

The easiest way to create the configuration files is to use the R package [RSTooDs](https://github.com/STooDs-tools/RSTooDs), which provides an user-friendly interface to manipulate data, specify a probabilistic model and estimate it: see the vignettes provided with this package.

It is also possible to create these configuration files by hand, or by implementing a user interface in a language other than R. The folder `examples` contains typical configuration files for a few simple case studies. A more complete documentation can be found in folder `doc`.

