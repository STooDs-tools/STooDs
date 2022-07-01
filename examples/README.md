## How to run these examples?
1. Get the latest STooDs executable for your system (Windows or Linux) [here](https://github.com/STooDs-tools/RSTooDs/tree/main/inst/bin) and save it into the current 'examples' folder.
2. Run STooDs with the workspace argument '-wk' pointing to the example folder of your choice; for instance `./STooDs -wk "Ex1_iid/"` (Linux) or `STooDs -wk "Ex1_iid/"` (Windows).

## Quick description of the examples

### Example 1: simplest iid model

**Data**: average spring streamflow (in $m^3.s^{-1}$) of the Barnard River measured at the [Barry station](http://www.bom.gov.au/water/hrs/#id=208009), 1951-2014. 

**Model**: data $y_t$ are assumed to be independant and identically distributed realizations from a log-normal distribution.
A flat prior distribution is used for parameter $\mu$,
a lognormal prior for the parameter $\sigma$.

$y_t \sim LN(\mu,\sigma)$

$\mu \sim U[-\infty;+\infty]$

$\sigma \sim LN(1,2)$

### Example 2: effect of a covariate

**Data**: same as in example 1 + average spring values of the standardized [nino3.4](https://psl.noaa.gov/gcos_wgsp/Timeseries/Nino34/) index. High values of this index denote [El Niño](https://en.wikipedia.org/wiki/El_Niño) episodes, low values denote [La Niña](https://en.wikipedia.org/wiki/La_Niña) episodes.

**Model**: data $y_t$ are assumed to be independant realizations from a log-normal distribution 
whose parameter $\mu$ varies as a function of 
the nino3.4 index $x_t$:

$y_t \sim LN(\mu_0+\mu_1*x_t,\sigma)$

$\mu_0 \sim U[-\infty;+\infty]$

$\mu_1 \sim U[-\infty;+\infty]$

$\sigma \sim LN(1,2)$

### Example 3: defining a (spatial) process

**Data**: average spring streamflow (in $m^3.s^{-1}$) for 21 rivers located in Northern New South Wales and Southern Queensland, Australia. Data can be retrieved on the Bureau Of Meteorology [website](http://www.bom.gov.au/water/hrs/). Average spring values of the standardized nino3.4 index are also used as in the previous example.

**Model**: the model can be viewed as a duplication of the model of example 2 at each site: data $y_{s,t}$ 
at site $s$ 
and time $t$ 
are assumed to be realizations from log-normal distributions, 
parameter $\mu$ varies as a function of 
the nino3.4 index $x_t$. 
The parameters $\mu_{0,s}$, 
$\mu_{1,s}$ and 
$\sigma_s$ are all site-specific, but distinct prior distributions are used:

1. For the parameters $\mu_{0,s} \text{ and } \sigma_s$, the same prior distributions as in previous example 2 are used.
2. The parameters $\mu_{1,s}$, controling the effect of the El Niño, are assumed to be realizations from a Gaussian distribution with _unknown_ parameters *m* and *s*. Parameter *m* quantifies the overall El Niño effect in the region,  parameter *s* quantifies the spatial variability of this effect.

$y_{s,t} \sim LN(\mu_{0,s}+\mu_{1,s}*x_t,\sigma_s)$

$\mu_{0,s} \sim U[-\infty;+\infty]$

$\mu_{1,s} \sim N(m,s)$

$\sigma_s \sim LN(1,2)$

