# Metapreg 
 A routing for meta-analysis, meta-regression and network meta-analysis of proportions from binomial data using generalized linear models in the frequentist and Bayesian approach.
 The program fits a fixed, random-effects, or a mixed-effects model assuming binomial, common-rho beta-binomial or binomial-normal distribution with a logit, loglog or the cloglog link to the data.

The command requires Stata 14.1 or later versions. Stata 16.1 or later versions is required to perform the bayesian analysis.

To install the from SSC, type
```
ssc install metapreg
```

To install the development version directly, type
```
net install metapreg, from("https://raw.githubusercontent.com/VNyaga/Metapreg/master/Build/")
```
## Intercept-only model 1.1 
### Summary by triage group
The dataset used in examples 1.1-1.2 was used previously to produce Figure 1 in Marc Arbyn et al.(2009). First, load the dataset:
```
use http://fmwww.bc.edu/repec/bocode/a/arbyn2009jcellmolmedfig1.dta
```
Then fit a mixed model to the data. 
```
metapreg num denom,  ///
	studyid(study) ///
	by(tgroup) sumtable(abs) noitable /// 
	cimethod(exact) ///
	label(namevar=author, yearvar=year) catpplot nofplot ///
	xlab(.25, 0.5, .75, 1) ///
	subti(Atypical cervical cytology, size(4)) ///
	texts(1.5) smooth gof
```

The options **by(tgroup)** requests to report summary estimates grouped by triage group, **catpplot** for a caterpilar plot i.e the forest plot sorted by the magnitude of the empirical proportions, and **smooth** lays the study-specific model-based estimates over the emprical estimates in the plot.

The graph is 
![](11.png)
