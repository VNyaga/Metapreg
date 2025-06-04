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
## Intercept-only model
### 1.1 Summary by triage group
The dataset used in examples 1.1-1.2 was used previously to produce Figure 1 in Marc Arbyn et al.(2009). First, load the dataset:
```
use http://fmwww.bc.edu/repec/bocode/a/arbyn2009jcellmolmedfig1.dta
```
Then fit a mixed model to the data. The options **by(tgroup)** requests to report summary estimates grouped by triage group, **catpplot** for a caterpilar plot i.e the forest plot sorted by the magnitude of the empirical proportions, and **smooth** lays the study-specific model-based estimates over the emprical estimates in the plot.

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
This produces the following Stata results:

```
    num ~ binomial(p, denom)
    logit(p) = mu + study
    study ~ N(0, tau2)


    Number of observations = 32
    Number of studies = 32


Goodness of Fit Criterion

        |      AIC       BIC  
--------+---------------------
  Value |   345.34    348.27  


Click to show the raw estimates

****************************************************************************************

Test of heterogeneity - LR Test: RE model vs FE model

-----------------------------------------------------------------------
           |       DF       Chisq           p        tau2       I2tau  
-----------+-----------------------------------------------------------
     Model |        1      591.57        0.00        0.07       67.34  
-----------------------------------------------------------------------
NOTE: H0: Between-study variance(s) = 0  vs. H1: Between-study variance(s) > 0

****************************************************************************************
            Conditional summary: Proportion
****************************************************************************************

-------------------------------------------------------------------------------------------
             Parameter |       Mean  SE(logit)  t(logit)      P>|t|    t Lower    t Upper  
-----------------------+-------------------------------------------------------------------
tgroup                 |                                                                   
                 ASCUS |       0.40       0.01    -29.02       0.00       0.39       0.41  
                ASC-US |       0.40       0.01    -29.02       0.00       0.39       0.41  
BORDERLINE DYSKARYOSIS |       0.40       0.01    -29.02       0.00       0.39       0.41  
                       |                                                                   
               Overall |       0.40       0.01    -29.02       0.00       0.39       0.41  
-------------------------------------------------------------------------------------------
NOTE: H0: Est = 0.5 vs. H1: Est != 0.5

****************************************************************************************

Population-averaged estimates: Proportion 

--------------------------------------------------------------------------------------------
             Parameter |     Mean        SE    Median     Lower     Upper      Sample size  
-----------------------+--------------------------------------------------------------------
tgroup                 |                                                                    
                 ASCUS |     0.42      0.01      0.42      0.39      0.45              800  
                ASC-US |     0.41      0.02      0.41      0.36      0.46              800  
BORDERLINE DYSKARYOSIS |     0.42      0.03      0.42      0.37      0.48              800  
                       |                                                                    
               Overall |     0.42      0.01      0.42      0.40      0.44              800  
--------------------------------------------------------------------------------------------
NOTE: % centiles obtained from 800 simulations of the posterior distribution
```
![Figure1.1](/Build/11.png)

### 1.2 Different models by triage group
With the **by(tgroup)** option in Example 1.1 the conditional estimates in each group are identical.  To fit different models and obtain seperate tables and graphs for each group, use instead the by prefix instead i.e **bysort tgroup:**  or **by tgroup:** if **tgroup** is already sorted. The option **rc0** ensures that the program runs in all groups even when there could be errors encountered in one of the sub-group analysis, without the option, the program stops running when the first error occurs in one of the groups.


