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
![Figure1.1](/Markdown/11.png)

### 1.2 Different models by triage group
With the **by(tgroup)** option in Example 1.1 the conditional estimates in each group are identical. To fit different models and obtain seperate tables and graphs for each group, use instead the by prefix instead i.e **bysort tgroup:**  or **by tgroup:** if **tgroup** is already sorted. 

```
bys tgroup, rc0: metapreg num denom,  ///
	studyid(study) sumtable(abs) noitable  ///
	cimethod(wilson)  ///
	label(namevar=author, yearvar=year) catpplot nofplot  ///
	xlab(.25, 0.5, .75, 1)  ///
	subti(Atypical cervical cytology, size(4))  ///
	texts(1.5) smooth 
```
The option **rc0** ensures that the program runs in all groups even when there could be errors encountered in one of the sub-group analysis, without the option, the program stops running when the first error occurs in one of the groups.

This produces the following Stata results:
```
-> tgroup = ASCUS

**************************************************************************
    num ~ binomial(p, denom)
    logit(p) = mu + study
    study ~ N(0, tau2)


    Number of observations = 20
    Number of studies = 20


Click to show the raw estimates

****************************************************************************************

Test of heterogeneity - LR Test: RE model vs FE model

-----------------------------------------------------------------------
           |       DF       Chisq           p        tau2       I2tau  
-----------+-----------------------------------------------------------
     Model |        1      493.24        0.00        0.08       74.08  
-----------------------------------------------------------------------
NOTE: H0: Between-study variance(s) = 0  vs. H1: Between-study variance(s) > 0

****************************************************************************************
            Conditional summary: Proportion
****************************************************************************************

-------------------------------------------------------------------------------
 Parameter |       Mean  SE(logit)  t(logit)      P>|t|    t Lower    t Upper  
-----------+-------------------------------------------------------------------
   Overall |       0.42       0.02    -19.14       0.00       0.41       0.42  
-------------------------------------------------------------------------------
NOTE: H0: Est = 0.5 vs. H1: Est != 0.5

****************************************************************************************

Population-averaged estimates: Proportion 

--------------------------------------------------------------------------------
 Parameter |     Mean        SE    Median     Lower     Upper      Sample size  
-----------+--------------------------------------------------------------------
   Overall |     0.43      0.02      0.43      0.40      0.46              800  
--------------------------------------------------------------------------------
NOTE: % centiles obtained from 800 simulations of the posterior distribution

NOTE: caterpillar plot name -> abscatpplot1

```
![Figure1.2.1](/Markdown/121.png)
```
-> tgroup = ASC-US

**************************************************************************
    num ~ binomial(p, denom)
    logit(p) = mu + study
    study ~ N(0, tau2)


    Number of observations = 6
    Number of studies = 6


Click to show the raw estimates

****************************************************************************************

Test of heterogeneity - LR Test: RE model vs FE model

-----------------------------------------------------------------------
           |       DF       Chisq           p        tau2       I2tau  
-----------+-----------------------------------------------------------
     Model |        1       36.90        0.00        0.05       68.22  
-----------------------------------------------------------------------
NOTE: H0: Between-study variance(s) = 0  vs. H1: Between-study variance(s) > 0

****************************************************************************************
            Conditional summary: Proportion
****************************************************************************************

-------------------------------------------------------------------------------
 Parameter |       Mean  SE(logit)  t(logit)      P>|t|    t Lower    t Upper  
-----------+-------------------------------------------------------------------
   Overall |       0.37       0.03    -15.72       0.00       0.35       0.39  
-------------------------------------------------------------------------------
NOTE: H0: Est = 0.5 vs. H1: Est != 0.5

****************************************************************************************

Population-averaged estimates: Proportion 

--------------------------------------------------------------------------------
 Parameter |     Mean        SE    Median     Lower     Upper      Sample size  
-----------+--------------------------------------------------------------------
   Overall |     0.41      0.02      0.41      0.36      0.46              800  
--------------------------------------------------------------------------------
NOTE: % centiles obtained from 800 simulations of the posterior distribution

NOTE: caterpillar plot name -> abscatpplot2
```
![Figure1.2.2](/Markdown/122.png)
```
-> tgroup = BORDERLINE DYSKARYOSIS

**************************************************************************
    num ~ binomial(p, denom)
    logit(p) = mu + study
    study ~ N(0, tau2)


    Number of observations = 6
    Number of studies = 6


Click to show the raw estimates

****************************************************************************************

Test of heterogeneity - LR Test: RE model vs FE model

-----------------------------------------------------------------------
           |       DF       Chisq           p        tau2       I2tau  
-----------+-----------------------------------------------------------
     Model |        1        0.14        0.35        0.02       27.22  
-----------------------------------------------------------------------
NOTE: H0: Between-study variance(s) = 0  vs. H1: Between-study variance(s) > 0

****************************************************************************************
            Conditional summary: Proportion
****************************************************************************************

-------------------------------------------------------------------------------
 Parameter |       Mean  SE(logit)  t(logit)      P>|t|    t Lower    t Upper  
-----------+-------------------------------------------------------------------
   Overall |       0.35       0.13     -4.84       0.01       0.27       0.43  
-------------------------------------------------------------------------------
NOTE: H0: Est = 0.5 vs. H1: Est != 0.5

****************************************************************************************

Population-averaged estimates: Proportion 

--------------------------------------------------------------------------------
 Parameter |     Mean        SE    Median     Lower     Upper      Sample size  
-----------+--------------------------------------------------------------------
   Overall |     0.40      0.03      0.40      0.32      0.46              800  
--------------------------------------------------------------------------------
NOTE: % centiles obtained from 800 simulations of the posterior distribution

NOTE: caterpillar plot name -> abscatpplot3
```
![Figure1.2.3](/Markdown/123.png)

### 1.3 Stratified analysis by triage group
The use of by prefix **bysort tgroup, rc0:** in Example 1.2 produced three seperare plots and tables for each triage group. To combine the plots and the consolidate the tables into one use the options **stratify by(tgroup)** instead. 

```
metapreg num denom,  ///
	by(tgroup) stratify ///
	studyid(study) sumtable(abs) noitable  ///
	cimethod(wilson)  ///
	label(namevar=author, yearvar=year) catpplot nofplot  ///
	xlab(.25, 0.5, .75, 1)  ///
	subti(Atypical cervical cytology, size(4))  ///
	texts(1.5) smooth 
```
This produces the following Stata results:
```
*********************************** Model for :tgroup = ASCUS ***************************************
    num ~ binomial(p, denom)
    logit(p) = mu + study
    study ~ N(0, tau2)


    Number of observations = 20
    Number of studies = 20


Click to show the raw estimates

*********************************** Model for :tgroup = ASC-US ***************************************
    num ~ binomial(p, denom)
    logit(p) = mu + study
    study ~ N(0, tau2)


    Number of observations = 6
    Number of studies = 6


Click to show the raw estimates

*********************************** Model for :tgroup = BORDERLINE DYSKARYOSIS ***************************************
    num ~ binomial(p, denom)
    logit(p) = mu + study
    study ~ N(0, tau2)


    Number of observations = 6
    Number of studies = 6


Click to show the raw estimates

*********************************** Model for : all studies ***************************************
    num ~ binomial(p, denom)
    logit(p) = mu + study
    study ~ N(0, tau2)


    Number of observations = 32
    Number of studies = 32


Click to show the raw estimates

****************************************************************************************

Test of heterogeneity - LR Test: RE model vs FE model

------------------------------------------------------------------------------------------
                              |       DF       Chisq           p        tau2       I2tau  
------------------------------+-----------------------------------------------------------
                 tgroup|ASCUS |        1      493.24        0.00        0.08       74.08  
                tgroup|ASC-US |        1       36.90        0.00        0.05       68.22  
tgroup|BORDERLINE DYSKARYOSIS |        1       32.37        0.00        0.17       75.46  
                  All studies |        1      591.57        0.00        0.07       67.34  
------------------------------------------------------------------------------------------
NOTE: H0: Between-study variance(s) = 0  vs. H1: Between-study variance(s) &gt; 0

****************************************************************************************
            Conditional summary: Proportion
****************************************************************************************

--------------------------------------------------------------------------------------------------
                    Parameter |       Mean  SE(logit)  t(logit)      P&gt;|t|    t Lower    t Upper  
------------------------------+-------------------------------------------------------------------
tgroup|ASCUS                  |                                                                   
                      Overall |       0.42       0.02    -19.14       0.00       0.41       0.42  
tgroup|ASC-US                 |                                                                   
                      Overall |       0.37       0.03    -15.72       0.00       0.35       0.39  
tgroup|BORDERLINE DYSKARYOSIS |                                                                   
                      Overall |       0.40       0.07     -6.35       0.00       0.35       0.44  
All studies                   |                                                                   
                      Overall |       0.40       0.01    -29.02       0.00       0.39       0.41  
--------------------------------------------------------------------------------------------------
NOTE: H0: Est = 0.5 vs. H1: Est != 0.5

****************************************************************************************

Population-averaged estimates: Proportion 

---------------------------------------------------------------------------------------------------
                    Parameter |     Mean        SE    Median     Lower     Upper      Sample size  
------------------------------+--------------------------------------------------------------------
tgroup|ASCUS                  |                                                                    
                      Overall |     0.43      0.02      0.43      0.40      0.46              800  
tgroup|ASC-US                 |                                                                    
                      Overall |     0.40      0.02      0.40      0.35      0.45              800  
tgroup|BORDERLINE DYSKARYOSIS |                                                                    
                      Overall |     0.42      0.04      0.43      0.34      0.51              800  
All studies                   |                                                                    
                      Overall |     0.42      0.01      0.42      0.40      0.44              800  
---------------------------------------------------------------------------------------------------
NOTE: % centiles obtained from 800 simulations of the posterior distribution

```
![Figure1.3](/Markdown/13.png)

### 1.4 Proportions near 0
The dataset used in this example produced the top-left graph in figure two in Ioanna Tsoumpou et al. (2009). First, load the dataset:
```
use http://fmwww.bc.edu/repec/bocode/t/tsoumpou2009cancertreatrevfig2WNL.dta
```
Logistic regression correctly handles the extreme cases appropriately without need for transformation or continuity correction.
```
metapreg p16p p16tot, ///
                studyid(study) ///
                label(namevar=author, yearvar=year) catpplot nofplot noitable ///
                sortby(year author) ///
                xlab(0, .1, .2, 0.3, 0.4, 0.5) ///
                xline(0, lcolor(black)) ///
                ti(Positivity of p16 immunostaining, size(4) color(blue)) ///
                subti("Cytology = WNL", size(4) color(blue)) ///
                pointopt(msymbol(X) msize(2)) ///
                texts(1.5) smooth gof
```
This produces the following Stata results:
```
 p16p ~ binomial(p, p16tot)
    logit(p) = mu + study
    study ~ N(0, tau2)


    Number of observations = 19
    Number of studies = 19


Goodness of Fit Criterion

        |      AIC       BIC  
--------+---------------------
  Value |   110.70    112.59  


Click to show the raw estimates

****************************************************************************************

Test of heterogeneity - LR Test: RE model vs FE model

-----------------------------------------------------------------------
           |       DF       Chisq           p        tau2       I2tau  
-----------+-----------------------------------------------------------
     Model |        1      212.06        0.00        2.42       46.96  
-----------------------------------------------------------------------
NOTE: H0: Between-study variance(s) = 0  vs. H1: Between-study variance(s) &gt; 0

****************************************************************************************
            Conditional summary: Proportion
****************************************************************************************

-------------------------------------------------------------------------------
 Parameter |       Mean  SE(logit)  t(logit)      P&gt;|t|    t Lower    t Upper  
-----------+-------------------------------------------------------------------
   Overall |       0.07       0.43     -6.12       0.00       0.03       0.15  
-------------------------------------------------------------------------------
NOTE: H0: Est = 0.5 vs. H1: Est != 0.5

****************************************************************************************

Population-averaged estimates: Proportion 

--------------------------------------------------------------------------------
 Parameter |     Mean        SE    Median     Lower     Upper      Sample size  
-----------+--------------------------------------------------------------------
   Overall |     0.12      0.06      0.17      0.08      0.29              800  
--------------------------------------------------------------------------------
NOTE: % centiles obtained from 800 simulations of the posterior distribution
```
![Figure1.4](/Markdown/14.png)

### 1.5 Proportions near 0 - Loglog link
The loglog regression is an extension of the logistic regression model and is particularly useful when event probabilities are very close to zero like in example 1.4.
To use the loglog link, add the option **link(loglog)**
```
metapreg p16p p16tot, link(loglog) ///
	studyid(study) ///
	label(namevar=author, yearvar=year) catpplot nofplot noitable ///
	sortby(year author) sumtable(abs) ///
	xlab(0, .1, .2, 0.3, 0.4, 0.5) ///
	xline(0, lcolor(black)) ///
	ti(Positivity of p16 immunostaining, size(4) color(blue)) ///
	subti("Cytology = WNL", size(4) color(blue)) ///
	pointopt(msymbol(X) msize(2)) ///
	texts(1.5) smooth gof
```
This produces the following Stata results. Compared to the logit link, the model fit with the loglog link is slightly better as indicated by the lower BIC.
```
p16p ~ binomial(p, p16tot)
    loglog(p) = mu + study
    study ~ N(0, tau2)


    Number of observations = 19
    Number of studies = 19


Goodness of Fit Criterion

        |      AIC       BIC  
--------+---------------------
  Value |   108.93    110.81  


Click to show the raw estimates

****************************************************************************************

Test of heterogeneity - LR Test: RE model vs FE model

-----------------------------------------------------------------------
           |       DF       Chisq           p        tau2       I2tau  
-----------+-----------------------------------------------------------
     Model |        1      213.84        0.00        0.42              
-----------------------------------------------------------------------
NOTE: H0: Between-study variance(s) = 0  vs. H1: Between-study variance(s) &gt; 0

****************************************************************************************
            Conditional summary: Proportion
****************************************************************************************

-------------------------------------------------------------------------------
 Parameter |       Mean  SE(log~g)  t(loglog)      P&gt;|t|    t Lower    t Upper  
-----------+-------------------------------------------------------------------
   Overall |       0.08       0.17      5.38       0.00       0.03       0.17  
-------------------------------------------------------------------------------
NOTE: H0: Est = 0.5 vs. H1: Est != 0.5

****************************************************************************************

Population-averaged estimates: Proportion 

--------------------------------------------------------------------------------
 Parameter |     Mean        SE    Median     Lower     Upper      Sample size  
-----------+--------------------------------------------------------------------
   Overall |     0.12      0.05      0.16      0.07      0.26              787  
--------------------------------------------------------------------------------
NOTE: % centiles obtained from 800 simulations of the posterior distribution

```
![Figure1.5](/Markdown/15.png)

## Meta-regression
### 2.1 Independent studies
The use of **stratify by(tgroup)** in Example 1.3 only allows informal testing of differences between the sub-groups. The formal testing is perfomed by fitting a logistic regression with triage used as a categorical covariate.

First, load the dataset:
```
 use http://fmwww.bc.edu/repec/bocode/a/arbyn2009jcellmolmedfig1.dta
```
The option **sumtable(rd rr)** requests for summary tables of the risk differences and risk ratios, while **summaryonly** directs the program to present only the summary diamonds in the graph.
```
metapreg num denom tgroup, ///
	studyid(study) ///
	sumtable(rd rr) catpplot nofplot noitable ///
	cimethod(exact) ///
	label(namevar=author, yearvar=year) ///
	xlab(.25, 0.5, .75) ///
	subti(Atypical cervical cytology, size(4)) ///
	texts(1.5) summaryonly
```
This produces the following Stata results:
```
num ~ binomial(p, denom)
    logit(p) = mu + tgroup + study
    study ~ N(0, tau2)

    Base levels

        Variable -- Base Level
        tgroup -- ASCUS


    Number of observations = 32
    Number of studies = 32


Click to show the raw estimates

****************************************************************************************

Test of heterogeneity - LR Test: RE model vs FE model

-----------------------------------------------------------
           |       DF       Chisq           p        tau2  
-----------+-----------------------------------------------
     Model |        1      561.75        0.00        0.08  
-----------------------------------------------------------
NOTE: H0: Between-study variance(s) = 0  vs. H1: Between-study variance(s) &gt; 0

****************************************************************************************
            Conditional summary: Proportion Difference
****************************************************************************************

-------------------------------------------------------------------------------------------
             Parameter |       Mean         SE         t      P&gt;|t|    t Lower    t Upper  
-----------------------+-------------------------------------------------------------------
tgroup                 |                                                                   
                 ASCUS |       0.00       0.00                            0.00       0.00  
                ASC-US |       0.01       0.01      1.84       0.08      -0.00       0.03  
BORDERLINE DYSKARYOSIS |      -0.04       0.01     -4.32       0.00      -0.06      -0.02  
-------------------------------------------------------------------------------------------

****************************************************************************************

Wald-type test for nonlinear hypothesis

    H0: All RD equal vs. H1: Some RD different

-------------------------------------------
 Parameter |     chi2        df         p  
-----------+-------------------------------
    tgroup |    27.97         2      0.00  
-------------------------------------------

****************************************************************************************

Population-averaged estimates: Proportion Difference 

--------------------------------------------------------------------------------------------
             Parameter |     Mean        SE    Median     Lower     Upper      Sample size  
-----------------------+--------------------------------------------------------------------
tgroup                 |                                                                    
                 ASCUS |     0.00      0.00      0.00      0.00      0.00                   
                ASC-US |     0.01      0.03      0.01     -0.05      0.08              800  
BORDERLINE DYSKARYOSIS |    -0.01      0.03     -0.01     -0.08      0.05              800  
--------------------------------------------------------------------------------------------
NOTE: % centiles obtained from 800 simulations of the posterior distribution

****************************************************************************************
            Conditional summary: Proportion Ratio
****************************************************************************************

-------------------------------------------------------------------------------------------
             Parameter |       Mean    SE(lrr)    t(lrr)      P&gt;|t|    t Lower    t Upper  
-----------------------+-------------------------------------------------------------------
tgroup                 |                                                                   
                 ASCUS |       1.00       0.00                            1.00       1.00  
                ASC-US |       0.96       0.02     -1.82       0.08       0.93       1.00  
BORDERLINE DYSKARYOSIS |       1.10       0.02      4.43       0.00       1.05       1.14  
-------------------------------------------------------------------------------------------
NOTE: H0: Est = 1 vs. H1: Est != 1

****************************************************************************************

Wald-type test for nonlinear hypothesis

    H0: All (log)RR equal vs. H1: Some (log)RR different

-------------------------------------------
 Parameter |     chi2        df         p  
-----------+-------------------------------
    tgroup |    29.34         2      0.00  
-------------------------------------------

****************************************************************************************

Population-averaged estimates: Proportion Ratio 

--------------------------------------------------------------------------------------------
             Parameter |     Mean        SE    Median     Lower     Upper      Sample size  
-----------------------+--------------------------------------------------------------------
tgroup                 |                                                                    
                 ASCUS |     1.00      0.00      1.00      1.00      1.00                   
                ASC-US |     0.97      0.08      0.97      0.82      1.12              800  
BORDERLINE DYSKARYOSIS |     1.03      0.08      1.02      0.88      1.18              800  
--------------------------------------------------------------------------------------------
NOTE: % centiles obtained from 800 simulations of the posterior distribution
```
![Figure2.1](/Markdown/21.png)

The results are interpreted as follows: conditional on the random-effect, there are differences between the groups, however, these differences disappear at the population level.

### 2.2 Comparative studies
#### Frequentist vs Bayesian estimation
In the frequentist approach, maximum likelihood estimation is used to obtain the model parameters. Maximum likelihood algorithms often struggle with convergence when event probabilities approach zero, as the likelihood function becomes flat, leading to large or undefined estimates and standard errors. They also struggle with estimation of between-study variance components when only a few studies are available. 

The Bayesian approach incoorporates weakly informative priors based on plausible assumptions to overcome non-convergence issues in extreme sparsity or zero-event scenarios. In such cases, the prior distribution provides necessary information to stabilize estimates, yielding more reasonable and statistically defensible results even with limited data. The inverse-gamma (scale and shape = 0.01) is used as the default prior for the between-study variance. This prior is conditionally conjugate for the variance components thereby enhancing computational efficiency and numerical stability.

Bayesian computations take increasingly longer time with increasing number of studies. However, Bayesian estimation is not always essential for large datasets, where frequentist estimates can be sufficient.  

### 2.2.1 Non-convergence issues in frequentist estimation when proportions near 0
Hemkens et al. (2016) investigated the risk of fatal stroke after long-term colchicine use.  Three out of four studies contained double-zero events. First, load the dataset:
```
use https://github.com/VNyaga/Metapreg/blob/master/Markdown/hemkens2016analysis110.dta?raw=1
```
Then sort the data in the right order. To correctly use the option **design(comparative)** there should be two rows of data per each **studyid**. The first row has the index/treatment data and the second row has the control data.

```
gsort study -treatment
```
The required varlist has the form **n N bicat** where **bicat** is the first covariate which should be a string or labelled numeric variable with two levels.
The option **cov(independent)** in *design(comparative, cov(independent))* requests for a mixed-effects logistic regression model assuming event probabilities (logit scale) in the control groups and treatment effects (log ORs) are independently normally distributed.
```
metapreg event total treatment,  ///
	studyid(study)  ///
	design(comparative, cov(independent))  ///
	smooth gof catpplot nofplot noitable ///
	outplot(rr) xline(1) sumstat(Risk Ratio)  ///
	xlab(0, 1, 30) logscale  ///
	texts(2.35) astext(80)   
```
This produces the following Stata results:
```
events ~ binomial(p, total)
    logit(p) = mu + treatment + treatment.study + study
    study ~ N(0, tau2)
    treatment.study ~ N(0, sigma2)

    Base levels

        Variable -- Base Level
        treatment -- Control


    Number of observations = 8
    Number of studies = 4


Goodness of Fit Criterion

        |      AIC       BIC  
--------+---------------------
  Value |     6.89      6.97  


Click to show the raw estimates

****************************************************************************************

Test of heterogeneity - LR Test: RE model vs FE model

-----------------------------------------------------------------------------------------------
           |       DF       Chisq           p        tau2      sigma2       I2tau     I2sigma  
-----------+-----------------------------------------------------------------------------------
     Model |        2                                0.00        0.00       98.18        1.82  
-----------------------------------------------------------------------------------------------
NOTE: H0: Between-study variance(s) = 0  vs. H1: Between-study variance(s) &gt; 0

****************************************************************************************
            Conditional summary: Proportion Difference
****************************************************************************************

-------------------------------------------------------------------------------
 Parameter |       Mean         SE         t      P&gt;|t|    t Lower    t Upper  
-----------+-------------------------------------------------------------------
   Overall |      -0.00       0.00     -1.00       0.42      -0.01       0.01  
-------------------------------------------------------------------------------

****************************************************************************************

Population-averaged estimates: Proportion Difference 

--------------------------------------------------------------------------------
 Parameter |     Mean        SE    Median     Lower     Upper      Sample size  
-----------+--------------------------------------------------------------------
   Overall |    -0.00      0.00     -0.00     -0.01     -0.00              800  
--------------------------------------------------------------------------------
NOTE: % centiles obtained from 800 simulations of the posterior distribution

****************************************************************************************
            Conditional summary: Proportion Ratio
****************************************************************************************

-------------------------------------------------------------------------------
 Parameter |       Mean    SE(lrr)    t(lrr)      P>|t|    t Lower    t Upper  
-----------+-------------------------------------------------------------------
   Overall |  143434.42       1.00     11.89       0.01    1949.92   1.06e+07  
-------------------------------------------------------------------------------
NOTE: H0: Est = 1 vs. H1: Est != 1

****************************************************************************************

Population-averaged estimates: Proportion Ratio 

--------------------------------------------------------------------------------
 Parameter |     Mean        SE    Median     Lower     Upper      Sample size  
-----------+--------------------------------------------------------------------
   Overall |  1.4e+05   2.5e+05   1.4e+05  19408.17   8.7e+05              800  
--------------------------------------------------------------------------------
NOTE: % centiles obtained from 800 simulations of the posterior distribution
```
![Figure2.2.1f](/Markdown/221f.png)
       
The estimation results of the fitted model are stored as **metapreg_modest**. These can be requested:
```
estimates replay metapreg_modest
```
This replays the following results:
```

Mixed-effects logistic regression               Number of obs     =          8
Binomial variable: total
Group variable: study                           Number of groups  =          4

                                                Obs per group:
                                                              min =          2
                                                              avg =        2.0
                                                              max =          2

Integration method: mvaghermite                 Integration pts.  =          7

                                                Wald chi2(1)      =     140.74
Log likelihood = -2.4437626                     Prob > chi2       =     0.0000
 ( 1)  [events]_cons = 0
----------------------------------------------------------------------------------
          events | Coefficient  Std. err.      z    P>|z|     [95% conf. interval]
-----------------+----------------------------------------------------------------
              mu |  -18.03693          .        .       .            .           .
                 |
       treatment |
     Colchicine  |   11.87574   1.001047    11.86   0.000     9.913725    13.83776
           _cons |          0  (omitted)
-----------------+----------------------------------------------------------------
study            |
 var(2.treatment)|   1.41e-31   2.32e-15                             .           .
          var(mu)|   7.62e-30          .                             .           .
----------------------------------------------------------------------------------
Warning: Convergence not achieved.
```
The model does not converge, the point estimates and standard errors are excessively large or undefined. The point estimates for the between-study variance components are very close to zero. 

To switch to Bayesian estimation, add the options **inference(bayesian) bwd(path)** where path is a directory path specifying the location where the simulation results containing MCMC samples should be saved. 
```
global wdir "C:\DATA\WIV\Projects\Stata\Metapreg\mcmcresults"

metapreg event total treatment,  ///
	studyid(study)  ///
	design(comparative, cov(independent))  ///
	inference(bayesian) bwd($wdir)  ///
	smooth gof sumtable(rd rr) catpplot nofplot noitable  ///
	outplot(rr) xline(1) sumstat(Risk Ratio)  ///
	xlab(0, 1, 30) logscale  ///
	texts(2.35) astext(80)
```
This procduces the following Stata results:
```
    events ~ binomial(p, total)
    logit(p) = mu + treatment + treatment.study + study
    study ~ N(0, tau2)
    treatment.study ~ N(0, sigma2)

    Base levels

        Variable -- Base Level
        treatment -- Control


    Number of observations = 8
    Number of studies = 4


Goodness of Fit Criterion

        |  Avg DIC      Avg log(ML)  
--------+----------------------------
  Value |     7.77            23.82  


Click to show the raw estimates

****************************************************************************************

Test of heterogeneity - LR Test: RE model vs FE model

-----------------------------------------------------------------------------------------------
           | Delta ML     log(BF)    Post P~b        tau2      sigma2       I2tau     I2sigma  
-----------+-----------------------------------------------------------------------------------
     Model |      -31       31.02        1.00        0.81        1.17       41.02       58.98  
-----------------------------------------------------------------------------------------------
NOTE: H0: Between-study variance(s) = 0  vs. H1: Between-study variance(s) > 0

****************************************************************************************
            Conditional summary: Proportion Difference
****************************************************************************************

-------------------------------------------------------------------------------
 Parameter |       Mean         SD      MCSE     Median  eti Lower  eti Upper  
-----------+-------------------------------------------------------------------
   Overall |      -0.02       0.11      0.00      -0.00      -0.37       0.02  
-------------------------------------------------------------------------------

****************************************************************************************

Population-averaged estimates: Proportion Difference 

--------------------------------------------------------------------------------
 Parameter |     Mean        SE    Median     Lower     Upper      Sample size  
-----------+--------------------------------------------------------------------
   Overall |    -0.00      0.00     -0.00     -0.01      0.00             7500  
--------------------------------------------------------------------------------
NOTE: % centiles obtained from 7500 simulations of the posterior distribution

****************************************************************************************
            Conditional summary: Proportion Ratio
****************************************************************************************

-------------------------------------------------------------------------------
 Parameter |       Mean         SD      MCSE     Median  eti Lower  eti Upper  
-----------+-------------------------------------------------------------------
   Overall |       5.34      21.61      0.91       1.18       0.01      32.79  
-------------------------------------------------------------------------------

****************************************************************************************

Population-averaged estimates: Proportion Ratio 

--------------------------------------------------------------------------------
 Parameter |     Mean        SE    Median     Lower     Upper      Sample size  
-----------+--------------------------------------------------------------------
   Overall |  6.7e+05   5.5e+07      2.52      0.09    482.71             7500  
--------------------------------------------------------------------------------
NOTE: % centiles obtained from 7500 simulations of the posterior distribution

        | The bayesian estimation commands saved a dataset C:\DATA\WIV\Projects\Stata\Metapreg\mcmcresults\metapreg_bayesest.dta
        | C:\DATA\WIV\Projects\Stata\Metapreg\mcmcresults\metapreg_bayesest.dta
        | containing the MCMC samples of the parameters to the disk.
        | It is your responsibility to erase the dataset
        | after it is no longer needed.
        Click to erase the dataset

```
The Bayesian estimate of the between-study variance components are 0.81 (variance of control group event probabilities in the logit scale) and 1.17 (treatment-effects in log OR scale).  

![Figure2.2.1b](/Markdown/221b.png)

The population-averaged summary RR is 2.52 (0.09, 483.71), which is smaller and closer to 1 than the frequentist estimate 1.4e+05(19408.17, 8.7e+05).

