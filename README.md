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

        | The bayesian estimation commands saved a dataset 
        | C:\DATA\WIV\Projects\Stata\Metapreg\mcmcresults\metapreg_bayesest.dta
        | containing the MCMC samples of the parameters to the disk.
        | It is your responsibility to erase the dataset
        | after it is no longer needed.
        Click to erase the dataset

```
The Bayesian estimate of the between-study variance components are 0.81 (variance of control group event probabilities in the logit scale) and 1.17 (treatment-effects in log OR scale).  

![Figure2.2.1b](/Markdown/221b.png)

The population-averaged summary RR is 2.52 (0.09, 483.71), which is smaller and closer to 1 than the frequentist estimate 1.4e+05(19408.17, 8.7e+05).

### 2.2.2 Frequentist approach under-estimates the between-study variance in meta-analysis of three studies
Bender et al. (2018) conducted a meta-analysis of three studies evaluating the risk of fever following sipuleucel-T therapy in asymptomatic or minimally symptomatic metastatic castrate-resistant prostate cancer. First, load the dataset:
```
use "https://github.com/VNyaga/Metapreg/blob/master/Build/bender2018fig2.dta?raw=1"
```
Fit a mixed-effects logistic regression model assuming event probabilities (logit scale) in the control groups and treatment effects (log ORs) are independently normally distributed i.e.
**design(comparative, cov(independent))**
```
metapreg event total treatment,  ///
	studyid(study) design(comparative, cov(independent))  ///
	smooth gof catpplot nofplot noitable cimethod(,wald)  ///
	outplot(rr) xline(1) sumstat(Risk Ratio)  ///
	xlab(1, 5, 30) logscale  ///
	texts(2) astext(70)
```
This produces the following Stata results:
```
event ~ binomial(p, total)
    logit(p) = mu + treatment + treatment.study + study
    study ~ N(0, tau2)
    treatment.study ~ N(0, sigma2)

    Base levels

        Variable -- Base Level
        treatment -- Control


    Number of observations = 6
    Number of studies = 3


Goodness of Fit Criterion

        |      AIC       BIC  
--------+---------------------
  Value |    30.41     29.99  


Click to show the raw estimates

****************************************************************************************

Test of heterogeneity - LR Test: RE model vs FE model

-----------------------------------------------------------------------------------------------
           |       DF       Chisq           p        tau2      sigma2       I2tau     I2sigma  
-----------+-----------------------------------------------------------------------------------
     Model |        2        3.58                    0.00        0.00       46.98       53.02  
-----------------------------------------------------------------------------------------------
NOTE: H0: Between-study variance(s) = 0  vs. H1: Between-study variance(s) > 0

****************************************************************************************
            Conditional summary: Proportion Difference
****************************************************************************************

-------------------------------------------------------------------------------
 Parameter |       Mean         SE         z      P>|z|    z Lower    z Upper  
-----------+-------------------------------------------------------------------
   Overall |      -0.19       0.03     -6.39       0.00      -0.24      -0.13  
-------------------------------------------------------------------------------

****************************************************************************************

Population-averaged estimates: Proportion Difference 

--------------------------------------------------------------------------------
 Parameter |     Mean        SE    Median     Lower     Upper      Sample size  
-----------+--------------------------------------------------------------------
   Overall |    -0.19      0.03     -0.19     -0.24     -0.13              800  
--------------------------------------------------------------------------------
NOTE: % centiles obtained from 800 simulations of the posterior distribution

****************************************************************************************
            Conditional summary: Proportion Ratio
****************************************************************************************

-------------------------------------------------------------------------------
 Parameter |       Mean    SE(lrr)    z(lrr)      P>|z|    z Lower    z Upper  
-----------+-------------------------------------------------------------------
   Overall |       2.62       0.19      5.05       0.00       1.80       3.81  
-------------------------------------------------------------------------------
NOTE: H0: Est = 1 vs. H1: Est != 1

****************************************************************************************

Population-averaged estimates: Proportion Ratio 

--------------------------------------------------------------------------------
 Parameter |     Mean        SE    Median     Lower     Upper      Sample size  
-----------+--------------------------------------------------------------------
   Overall |     2.62      0.50      2.63      1.83      3.71              800  
--------------------------------------------------------------------------------
NOTE: % centiles obtained from 800 simulations of the posterior distribution

```
The frequentist estimates for the between-study variance components are both very close to zero reducing the model to fixed-effects logistic regression.
![Figure 2.2.2f](/Markdown/222f.png)
The population-averaged summary RR is 2.63 (1.83, 3.71).

Switch to Bayesian estimation by adding the options **inference(bayesian) and bwd($wdir)**. The global macro **wdir** was declared earlier.
```
metapreg event total treatment,  inference(bayesian) bwd($wdir) ///
	studyid(study) design(comparative, cov(independent))  ///
	smooth gof catpplot nofplot noitable cimethod(,wald)  ///
	outplot(rr) xline(1) sumstat(Risk Ratio)  ///
	xlab(1, 5, 30) logscale  ///
	texts(2) astext(70)
```
This produces the following Stata results:
```
event ~ binomial(p, total)
    logit(p) = mu + treatment + treatment.study + study
    study ~ N(0, tau2)
    treatment.study ~ N(0, sigma2)

    Base levels

        Variable -- Base Level
        treatment -- Control


    Number of observations = 6
    Number of studies = 3


Goodness of Fit Criterion

        |  Avg DIC      Avg log(ML)  
--------+----------------------------
  Value |    36.24           -15.03  


Click to show the raw estimates

****************************************************************************************

Test of heterogeneity - LR Test: RE model vs FE model

-----------------------------------------------------------------------------------------------
           | Delta ML     log(BF)    Post P~b        tau2      sigma2       I2tau     I2sigma  
-----------+-----------------------------------------------------------------------------------
     Model |       -6        6.45        1.00        0.05        0.06       44.31       55.69  
-----------------------------------------------------------------------------------------------
NOTE: H0: Between-study variance(s) = 0  vs. H1: Between-study variance(s) > 0

****************************************************************************************
            Conditional summary: Proportion Difference
****************************************************************************************

-------------------------------------------------------------------------------
 Parameter |       Mean         SD      MCSE     Median  eti Lower  eti Upper  
-----------+-------------------------------------------------------------------
   Overall |      -0.19       0.07      0.00      -0.19      -0.34      -0.05  
-------------------------------------------------------------------------------

****************************************************************************************

Population-averaged estimates: Proportion Difference 

--------------------------------------------------------------------------------
 Parameter |     Mean        SE    Median     Lower     Upper      Sample size  
-----------+--------------------------------------------------------------------
   Overall |    -0.19      0.03     -0.19     -0.25     -0.13             7500  
--------------------------------------------------------------------------------
NOTE: % centiles obtained from 7500 simulations of the posterior distribution

****************************************************************************************
            Conditional summary: Proportion Ratio
****************************************************************************************

-------------------------------------------------------------------------------
 Parameter |       Mean         SD      MCSE     Median  eti Lower  eti Upper  
-----------+-------------------------------------------------------------------
   Overall |       2.76       0.82      0.04       2.67       1.41       4.70  
-------------------------------------------------------------------------------

****************************************************************************************

Population-averaged estimates: Proportion Ratio 

--------------------------------------------------------------------------------
 Parameter |     Mean        SE    Median     Lower     Upper      Sample size  
-----------+--------------------------------------------------------------------
   Overall |     2.84      0.65      2.75      1.90      4.35             7500  
--------------------------------------------------------------------------------
NOTE: % centiles obtained from 7500 simulations of the posterior distribution

```
The Bayesian estimate of the between-study variance components are 0.05 (variance of control group event probabilities in the logit scale) and 0.06 (treatment-effects in log OR scale). 
![Figure 2.2.2b](/Markdown/222b.png)
The population-averaged summary RR is 2.75 (1.90, 4.35). The credible interval is wider than the frequentist confidence interval. 

### 2.2.3 Frequentist approach under-estimates the between-study variance in meta-analysis of sparse studies
Hemkens et al. (2016) investigated the risk of cardiovascular mortality after long-term colchicine use.  The meta-analysis included seven studies, two with double-zero events and three with single-zero events. First, load the dataset and sort it accordingly:
```
use "https://github.com/VNyaga/Metapreg/blob/master/Build/hemkens2016analysis18.dta?raw=1"

gsort study -treatment
```
Fit a mixed-effects logistic regression model assuming event probabilities (logit scale) in the control groups are homogeneous and treatment effects (log ORs) are normally distributed i.e.
**design(comparative, cov(commonint))**
```
metapreg event total treatment, ///
	studyid(study) ///
	design(comparative, cov(commonint)) ///
	smooth gof catpplot nofplot sumtable(rd rr) noitable ///
	outplot(rr) xline(1) sumstat(Risk Ratio) ///
	xlab(0.01, 1, 100) logscale ///
	texts(1.75) astext(60)
```
This produces the following Stata results:
```
 events ~ binomial(p, total)
    logit(p) = mu + treatment + treatment.study + study
    treatment.study ~ N(0, sigma2)

    Base levels

        Variable -- Base Level
        treatment -- Control


    Number of observations = 14
    Number of studies = 7


Goodness of Fit Criterion

        |      AIC       BIC  
--------+---------------------
  Value |    34.46     35.73  


Click to show the raw estimates

****************************************************************************************

Test of heterogeneity - LR Test: RE model vs FE model

-----------------------------------------------------------------------------------------------
           |       DF       Chisq           p        tau2      sigma2       I2tau     I2sigma  
-----------+-----------------------------------------------------------------------------------
     Model |        1        0.00                    0.00        0.00        0.00      100.00  
-----------------------------------------------------------------------------------------------
NOTE: H0: Between-study variance(s) = 0  vs. H1: Between-study variance(s) > 0

****************************************************************************************
            Conditional summary: Proportion Difference
****************************************************************************************

-------------------------------------------------------------------------------
 Parameter |       Mean         SE         t      P>|t|    t Lower    t Upper  
-----------+-------------------------------------------------------------------
   Overall |       0.02       0.01      2.78       0.02       0.00       0.04  
-------------------------------------------------------------------------------

****************************************************************************************

Population-averaged estimates: Proportion Difference 

--------------------------------------------------------------------------------
 Parameter |     Mean        SE    Median     Lower     Upper      Sample size  
-----------+--------------------------------------------------------------------
   Overall |     0.02      0.01      0.02      0.01      0.04              800  
--------------------------------------------------------------------------------
NOTE: % centiles obtained from 800 simulations of the posterior distribution

****************************************************************************************
            Conditional summary: Proportion Ratio
****************************************************************************************

-------------------------------------------------------------------------------
 Parameter |       Mean    SE(lrr)    t(lrr)      P>|t|    t Lower    t Upper  
-----------+-------------------------------------------------------------------
   Overall |       0.20       0.63     -2.54       0.03       0.05       0.84  
-------------------------------------------------------------------------------
NOTE: H0: Est = 1 vs. H1: Est != 1

****************************************************************************************

Population-averaged estimates: Proportion Ratio 

--------------------------------------------------------------------------------
 Parameter |     Mean        SE    Median     Lower     Upper      Sample size  
-----------+--------------------------------------------------------------------
   Overall |     0.20      0.16      0.20      0.06      0.63              800  
--------------------------------------------------------------------------------
NOTE: % centiles obtained from 800 simulations of the posterior distribution
```
The frequentist estimate of the treatment effect between-study variance is very close zero (4.69e-33) reducing the model to fixed-effects logistic regression. 
![Figure 2.2.3f](/Markdown/223f.png)
The population-averaged summary RR is 0.20 (0.06, 0.63).

To switch to Bayesian estimation, add the options **inference(bayesian) and bwd($wdir)**. 
```
 metapreg event total treatment, ///
	studyid(study) ///
	design(comparative, cov(commonint)) ///
	inference(bayesian) bwd($wdir) ///
	smooth gof catpplot nofplot sumtable(rd rr) noitable ///
	outplot(rr) sumstat(Risk Ratio) ///
	xlab(0.01, 1, 100) logscale ///
	xline(1) texts(1.75) astext(60)
```
This produces the following Stata results:
```
events ~ binomial(p, total)
    logit(p) = mu + treatment + treatment.study + study
    treatment.study ~ N(0, sigma2)

    Base levels

        Variable -- Base Level
        treatment -- Control


    Number of observations = 14
    Number of studies = 7


Goodness of Fit Criterion

        |  Avg DIC      Avg log(ML)  
--------+----------------------------
  Value |    34.67            -4.66  


Click to show the raw estimates

****************************************************************************************

Test of heterogeneity - LR Test: RE model vs FE model

-----------------------------------------------------------------------------------------------
           | Delta ML     log(BF)    Post P~b        tau2      sigma2       I2tau     I2sigma  
-----------+-----------------------------------------------------------------------------------
     Model |      -15       15.44        1.00        0.00        0.15        0.00      100.00  
-----------------------------------------------------------------------------------------------
NOTE: H0: Between-study variance(s) = 0  vs. H1: Between-study variance(s) > 0

****************************************************************************************
            Conditional summary: Proportion Difference
****************************************************************************************

-------------------------------------------------------------------------------
 Parameter |       Mean         SD      MCSE     Median  eti Lower  eti Upper  
-----------+-------------------------------------------------------------------
   Overall |       0.02       0.01      0.00       0.02       0.01       0.04  
-------------------------------------------------------------------------------

****************************************************************************************

Population-averaged estimates: Proportion Difference 

--------------------------------------------------------------------------------
 Parameter |     Mean        SE    Median     Lower     Upper      Sample size  
-----------+--------------------------------------------------------------------
   Overall |     0.02      0.01      0.02      0.01      0.03             7500  
--------------------------------------------------------------------------------
NOTE: % centiles obtained from 7500 simulations of the posterior distribution

****************************************************************************************
            Conditional summary: Proportion Ratio
****************************************************************************************

-------------------------------------------------------------------------------
 Parameter |       Mean         SD      MCSE     Median  eti Lower  eti Upper  
-----------+-------------------------------------------------------------------
   Overall |       0.24       0.17      0.01       0.19       0.04       0.68  
-------------------------------------------------------------------------------

****************************************************************************************

Population-averaged estimates: Proportion Ratio 

--------------------------------------------------------------------------------
 Parameter |     Mean        SE    Median     Lower     Upper      Sample size  
-----------+--------------------------------------------------------------------
   Overall |     0.26      0.17      0.22      0.05      0.71             7500  
--------------------------------------------------------------------------------
NOTE: % centiles obtained from 7500 simulations of the posterior distribution

        | The bayesian estimation commands saved a dataset C:\DATA\WIV\Projects\Stata\Metapreg\mcmcresults\metapreg_bayesest.dta
        | C:\DATA\WIV\Projects\Stata\Metapreg\mcmcresults\metapreg_bayesest.dta
        | containing the MCMC samples of the parameters to the disk.
        | It is your responsibility to erase the dataset
        | after it is no longer needed.
        Click to erase the dataset
```
The Bayesian estimate of the treatment effect between-study variance is 0.15. 

![Figure 2.2.3b](/Markdown/223b.png)

The population-averaged summary RR was 0.22 [0.05, 0.71]. The credible interval is wider than the frequentist confidence interval 0.20 (0.06, 0.63). 

### 2.2.4 Interaction between two categorical covariates
This dataset, based on a Cochrane review, was used to evaluate the impact of missing data on clinical outcomes. The outcome of interest is clinical improvement and risk ratios (RRs) larger than 1 favour haloperidol over placebo in treatment of schizophrenia. Chaimani et al. (2014) fitted the fixed effect as well as random effects models for log RR using **metan**. They also conducted a subgroup analysis where studies were classified according to thr presenence or absence of missing outcome data in bith arms. 

Chaimani et al. (2014) states
>A common misconception about subgroup analyses is that results differ between subgroups when the summary effect for one subgroup is statistically significant and not for the other. However, inference on subgroup differences should be based on an interaction test (i.e. the test for subgroup differences ) that compares statistically the two subgroup means accounting for their uncertainty. Differences between subgroups can be also identified visually by looking at the overlap of the CIs in their summary estimates.

First, load the dataset and sort it so that for each study, the first row contains data from the placebo arm. For graphical aesthetics, rename the **response** variable to a shorter name **n**:

```
 use "http://fmwww.bc.edu/repec/bocode/s/schizo.dta", clear

 gen n = response
 
 gen studyid = firstauthor + " " + string(year)
 
 gsort studyid -arm
```
 
regularization is explicitly introduced via penalization methods like ridge, Lasso and penalized ML logistic regression 

### 2.2.4 Interaction between categorical and continous covariates
The data used in this example were presented in table IV of Berkey et al. (1995). Sharp (1998) used the data to demonstrate meta-analysis of odds-ratios with the **meta** command and investigaged the effect of latitude on BCG vaccination with the **metareg** command. The **meta** and **metareg** commands work on a dataset containing an estimate of the effect and its variability or confidence interval, for each study.

Their analysis suggested that BCG vaccination was more effective at higher absolute latitude. First, load the dataset:
```
use "http://fmwww.bc.edu/repec/bocode/b/bcg.dta", clear
```
Fit a logistic regression model with *bcg*, a categorical variable for the arm and *lat*, a continous variable with absolute latitude.
Activated by the option **interaction**, introducing an interaction term in the model allows to assess whether the log OR for arm vary by absolute latitude.
```
metapreg cases_tb population bcg lat, ///
	studyid(study) model(mixed, intmethod(mv)) ///
	sortby(lat) ///
	design(comparative, cov(commonslope)) ///
	outplot(rr) sumtable(rd rr) noitable ///
	interaction ///
	xlab(0, 1, 2) ///
	xtick(0, 1, 2) ///
	rcols(cases_tb population) ///
	astext(80) ///
	texts(1.5) logscale smooth gof
```
This produces the following Stata results:
```

    cases_tb ~ binomial(p, population)
    logit(p) = mu + bcg + lat + lat*bcg + study
    study ~ N(0, tau2)

    Base levels

        Variable -- Base Level
        bcg -- No


    Number of observations = 26
    Number of studies = 13


Goodness of Fit Criterion

        |      AIC       BIC  
--------+---------------------
  Value |   247.50    253.79  


Click to show the raw estimates

****************************************************************************************

Test of heterogeneity - LR Test: RE model vs FE model

-----------------------------------------------------------------------------------------------
           |       DF       Chisq           p        tau2      sigma2       I2tau     I2sigma  
-----------+-----------------------------------------------------------------------------------
     Model |        1     2504.11        0.00        1.26        0.00                          
-----------------------------------------------------------------------------------------------
NOTE: H0: Between-study variance(s) = 0  vs. H1: Between-study variance(s) > 0

****************************************************************************************
            Conditional summary: Proportion Difference
****************************************************************************************

-------------------------------------------------------------------------------
 Parameter |       Mean         SE         t      P>|t|    t Lower    t Upper  
-----------+-------------------------------------------------------------------
   Overall |       0.02       0.01      2.42       0.02       0.00       0.03  
-------------------------------------------------------------------------------

****************************************************************************************

Population-averaged estimates: Proportion Difference 

--------------------------------------------------------------------------------
 Parameter |     Mean        SE    Median     Lower     Upper      Sample size  
-----------+--------------------------------------------------------------------
   Overall |     0.03      0.02      0.03      0.01      0.07              800  
--------------------------------------------------------------------------------
NOTE: % centiles obtained from 800 simulations of the posterior distribution

****************************************************************************************
            Conditional summary: Proportion Ratio
****************************************************************************************

-------------------------------------------------------------------------------
 Parameter |       Mean    SE(lrr)    t(lrr)      P>|t|    t Lower    t Upper  
-----------+-------------------------------------------------------------------
   Overall |       0.49       0.05    -14.90       0.00       0.45       0.54  
-------------------------------------------------------------------------------
NOTE: H0: Est = 1 vs. H1: Est != 1

****************************************************************************************

Population-averaged estimates: Proportion Ratio 

--------------------------------------------------------------------------------
 Parameter |     Mean        SE    Median     Lower     Upper      Sample size  
-----------+--------------------------------------------------------------------
   Overall |     0.55      0.03      0.56      0.51      0.62              800  
--------------------------------------------------------------------------------
NOTE: % centiles obtained from 800 simulations of the posterior distribution
```
![Figure 2.2.6](/Markdown/226.png)

Replay the estimation results:
```
estimates replay metapreg_modest
```
This producest the following results:
```
Mixed-effects logistic regression               Number of obs     =         26
Binomial variable: population
Group variable: study                           Number of groups  =         13

                                                Obs per group:
                                                              min =          2
                                                              avg =        2.0
                                                              max =          2

Integration points =   7                        Wald chi2(4)      =     453.80
Log likelihood = -118.75083                     Prob > chi2       =     0.0000

------------------------------------------------------------------------------
    cases_tb | Coefficient  Std. err.      z    P>|z|     [95% conf. interval]
-------------+----------------------------------------------------------------
          mu |  -6.601802   .8197557    -8.05   0.000    -8.208493    -4.99511
             |
         bcg |
        Yes  |   .3995109   .0822356     4.86   0.000      .238332    .5606897
         lat |   .0735971   .0226449     3.25   0.001     .0292138    .1179804
             |
   bcg#c.lat |
        Yes  |  -.0333459   .0027862   -11.97   0.000    -.0388067    -.027885
------------------------------------------------------------------------------

------------------------------------------------------------------------------
  Random-effects parameters  |   Estimate   Std. err.     [95% conf. interval]
-----------------------------+------------------------------------------------
study: Identity              |
                     var(mu) |   1.261433   .5078444      .5730283    2.776848
------------------------------------------------------------------------------
LR test vs. logistic model: chibar2(01) = 2504.11     Prob >= chibar2 = 0.0000
```
The interaction term from **metapreg** and the coefficient for lat using **metareg** as was done by Sharp (1998) match closely, i.e. -0.0333459(-0.0388067, -0.027885) vs. -0.0315724(-0.0436704, -0.0194744).


