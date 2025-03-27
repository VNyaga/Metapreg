/*

Refeference: Comparison of random-effects meta-analysis models for the relative risk in the case of rare events: A simulation study, 2020, 10.1002/bimj.201900379

Fix hexact then try it.
*/
*=======================HOUSE-KEEPING===============================
{
/*Install packages*/
/*
ssc install metapreg
help metapreg

ssc install metan
help metan

ssc install moremata

ssc install integrate 

ssc install rsource
help rsource
*/
//Summary 
{
	global summaryfile "C:/DATA/WIV/Projects/GitHub/Metapreg/Data/summary.txt"

	/*
	local s = "Dataset; Statistic; Metapreg; Max; Min"
	file open results using $summaryfile , text write replace //append 
	file write results "`s'" _n
	file close results
	*/
}
global review "hemkens2"

global rootdir "C:/DATA/WIV/Projects/GitHub/Metapreg"
*global rootdir "M:/1120_Cancer_Employee/BMH/Victoria/GitHub/Metapreg"

global graphs "C:/DATA/WIV/Projects/GitHub/Metapreg/Data/$review/Graphs"

mkdir "$rootdir/Data/$review"
*mkdir "C:/DATA/WIV/Projects/GitHub/Metapreg/Data/$review"
mkdir "$graphs"
global wdir "$rootdir/Data/$review"
*global wdir "C:/DATA/WIV/Projects/GitHub/Metapreg/Data/$review"

global observedresults "$wdir/$review.txt"
global simresults "$wdir/simresults.txt"

//Heads for the results file 
foreach name in observedresults /*simresults*/ {
	local s = "sim; model; stat; esthat; esthatlo; esthatup; tau2hat; tau2hatlo; tau2hatup; sigma2hat; sigma2hatlo; sigma2hatup; IC; truemu0; trueor;  truetausq; truesigmasq; k; nstudies; studysize"
	file open results using $`name',  text write replace //append or replace
	file write results "`s'" _n
	file close results
}

global scriptdir "$rootdir/Scripts"
*global Rterm_path `"C:\Program Files\R\R-4.4.0\bin\x64\Rterm.exe"'
global Rterm_path `"C:\DATA\Software\R-4.4.2\bin\x64\Rterm.exe"'

global Rterm_options `"--vanilla"'

do "C:\DATA\WIV\Projects\GitHub\Metapreg\Build\metapreg.ado"

}
*===========================ORIGINAL Long DATA===========================
{

frame create longdata
frame change longdata

import delimited "$rootdir/hemkens-cochrane-2016-analysis1-10.csv", varnames(1) clear

gsort study -treatment

rename events event

gen group = treatment
gsort study -group

set more off
set trace off
metapreg event total group,  model(mixed) smooth gof  ///
	/*inference(bayesian) bwd(C:\DATA\WIV\Projects\Stata\Metapreg\mcmcresults)*/  ///
	studyid(study) design(comparative, cov(independent)) sumtable(none) 	///
	xlab(0, 0.05, 0.10) ///
	lcols(event total) texts(2.35)  astext(80)  outplot(abs)   ///
	 graphregion(color(gray))  fxsize(100)	nooverall subline

*Diagnostics 
/*
estimates restore metapreg_modest

bayesgraph diagnostics {event:mu}

bayesgraph diagnostics {tausq}

bayesgraph diagnostics {sigmasq}
*/

metapreg event total group,  model(mixed) smooth gof  catpplot nofplot  ///
	/*inference(bayesian) bwd(C:\DATA\WIV\Projects\Stata\Metapreg\mcmcresults)*/  ///
	studyid(study) design(comparative, cov(independent)) sumstat(Risk Difference)	///
	xlab(-0.20, 0, 0.20)  ///
	texts(2.35)  astext(75)  outplot(rd)   ///
	xline(0) graphregion(color(gray))  fxsize(90)	
	
metapreg event total group,  model(mixed) smooth gof  catpplot nofplot  ///
	/*inference(bayesian) bwd(C:\DATA\WIV\Projects\Stata\Metapreg\mcmcresults)*/  ///
	studyid(study) design(comparative, cov(independent)) sumstat(Odds Ratio)	///
	xlab(0, 1, 30) logscale ///
	texts(2.35)  astext(80)  outplot(or rr)  xline(1) ///
		graphregion(color(gray)) fxsize(100)
	 
metapreg event total group,  model(mixed) smooth gof  catpplot nofplot  ///
	/*inference(bayesian) bwd(C:\DATA\WIV\Projects\Stata\Metapreg\mcmcresults)*/  ///
	studyid(study) design(comparative, cov(independent)) sumstat(Risk Ratio)	///
	xlab(0, 1, 30) logscale ///
	texts(2.35)  astext(80)  outplot(rr)  xline(1) ///
		graphregion(color(gray)) fxsize(100)
		
gr combine absfplot rrcatpplot rdcatpplot  orcatpplot , ///
		graphregion(color(gray)  margin(tiny))  cols(2) imargin(0 0 0 0) ///
		 name(stats, replace)	ysize(10) xsize(17.5)	 
	 		 
		 
graph export "$graphs\data-F.png", as(png) width(2000) height(1000) replace	


	
estimates restore metapreg_modest

bayesgraph diagnostics {event:2.group}

bayesgraph diagnostics {event:mu}


bayesgraph diagnostics {tausq}

bayesgraph diagnostics {sigmasq}
	
//----------------------HEXACT
{
	set more off
set trace off
metapreg event total group,  model(hexact)  ///
	studyid(study) design(comparative) 	///
	xlab(0, 1, 30) logscale  ///
	lcols(event total) texts(2)  astext(70) xline(1) ///
	 ysize(5) xsize(7.5)	


metapreg event total group,  model(crbetabin) smooth gof  ///
	/*inference(bayesian) bwd(C:\DATA\WIV\Projects\Stata\Metapreg\mcmcresults) */ ///
	studyid(study) design(comparative, cov(independent)) sumtable(none) 	///
	xlab(0, 0.1) ///
	lcols(event total) texts(2.60)  astext(80)  outplot(abs)   ///
	 graphregion(color(gray))  fxsize(110)	

metapreg event total group,  model(crbetabin) smooth gof  catpplot nofplot  ///
	/*inference(bayesian) bwd(C:\DATA\WIV\Projects\Stata\Metapreg\mcmcresults) */ ///
	studyid(study) design(comparative, cov(independent)) nowt sumstat(Proportion Difference)	///
	xlab(-0.25, 0, 0.25)  ///
	texts(2.60)  astext(70)  outplot(rd)   ///
	xline(0) graphregion(color(gray))  fxsize(110)	
	
metapreg event total group,  model(crbetabin) smooth gof  catpplot nofplot  ///
	/*inference(bayesian) bwd(C:\DATA\WIV\Projects\Stata\Metapreg\mcmcresults) */ ///
	studyid(study) design(comparative, cov(independent)) nowt popstat(median)	///
	xlab(0, 1, 30) logscale ///
	texts(2.60)  astext(80)  outplot(rr or)  xline(1) ///
		graphregion(color(gray)) fxsize(90)
	 
gr combine absfplot rrcatpplot rdcatpplot  orcatpplot , ///
		graphregion(color(gray)  margin(tiny))  cols(2) imargin(0 0 0 0) ///
		 name(stats, replace)	ysize(10) xsize(20)	 	 
}	

	
				
}
*===========================Fit all models===========================
{
frame create longdata
frame change longdata

import delimited "$rootdir/hemkens-cochrane-2016-analysis1-10.csv", varnames(1) clear

gsort study treatment

rename events event
rename study studyid

drop riskgroup


//To wide format
gen code = 0 if treatment == "Control"
replace code = 1 if treatment != "Control"

reshape wide event total treatment, i(study ) j(code)

keep studyid event* total*
	 
//Fit all the models
global outfile = "$observedresults"
do "$scriptdir\mainfit.do"

}

*========================Process observed data
{
//Read the observed results
frame create observedresults
frame change observedresults
import delimited "$observedresults", varnames(1) clear

do "$scriptdir\process.do"

gen ciwidth = esthatup - esthatlo


sort package inference AIC BIC DIC  Dist Design Covariance Link Slope Weighting Sigmethod CImethod loglest sumethod param Prior

format AIC BIC DIC logBF  sigma2hat tau2hat   %10.2f
tostring ciwidth, gen(WidthCI) format(%4.2f) force

format ciwidth  %10.2f
do "C:\DATA\WIV\Projects\Stata\Blobbogram\blobbogram.ado"
//Generate plots
mat optimalor = (1.2846097, .00994537, 61.374548, .0001214, 24.022034)
mat optimalabs = (0.00, 0.00, 0.01)
global text2add "metapreg:tau2=0.81, sigma2=1.67"
global text2addright "1.28 eti:(0.01, 61.37)  hpd:(0.00, 24.02)"
global obsrange "0, 0.1"

global orrange "0.0001, 1, 1000"
global orline "1"

global simorrange "0.01, 1, 50"
global allrrange "0.01, 0.5, 15"

global rrrange "0, 1, 5"
global rrline "1"
mat optimalrr = (1.39, 0.5, 2)

global rdrange "-0.02, 0, 0.02"
global rdline "0"
mat optimalrd = (-0.0022, -0.0066, 0.0021)


duplicates drop stat Weighting CImethod if package == "meta" & slope=="common", force

//Relative ratio
gen select = 1 if package =="metapreg" & link == "logit" & ///
	(Covariance == "IND" & Design == "comparative") & Dist == "BN"  & ///
	((stat == "median-orout" & Inf == "B")|(stat == "mean-orout" & Inf == "F") )
	
sort select inference

gen bselect = esthat[1]  //bayesian
gen fselect = esthat[3]  //frequentist

gen bratio = esthat/bselect 
tostring bratio, replace format(%4.2f) force

gen fratio = esthat/fselect 
tostring fratio, replace format(%4.2f) force

frame copy observedresults cleanresults, replace
frame change cleanresults
do "$scriptdir\Modelselection.do"
//Step 1: Go to modelselection.do

frame copy observedresults cleanresults, replace

//Go to plotobsdata
do "$scriptdir\plotobsdata.do"

frame copy observedresults cleanresults, replace
do "$scriptdir\widen.do"

//Step 2: Pick best model
gen conditional = 0
replace conditional = 1 if package == "metapreg" & link == "logit" & ///
						(Covariance == "IND" & Design == "comparative") & Dist == "BN"  

//Step 3: Add BB
gen pick = conditional
replace pick = 1 if strpos(Dist, "BB") != 0 & link == "logit" 

//Step 3.1: binomial
replace pick = 1 if inlist(Dist, "B2") & link == "logit" 

//Step 4: Get stats from the best model and BB estimates ==> run StatsBestPick.do


//Step 5: Pick similar/alternative conditional models
replace conditional = 1 if  package == "metastan" & inlist(Dist, "BN")  
replace conditional = 1 if package== "bayesmeta" & inlist(Dist,  "NN")  
drop if package == "metabma" | package == "metasem" | package == "randmeta"
replace conditional = 1 if inlist(Sigmethod, "mp", "sj") & Weighting == "SSW" & package == "meta" & Env == "R"
replace conditional = 0 if slope == "common"
replace conditional = 1 if (Dist == "QN" | Sigmethod == "mp" | Sigmethod == "pl") & package == "metan" 
replace conditional = 1 if Sigmethod == "sj" & package == "meta" & Env == "Stata"
replace conditional = 1 if  package == "metaplus"		 
replace conditional = 1 if package == "metafor" & Dist == "BN" & (Covariance == "IND" )

//Step 6: Plot all conditionals --> AllCond.do
do "$scriptdir\AllCond.do"	

}