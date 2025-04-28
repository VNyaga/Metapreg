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

//Specify directories
global review "hemkens2016b"
global rootdir "C:/DATA/WIV/Projects/GitHub/Metapreg/RSM2025"
global graphs "$rootdir/$review/Graphs"

global wdir "$rootdir/$review"
global scriptdir "$rootdir/Scripts"

global observedresults "$wdir/$review.txt"
global simresults "$wdir/simresults.txt"

//Make folders
mkdir "$rootdir/$review"
mkdir "$graphs"

//Put headers on results files 
foreach name in observedresults simresults  {
	local s = "sim; model; stat; esthat; esthatlo; esthatup; tau2hat; tau2hatlo; tau2hatup; sigma2hat; sigma2hatlo; sigma2hatup; IC; truemu0; trueor;  truetausq; truesigmasq; k; nstudies; studysize"
	file open results using $`name',  text write replace 
	file write results "`s'" _n
	file close results
}

//Specify where R is installed
global Rterm_path `"C:\DATA\Software\R-4.4.2\bin\x64\Rterm.exe"'
global Rterm_options `"--vanilla"'

//Run metapreg if not installed
do "$scriptdir/metapreg.ado"

}
*===========================ORIGINAL Long DATA===========================
{
/*

frame create longdata
frame change longdata

import excel using "$rootdir/Data/metadata.xlsx",sheet("hemkens2016-analysis1-10") ///
      firstrow clear
	  
gsort study -treatment

rename events event

set more off
set trace off
metapreg event total treatment,  model(mixed) smooth gof  ///
	inference(bayesian) bwd($wdir)  ///
	studyid(study) design(comparative, cov(independent)) sumtable(none) 	///
	xlab(0, 0.05, 0.10) ///
	lcols(event total) texts(2.35)  astext(80)  outplot(abs)   ///
	 graphregion(color(gray))  fxsize(100)	nooverall subline

metapreg event total treatment,  model(mixed) smooth gof  catpplot nofplot  ///
	inference(bayesian) bwd($wdir) ///
	studyid(study) design(comparative, cov(independent)) sumstat(Risk Difference)	///
	xlab(-0.20, 0, 0.20)  ///
	texts(2.35)  astext(75)  outplot(rd)   ///
	xline(0) graphregion(color(gray))  fxsize(90)	
	
metapreg event total treatment,  model(mixed) smooth gof  catpplot nofplot  ///
	inference(bayesian) bwd($wdir) ///
	studyid(study) design(comparative, cov(independent)) sumstat(Odds Ratio)	///
	xlab(0, 1, 30) logscale ///
	texts(2.35)  astext(80)  outplot(or rr)  xline(1) ///
		graphregion(color(gray)) fxsize(100)
	 
metapreg event total treatment,  model(mixed) smooth gof  catpplot nofplot  ///
	inference(bayesian) bwd($wdir) ///
	studyid(study) design(comparative, cov(independent)) sumstat(Risk Ratio)	///
	xlab(0, 1, 30) logscale ///
	texts(2.35)  astext(80)  outplot(rr)  xline(1) ///
		graphregion(color(gray)) fxsize(100)
		
gr combine absfplot rrcatpplot rdcatpplot  orcatpplot , ///
		graphregion(color(gray)  margin(tiny))  cols(2) imargin(0 0 0 0) ///
		 name(stats, replace)	ysize(10) xsize(17.5)	 
	 		 
		 
graph export "$graphs\data-B.png", as(png) width(2000) height(1000) replace	
				
*/
}
*===========================Fit all models===========================
{
/*
	frame create widedata
	frame change widedata

	import excel using "$rootdir/Data/metadata.xlsx",sheet("hemkens2016-analysis1-10") ///
		firstrow clear
		  
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

*/
}

*========================Process observed data
{
//Read the observed results
cap frame drop observedresults
frame create observedresults
frame change observedresults
import delimited "$observedresults", varnames(1) clear

do "$scriptdir\process.do"

gen ciwidth = esthatup - esthatlo

sort package inference AIC BIC DIC  Dist Design Covariance Link Slope Weighting Sigmethod CImethod loglest sumethod param Prior

format AIC BIC DIC logBF  sigma2hat tau2hat   %10.2f
tostring ciwidth, gen(WidthCI) format(%4.2f) force

format ciwidth  %10.2f
do "$scriptdir\blobbogram.ado"

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