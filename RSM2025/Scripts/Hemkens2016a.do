/*
References
Colchicine for prevention of cardiovascular events, 2016, https://doi.org/10.1002/14651858.CD011047.pub2 analysis 1.8
Comparison of random-effects meta-analysis models for the relative risk in the case of rare events: A simulation study, 2020, 10.1002/bimj.201900379

*/

{
//=================================Install packages
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
global review "hemkens2016a"
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
*===========================Real Long Data===========================
{
/*
frame create longdata
frame change longdata

import excel using "$rootdir/Data/metadata.xlsx",sheet("hemkens2016-analysis1-8") ///
      firstrow clear
	  
gsort study -treatment

rename events event

metapreg event total treatment,  model(mixed) smooth gof  ///
	inference(bayesian) bwd($wdir)  ///
	studyid(study) design(comparative, cov(commonint)) sumtable(all) sumstat(Proportion)	///
	xlab(0, 0.05, 0.20)  ///
	lcols(event total) texts(2.5)  astext(75)  outplot(abs)   ///
	 graphregion(color(gray))  fxsize(90)	nooverall subline


metapreg event total treatment,  model(mixed) smooth gof  catpplot nofplot  ///
	/*inference(bayesian) bwd($wdir)*/  ///
	studyid(study) design(comparative, cov(commonint))  sumstat(Risk Difference)	///
	xlab(-0.15, 0, 0.15)    ///
	texts(2.5)  astext(75)  outplot(rd)   ///
	xline(0) graphregion(color(gray))  fxsize(90)	
	
metapreg event total treatment,  model(mixed) smooth gof  catpplot nofplot  ///
	/*inference(bayesian) bwd($wdir)*/  ///
	studyid(study) design(comparative, cov(commonint)) sumstat(Risk Ratio)	///
	xlab(0, 1, 10) logscale  ///
	texts(2.5)  astext(75)  outplot(rr)  xline(1) ///
		graphregion(color(gray)) fxsize(80)
		
metapreg event total treatment,  model(mixed) smooth gof  catpplot nofplot  ///
	/*(bayesian) bwd($wdir)*/  ///
	studyid(study) design(comparative, cov(commonint)) sumstat(Odds Ratio)	///
	xlab(0, 1, 10) logscale  ///
	texts(2.5)  astext(75)  outplot(or)  xline(1) ///
		graphregion(color(gray)) fxsize(80)	
		
gr combine absfplot rrcatpplot rdcatpplot  orcatpplot , ///
		graphregion(color(gray)  margin(tiny))  cols(2) imargin(0 0 0 0) ///
		 name(stats, replace)	ysize(10) xsize(17.5)	 
		 
		 
graph export "$graphs\data-F-BN-CI.png", as(png) width(2000) height(1000) replace	
		

*To compare with metafor
metapreg event total treatment if !inlist(study, "Yurdakul 2001", "Parise 1995"),  model(mixed) smooth gof  catpplot nofplot  ///
	studyid(study) design(comparative, cov(independent)) sumstat(Odds Ratio)	///
	xlab(0, 1, 10) logscale  ///
	texts(2.5)  astext(75)  outplot(or)  xline(1) ///
		graphregion(color(gray)) fxsize(80)				
	 
*/					
}
*===========================Fit all models===========================
{
/*
	frame create widedata
	frame change widedata

	import excel using "$rootdir/Data/metadata.xlsx",sheet("hemkens2016-analysis1-8") ///
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
*========================Process results from real data
{
//Read the observed results
cap frame drop observedresults
frame create observedresults
frame change observedresults
import delimited "$observedresults", varnames(1) clear

do "$scriptdir\process.do"

gen ciwidth = esthatup - esthatlo

sort package inference AIC BIC DIC  Dist Design Covariance Link Slope Weighting Sigmethod CImethod loglest sumethod param Prior

format AIC BIC DIC logBF  sigma2hat tau2hat ciwidth  %10.2f

tostring ciwidth, gen(WidthCI) format(%4.2f) force

do "$scriptdir\blobbogram.ado"
//Generate plots

//BN-CI
mat optimalor = (.1907646, .0363202, .67332775,  .0052328, .54379046)
global text2add "metapreg:tau2=0, sigma2=0.15"
global text2addright "0.19 eti:(0.04, 0.67)  hpd:(0.01, 0.54)"

global simorrange "0.00, 0.15, 2"
global allrrange "0.01, 0.15, 1.5"

global rrrange "0, 1, 5"
global rrline "1"
mat optimalrr = (0.20, 0.06, 0.69)


global orrange "0.01, 0.10, 2"
global orline "1"


global rdrange "-.1, 0, 0.1"
global rdline "0"
mat optimalrd = (0.0205, 0.0058, 0.0351)


duplicates drop stat Weighting CImethod if package == "meta" & slope=="common", force

//Relative ratio
gen select = 1 if package =="metapreg" & link == "logit" & ///
	(Covariance == "CI" & Design == "comparative") & Dist == "BN"  & ///
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
do "$scriptdir\Modelselection.do"    //Step 1: Go to modelselection.do

frame copy observedresults cleanresults, replace


do "$scriptdir\plotobsdata.do"   //Go to plotobsdata

frame copy observedresults cleanresults, replace
do "$scriptdir\widen.do"

//Step 2: Pick best model

gen conditional = 0
replace conditional = 1 if package == "metapreg" & link == "logit" & /// 
						(Covariance == "CI" & Design == "comparative") & Dist == "BN"  
*replace conditional = 1 if package == "metapreg" & link == "logit" & Covariance == "FI"  & Dist == "BN"  & inference == "bayesian"

//Step 3: Add BB
gen pick = conditional
replace pick = 1 if strpos(Dist, "BB") != 0  & link == "logit"
drop if strpos(stat, "pop") != 0 & strpos(model, "betabin") != 0

//Step 3.1: binomial
replace pick = 1 if Dist=="B2" & link == "logit" 

//Step 5: Pick similar/alternative conditional models
replace conditional = 1 if  package == "metastan" & inlist(Dist, "BN")  
replace conditional = 1 if package== "bayesmeta" & inlist(Dist,  "NN")  
replace conditional = 1 if inlist(Sigmethod, "mp", "sj") & Weighting == "SSW" & package == "meta" & Env == "R"
replace conditional = 0 if slope == "common"
replace conditional = 1 if (Dist == "QN" | Sigmethod == "mp" | Sigmethod == "pl") & package == "metan" 
replace conditional = 1 if Sigmethod == "sj" & package == "meta" & Env == "Stata"
replace conditional = 1 if  package == "metaplus"		 
replace conditional = 1 if package == "metafor" & Dist == "BN" & (Covariance == "IND" )

//Step 6: Plot all conditionals --> AllCond.do
do "$scriptdir\AllCond.do"	

}

*===========================Replicate Best Fit ===========================
{
/*
	//Parameters from the second `best' bayesian fit - BN-CI
	mat truemeans = (-3.661096,  -1.656715)  //medians
	mat truevarcov = ( 0 ,0 \ 0,  .1507779)  //medians
	
	import excel using "$rootdir/Data/metadata.xlsx",sheet("hemkens2016-analysis1-8") ///
      firstrow clear

	gsort study treatment

	rename events event
	rename study studyid

	//To wide format
	gen code = 0 if treatment == "Control"
	replace code = 1 if treatment != "Control"

	reshape wide event total treatment, i(study ) j(code)

	keep studyid event* total*
	 
	rename event0 obsevent0
	rename event1 obsevent1

	//Replicate and Fit all the models
	global outfile = "$simresults"
	do "$scriptdir/replicatefit.do"
	
	global brun 1
	global frun 1
	global run 1
	global test "no"
	
	forvalues r = 1(1)100 {
		replicatefit, truevarcov(truevarcov) truemeans(truemeans) seed(`r') 		
	}
	
*/
}
*===========================Process Simulated results===========================
cap frame drop simresults
frame create simresults
frame change simresults

//Read the simulated results
import delimited "$simresults", varnames(1) clear

//Process the variables
do "$scriptdir\process.do"
save "$wdir/simulatedresults.dta", replace

 *=========Run R script
rsource using "$scriptdir/hemkens2016a.R"  , noloutput

 

