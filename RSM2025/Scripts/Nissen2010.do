/*
Rare events meta-analysis using the Bayesian beta-binomial mode, 2023, 10.1002/jrsm.1662
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
global review "nissen2010"
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
*/
}
*===========================ORIGINAL Wide DATA===========================
{
/*

frame create widedata
frame change widedata

import excel using "$rootdir/Data/metadata.xlsx",sheet("nissen2010") ///
      firstrow clear

//To long format
reshape long event total, i(trial) j(class)

gen treatment = "Control" if class==0
replace treatment = "Treatment" if class==1	

rename trial study

set more off

metapreg event total treatment,  model(mixed) smooth gof progress  ///
	/*inference(bayesian)  bwd($wdir)*/ ///
	studyid(study) design(comparative, cov(commonslope)) sumtable(none) 	///
	xlab(0, 0.05, 0.5)  ///
	texts(0.75)  astext(75)  outplot(abs)   ///
	nooverall subline fxsize(100) fysize(200) 
	
graph export "$graphs\data-P-F.png", as(png) width(2000) height(1750) replace	

metapreg event total treatment,  model(mixed) smooth gof  catpplot nofplot  ///
	inference(bayesian)  bwd($wdir) ///
	studyid(study) design(comparative, cov(commonslope)) 	///
	xlab(-.15, 0, .15)  ///
	lcols(event total) texts(1.25)  astext(70)  outplot(rd) sumstat(Risk Difference) ///
	xline(0)   ysize(17.5) xsize(20)
	
graph export "$graphs\data-RD-F.png", as(png) width(2000) height(1750) replace	


metapreg event total treatment,  model(mixed ) smooth gof  catpplot nofplot  ///
	inference(bayesian)  bwd($wdir) ///
	studyid(study) design(comparative, cov(commonslope)) sumtable(all)	///
	xlab(0.01, 1,  30) logscale  ///
	lcols(event total) texts(1.25)   astext(70)  outplot(or)  xline(1)  ///
	 ysize(17.5) xsize(20)

graph export "$graphs\data-OR-F.png", as(png) width(2000) height(1750) replace	
	 
metapreg event total treatment,  model(mixed) smooth gof  catpplot nofplot  ///
	inference(bayesian)  bwd($wdir)  ///
	studyid(study) design(comparative, cov(commonslope)) sumtable(all)	///
	xlab(0.01, 1,  30) logscale sumstat(Risk Ratio) ///
	lcols(event total) texts(1.25)  astext(70)  outplot(rr)  xline(1)  ///
	ysize(17.5) xsize(20)
	 
graph export "$graphs\data-RR-F.png", as(png) width(2000) height(1750) replace	
	
*/
}
*===========================Fit all models===========================
{
/*
frame create widedata
frame change widedata

import excel using "$rootdir/Data/metadata.xlsx", sheet("nissen2010") ///
      firstrow clear
 
rename trial studyid
	 
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

sort package inference AIC BIC DIC Dist Design Covariance Link Slope Weighting Sigmethod CImethod loglest sumethod param Prior

format AIC BIC DIC logBF sigma2hat tau2hat ciwidth  %10.2f

tostring ciwidth, gen(WidthCI) format(%4.2f) force

do "$scriptdir\blobbogram.ado"

//Generate plots
mat optimalor = (1.234775, .97741935, 1.5654786, .96336166, 1.5480207)
mat optimalabs = (0.06, 0.04, 0.07)
global text2add "metapreg: tau2=0.90, sigma2=0.00"
global text2addright "1.23 eti:(0.98, 1.57)  hpd:(0.96, 1.54)"

global obsrange "0, 0.5"

global obsrange "0, 0.1"

global simorrange "0.75, 1.25, 2"
global allrrange "0.15, 1, 3"
global orrange "0.75, 1.25, 1.75"
global orline "1"

global rrrange "0, 1, 5"
global rrline "1"
mat optimalrr = (1.49, 1.14, 1.95)

global rdrange "-.1, 0, 0.1"
global rdline "0"
mat optimalrd = (-0.028036, -0.05083, -0.00633)

duplicates drop stat Weighting CImethod if package == "meta" & slope=="common", force

//Relative ratio
gen select = 1 if package =="metapreg" & link == "logit" & ///
				(Covariance == "CS" & Design == "comparative") & Dist == "BN"  & ///
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
							(Covariance == "CS" & Design == "comparative") & Dist == "BN"  


//Step 3: Add BB
gen pick = conditional
replace pick = 1 if strpos(Dist, "BB") != 0 & link == "logit" /*& CImethod == "z"*/

//Step 4: Get stats from the best model and BB estimates ==> run StatsBestPick.do
replace conditional = 1 if  package == "metastan" & inlist(Dist, "BN")  
replace conditional = 1 if package== "bayesmeta" & inlist(Dist,  "NN")  
replace conditional = 1 if inlist(Sigmethod, "mp", "sj") & Weighting == "SSW" & package == "meta" & Env == "R"
replace conditional = 0 if slope == "common"
replace conditional = 1 if (Dist == "QN" | Sigmethod == "mp" | Sigmethod == "pl") & package == "metan" 
replace conditional = 1 if Sigmethod == "sj" & package == "meta" & Env == "Stata"
replace conditional = 1 if  package == "metaplus"
	 
replace conditional = 1 if package == "metafor" & Dist == "BN" & Slope == "C"

//Step 5: Plot all conditionals --> AllCond.do
do "$scriptdir\AllCond.do"	

}
*===========================Replicate Best Model===========================
{
/*

	//Parameters from the best bayesian fit; highest* logBF - cs; good mixing
	mat truemeans = (-5.820007 , .2108263)  //median
	mat truevarcov = (.9047816, 0 \ 0,  0) //median
	
	import excel using "$rootdir/Data/metadata.xlsx", sheet("nissen2010") ///
      firstrow clear

	rename trial studyid
	 
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
*drop package- v10
do "$scriptdir\process.do"
save "$wdir/simulatedresults.dta", replace

 *=========Run R script
rsource using "$scriptdir/nissen2010.R"  , noloutput

