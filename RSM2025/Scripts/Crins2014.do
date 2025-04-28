/*
Interleukin-2 receptor antagonists for pediatric liver transplant recipients: A systematic review and meta-analysis of controlled studies, 2014, https://doi.org/10.1111/petr.12362 - Figure 2
Meta-analysis of few small studies in orphan diseases, 2017, 10.1002/jrsm.1217 - Figure 1, plot on top.


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

global review "crins2014"
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

import excel using "$rootdir/Data/metadata.xlsx",sheet("crins2014") ///
      firstrow clear

gsort study group

rename group treatment
rename events event

metapreg event total treatment,  model(mixed) smooth gof  ///
	inference(bayesian) bwd($wdir) ///
	studyid(study) design(comparative, cov(independent)) sumtable(all) sumstat(Proportion)	///
	xlab(0, 0.5, 1)  ///
	lcols(event total) texts(2)  astext(75)  outplot(abs)   ///
	 graphregion(color(gray))  fxsize(90)	nooverall subline


metapreg event total treatment,  model(mixed) smooth gof  catpplot nofplot  ///
	inference(bayesian) bwd($wdir)  ///
	studyid(study) design(comparative, cov(independent))  sumstat(Risk Difference)	///
	xlab(-0.15, 0, 0.75)    ///
	texts(2)  astext(75)  outplot(rd)   ///
	xline(0) graphregion(color(gray))  fxsize(90)	
	
metapreg event total treatment,  model(mixed) smooth gof  catpplot nofplot  ///
	inference(bayesian) bwd($wdir)  ///
	studyid(study) design(comparative, cov(independent)) sumstat(Risk Ratio)	///
	xlab(0.01, 1, 2) logscale  ///
	texts(2)  astext(75)  outplot(rr)  xline(1) ///
		graphregion(color(gray)) fxsize(80)
		
metapreg event total treatment,  model(mixed) smooth gof  catpplot nofplot  ///
	inference(bayesian) bwd($wdir)  ///
	studyid(study) design(comparative, cov(independent)) sumstat(Odds Ratio)	///
	xlab(0.01, 1, 2) logscale  ///
	texts(2)  astext(75)  outplot(or)  xline(1) ///
		graphregion(color(gray)) fxsize(80)	

gr combine absfplot rrcatpplot rdcatpplot  orcatpplot , ///
		graphregion(color(gray)  margin(tiny))  cols(2) imargin(0 0 0 0) ///
		 name(stats, replace)	ysize(10) xsize(17.5)	 
		 
		 
graph export "$graphs\data-F-IND.png", as(png) width(2000) height(1000) replace	

*/
}
*===========================Fit all models===========================
{
/*
	cap frame drop widedata
	frame create widedata
	frame change widedata

	import excel using "$rootdir/Data/metadata.xlsx",sheet("crins2014") ///
      firstrow clear
	  
	gsort study group
	
	rename study studyid
	rename group treatment
	rename events event

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

qui do "$scriptdir\process.do"

gen ciwidth = esthatup - esthatlo

sort package inference AIC BIC DIC Dist Design Covariance Link Slope Weighting Sigmethod CImethod loglest sumethod param Prior

format AIC BIC DIC logBF sigma2hat tau2hat ciwidth  %10.2f

tostring ciwidth, gen(WidthCI) format(%4.2f) force

qui do "C:\DATA\WIV\Projects\Stata\Blobbogram\blobbogram.ado"
//Generate plots
mat optimalor = (.1915889, .0782588, .4200697, .0573697, .3689821)
global orbest ".20837106"

global poprrbest ".36775263"

mat optimalabs = (0.47227, 0.24604, 0.78707)
global text2add "metapreg:tau2=1.83, sigma2=0.21"
global text2addright "0.19 eti:(0.08, 0.42)  hpd:(0.06, 0.37)"
global obsrange "0.25, 0.75"

global simorrange "0.01, 0.15, 1"
global allrrange "0.02, 0.25, 1"

global orrange "0.05, 0.5, 1"
global orline "1"

global rrrange "0, 1"
global rrline "1"
mat optimalrr = (0.3249, 0.1483, 0.521147)

global rdrange "-.1, 0, 0.1"
global rdline "0"
mat optimalrd = (0.31, 0.17, 0.44)

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

//Step 1: Go to modelselection.do
do "$scriptdir\Modelselection.do"

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
drop if strpos(stat, "pop") != 0 & strpos(model, "betabin") != 0

//Step 4: binomial
replace pick = 1 if Dist=="B2" & link == "logit"  //exact

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
*===========================Replicate Best Model===========================
{
/*
	
//Parameters from the best bayesian fit; highest logBF - logit independent igamma prior
	mat truemeans = ( -.1258962, -1.649212)  //medians
	mat truevarcov = (1.824899, 0 \ 0,  .2331431)  //medians
	
	import delimited "C:\DATA\WIV\Projects\Stata\Metapreg\Data\crins2014.csv", varnames(1) clear

	gsort author group

	gen studyid = author + " " + string(year)
	drop year author

	rename group treatment
	rename events event

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
rsource using "$scriptdir/crins2014.R"  , noloutput
