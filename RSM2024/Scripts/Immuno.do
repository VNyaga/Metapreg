/*
Meta-analysis of few small studies in orphan diseases, 2017, 10.1002/jrsm.1217


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
global review "immuno"

global rootdir "C:/DATA/WIV/Projects/GitHub/Metapreg"
*global rootdir "M:/1120_Cancer_Employee/BMH/Victoria/GitHub/Metapreg"

global graphs "C:/DATA/WIV/Projects/GitHub/Metapreg/Data/$review/Graphs"

mkdir "C:/DATA/WIV/Projects/GitHub/Metapreg/Data/$review"
mkdir "$graphs"

global wdir "C:/DATA/WIV/Projects/GitHub/Metapreg/Data/$review"

global observedresults "$wdir/$review.txt"
global simresults "$wdir/simresults.txt"

//Heads for the results file 
foreach name in observedresults /*simresults*/ {
	local s = "sim; model; stat; esthat; esthatlo; esthatup; tau2hat; tau2hatlo; tau2hatup; sigma2hat; sigma2hatlo; sigma2hatup; IC; truemu0; trueor;  truetausq; truesigmasq; k; nstudies; studysize"
	file open results using $`name',  text write replace //append or replace
	file write results "`s'" _n
	file close results
}

global scriptdir "C:/DATA/WIV/Projects/GitHub/Metapreg/Scripts"
global Rterm_path `"C:\DATA\Software\R-4.4.2\bin\x64\Rterm.exe"'
global Rterm_options `"--vanilla"'

do "C:\DATA\WIV\Projects\GitHub\Metapreg\Build\metapreg.ado"
	
}
*===========================ORIGINAL Long DATA===========================
{
frame create longdata
frame change longdata
import delimited "C:\DATA\WIV\Projects\Stata\Metapreg\Data\immuno.csv", varnames(1) clear

gsort author group

rename author study
rename group treatment
rename events event

gen T = "T" if treatment ! = "Control"
replace T = "C" if treatment == "Control"
	

gen group = treatment
gsort study group

metapreg event total treatment,  model(mixed) smooth gof  ///
	/*inference(bayesian) bwd(C:\DATA\WIV\Projects\Stata\Metapreg\mcmcresults)*/ ///
	studyid(study) design(comparative, cov(independent)) sumtable(all) sumstat(Proportion)	///
	xlab(0, 0.5, 1)  ///
	lcols(event total) texts(2)  astext(75)  outplot(abs)   ///
	 graphregion(color(gray))  fxsize(90)	nooverall subline


metapreg event total treatment,  model(mixed) smooth gof  catpplot nofplot  ///
	/*inference(bayesian) bwd(C:\DATA\WIV\Projects\Stata\Metapreg\mcmcresults)*/  ///
	studyid(study) design(comparative, cov(independent))  sumstat(Risk Difference)	///
	xlab(-0.15, 0, 0.75)    ///
	texts(2)  astext(75)  outplot(rd)   ///
	xline(0) graphregion(color(gray))  fxsize(90)	
	
metapreg event total treatment,  model(mixed) smooth gof  catpplot nofplot  ///
	/*inference(bayesian) bwd(C:\DATA\WIV\Projects\Stata\Metapreg\mcmcresults)*/  ///
	studyid(study) design(comparative, cov(independent)) sumstat(Risk Ratio)	///
	xlab(0.01, 1, 2) logscale  ///
	texts(2)  astext(75)  outplot(rr)  xline(1) ///
		graphregion(color(gray)) fxsize(80)
		
metapreg event total treatment,  model(mixed) smooth gof  catpplot nofplot  ///
	/*inference(bayesian) bwd(C:\DATA\WIV\Projects\Stata\Metapreg\mcmcresults) */ ///
	studyid(study) design(comparative, cov(independent)) sumstat(Odds Ratio)	///
	xlab(0.01, 1, 2) logscale  ///
	texts(2)  astext(75)  outplot(or)  xline(1) ///
		graphregion(color(gray)) fxsize(80)	

gr combine absfplot rrcatpplot rdcatpplot  orcatpplot , ///
		graphregion(color(gray)  margin(tiny))  cols(2) imargin(0 0 0 0) ///
		 name(stats, replace)	ysize(10) xsize(17.5)	 
		 
		 
graph export "$graphs\data-F-IND.png", as(png) width(2000) height(1000) replace	

		 
				
}
*===========================Fit all models===========================
{
	cap frame drop widedata
	frame create widedata
	frame change widedata

	import delimited "C:\DATA\WIV\Projects\Stata\Metapreg\Data\immuno.csv", varnames(1) clear

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

	//Fit all the models
	global outfile = "$observedresults"
	
	do "$scriptdir\mainfit.do"
	
	/*
	drop if  strpos(model, "bayesian") != 0 & strpos(model, "metapreg") != 0 
	export delimited using "$observedresults",   delimit(";") replace

	*/
}

*========================Process observed data
{
//Read the observed results
frame create observedresults
frame change observedresults
import delimited "$observedresults", varnames(1) clear

qui do "$scriptdir\process.do"

gen ciwidth = esthatup - esthatlo

sort package inference AIC BIC DIC Dist Design Covariance Link Slope Weighting Sigmethod CImethod loglest sumethod param Prior

format AIC BIC DIC logBF sigma2hat tau2hat ciwidth  %10.2f

tostring ciwidth, gen(WidthCI) format(%4.2f) force

/*
replace esthat = ln(esthat) if stat=="mean-orout"
replace esthatlo = ln(esthatlo) if stat=="mean-orout"
replace esthatup = ln(esthatup) if stat=="mean-orout"
*/
qui do "C:\DATA\WIV\Projects\Stata\Blobbogram\blobbogram.ado"
//Generate plots
mat optimalor = (.1915889, .0782588, .4200697, .0573697, .3689821)
global orbest ".20837106"

global poprrbest ".36775263"

*mat optimalor = (0.2156, 0.05056, 0.3497)
*mat optimalor = (-1.64, -2.31, -0.97)  //log scale
mat optimalabs = (0.47227, 0.24604, 0.78707)
global text2add "metapreg:tau2=1.83, sigma2=0.21"
global text2addright "0.19 eti:(0.08, 0.42)  hpd:(0.06, 0.37)"
global obsrange "0.25, 0.75"

*global orrange "-3, 0, 1"
*global orline "0"

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

//Step 3.1: binomial
replace pick = 1 if Dist=="B2" & link == "logit"  //exact
*replace pick = 1 if Dist=="B1" & Covariance=="CI" & link == "logit" 

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
*===========================Replicate Best Model===========================
{
	/*Parameters from the best frequentist fit - zeroslop loglog*/
	mat truemeans = (-.3447732 , .9553016)  //0.49
	mat truevarcov = (.4262142, 0 \  0,  0) 
	
	/*Parameters from the best bayesian fit; small DIC - logit unstructured*/
	mat truemeans = (-.1154955 , -1.698386)  //0.41
	mat truevarcov = (1.220993, .4402483 \ .4402483,  .813241) 
	
	/*Parameters from the best bayesian fit; highest logBF - logit independent jeffreys prior*/
	mat truemeans = ( -.1067554 , -1.654238)  //0.41
	mat truevarcov = (1.997152, 0 \ 0,  .0631955) 
	
/*Parameters from the best bayesian fit; highest logBF - logit independent igamma prior*/
	mat truemeans = ( -.1263567, -1.665534)  //0.41
	mat truevarcov = (1.824899, 0 \ 0,  .2331431) 
	
/*Parameters from the best bayesian fit; highest logBF - logit independent igamma prior*/
	mat truemeans = ( -.1258962, -1.649212)  //medians
	mat truevarcov = (1.824899, 0 \ 0,  .2331431)  //medians
	
	
	import delimited "C:\DATA\WIV\Projects\Stata\Metapreg\Data\immuno.csv", varnames(1) clear

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
		replicatefit, ///
			truevarcov(truevarcov) truemeans(truemeans) seed(`r') link(loglog)
	}
}
*===========================Process Simulated results===========================
{
forvalues loop=1/4 {
	import delimited "$wdir/simresults`loop'.txt", varnames(1) clear

	if `loop' > 1 {
		append using "$wdir/simresults.dta", force
	}
	
	save "$wdir/simresults.dta", replace
}
export delimited using "$simresults",   delimit(";") replace

import delimited "$simresults", varnames(1) clear
}

frame create simresults
frame change simresults
//Read the simulated results
import delimited "$simresults", varnames(1) clear

//Process the variables
*drop package- v10
do "$scriptdir\process.do"
save "$wdir/simulatedresults.dta", replace
