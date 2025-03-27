/*
Diuretics
Reference:
Beta-binomial model for meta-analysis of odds ratios, 2017, 10.1002/sim.7233

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
global review "diuretics"
global rootdir "C:/DATA/WIV/Projects/GitHub/Metapreg"
*global rootdir "M:/1120_Cancer_Employee/BMH/Victoria/GitHub/Metapreg"

global graphs "C:/DATA/WIV/Projects/GitHub/Metapreg/Data/$review/Graphs"

mkdir "$rootdir/Data/$review"
mkdir "$graphs"

global wdir "$rootdir/Data/$review"

global observedresults "$wdir/$review.txt"
global simresults "$wdir/simresults.txt"

global scriptdir "C:/DATA/WIV/Projects/GitHub/Metapreg/Scripts"

//Heads for the results file 
foreach name in observedresults /*simresults*/ {
	local s = "sim; model; stat; esthat; esthatlo; esthatup; tau2hat; tau2hatlo; tau2hatup; sigma2hat; sigma2hatlo; sigma2hatup; IC; truemu0; trueor;  truetausq; truesigmasq; k; nstudies; studysize"
	file open results using $`name',  text write replace //append or replace
	file write results "`s'" _n
	file close results
}

global Rterm_path `"C:\DATA\Software\R-4.4.2\bin\x64\Rterm.exe"'
global Rterm_options `"--vanilla"'

do "C:\DATA\WIV\Projects\GitHub\Metapreg\Build\metapreg.ado"	

}
*===========================ORIGINAL Long DATA===========================
{
frame create obsdata
frame change obsdata

import delimited "$rootdir/$review.csv", varnames(1) clear

gsort trial group

*gen study = "Study " + string(trial)
rename group treatment
rename events event

gen T = "T" if treatment ! = "Control"
replace T = "C" if treatment == "Control"


metapreg event total treatment,  model(mixed) smooth gof  ///
	/*inference(bayesian) bwd(C:\DATA\WIV\Projects\Stata\Metapreg\mcmcresults) */ ///
	studyid(study) design(comparative, cov(independent)) sumtable(all) sumstat(Proportion)	///
	xlab(0, 0.30, 0.75)  ///
	lcols(event total) texts(2)  astext(75)  outplot(abs)   ///
	 graphregion(color(gray))  fxsize(90)	nooverall subline


metapreg event total treatment,  model(mixed) smooth gof  catpplot nofplot  ///
	/*inference(bayesian) bwd(C:\DATA\WIV\Projects\Stata\Metapreg\mcmcresults) */ ///
	studyid(study) design(comparative, cov(independent))  sumstat(Risk Difference)	///
	xlab(-0.15, 0, 0.15)    ///
	texts(2)  astext(75)  outplot(rd)   ///
	xline(0) graphregion(color(gray))  fxsize(90)	
	
metapreg event total treatment,  model(mixed) smooth gof  catpplot nofplot  ///
	/*inference(bayesian) bwd(C:\DATA\WIV\Projects\Stata\Metapreg\mcmcresults) */ ///
	studyid(study) design(comparative, cov(independent)) sumstat(Risk Ratio)	///
	xlab(0.15, 1, 2) logscale  ///
	texts(2)  astext(75)  outplot(rr)  xline(1) ///
		graphregion(color(gray)) fxsize(80)
		
metapreg event total treatment,  model(mixed) smooth gof  catpplot nofplot  ///
	/*inference(bayesian) bwd(C:\DATA\WIV\Projects\Stata\Metapreg\mcmcresults)*/  ///
	studyid(study) design(comparative, cov(independent)) sumstat(Odds Ratio)	///
	xlab(0.15, 1, 2) logscale  ///
	texts(2)  astext(75)  outplot(or)  xline(1) ///
		graphregion(color(gray)) fxsize(80)	

gr combine absfplot rrcatpplot rdcatpplot  orcatpplot , ///
		graphregion(color(gray)  margin(tiny))  cols(2) imargin(0 0 0 0) ///
		 name(stats, replace)	ysize(10) xsize(17.5)	 
		 
		 
graph export "$graphs\data-F-IND.png", as(png) width(2000) height(1000) replace	

	
					
}
*===========================Fit all models===========================
{
	frame create widedata
	frame change widedata
	import delimited "$rootdir/$review.csv", varnames(1) clear

	gsort trial group

	gen studyid = "Study " + string(trial)
	rename group treatment
	rename events event

	//To wide format
	gen code = 0 if treatment == "Control"
	replace code = 1 if treatment != "Control"


	reshape wide event total treatment, i(studyid ) j(code)
	
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

*export delimited using "$observedresults",  delimit(";") replace

qui do "$scriptdir\process.do"

gen ciwidth = esthatup - esthatlo

sort package inference AIC BIC DIC Dist Design Covariance Link Slope Weighting Sigmethod CImethod loglest sumethod param Prior

format AIC BIC DIC logBF sigma2hat tau2hat ciwidth  %10.2f

tostring ciwidth, gen(WidthCI) format(%4.2f) force

do "C:\DATA\WIV\Projects\Stata\Blobbogram\blobbogram.ado"
//Generate plots
mat optimalor = (.5972605, .3652776, 1.079219, .3305971, 1.002854)
mat optimalabs = (0.14, 0.07, 0.26)
global text2add "metapreg:tau2=1.90, sigma2=0.36"
global text2addright "0.60 eti:(0.37, 1.08)  hpd:(0.33, 1.00)"

global obsrange "0, 0.5"

global orrange "0.25, 1, 1.5"
global orline "1"

global simorrange "0.25, 1, 2"
global allrrange "0.25, 1, 2"

global rrrange "0, 1, 5"
global rrline "1"
mat optimalrr = (0.63, 0.44, 0.92)

global rdrange "-0.1, 0, 0.1"
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

//Step 3.1: binomial
*replace pick = 1 if Dist=="B2" & link == "logit"  //exact not feasible
*replace pick = 1 if Dist=="B1" & Covariance=="CI" & link == "logit" //

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
/*Parameters from the best fit - abnetwork logit*/
mat truemeans = (-1.830137,  -2.346125)
mat truevarcov = (1.290397 , 0 \  0,  .1334057) 

/*Parameters from the best fit: comparative, independent, logit*/
mat truemeans = (-1.827625,  -.5185274)
mat truevarcov = ( 1.351328 , 0 \  0,  .2297258 ) 

/*Parameters from the best frequentist fit: comparative, correlated, logit*/
mat truemeans = (-1.829863,  -.5163397)
mat truevarcov = (1.411039 ,  -.1202088  \   -.1202088 ,  .2642035 )

/*Parameters Bayesian fit: comparative, correlated, logit; jeffreys*/
mat truemeans = (-1.82953,  -.5201875)
mat truevarcov = (1.796531 ,  -.0582516   \   -.0582516  ,  .5875037)

/*Parameters Bayesian fit, highest logBF: comparative, independent, logit; jeffreys*/
mat truemeans = (-1.835723,  -.5061944)
mat truevarcov = (1.73345 ,  0  \   0  ,  .2686132)

/*Parameters Bayesian fit, highest logBF: comparative, independent, logit; igamma(0.01, 0.01)*/
/*Trace plots indicate convergence and good mixing*/
mat truemeans = (-1.815705,  -.5035812) //medians
mat truevarcov = (1.873105 ,  0  \   0  ,  .3493668)  //median

/*Parameters Bayesian fit, highest logBF: comparative, commonslope, logit; igamma(0.01, 0.01)*/
/*Trace plots indicate convergence and good mixing*/
*mat truemeans = (-1.838223,  -.4088812) //medians
*mat truevarcov = (1.678856 ,  0  \   0  ,  0)  //median

	
	import delimited "C:\DATA\WIV\Projects\Stata\Metapreg\Data\diuretics.csv", varnames(1) clear

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
			truevarcov(truevarcov) truemeans(truemeans) seed(`r') 
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
*===========================Process Simulated results ABS===========================
