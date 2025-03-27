/*
Reference
Comparison of random-effects meta-analysis models for the relative risk in the case of rare events: A simulation study, 2020, 10.1002/bimj.201900379


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
global review "hemkens1"

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
foreach name in observedresults /*simresults*/  {
	local s = "sim; model; stat; esthat; esthatlo; esthatup; tau2hat; tau2hatlo; tau2hatup; sigma2hat; sigma2hatlo; sigma2hatup; IC; truemu0; trueor;  truetausq; truesigmasq; k; nstudies; studysize"
	file open results using $`name',  text write replace //append or replace
	file write results "`s'" _n
	file close results
}

global scriptdir "C:/DATA/WIV/Projects/GitHub/Metapreg/Scripts"
global Rterm_path `"C:\DATA\Software\R-4.4.2\bin\x64\Rterm.exe"'
global Rterm_options `"--vanilla"'

do "C:\DATA\WIV\Projects\GitHub\Metapreg\Build\metapreg.ado"

/*
global scriptdir "$rootdir/Scripts"
global Rterm_path `"C:\Program Files\R\R-4.4.0\bin\x64\Rterm.exe"'
global Rterm_options `"--vanilla"'

do "$scriptdir\metapreg.ado"
*/

}
*===========================ORIGINAL Long DATA===========================
{

frame create longdata
frame change longdata

import delimited "$rootdir/hemkens-cochrane-2016-analysis1-8.csv", varnames(1) clear

gsort study -treatment

rename events event

gen T = "T" if treatment ! = "Control"
replace T = "C" if treatment == "Control"

global bwd 	"C:\DATA\WIV\Projects\GitHub\Metapreg\Data\hemkens1"

gen group = treatment
gsort study -group

*encode group, gen(category)

*elasticnet poisson event category, exposure(total) alpha(0)

*firthlogit event category

*penlogit event category, binomial(total) or

metapreg event total group,  model(mixed) smooth gof  ///
	/*inference(bayesian) bwd(C:\DATA\WIV\Projects\Stata\Metapreg\mcmcresults) */ ///
	studyid(study) design(comparative, cov(independent)) sumtable(all) sumstat(Proportion)	///
	xlab(0, 0.05, 0.20)  ///
	lcols(event total) texts(2.5)  astext(75)  outplot(abs)   ///
	 graphregion(color(gray))  fxsize(90)	nooverall subline


metapreg event total group,  model(mixed) smooth gof  catpplot nofplot  ///
	/*inference(bayesian) bwd(C:\DATA\WIV\Projects\Stata\Metapreg\mcmcresults)*/  ///
	studyid(study) design(comparative, cov(independent))  sumstat(Risk Difference)	///
	xlab(-0.15, 0, 0.15)    ///
	texts(2.5)  astext(75)  outplot(rd)   ///
	xline(0) graphregion(color(gray))  fxsize(90)	
	
metapreg event total group,  model(mixed) smooth gof  catpplot nofplot  ///
	/*inference(bayesian) bwd(C:\DATA\WIV\Projects\Stata\Metapreg\mcmcresults) */ ///
	studyid(study) design(comparative, cov(independent)) sumstat(Risk Ratio)	///
	xlab(0, 1, 10) logscale  ///
	texts(2.5)  astext(75)  outplot(rr)  xline(1) ///
		graphregion(color(gray)) fxsize(80)
		
metapreg event total group,  model(mixed) smooth gof  catpplot nofplot  ///
	/*inference(bayesian) bwd(C:\DATA\WIV\Projects\Stata\Metapreg\mcmcresults)*/  ///
	studyid(study) design(comparative, cov(independent)) sumstat(Odds Ratio)	///
	xlab(0, 1, 10) logscale  ///
	texts(2.5)  astext(75)  outplot(or)  xline(1) ///
		graphregion(color(gray)) fxsize(80)	
		
gr combine absfplot rrcatpplot rdcatpplot  orcatpplot , ///
		graphregion(color(gray)  margin(tiny))  cols(2) imargin(0 0 0 0) ///
		 name(stats, replace)	ysize(10) xsize(17.5)	 
		 
		 
graph export "$graphs\data-F.png", as(png) width(2000) height(1000) replace	
		

*To compare with metafor
metapreg event total group if !inlist(study, "Yurdakul 2001", "Parise 1995"),  model(mixed) smooth gof  catpplot nofplot  ///
	/*inference(bayesian) bwd(C:\DATA\WIV\Projects\Stata\Metapreg\mcmcresults)*/  ///
	studyid(study) design(comparative, cov(independent)) sumstat(Odds Ratio)	///
	xlab(0, 1, 10) logscale  ///
	texts(2.5)  astext(75)  outplot(or)  xline(1) ///
		graphregion(color(gray)) fxsize(80)				
	 

	
	
//----------------------HEXACT
{
set more off
set trace off
metapreg event total group,  model(hexact) smooth gof  catpplot nofplot  ///
	studyid(study) design(comparative, cov(independent)) nowt 	///
	xlab(0, 1, 10) logscale ///
	texts(2.75)  astext(75)  outplot(rr or)  xline(1) ///
		graphregion(color(gray)) fxsize(90)	 
		
	
metapreg event total group,  model(fixed) smooth gof  catpplot nofplot  ///
	studyid(study) design(comparative, cov(independent)) nowt 	///
	xlab(0, 1, 10) logscale ///
	texts(2.75)  astext(75)  outplot(rr or)  xline(1) ///
		graphregion(color(gray)) fxsize(90)
}	
					
}
*===========================Fit all models===========================
{
	frame create widedata
	frame change widedata

	import delimited "$rootdir/hemkens-cochrane-2016-analysis1-8.csv", varnames(1) clear

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

/*
drop in 643/644
*/

do "$scriptdir\process.do"

gen ciwidth = esthatup - esthatlo


sort package inference AIC BIC DIC  Dist Design Covariance Link Slope Weighting Sigmethod CImethod loglest sumethod param Prior

format AIC BIC DIC logBF  sigma2hat tau2hat ciwidth  %10.2f

tostring ciwidth, gen(WidthCI) format(%4.2f) force

do "C:\DATA\WIV\Projects\Stata\Blobbogram\blobbogram.ado"
//Generate plots
mat optimalor = (.1686228, .0305341, .67442966,  .0028639, .52713314)
mat optimalabs = (0.03, 0.02, 0.04)
global text2add "metapreg:tau2=0.12, sigma2=0.22"
global text2addright "0.17 eti:(0.03, 0.67)  hpd:(0.00, 0.53)"


global simorrange "0., 0.15, 2"
global allrrange "0.01, 0.15, 1.5"

global rrrange "0, 1, 5"
global rrline "1"
mat optimalrr = (0.20, 0.06, 0.69)


global orrange "0.001, 0.10, 2"
global orline "1"


global rdrange "-.1, 0, 0.1"
global rdline "0"
mat optimalrd = (0.0205, 0.0058, 0.0351)

global orbest ".2113702"

global poprrbest ".2621639"

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
do "$scriptdir\Modelselection.do"    //Step 1: Go to modelselection.do

frame copy observedresults cleanresults, replace


do "$scriptdir\plotobsdata.do"   //Go to plotobsdata

frame copy observedresults cleanresults, replace
do "$scriptdir\widen.do"

//Step 2: Pick best model

gen conditional = 0
replace conditional = 1 if package == "metapreg" & link == "logit" & /// 
						(Covariance == "IND" & Design == "comparative") & Dist == "BN"  
*replace conditional = 1 if package == "metapreg" & link == "logit" & Covariance == "FI"  & Dist == "BN"  & inference == "bayesian"

//Step 3: Add BB
gen pick = conditional
replace pick = 1 if strpos(Dist, "BB") != 0  & link == "logit"
drop if strpos(stat, "pop") != 0 & strpos(model, "betabin") != 0


//Step 3.1: binomial
replace pick = 1 if Dist=="B2" & link == "logit" 
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
	/*Parameters from the best frequentist -logit- fit - fixed*/
	mat truemeans = (-3.637586, -1.631989)
	mat truevarcov = (0, 0 \  0,  0) 

	//Parameters from the best bayesian fit - unstructured
	mat truemeans = (-3.866031,  -1.886554)
	mat truevarcov = (.5418781, -.2858268 \  -.2858268,  .8584113 ) 
	
	//Parameters from the `best' bayesian fit - independent gamma
	mat truemeans = (-3.679317,  -1.989575)
	mat truevarcov = (.0882398,0 \ 0,  .1992353) 
	
		//Parameters from the `best' bayesian fit - independent gamma, poisson
	mat truemeans = (-3.722558,  -1.807533)  //medians
	mat truevarcov = ( .0744301 ,0 \ 0,  .1978213)  //medians
	
	//Parameters from the `best' bayesian fit - independent gamma, binomial-normal
	mat truemeans = (-3.695966,  -1.780091)  //medians
	mat truevarcov = ( .1177891 ,0 \ 0,  .2163384)  //medians
	
	import delimited "$rootdir/hemkens-cochrane-2016-analysis1-8.csv", varnames(1) clear

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
	forvalues r = 1(1)100 {
		replicatefit, ///
			truevarcov(truevarcov) truemeans(truemeans) seed(`r') 		
	}
}
*===========================Process Simulated results===========================
{
forvalues loop=1/4 {
	import delimited "$wdir/simresults3.txt", varnames(1) clear

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

 
 

