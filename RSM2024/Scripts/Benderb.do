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

global review "bender2018b"
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
global Rterm_path `"C:\DATA\Software\R-4.4.2\bin\x64\Rterm.exe"'

*global Rterm_path `"C:\Program Files\R\R-4.4.0\bin\x64\Rterm.exe"'
global Rterm_options `"--vanilla"'

*do "$scriptdir\metapreg.ado"
do "C:\DATA\WIV\Projects\GitHub\Metapreg\Build\metapreg.ado"	

}
*===========================ORIGINAL Long DATA===========================
{
frame create longdata
frame change longdata

import excel using "$rootdir/Data/metadat.xlsx",sheet("Bender2018b") ///
      firstrow clear

//To long format
reshape long event total, i(Study) j(class)

gen treatment = "Control" if class==0
replace treatment = "Treatment" if class==1	

gen T = "T" if treatment ! = "Control"
replace T = "C" if treatment == "Control"

rename Study study
gen group = treatment

gsort study group
set more off
set trace off

#set rmsg on 

metapreg event total group,   smooth gof model(mixed) ///
	/*inference(bayesian) bwd(C:\DATA\WIV\Projects\Stata\Metapreg\mcmcresults)*/  ///
	studyid(study) design(comparative, cov(independent)) sumtable(all) sumstat(Proportion)	///
	xlab(0, 0.25, 0.5)  ///
	lcols(event total) texts(2.5)  astext(80)  outplot(abs)   ///
	 graphregion(color(gray))  fxsize(90)	nooverall subline

#F = 2.04
#B = 41.77

metapreg event total group,   smooth gof  catpplot nofplot  ///
	/*inference(bayesian) bwd(C:\DATA\WIV\Projects\Stata\Metapreg\mcmcresults)*/  ///
	studyid(study) design(comparative, cov(independent)) sumstat(Risk Difference)	///
	xlab(-0.5, -.25, 0)  ///
	texts(2.5)  astext(70)  outplot(rd)   ///
	xline(0) graphregion(color(gray))  fxsize(90)	

/*	
metapreg event total group,   smooth gof  catpplot nofplot  ///
	inference(bayesian) bwd(C:\DATA\WIV\Projects\Stata\Metapreg\mcmcresults)  ///
	studyid(study) design(comparative, cov(independent))  	///
	xlab(1, 10, 50) logscale  ///
	texts(2.25)  astext(70)  outplot(rr or)  xline(1) ///
		graphregion(color(gray)) fxsize(80)
*/	
	
metapreg event total group,   smooth gof  catpplot nofplot  ///
	/*inference(bayesian) bwd(C:\DATA\WIV\Projects\Stata\Metapreg\mcmcresults)*/  ///
	studyid(study) design(comparative, cov(independent)) sumstat(Risk Ratio) 	///
	xlab(1, 10, 50) logscale  ///
	texts(2.5)  astext(70)  outplot(rr)  xline(1) ///
		graphregion(color(gray)) fxsize(80)	

metapreg event total group,   smooth gof  catpplot nofplot  ///
	/*inference(bayesian) bwd(C:\DATA\WIV\Projects\Stata\Metapreg\mcmcresults)*/  ///
	studyid(study) design(comparative, cov(independent)) sumstat(Odds Ratio) 	///
	xlab(1, 10, 50) logscale  ///
	texts(2.5)  astext(70)  outplot(or)  xline(1) ///
		graphregion(color(gray)) fxsize(80)			
	 
gr combine absfplot rrcatpplot rdcatpplot  orcatpplot , ///
		graphregion(color(gray)  margin(zero))  cols(2) imargin(0 0 0 0) ///
		 name(stats, replace)	ysize(10) xsize(17.5)	 
		 

graph export "$graphs\data-F.png", as(png) width(2000) height(1000) replace	

		 
}
*===========================Fit all models===========================
{
	frame create widedata
	frame change widedata
	import excel using "$rootdir/Data/metadat.xlsx",sheet("Bender2018b") ///
		  firstrow clear

	rename Study studyid
		 
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

*export delimited using "$observedresults",   delimit(";") replace

do "$scriptdir\process.do"

gen ciwidth = esthatup - esthatlo

sort Env package inference AIC BIC DIC Dist Design Covariance Link Slope Weighting Sigmethod CImethod loglest sumethod param Prior

format AIC BIC DIC logBF sigma2hat tau2hat ciwidth  %10.2f

tostring ciwidth, gen(WidthCI) format(%4.2f) force

do "C:\DATA\WIV\Projects\Stata\Blobbogram\blobbogram.ado"
//Generate plots
mat optimalor = (3.418909, 1.584143, 7.1570483, 1.293192, 6.5276114)   //eti
mat optimalabs = (0.06, 0.04, 0.07)
global text2add "metapreg:tau2=0.05, sigma2=0.06"
global text2addright "3.42 eti:(1.58, 7.16)  hpd:(1.29, 6.53)"
global obsrange "0, 0.5"

global obsrange "0, 0.1"

global simorrange "0.15, 1, 20"

global orrange "0.5, 1, 3, 20"
global orline "1"

global allrrange "1, 3, 10"

global rrrange "0, 1, 5"
global rrline "1"
mat optimalrr = (1.49, 1.14, 1.95)

global rdrange "-.1, 0, 0.1"
global rdline "0"
mat optimalrd = (-0.028036, -0.05083, -0.00633)


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
replace conditional = 1 if package == "metapreg" & link == "logit"  & ///
						(Covariance == "IND"  & Design == "comparative") & Dist == "BN"  

//Step 3: Add BB
gen pick = conditional
replace pick = 1 if strpos(Dist, "BB") != 0 & link == "logit" 

//Step 3.1: binomial
replace pick = 1 if inlist(Dist, "B2") & link == "logit" 
*replace pick = 1 if Dist=="B1" & Covariance=="CI" & link == "logit" 
*br if package=="metapreg" &CImethod =="exact"

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

	//Parameters from the best bayesian fit; highest* logBF - independent; good mixing
	mat truemeans = (-2.073778 , 1.233008)  //median
	mat truevarcov = (.0532203, 0 \ 0,  .0667226) //meadian
	
	import excel using "$rootdir/Data/metadat.xlsx",sheet("Bender2018b") ///
      firstrow clear

	rename study studyid
	 
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
