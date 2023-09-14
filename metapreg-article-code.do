
*Codesnippet 1
use "https://github.com/VNyaga/Metapreg/blob/master/dolmanrcog2014.dta?raw=true", clear

*Codesnippet 2 - hexact
set more off
metapreg cured treated, model(hexact) ///
	studyid(study) cimethod(exact, exact) ///
	sumtable(all) graphregion(color(white)) ///
	texts(2) astext(80)  ///
	xlab(0.5, 0.75, 1) xtick(0.5, 0.75, 1) ///
	 gof  smooth ///
	name(hexact, replace)  ysize(10) xsize(20) subti("CE")
		
	
*Codesnippet 3 - random
set more off
metapreg cured treated, model(random)  ///
	studyid(study) cimethod(exact)  ///
	 sumtable(all) graphregion(color(white)) ///
	texts(2) astext(80)  ///
	xlab(0.5, 0.75, 1) xtick(0.5, 0.75, 1) ///
	smooth gof  lcols(cured treated)  ///
	name(random, replace) ysize(10) xsize(10) subti("RE") 

*Figure 1
gr combine hexact  random , graphregion(color(white)) rows(1)   imargin(0 0 0 0) ysize(10) xsize(20) name(ex1, replace)

graph save Fig1.gph, replace	
graph export "$dir\Fig1.tif", as(tif) replace
graph export "$dir\Fig1.pdf", as(pdf) replace
graph export "$dir\Fig1.eps", as(eps) replace
graph export "$dir\Fig1.png", as(png) replace
 
*Codesnippet 4
use "https://github.com/VNyaga/Metapreg/blob/master/margins.dta?raw=true", clear

set more off
metapreg tpos n treatment, studyid(study) sumtable(abs rr) power(2) /// 
	graphregion(color(white)) ///
	xlab(0, 20, 40) xtick(0, 20, 40) ///
	sumstat(+ve Margins (%)) nowt ///
	texts(2)  gof summaryonly name(re, replace) ysize(10) xsize(10) subti("ME") 
	
set more off
metapreg tpos n treatment, model(fixed) studyid(study) sumtable(abs rr) power(2) /// 
	graphregion(color(white))  ///
	xlab(0, 20, 40) xtick(0, 20, 40) ///
	sumstat(+ve Margins (%)) ///
	texts(2)  gof summaryonly name(ce, replace) ysize(10) xsize(10) subti("FE") nowt

gr combine re  ce, graphregion(color(white)) rows(1)   imargin(0 0 0 0) ysize(10) xsize(20) name(cere, replace)

*Figure 2
graph save Fig2.gph, replace	
graph export "$dir\Fig2.tif", as(tif) replace
graph export "$dir\Fig2.pdf", as(pdf) replace
graph export "$dir\Fig2.eps", as(eps) replace
graph export "$dir\Fig2.png", as(png) replace

//stratified analysis
set more off
metapreg tpos n, studyid(study) sumtable(none) power(2) /// 
	graphregion(color(white)) by(treatment) stratify ///
	xlab(0, 50, 100) xtick(0, 50, 100) nowt sumstat(+ve Margins (%)) ///
	texts(2) summaryonly name(re, replace) ysize(10) xsize(20) subti("ME") 


**Codesnippet 7
use "https://github.com/VNyaga/Metapreg/blob/master/ivmg_acute_myo_infarct_wide.dta?raw=true", clear

gen noeventc = totalc - eventc
gen noeventt = totalt - eventt
sort study 

//CE - IV
set more off
metan eventt noeventt eventc noeventc, ///
	or model(fe) study(study)  ///
	forestplot(xsize(10) ysize(10) ///
	textsize(120) astext(40) ///
	graphregion(color(white)) ///
	subti("CE") name(ce, replace))

//RE - IV
set more off
metan eventt noeventt eventc noeventc, ///
	or model(re) study(study)  ///
	forestplot(xsize(10) ysize(10) ///
	textsize(120) astext(40) ///
	graphregion(color(white)) ///
	subti("RE") name(re, replace))

//IVHET - IV
set more off
metan eventt noeventt eventc noeventc, ///
	or model(ivhet) study(study)  ///
	forestplot(graphregion(color(white)) name(ivhet, replace))


*Figure 3
gr combine ce re, graphregion(color(white)) cols(2)  imargin(0 0 0 0) ysize(10) xsize(20) name(rr, replace)
graph save Fig3.gph, replace	
graph export "$dir\Fig3.tif", as(tif) replace
graph export "$dir\Fig3.pdf", as(pdf) replace
graph export "$dir\Fig3.eps", as(eps) replace
graph export "$dir\Fig3.png", as(png) replace


**Codesnippet 8
use "https://github.com/VNyaga/Metapreg/blob/master/ivmg_acute_myo_infarct_long.dta?raw=true", clear
sort study group 

//Fixed baselines + common OR 
set more off
metapreg event total group  study, model(fixed) ///
studyid(study) sumtable(all) design(comparative) outplot(rr) ///
xlab(0, 1, 3, 5) xtick(0, 1, 3, 5)  ///
texts(1.75) logscale sumstat(Relative risk of death) ///
lcols(event total) astext(80) gof boxopts(mcolor(green)) ///
graphregion(color(white)) smooth subtit("FE") name(fixedrr, replace) 


//random baselines + common OR 
set more off
metapreg event total group, model(mixed)  ///
studyid(study) sumtable(all) design(comparative) outplot(rr) ///
xlab(0, 1, 3, 5) xtick(0, 1, 3, 5)   ///
texts(1.75) logscale sumstat(Relative risk of death) ///
astext(80) gof boxopts(mcolor(green)) ///
graphregion(color(white))  smooth  subtit("ME") name(mixedrr, replace) 


estimates replay metapreg_modest , or cformat("%4.2f")

*Figure 4
gr combine fixedrr mixedrr, graphregion(color(white)) cols(2)  imargin(0 .5 0 0)ysize(10) xsize(20) name(rr, replace)
graph save Fig4.gph, replace	
graph export "$dir\Fig4.tif", as(tif) replace
graph export "$dir\Fig4.pdf", as(pdf) replace
graph export "$dir\Fig4.eps", as(eps) replace
graph export "$dir\Fig4.png", as(png) replace


//Fixed baselines + varying OR == Saturated model
set more off
metapreg event total group study,  ///
model(fixed) design(comparative) interaction outplot(rr) ///
studyid(study) sumtable(all) nograph noitable gof

//random baselines + random OR 
set more off
metapreg event total group, model(mixed)  ///
studyid(study) sumtable(all) design(comparative, cov(independent)) outplot(rr) ///
xlab(0, 1, 3) xtick(0, 1, 3)   ///
texts(1.75) logscale sumstat(Relative risk of death) ///
astext(80) lcols(event total) gof boxopts(mcolor(green)) ///
graphregion(color(white))  smooth  subtit("Independent covariance") name(mixedindrr, replace)  

estimates replay metapreg_modest , or cformat("%4.2f")


//random baselines + random OR ; correlated
set more off
metapreg event total group, model(mixed)  ///
studyid(study) sumtable(all) design(comparative, cov(unstructured)) outplot(rr) ///
xlab(0, 1, 3) xtick(0, 1, 3)  ///
texts(1.75) logscale sumstat(Relative risk of death) ///
astext(80) gof boxopts(mcolor(green))  ///
graphregion(color(white))  smooth  subtit("Unstructured covariance") name(mixedcorr, replace) 

estimates replay metapreg_modest , or cformat("%4.2f")

*Figure 5
gr combine mixedindrr mixedcorr, graphregion(color(white)) cols(2)  imargin(0 .5 0 0)ysize(10) xsize(20) name(rr, replace)

graph save Fig5.gph, replace	
graph export "$dir\Fig5.tif", as(tif) replace
graph export "$dir\Fig5.pdf", as(pdf) replace
graph export "$dir\Fig5.eps", as(eps) replace
graph export "$dir\Fig5.png", as(png) replace


set more off
metapreg event total group, model(mixed)  ///
studyid(study) sumtable(all) design(comparative, cov(unstructured)) outplot(abs) ///
xlab(0, .2, .4) xtick(0, .2, .4) ///
texts(1.35)  sumstat(Risk of death) ///
astext(80) lcols(event total) gof boxopts(mcolor(green))  ///
graphregion(color(white))  smooth  subtit("Unstructured covariance") name(mixedcorrabs, replace)  
 
estimates replay metapreg_modest , or cformat("%4.2f")


//Dealing with discrepant studies
gen discrepant = "No"

replace discrepant = "Yes" if study=="ISIS-4 1995b" |study=="MAGIC 2000" | study=="ISIS-4 1995a"


//random baselines + random OR ; correlated + interaction
set more off
metapreg event total group discrepant, model(mixed)  ///
studyid(study) sumtable(all) design(comparative, cov(unstructured)) outplot(rr) ///
xlab(0, 1, 3) xtick(0, 1, 3)   ///
texts(1.75) logscale sumstat(Relative risk of death) ///
astext(80) gof boxopts(mcolor(green))  interaction ///
graphregion(color(white))  smooth  subtit("ME") name(mixedcorr, replace) 

estimates replay metapreg_modest , or cformat("%4.2f")

//random baselines + interaction
sort discrepant study group
set more off
metapreg event total group discrepant, model(mixed)  ///
studyid(study) sumtable(all) design(comparative) outplot(rr) ///
xlab(0, 1, 3) xtick(0, 1, 3)  sortby(study) ///
texts(1.75) logscale sumstat(Relative risk of death) ///
astext(80) gof boxopts(mcolor(green)) lcols(event total)  interaction ///
graphregion(color(white)) nooverall smooth  subtit("ME") name(mixedrr, replace)

estimates replay metapreg_modest , or cformat("%4.2f")


//random baselines + interaction
gsort -discrepant study group
set more off
metapreg event total group discrepant, model(mixed)  ///
studyid(study) sumtable(all) design(comparative) outplot(abs) ///
xlab(0, .2, .4) xtick(0, .2, .4) sortby(study) ///  
texts(1.5) sumstat(Risk of death) ///
astext(80) gof boxopts(mcolor(green)) interaction ///
graphregion(color(white))  smooth  subtit("Working model") name(mixedabs, replace)

gr combine mixedcorrabs mixedabs, graphregion(color(white)) cols(2)  imargin(0 .5 0 0)ysize(10) xsize(20) name(rr, replace)
*Figure 6
graph save Fig6.gph, replace	
graph export "$dir\Fig6.tif", as(tif) replace
graph export "$dir\Fig6.pdf", as(pdf) replace
graph export "$dir\Fig6.eps", as(eps) replace
graph export "$dir\Fig6.png", as(png) replace


//constant baselines + interaction
encode study, gen(sid)
encode group, gen(treatment)

gen largertrial = 0
replace largertrial = 1 if discrepant == "Yes"
gen rest =1 - largertrial  //rest of studies
gen x = 1  //constant 

*Fit the model
set more off
binreg event x#c.largertrial i.sid#c.largertrial c.largertrial#i.treatment  ///
			 x#c.rest  	   i.sid#c.rest  	 c.rest#i.treatment  ///
	, noconstant n(total) ml
estat ic	

*Obtain the log odds
margins treatment, over(largertrial) predict(xb) post

*obtain the conditional OR
nlcom (largertrial: exp(_b[1.largertrial#2.treatment] - _b[1.largertrial#1.treatment])) ///
	  (rest: exp(_b[0.largertrial#2.treatment] - _b[0.largertrial#1.treatment])), ///
	   cformat("%4.2f")

*obtain the conditional RR
nlcom (largertrial: exp(ln(invlogit(_b[1.largertrial#2.treatment])) - ln(invlogit(_b[1.largertrial#1.treatment])))) ///
	  (rest: exp(ln(invlogit(_b[0.largertrial#2.treatment])) - ln(invlogit(_b[0.largertrial#1.treatment])))) ,  ///
	  cformat("%4.2f")



**Codesnippet 6
use "https://github.com/VNyaga/Metapreg/blob/master/Build/bcg.dta?raw=true", clear

//Random 
set more off
metapreg cases_tb population bcg lat, ///
	model(random, laplace) ///
	studyid(study)  sortby(lat) ///
	sumtable(all) design(comparative) ///
	outplot(rr) interaction  ///
	graphregion(color(white)) ///
	xlab(0.1, 1, 2) xtick(0.1, 1, 2) ///
	lcols(lat) rcols(cases_tb population) ///
	astext(80) texts(2) logscale ///
	sumstat(Relative risk of TB)  gof ///
	xline(1, lcolor(black)) smooth ysize(10) xsize(20) name(random, replace)

*Figure 7
graph save Fig7.gph, replace	
graph export "$dir\Fig7.tif", as(tif) replace
graph export "$dir\Fig7.pdf", as(pdf) replace
graph export "$dir\Fig7.eps", as(eps) replace
graph export "$dir\Fig7.png", as(png) replace
	
	
**Codesnippet 10 - contrast based network
use "https://github.com/VNyaga/Metapreg/blob/master/Build/matched.dta?raw=true", clear

*Figure 8
netplot comparator index, label arrows type(circle) 

graph save Fig8.gph, replace	
graph export "$dir\Fig8.tif", as(tif) replace
graph export "$dir\Fig8.pdf", as(pdf) replace
graph export "$dir\Fig8.eps", as(eps) replace
graph export "$dir\Fig8.png", as(png) replace

 //random
	set more off
	metapreg a b c d index comparator,  l(90) /// 
    studyid(study) ///
    model(random)  /// 
    sumtable(all) ///
    design(mcbnetwork)  ///
	by(index) ///
    outplot(rr) ///
    graphregion(color(white)) /// 
    xlab(0.9, 1, 1.1) /// 
    xtick(0.9, 1, 1.1)  ///
    astext(80) texts(1.75) logscale ///
	sumstat(Rel. Sensitivity)	///
	xline(.9, lcolor(black))  gof smooth subtit("ME") name(mixed, replace)

//fixed
set more off
metapreg a b c d index comparator,  l(90) /// 
    studyid(study) ///
    model(fixed)  /// 
    sumtable(all) ///
    design(mcbnetwork)  ///
	by(index) ///
    outplot(rr) ///
    graphregion(color(white)) /// 
    xlab(0.9, 1, 1.1) /// 
    xtick(0.9, 1, 1.1)  ///
    lcols(a b c d comparator) /// 
    astext(80) texts(1.75) logscale ///
	sumstat(Relative Sensitivity)	///
	xline(.9, lcolor(black))  gof smooth subtit("FE") name(fixed, replace)
	

*Figure 9
gr combine fixed mixed, graphregion(color(white)) cols(2)  imargin(0 0 0 0) ysize(10) xsize(20) name(cbnet, replace)

graph save Fig9.gph, replace	
graph export "$dir\Fig9.tif", as(tif) replace
graph export "$dir\Fig9.pdf", as(pdf) replace
graph export "$dir\Fig9.eps", as(eps) replace
graph export "$dir\Fig9.png", as(png) replace


**Codesnippet 12 - stratified analysis
use "https://github.com/VNyaga/Metapreg/blob/master/maniacefficacy_expanded.dta?raw=true", clear

*Network plots Figure 10
netplot drug treatment if drug!=treatment & PLA, label arrows type(circle)  
graph rename PLA, replace

netplot drug treatment if drug!=treatment, label arrows type(circle)  
graph rename ALL, replace

gr combine ALL PLA, graphregion(color(white)) cols(2)  imargin(0 0 3 3) ysize(10) xsize(20) name(network, replace)

graph save Fig10.gph, replace	
graph export "$dir\Fig10.tif", as(tif) replace
graph export "$dir\Fig10.pdf", as(pdf) replace
graph export "$dir\Fig10.eps", as(eps) replace
graph export "$dir\Fig10.png", as(png) replace

**Codesnippet 11 - stratified analysis
use "https://github.com/VNyaga/Metapreg/blob/master/maniacefficacy.dta?raw=true", clear

set more off
metapreg event total, stratify by(drug) ///
	studyid(study) model(random) sumtable(all) ///
	graphregion(color(white))  ///
	xlab(0, .5, 1) xtick(0, .5, 1) astext(70) texts(2) xline(0.5) ///
	sumstat(Response rate) summaryonly nooverall  subtit("Stratified RE")  name(abs, replace)

/*
**Reshaping code
egen sortorder = group(drug)

gen placebo = 0
replace placebo = 1 if drug =="PLA"
gsort study -placebo
bys study: egen T = seq()

bys study: egen nt = count(drug)
count
global nobs = r(N)

gen nexpand = 1
replace nexpand = 2 if placebo & nt == 3
expand nexpand

replace T = 4 if nt == 3 & T == 3

replace T = 3 if _n>$nobs

gsort study -T

bys study : egen PLA = max(placebo) //studies with placebo
egen O = seq(),f(1) t(2) b(1)
gen treatment = drug if O == 1
replace treatment = treatment[_n - 1] if O == 2

bys treatment: egen Nstudies = count(study)
replace Nstudies = 0.5*Nstudies

gsort treatment study -placebo

gen group = "a-ARI" if treatment=="ARI"
replace group = "b-PLA" if treatment=="PLA"
replace group = "c-HAL" if treatment=="HAL"
replace group = "d-QUE" if treatment=="QUE"
replace group = "e-LITH" if treatment=="LITH"
replace group = "f-ZIP" if treatment=="ZIP"
replace group = "g-OLA" if treatment=="OLA"
replace group = "h-DIV" if treatment=="DIV"
replace group = "i-RIS" if treatment=="RIS"
replace group = "j-CARB" if treatment=="CARB"
replace group = "k-LAM" if treatment=="LAM"
replace group = "l-PAL" if treatment=="PAL"
replace group = "m-TOP" if treatment=="TOP"
replace group = "n-ASE" if treatment=="ASE"

gsort PLA group study -placebo
*/

**Codesnippet 12 - stratified analysis
use "https://github.com/VNyaga/Metapreg/blob/master/maniacefficacy_expanded.dta?raw=true", clear

gsort PLA group study -placebo

set more off
metapreg event total drug if PLA,  ///
	studyid(study) nomc  ///
	model(random) design(comparative) sumtable(all) ///
	outplot(rr) stratify by(group) sortby(study) ///
	graphregion(color(white) margin(zero)) ///
	xlab(.5, 1, 2, 2.5) xtick(.5, 1, 2, 2.5) ///
	sumstat(Relative response rate)  ///
	texts(1.75) logscale xline(1)  summaryonly  name(rrpla, replace) subtit("Stratified ME")


**Codesnippet 10 - arm based netwrok meta-analysis
use "https://github.com/VNyaga/Metapreg/blob/master/maniacefficacy.dta?raw=true", clear

set more off

gsort study drug 
metapreg event total drug, studyid(study) ///
	sumtable(all) outplot(rr)  ///
	design(abnetwork, baselevel(PLA)) ///
	graphregion(color(white)) ///
	xlab(.5, 1, 1.5, 2) xtick(.5, 1, 1.5, 2) ///
	sumstat(Relative response rate) ///
	texts(2) logscale xline(1) nooverall  name(abnetrr, replace) subtit("AB Network") 

*Figure 12
gr combine rrpla abnetrr, graphregion(color(white)) cols(2)  imargin(0 0 0 0) ysize(10) xsize(20) name(abnet, replace)

graph save Fig12.gph, replace	
graph export "$dir\Fig12.tif", as(tif) replace
graph export "$dir\Fig12.pdf", as(pdf) replace
graph export "$dir\Fig12.eps", as(eps) replace
graph export "$dir\Fig12.png", as(png) replace

gsort study drug 
set more off
metapreg event total drug, studyid(study) ///
	sumtable(all) outplot(abs)  ///
	design(abnetwork, baselevel(PLA)) ///
	graphregion(color(white)) ///
	xlab(0, .5, 1) xtick(0, .5, 1)  texts(2) xline(0.5) ///
	sumstat(Response rate) summaryonly nooverall  name(abnetabs, replace) subtit("AB Network")  smooth nowt/* ysize(10) xsize(10)*/


*Figure 11
gr combine abs abnetabs, graphregion(color(white)) cols(2)  imargin(0 0 0 0) ysize(10) xsize(20) name(abnet, replace)

graph save Fig11.gph, replace	
graph export "$dir\Fig11.tif", as(tif) replace
graph export "$dir\Fig11.pdf", as(pdf) replace
graph export "$dir\Fig11.eps", as(eps) replace
graph export "$dir\Fig11.png", as(png) replace
