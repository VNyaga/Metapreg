
*Codesnippet 1
use "https://github.com/VNyaga/Metapreg/blob/master/dolmanrcog2014.dta?raw=true", clear

set more off
metapreg cured treated, model(random, sformat(%8.4f) ) ///
	studyid(study) cimethod(wilson) dp(4) ///
	nograph sumstat(ES) sumtable(abs)

set more off
 metaprop_one cured treated, ///
 lcols(study) random groupid(study) nograph logit dp(4)
 
*Codesnippet 2
use "https://github.com/VNyaga/Metapreg/blob/master/margins.dta?raw=true", clear

metapreg tpos n, studyid(study) sumtable(abs) dp(2) power(2) /// 
	graphregion(color(white)) by(treatment) stratify ///
	xlab(0, 25, 50, 75, 100) xtick(0, 25, 50, 75, 100) ///
	sumstat(Positive Margins (%)) ///
	texts(1.2)

**Figure 1
set more off
metapreg tpos n if treatment=="CKC", studyid(study) sumtable(abs) dp(2) power(2) /// 
	graphregion(color(white)) ///
	xlab(0, 25, 50, 75, 100) xtick(0, 25, 50, 75, 100) ///
	sumstat(+ve Margins) tit("CKC") fysize(35) ///
	texts(1.65) astext(80) graphsave("f1.gph")

*Graph 2 	
set more off
metapreg tpos n if treatment=="LC", studyid(study) sumtable(abs) dp(2) power(2) /// 
	graphregion(color(white)) ///
	xlab(0, 25, 50, 75, 100) xtick(0, 25, 50, 75, 100) ///
	sumstat(+ve Margins) tit("LC") ///
	texts(1.75) astext(80) graphsave("f2.gph")

set more off
metapreg tpos n if treatment=="LLETZ", studyid(study) sumtable(abs) dp(2) power(2) /// 
	graphregion(color(white)) ///
	xlab(0, 25, 50, 75, 100) xtick(0, 25, 50, 75, 100) ///
	sumstat(+ve Margins) tit("LLETZ")	 ///
	texts(1.5) astext(80) graphsave("f3.gph")

set more off
metapreg tpos n if treatment=="Mixed", studyid(study) sumtable(abs) dp(2) power(2) /// 
	graphregion(color(white)) ///
	xlab(0, 25, 50, 75, 100) xtick(0, 25, 50, 75, 100) ///
	sumstat(+ve Margins) tit("Mixed") ///
	texts(1.75) astext(80) graphsave("f4.gph")


gr combine f1.gph f3.gph, graphregion(color(white)) cols(1)  imargin(0 0 0 0)
graph save 13.gph, replace

gr combine f2.gph f4.gph, graphregion(color(white)) cols(1)  imargin(0 0 0 0)	
graph save 24.gph, replace

gr combine 13.gph 24.gph, cols(2) graphregion(color(white))  iscale(1)

**Codesnippet 3
set more off
metapreg tpos n treatment, studyid(study) sumtable(abs rr) dp(4)  /// 
	graphregion(color(white)) ///
	xlab(0, .25, .50, .75, 1) xtick(0, .25, .50, .75, 1) ///
	sumstat(+ve Margins (%)) ///
	texts(1.2)


**Codesnippet 4
use "https://github.com/VNyaga/Metapreg/blob/master/Build/bcg.dta?raw=true", clear

//Random - Adjusted
set more off
metapreg cases_tb population bcg lat, ///
	studyid(study) sortby(lat) ///
	sumtable(all) design(comparative) ///
	outplot(rr) interaction ///
	graphregion(color(white)) ///
	xlab(0.1, 1, 2) xtick(0.1, 1, 2) ///
	lcols(lat) rcols(cases_tb population) ///
	astext(80) texts(1.5) logscale ///
	sumstat(TB risk ratio)  ///
	xline(1, lcolor(black)) gof
	

//Random - Unadjusted
set more off
metapreg cases_tb population bcg, ///
	studyid(study) sortby(lat) ///
	sumtable(all) design(comparative) ///
	outplot(rr) ///
	graphregion(color(white)) ///
	xlab(0.1, 1, 2) xtick(0.1, 1, 2) ///
	lcols(lat) rcols(cases_tb population) ///
	astext(80) texts(1.5) logscale ///
	sumstat(TB risk ratio)  ///
	xline(1, lcolor(black)) gof
	

//Fixed - Adjusted
metapreg cases_tb population bcg lat, model(fixed) ///
	studyid(study) sortby(lat) ///
	sumtable(all) design(comparative) ///
	outplot(rr) interaction ///
	graphregion(color(white)) ///
	xlab(0.1, 1, 2) xtick(0.1, 1, 2) ///
	lcols(lat) rcols(cases_tb population) ///
	astext(80) texts(1.5) logscale ///
	sumstat(TB risk ratio)  ///
	xline(1, lcolor(black))	gof

//Fixed - Unadjusted
set more off
metapreg cases_tb population bcg, model(fixed) ///
	studyid(study) sortby(lat) ///
	sumtable(all) design(comparative) ///
	outplot(rr) ///
	graphregion(color(white)) ///
	xlab(0.1, 1, 2) xtick(0.1, 1, 2) ///
	lcols(lat) rcols(cases_tb population) ///
	astext(80) texts(1.5) logscale ///
	sumstat(TB risk ratio)  ///
	xline(1, lcolor(black)) gof
	

**Codesnippet 5
use "https://github.com/VNyaga/Metapreg/blob/master/ivmg_acute_myo_infarct_wide.dta?raw=true", clear
gen noeventc = totalc - eventc
gen noeventt = totalt - eventt
sort at study 
metan eventt noeventt eventc noeventc, rr model(fe) study(study) by(at)
metan eventt noeventt eventc noeventc, rr model(re) study(study) by(at)
metan eventt noeventt eventc noeventc, rr model(ivhet) study(study) by(at)


drop _*
reshape long event noevent total, i(study) j(group) string

sort at year study group 
set more off

//Fixed
metapreg event total group at, model(fixed) ///
studyid(study) sumtable(all) design(comparative) outplot(rr) ///
xlab(0, 1, 3) xtick(0, 1, 3) sortby(year) ///
texts(1.5) logscale sumstat(RR of death) ///
rcols(event total) astext(80) gof ///
graphregion(color(white)) 


//Random
set more off
metapreg event total group at,  ///
studyid(study) sumtable(all) design(comparative) outplot(rr) ///
xlab(0, 1, 3) xtick(0, 1, 3) sortby(year) ///
texts(2) logscale sumstat(RR of death) ///
astext(80)	///
rcols(event total) ///
graphregion(color(white)) gof

	
**Codesnippet 6 - contrast based network
use "https://github.com/VNyaga/Metapreg/blob/master/Build/matched.dta?raw=true", clear

netplot comparator index, label arrows type(circle)

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
    astext(80) texts(1.5) logscale ///
	sumstat(Rel. Sensitivity)	///
	xline(.9, lcolor(black)) dp(3) gof
	

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
    lcols(a b c d comparator) /// 
    astext(80) texts(1.5) logscale ///
	sumstat(Rel. Sensitivity)	///
	xline(.9, lcolor(black)) dp(3) gof
	

**Codesnippet 7 - stratified analysis
use "https://github.com/VNyaga/Metapreg/blob/master/maniacefficacy.dta?raw=true", clear

set more off
metapreg event total, stratify by(drug) ///
	studyid(study) model(random) sumtable(all) ///
	graphregion(color(white)) ///
	xlab(0, .5, 1) xtick(0, .5, 1)  texts(2.25) xline(0.5) ///
	sumstat(Response Rate) summaryonly nooverall

**Codesnippet 8
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

*Network plot 
netplot drug treatment if drug!=treatment, label arrows type(circle)  

**Codesnippet 9 - Comparative stratified analysis
set more off
metapreg event total drug if PLA, ///
	studyid(study) nomc  ///
	model(random) design(comparative) sumtable(all) ///
	outplot(rr) stratify by(treatment) ///
	graphregion(color(white) margin(zero)) ///
	xlab(.5, 1, 2, 2.5) xtick(.5, 1, 2, 2.5) ///
	sumstat(Response Rate Ratio)  ///
	texts(1.5) logscale xline(1) astext(40) xsize(7) ysize(10)

**Codesnippet 10 - arm based netwrok meta-analysis
use "https://github.com/VNyaga/Metapreg/blob/master/maniacefficacy.dta?raw=true", clear

sort drug study
set more off
metapreg event total drug, studyid(study) ///
	sumtable(all) outplot(rr) ///
	design(abnetwork, baselevel(PLA)) ///
	graphregion(color(white)) ///
	xlab(.5, 1, 1.5, 2) xtick(.5, 1, 1.5, 2) ///
	sumstat(Response Rate Ratio) ///
	texts(1.5) logscale xline(1)


	
