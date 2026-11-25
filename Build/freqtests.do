set more off
	
cap log close
log using "C:\DATA\WIV\Projects\Stata\Metapreg\Logs\freqtests.txt", text replace

global inference 
global wdir 


//1.1
	use "http://fmwww.bc.edu/repec/bocode/a/arbyn2009jcellmolmedfig1.dta", clear
	
	metapreg num denom,   ///
		studyid(study) ///
		by(tgroup) ///
		inference($inference) bwd($wdir)  /// 
		cimethod(exact) ///
		label(namevar=author, yearvar=year) catpplot  ///
		xlab(.25, 0.5, .75, 1) ///
		subti(Atypical cervical cytology, size(4)) ///
		texts(1.25) smooth gof	


//1.2

	use "http://fmwww.bc.edu/repec/bocode/a/arbyn2009jcellmolmedfig1.dta", clear

	bys tgroup, rc0: metapreg num denom,  ///
		inference($inference) bwd($wdir) /// 
		studyid(study) ///
		cimethod(wilson) ///
		label(namevar=author, yearvar=year) catpplot ///
		xlab(0, .25, 0.5, .75, 1)  ///
		subti(Atypical cervical cytology, size(4)) ///
		graphregion(color(white)) /// 
		texts(1.5) smooth


//1.3

	use "http://fmwww.bc.edu/repec/bocode/t/tsoumpou2009cancertreatrevfig2WNL.dta", clear

	metapreg p16p p16tot, ///		
		studyid(study)  ///
		label(namevar=author, yearvar=year) catpplot  ///
		inference($inference) bwd($wdir)  /// 
		sortby(year author) ///
		xlab(0, .1, .2, 0.3, 0.4, 0.5) ///
		xline(0, lcolor(black)) ///
		ti(Positivity of p16 immunostaining, size(4) color(blue)) ///
		subti(Cytology = WNL, size(4) color(blue)) ///
		pointopt(msymbol(X)msize(2)) /// 
		texts(1.5) smooth gof


//1.4

	use "http://fmwww.bc.edu/repec/bocode/t/tsoumpou2009cancertreatrevfig2WNL.dta", clear
		
	metapreg p16p p16tot, link(loglog) ///		
		studyid(study) sumtable(all) ///
		label(namevar=author, yearvar=year) catpplot ///
		sortby(year author)  ///
		xlab(0, .1, .2, 0.3, 0.4, 0.5) ///
		xline(0, lcolor(black)) ///
		ti(Positivity of p16 immunostaining, size(4) color(blue)) ///
		subti(Cytology = WNL, size(4) color(blue)) ///
		pointopt(msymbol(X)msize(2)) /// 
		texts(1.5) smooth gof

//2.1

	use "http://fmwww.bc.edu/repec/bocode/a/arbyn2009jcellmolmedfig1.dta", clear

	metapreg num denom tgroup,  ///
		inference($inference) bwd($wdir) /// 
		studyid(study) ///
		cimethod(exact) ///
		sumtable(all)  ///
		label(namevar=author, yearvar=year)  ///
		xlab(.25, 0.5, .75) ///
		subti(Atypical cervical cytology, size(4)) ///
		texts(1.5)  summaryonly 

//2.2.1

	use "https://github.com/VNyaga/Metapreg/blob/master/Build/hemkens2016analysis110.dta?raw=1", clear

	gsort study -treatment

	metapreg events total treatment,   ///
	studyid(study)  design(comparative, cov(independent))   ///
	inference($inference) bwd($wdir) ///
	smooth gof  catpplot nofplot   ///
	outplot(rr)  xline(1) sumstat(Risk Ratio)  ///	
	xlab(0, 1, 30) logscale  ///
	texts(2.35)  astext(80) 



//2.2.2 Three studies -  bender 2018 fig2

	use "https://github.com/VNyaga/Metapreg/blob/master/Build/bender2018fig2.dta?raw=1", clear

	metapreg event total treatment,  ///
	studyid(study) design(comparative, cov(independent))  ///
	smooth gof  catpplot nofplot cimethod(,wald)  ///
	inference($inference) bwd($wdir) ///
	outplot(rr)  xline(1) sumstat(Risk Ratio)  	///
	xlab(1, 5, 30) logscale  ///
	texts(2) astext(70) 

//2.2.3 - Extreme data sparsity

	use "https://github.com/VNyaga/Metapreg/blob/master/Build/hemkens2016analysis18.dta?raw=1", clear

	gsort study -treatment	
		
	metapreg event total treatment,   ///
		studyid(study) ///
		design(comparative, cov(commonint))  ///
		inference($inference) bwd($wdir) ///
		smooth gof  catpplot nofplot  ///
		outplot(rr) xline(1) sumstat(Risk Ratio)	///
		xlab(0.01, 1, 100) logscale  ///
		texts(1.75) astext(60)

//2.2.4
	use "http://fmwww.bc.edu/repec/bocode/b/bcg.dta", clear
	
	
	#delimit ;
	metapreg cases_tb population bcg, 
		studyid(study) model(mixed, intmethod(mv))
		design(comparative, cov(commonslope))	
		inference($inference) bwd($wdir) 
		outplot(rr) 
		sumstat(Risk ratio) 
		xlab(0, 1, 2) 
		xtick(0, 1, 2)  
		rcols(cases_tb population) 
		astext(80) 
		texts(1.5) logscale smooth gof
		;
	#delimit cr	

//2.2.5
	use "http://fmwww.bc.edu/repec/bocode/b/bcg.dta", clear
	
	metapreg cases_tb population lat,  /// 
		inference($inference) bwd($wdir) /// 
		studyid(study) model(mixed, intmethod(mv)) ///
		sortby(lat) by(bcg) ///
		xlab(0, 0.05, 0.15) /// 
		xtick(0, 0.05, 0.15)  /// 
		sumstat(Proportion) ///
		rcols(cases_tb population) /// 
		astext(80) /// 
		texts(1.5) smooth gof		

//2.2.6
	use "http://fmwww.bc.edu/repec/bocode/b/bcg.dta", clear
	
	metapreg cases_tb population bcg lat,  /// 
		inference($inference) bwd($wdir) /// 
		studyid(study) model(mixed, intmethod(mv)) ///
		sortby(lat) ///
		design(comparative, cov(commonslope))  ///
		outplot(rr) ///
		interaction ///
		xlab(0, 1, 2) /// 
		xtick(0, 1, 2)  /// 
		rcols(cases_tb population) /// 
		astext(80)  /// 
		texts(1.5) logscale smooth gof

//2.2.7  
	use "http://fmwww.bc.edu/repec/bocode/s/schizo.dta", clear

	set more off

	gsort firstauthor -arm

	metapreg response total arm missingdata, /// 
		studyid(firstauthor) ///
		sortby(year)  ///
		sumtable(all) ///
		design(comparative, cov(commonslope))  ///
		inference($inference) bwd($wdir)  /// 
		outplot(rr) ///
		interaction ///
		xlab(0, 0.5, 1.5) /// 
		xtick(0, .5, 1.5)  /// 
		sumstat(Rel Ratio) ///
		lcols(response total) /// 
		astext(70) texts(1.5)  /// 
		logscale smooth xsize(12) ysize(6) gof

//2.3.1

	use "https://github.com/VNyaga/Metapreg/blob/master/Build/repro.dta?raw=1", clear
	
	set more off

	metapreg pp pn np nn,  ///
		design(mpair, cov(commonslope)) ///
		inference($inference) bwd($wdir)  /// 
		studyid(paper) sumtable(all) ///
		stratify by(type) ///
		xlab(0.5, 1, 2) ///
		sumstat(Positivity Ratio) ///
		lcols(test)  ///
		boxopts(msize(0.1) mcolor(black)) pointopt(msymbol(none)) 	///
		astext(50)  logscale  xline(1) smooth 

//2.3.2
	use "http://fmwww.bc.edu/repec/bocode/m/matched.dta", clear
			
	metapreg a b c d index comparator,  /// 
		inference($inference) bwd($wdir)  /// 
		studyid(study)  ///
		model(fixed)  /// 
		design(mcbnetwork) ///
		by(comparator) ///
		xlab(0.9, 1, 1.1) /// 
		xtick(0.9, 1, 1.1)  /// 
		sumstat(Ratio) ///
		lcols(a b c d comparator index) /// 
		astext(80) /// 
		texts(1.5) logscale smooth gof