/*
CREATED:	8 Sep 2017
AUTHOR:		Victoria N Nyaga
PURPOSE: 	Generalized linear fixed, mixed & random effects modelling of binomial data.
VERSION: 	4.0.0
NOTES
1. Variable names and group names should not contain underscore(_)
2. Data should be sorted and no duplicates
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
UPDATES
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
DATE:						DETAILS:
24.08.2020
							grid: Grid lines between studies
							noveral in help file changed to noOVerall
							Print full matrix of rr when data is not repeated
							Print # of studies
							Correct the tick & axis position for sp
							Correct computation of the I2
							graphsave(filename) option included
03.09.2020					Change paired to comparative
14.19.2020					Correct way of counting the distinct groups in the meta-analysis
							Check to ensure variable names do no contain underscore.
12.02.2021					paired data: a b c d comparator index covariates, by(byvar)
							comparator, index, byvar need to be string
							Need to test more with more covariates!!!
09.07.2021					cimodel > cimethod
							Overal isq not showing with dp>2
01.02.2022					Change paired to matched
							paired data: n1 n2 N comparator index covariates, by(byvar)
							Subgroup analysis with superimposed graphs; stratify option
							stratify not an option for paired, matched or network
15.03.2022					design(independent|matched|paired|comparative|network, baselevel(string))
						    network data: n N Assignment covariates ....repeated measurements per study	
27.05.2022					Corrections on absoutp
							stratify + comparative + outplot(RR) 
11.02.2023					change network to abnetwork, paired to pcbnetwork, matched to mcbnetwork							
							Work on: hetout, stratify marginal results with 1 study
07.03.2023					change independent to general
21.03.2023					Compute weights from the maximized log likelihood
							Exact inference in few studies
18.04.2023					Use t-distribution for summaries
05.06.2023					Simulate posterior distributions
							smooth:Option to generate smooth estimates							
26.07.2023 					if version 16; use melogit instead of meqrlogit	
20.09.2023					Include outplot(OR)	
23.10.2023					nsims;how many times to simulate the posterior distributions	
31.10.2023					cloglog ink (might be better if p > .9)	loglog ink (might be better if p < .1)	
10.01.2024					Fit FE/Hexact if RE fails.
							If isq = ., suppress the text in the graph
01.02.2024					Introduce beta-binomial regression :cbbetabin - common beta beta-binomial, crbetabin - common rho beta-binomial	
05.02.2024					Introduce catterpillar plot	
							seperate forest and catterpillar plot options
26.02.2024					Introduce Bayesian in version >16.1
							inference(frequentist|bayesian)
							More options for rr ci's; katz, bailey, noether, asihn
28.03.2024					popstat(median|mean)
11.04.2024					cov(commonint|commonslope)
18.04.2024					Request for working directory if bayesian and inform the user to delete the datasets afterwards	
19.04.2024					Introduce outplot(lor|lrr)	
04.05.2024					aliasdesign; for abnetwork/general processing as comparative
24.05.2024					link(log); only for fixed|mixed and frequentist
							cov(freeint) - each study has its own intercept aka stratified intercepts
04.06.2024					outplot(rd)	
24.06.2024					Print all table unless otherwise
02.07.2024					Generate graphs in all scales	
25.10.2024					Regresss a b c d data when there is one test, i.e simplified mcbnetwork.matched-pair, introduce design(mpair)	
25.02.2025					Do not print marginal estimates of logits, p, rd and rr if cov(freeint) since they are irrelevant.	
28.04.2025					Work with encoded variables as well as string as categorical
09.05.2025					popstat(median|mean) ------> stat(median|mean)	
15.05.2025					ocimethod ------>scimethod summary CI method
20.05.2025					Remove prediction					
*/

/*++++++++++++++++++++++	METAPREG +++++++++++++++++++++++++++++++++++++++++++
						WRAPPER FUNCTION
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
cap program drop metapreg
program define metapreg, eclass sortpreserve byable(recall)

	version 14.1
	
	if _caller() >= 16 {
		version 16.1
	}

	#delimit ;
	syntax varlist(min=2) [if] [in], 
		STudyid(varname) [
		
		/*Model options*/
		Model(string) //model(random|mixed|fixed|hexact|cbbetabin|crbetabin, options)
		INFerence(string) //inference(FREQuentist|BAYesian)
		BWD(string asis) /*working directory to save bayesian estimates*/
		link(string) /*logit|cloglog|loglog|log*/
		DESign(string asis) //design(general|mpair|mcbnetwork|paired-pcbnetwork|comparative|network-abnetwork, baselevel(string) | cov(commonslope|commonint|freeint|inde|unstr))
		ALIASdesign(string) //aliasdesign(comparative)
		MC /*Model comparison - Saves time*/
		PROGress /*See the model fitting*/
		by(varname)
		STRatify  /*Stratified analysis, requires byvar()*/	
		ALphasort
		INTeraction
		SMooth nsims(integer 800) //max dim in stata ic
		
		/*Display options*/
		CImethod(string asis) //i=[wald, exact, score], s=[wald, exact, score, t]
		GOF //Goodness of fit
		DOWNload(string asis) 
		DP(integer 2)
		POwer(integer 0)		
		Level(integer 95) 
		LABEL(string) 
		noGRaph //Synonym with nofplot
		noFPlot
		CATPplot
		noSUBgroup 
		noOVerall 
		noITAble
		noWT
		outplot(string) //abs|rr|or|lor|lrr|rd
		SUMTable(string) //none|logit|abs|rr|or|all
		SUMMARYonly
		SUMStat(string)
		STAt(string) //median|mean ; median is default
		FOptions(string asis) /*Options specific to the forest plot*/
		COptions(string asis) /*Options specific to the catterpilar plot*/

		/*passthrough options that go to the forest plot*/
		noOVLine 
		noSTats 
		noBox
		DOUBLE 
		AStext(integer 50) 
		CIOpts(passthru) 
		DIAMopts(passthru) 
		OLineopts(passthru) 
		POINTopts(passthru) 
		BOXopts(passthru) 
		/*PREDciOpt(passthru)*/
		RCols(varlist) 
		/*PREDIction*/  //prediction
		SORtby(varlist) //varlist
		LCols(varlist) 		
		SUBLine
		TEXts(real 1.0) 
		XLAbel(string asis)
		PXlabel(string asis)  //proportions labels
		RXlabel(string asis)  //Ratios labels
		XLIne(passthru)	/*silent option*/	
		XTick(passthru)  
		graphsave(passthru)
		logscale		
		*] ;
	#delimit cr
	
	preserve
	
	marksample touse, strok 
	qui drop if !`touse'
	
	tempvar rid event nonevent total invtotal use id cid neolabel ///
			es se lci uci grptotal uniq mu use rid lpi upi obsid ///
			modeles modellci modeluci holder uniqstudyid clone clones strata numsid  ///
			or orlci oruci rr rrlci rruci rd rdlci rduci lor lorlci loruci lrr lrrlci lrruci abs abslci absuci
			
	tempname nltestRR nltestOR nltestRD mctest samtrix rawest rawesti logodds rrout orout absout logoddsi orouti ///
			rdout rdouti rrouti absouti  exactabsouti exactabsout absexact ///
			coefmat coefvar BVar WVar  omat isq2 bghet bshet lrtestp V dftestnl ptestnl lrtest matgof ///
			outr absoutp absoutpi hetout hetouti popabsout popabsouti poprrout poprdout poplrrout poprdouti poprrouti poplrrouti poporout ///
			poporouti poplorout poplorouti exactorouti  exactlorouti exactorout  exactlorout covmat covmati ///
			neorrout neoorout neoabsout neordout neorawest
	/*Check for mu/cons; its reserved*/
	qui ds
	local vlist = r(varlist)
	foreach v of local vlist {
		if "`v'" == "mu" {
			di in re "mu are a reserved variables name; drop or rename mu"
			exit _rc
		}
	}
	qui {		
		cap gen mu = 1
		cap gen _ESAMPLE = 0
		cap drop _WT
		gen _WT = .		
	}
	//Check link
	if "`link'" == "" {
		local link "logit"
	}
	else {
		if strpos("`link'", "cl") == 1 {
			local link "cloglog"
		}
		else if strpos("`link'", "logl") == 1 {
			local link "loglog"
		}
		else if strpos("`link'", "logi") == 1{
			local link "logit"
		}
		else if "`link'" == "log"{
			local link "log"
		}
		else  {
			di as error "The link `link' not allowed."
				exit
		}
	}
	
	if _by() {
		global by_index_ = _byindex()
		if "`graph'" == "" & "$by_index_" == "1" {
			cap graph drop _all
		}
	}
	else {
		global by_index_ 
	}

	//Mixed or random are synonym
	if 	"`model'" == "" {
		local model "random"
	}
	else {
		*tokenize "`model'", parse(",")
		gettoken model modelopts: model, parse(",")
		gettoken comma modelopts: modelopts //remove the comma
		*local model `1'
		*local modelopts "`3'"
	}
	if strpos("`model'", "f") == 1 {
		local model "fixed"
	}
	else if (strpos("`model'", "r") == 1) | (strpos("`model'", "m") == 1) {
		local model "random"
	}
	else if strpos("`model'", "h") == 1 {
		local model "hexact"
	}
	else if strpos("`model'", "cr") == 1 {
		local model "crbetabin"
	}
	else if strpos("`model'", "cb") == 1 {
		local model "cbbetabin"
	}
	else {
		di as error "Invalid option `model'"
		di as error "Specify either fixed, random, mixed, crbetabin or hexact"
		exit
	}
	
	if "`design'" == "" {
		local design = "general"
	}
	else {
		tokenize "`design'", parse(",")
		local design `1'
		
		//options
		local desopts "`3'"
		while "`desopts'" != "" {
			gettoken option desopts : desopts
			macro shift
			if strpos("`option'", "base") != 0 {
				local baselevel "`option'"
			}
			else if strpos("`option'", "cov") != 0 {
				local cov "`option'"
			}
			else {
				di as error "`option' not allowed in specifying the design()"
				exit
			}
		}
		if "`cov'" != "" {
			cap assert ("`design'" == "comparative" | "`design'" == "mpair")
			if _rc!=0 {
				di as error "The option `cov' only allowed in comparative meta-analysis"
				exit
			}
			if strpos("`cov'", "ind") !=0 {
				local cov "independent"
			}
			else if strpos("`cov'", "unst") !=0 {
				local cov "unstructured"
			}
			else if strpos("`cov'", "int") !=0 {
				if strpos("`cov'", "com") != 0 {
					local cov "commonint"
				}
				if strpos("`cov'", "fre") != 0 {
					local cov "freeint"
				}
			}
			else if strpos("`cov'", "slo") !=0 {
				local cov "commonslope"
			}
			else {
				di as error "`cov' not allowed in specifying the design()"
				exit
			}
		}
		else {
			//Default
			if ("`design'" == "comparative" | "`design'" == "mpair") & strpos("`model'", "random") {
				local cov "independent"
			}
		}
	}
	
	//depracated options
	if 	"`design'" == "paired" {
		di as res "Use of the option -design(paired)- is deprecated and replaced with -design(pcbnetwork)-"
		local design "pcbnetwork"
	}
	
	if "`design'" == "matched"  {
		di as res "Use of the option -design(matched)- is deprecated and replaced with -design(mcbnetwork)-"
		local design "mcbnetwork"
	}
	if "`design'" == "network" {
		di as res "Use of the option -design(network)- is deprecated and replaced with -design(abnetwork)-"
		local design "abnetwork"
	}
	
	if ("`design'" == "mcbnetwork") | ("`design'" == "pcbnetwork") {
		tempvar index byvar assignment idpair ipair
	}
	
	if ("`design'" == "mpair")  {
		tempvar byvar idpair ipair
	}
	
	local fopts `"`options'"'
	
	/*Check if variables exist*/
	foreach var of local varlist {
		cap confirm var `var'
		if _rc!=0  {
			di in re "Variable `var' not in the dataset"
			exit _rc
		}
	}
		
	//Default frequentist inference
	if "`inference'" == "" | strpos("`inference'", "freq") == 1 {
		local inference "frequentist"
	}
	else if strpos("`inference'", "bay") == 1 {
		cap assert _caller() >= 16
		if _rc != 0 {
			di as error "Bayesian modelling needs Stata 16.1 or later versions"
			exit
		}
		local inference "bayesian"
		local model = "bayes" + "`model'"
		
		//Request for the working directory
		if "`bwd'" == "" {
			di as re "{pmore} *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-{p_end}"
			di as re "{pmore}| The bayesian estimation commands saves several datasets  {p_end}" 
			di as re "{pmore}| containing the MCMC samples of the parameters to disk.{p_end}" 
			di as re "{pmore}| Please specify a directory path to save the datasets {p_end}"
			di as re "{pmore}| using the option bwd() and rerun the code.{p_end}"
			di as re "{pmore}| It is your responsibility to erase the datasets {p_end}" 
			di as re "{pmore}| (prefixed metapreg_bayes*) after they are no longer needed.{p_end}"
			di as re "{pmore} *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-{p_end}"
			exit
		}
	}
	else {
		di as error "Invalid inference(`inference') option"
	}
	
	//check frequentist modeloptions
	if "`model'" == "fixed" & strpos("`modelopts'", "ml") != 0 {
		di as error "Option ml not allowed in `modelopts'"
		exit
	}
	if "`model'" == "fixed" & strpos("`modelopts'", "irls") != 0 {
		di as error "Option irls not allowed in `modelopts'"
		exit
	}
	
	//Default bayesian options
	if "`inference'" == "bayesian" {
		metabayesoptscheck, `modelopts'
		local modelopts = r(modelopts)
		local nsims = r(mcmcsize)
		local refsampling = r(refsampling)
	}
	
		//Which population statistic to use on the forest plot/itable
	if "`stat'" != "" {
		if strpos("`stat'", "med") != 0 {
			local stat "Median"
		}
		else if strpos("`stat'", "mea") != 0 {
			local stat "Mean"
		}
		else {
			di as error "Invalid option stat(`stat')"
			di as error "Specify either Median or Mean"
			exit 
		}
	}
	else {
		local stat "Median"  //Default
	}
		
	//Check if crbetabin is installed
	if "`model'" == "crbetabin" {
		capture which betabin
		if _rc != 0 {
			di as res "The user-package betabin is required"
			di  `"{stata "search betabin": Click to search and install the package}"'
			exit
		}
	}
	
	//Avoid Incosistencies & Redundancies
	if "`stratify'" != "" & "`summaryonly'" != "" {
		local wt nowt
	}
		
	qui count
	if `=r(N)' < 2 {
		di as err "Insufficient data to perform meta-analysis"
		exit 
	}
	if `=r(N)' < 3 {
		if "`inference'" == "frequentist" {
			if "`model'" != "hexact" {
				local model hexact //If less than 3 studies, use exact
				di as res _n  "Note: Homo-exact model imposed whenever number of studies is less than 3."
				if "`modelopts'" != "" {
					local modelopts
					di as res _n  "Warning: Model options ignored."
					di as res _n  "Warning: Consider re-specifying options for the fixed-effects model should the model not converge."
				}
			}
		}
		else {
			if "`model'" == "bayesrandom" {
				local model bayesfixed //If less than 3 studies, use exact
				di as res _n  "Note: FE model imposed whenever number of studies is less than 3."
				if "`modelopts'" != "" {
					di as res _n  "Warning: Consider re-specifying options for the fixed-effects model should the model not converge."
				}
			}
		}
	}
	if `level'<1 {
		local level `level'*100
	}
	if `level'>99 | `level'<10 {
		local level 95
	}
	if `astext'>99 | `astext' <1 {
		local astext 50
	}

	//Number of studies in the analysis
	cap assert "`studyid'" != ""
	if _rc != 0 {
		di as err "The study identifier variable is not specified"
		di as err "Specify it with STUDYID(varname) "
		exit _rc
	}
	
	tokenize `varlist'
	if "`design'" == "general" | "`design'" == "comparative" | "`design'" == "abnetwork"  {
		local event = "`1'"
		local total = "`2'"
		
		gen `nonevent' = `total' - `event'
				
		forvalues num = 1/2 {
			cap confirm integer number `num'
			if _rc != 0 {
				di as error "`num' found where integer expected"
				exit
			}
		}
		
		cap assert `total' >= `event' if (`event' ~= .)
		if _rc != 0 {
			di as err "Order should be {n N}. Check your data."
			exit _rc
		}
		
		cap assert `total' > 0 if (`total' ~= .)
		if _rc != 0 {
			di as err "`total' has nonpositive values. Check your data. Ensure all `total' > 0"
			exit _rc
		}
		local depvars "`1' `2'" 
	
		macro shift 2
	}
	else if "`design'" == "mcbnetwork" {
		local a = "`1'"
		local b = "`2'"
		local c = "`3'"
		local d = "`4'"
		cap assert "`6'" != ""
		if _rc != 0 {
			di as err "mcbnetwork data requires atleast 6 variable"
			exit _rc
		}
		local depvars "`1' `2' `3' `4'"
		local Comparator = "`6'"
		local Index = "`5'"
		
		forvalues num = 1/4 {
			cap confirm integer number `num'
			if _rc != 0 {
				di as error "`num' found where integer expected"
				exit
			}
		}
		cap confirm string variable `5'
		if _rc != 0 {
			di as error "The first & second covariate in cbnetwork analysis should be a string"
			exit _rc
		}
		cap confirm string variable `6'
		if _rc != 0 {
			di as error "The first & second covariate in cbnetwork analysis should be a string"
			exit _rc
		}
		macro shift 6
	}
	else if "`design'" == "mpair" {
		local a = "`1'"
		local b = "`2'"
		local c = "`3'"
		local d = "`4'"
		cap assert "`4'" != ""
		if _rc != 0 {
			di as err "mpair data requires atleast 4 variable"
			exit _rc
		}
		local depvars "`1' `2' `3' `4'"
		
		forvalues num = 1/4 {
			cap confirm integer number `num'
			if _rc != 0 {
				di as error "`num' found where integer expected"
				exit
			}
		}
		macro shift 4
	}
	else if "`design'" == "pcbnetwork" {
		local event1 = "`1'"
		local event2 = "`2'"
		local Total = "`3'"
		cap assert "`5'" != ""
		if _rc != 0 {
			di as err "pcbnetwork data requires atleast 5 variable"
			exit _rc
		}
		local depvars "`1' `2' `3'"
		local Comparator = "`5'"
		local Index = "`4'"
		
		cap assert ((`Total' >= `event1') & (`Total' >= `event2')) if ((`event1' ~= .) & (`event2' ~= .))
		if _rc != 0 {
			di as err "Order should be {n1 n2 N}. Check your data."
			exit _rc
		}
		
		forvalues num = 1/3 {
			cap confirm integer number `num'
			if _rc != 0 {
				di as error "`num' found where integer expected"
				exit
			}
		}
		cap confirm string variable `4'
		if _rc != 0 {
			di as error "The first & second covariate in pcbnetwork analysis should be a string"
			exit _rc
		}
		cap confirm string variable `5'
		if _rc != 0 {
			di as error "The first & second covariate in pcbnetwork analysis should be a string"
			exit _rc
		}
		macro shift 5
	}
	
	local regressors "`*'"
	local p: word count `regressors'
	
	/*
	if "`model'" == "hexact" {
		cap assert `p' == 0	
		if _rc != 0 {
			di as error "Covariates not allowed in the `model' model. Specify model(fixed) or model(mixed)"
			exit _rc
		}
	}
	*/
	
	if "`design'" == "comparative" | "`design'" == "abnetwork" {
		*if "`model'" != "fixed" {
			cap assert `p' > 0
			if _rc != 0 {
				di as error "`design' analysis requires at least 1 covariate to be specified"
				exit _rc
			}
		/*}
		else {
			cap assert `p' > 1
			if _rc != 0 {
				di as error "`design' analysis requires at least 2 covariate to be specified"
				exit _rc
			}
		}*/
	}
	
	if "`design'" == "comparative" | "`design'" == "abnetwork"  {
		gettoken first confounders : regressors
		if "`first'" != "" {
			/*
			cap confirm string variable `first'
			if _rc != 0 {
				di as error "The first covariate in `design' analysis should be a string"
				exit _rc
			}
			*/
			
			cap confirm string variable `first'
		
			if _rc != 0 {
				cap label list `first' //see if labelled
				
				if _rc != 0 {
					di as error "The first covariate in `design' analysis should be a string"
					exit, _rc
				}
			}
		}
		if ("`by'" != "") & strpos("`outplot'", "abs" ) != 0 {
			cap assert ("`first'" != "`by'") 
			if _rc != 0 { 
					di as error "Remove the option by(`by') or specify a different by-variable"
					exit _rc
			}
		}
	}
	
	if "`outplot'" == "" {
		if ("`design'" == "mcbnetwork" | "`design'" == "pcbnetwork" | "`design'" == "mpair")  {
			local outplot = "rr"  //default
		}
		else {
			local outplot = "abs"
		}
	}
	else {
		if strpos("`outplot'", "rr") != 0 | strpos("`outplot'", "or") != 0 | strpos("`outplot'", "rd") != 0 {
			cap assert "`design'" != "general"
			if _rc != 0 {
				di as error "Option outplot(`outplot') only avaialable for comparative/mcbnetwork/pcbnetwork/abnetwork designs with first covariate as string"
				di as error "Specify the first string covariate and the appropriate design(comparative/mcbnetwork/pcbnetwork/abnetwork)"
				exit _rc
			}
		}
	}
	
	//
	if ("`model'" == "hexact") & strpos("`outplot'",  "or") == 0 & "`design'" == "comparative" { 
		local outplot = "or"  //or default
	}
		
	//smooth is redundant in betabin, we cannot recover the individual estimates
	if (("`model'" == "hexact"  & "`design'" == "comparative") | "`model'" == "crbetabin" ) & "`smooth'" != "" {
		local smooth
		di as res "The option -smooth- is ignored. The model-based study estimates are irrecoverable."
	}
	
	foreach metric of local outplot {	
		cap assert ("`metric'" == "lrr") | ("`metric'" == "lor") |("`metric'" == "rr") | ("`metric'" == "or") | ("`metric'" == "rd") | ("`metric'" == "abs") 
		if _rc != 0  {
			di as error "`metric' not allowed in outplot(`outplot')"
			exit _rc
		}
	}
	
	//Empty 
	local modeles
	local modellci
	local modeluci
	
	
	//generate variables for model-based estimates
	local metrics "p rr rd lrr or lor"
	
	foreach metric of local metrics {
		tempvar model`metric'  model`metric'lci  model`metric'uci 
		qui {
			cap gen `model`metric'' = .
			cap gen `model`metric'lci' = .
			cap gen `model`metric'uci' = .	
		}
			
		local modeles "`modeles' `model`metric''"
		local modellci "`modellci' `model`metric'lci'"
		local modeluci "`modeluci' `model`metric'uci'"
	}
		
	/*	
	if "`prediction'" != "" & "`outplot'" != "abs" {
		local prediction 
		di as res "NOTE: Predictions only computed for absolute measures. The option _prediction_ will be ignored"
	}*/
		
	//check no underscore in the variable names
	if strpos("`regressors'", "_") != 0  {
		di as error "Underscore is a reserved character and covariate(s) containing underscore(s) is(are) not allowed"
		di as error "Rename the covariate(s) to remove the underscore(s) character(s)"
		exit	
	}
	
	if `p' < 2 & "`interaction'" !="" & ("`design'" != "`mcbnetwork'" | "`design'" != "`pcbnetwork'" ) {
		di as error "Interactions allowed with atleast 2 covariates"
		exit
	}
	
	//=======================================================================================================================
	tempfile master
	qui save "`master'"
	
	if strpos("`model'", "bayes") == 1 {
		tempfile metapregbayesreps
	}
		
	*declare study labels for display
	if "`label'"!="" {
		tokenize "`label'", parse("=,")
		while "`1'"!="" {
			cap confirm var `3'
			if _rc!=0  {
				di as err "Variable `3' not defined"
				exit
			}
			local `1' "`3'"
			mac shift 4
		}
	}	
	qui {
		*put name/year variables into appropriate variable
		if "`namevar'"!="" {
			local lbnvl : value label `namevar'
			if "`lbnvl'"!=""  {
				quietly decode `namevar', gen(`neolabel')
			}
			else {
				gen str10 `neolabel'=""
				cap confirm string variable `namevar'
				if _rc==0 {
					replace `neolabel'=`namevar'
				}
				else if _rc==7 {
					replace `neolabel'=string(`namevar')
				}
			}
		}
		if "`namevar'"==""  {
			cap confirm numeric variable `studyid'
			if _rc != 0 {
				gen `neolabel' = `studyid'
			}
			if _rc == 0{
				gen `neolabel' = string(`studyid')
			}
		}
		if "`yearvar'"!="" {
			local yearvar "`yearvar'"
			cap confirm string variable `yearvar'
			if _rc==7 {
				local str "string"
			}
			if "`namevar'"=="" {
				replace `neolabel'=`str'(`yearvar')
			}
			else {
				replace `neolabel'=`neolabel'+" ("+`str'(`yearvar')+")"
			}
		}
	}
	if "`design'" == "mcbnetwork" | "`design'" =="pcbnetwork" | "`design'" == "mpair" {
		longsetup `varlist', rid(`rid') assignment(`assignment') event(`event') total(`total') idpair(`idpair') `design'
		if "`design'" == "mpair" {
			qui gen `ipair' = "Response1"
			qui replace `ipair' = "Response2" if `idpair'
		}
		else {
			qui gen `ipair' = "Yes"
			qui replace `ipair' = "No" if `idpair'
		}

		qui gen `nonevent' = `total' - `event'
	}
	else {
		qui gen `rid' = _n
	}
	
	//panelize data
	if "`model'" == "cbbetabin" {
		tempvar count
		qui drop `rid' mu
		longsetup `event' `nonevent', rid(`rid') idpair(mu) panelize event(`count')
		qui replace `event' = `count'
		qui drop `count'
	}
	
	//byvar
	/*
	if "`by'" != "" {		
		cap confirm string variable `by'
		if _rc != 0 {
			di as error "The by() variable should be a string"
			exit _rc
		}
		if strpos(`"`varlist'"', "`by'") == 0 {
			tempvar byvar
			my_ncod `byvar', oldvar(`by')
			qui drop `by'
			rename `byvar' `by'
		}
	}
	*/
	if "`by'" != "" {
		cap confirm string variable `by'
		
		if _rc != 0 {
			cap label list `by' //see if labelled
			
			if _rc != 0 {
				di as error "The by() variable should be a string or coded integer"
				exit, _rc
			}
		}
		else {
			//Check if the byvariable is part of varlist
			local found 0
			foreach v of local varlist {
				if "`v'" == "`by'" {
					local found = 1
					continue, break
				}
			}
			
			if !`found' {
			/*if strpos(`"`varlist'"', "`by'") == 0*/ 
				tempvar byvar
				my_ncod `byvar', oldvar(`by')
				drop `by'
				rename `byvar' `by'
			}
		}
	}
	
	buildregexpr `varlist', `interaction' `alphasort' `design' ipair(`ipair') `baselevel'  studyid(`studyid') model(`model')
	
	local catreg = r(catreg)
	local contreg = r(contreg)
	local basecode = r(basecode)
	if "`model'" == "cbbetabin" {
		local regexpression2 = r(regexpression)  
		local regexpression = r(regexpression2) //for the conventional logistic
	}
	else {
		local regexpression = r(regexpression)
	}
	
	if "`interaction'" != "" { 
		local varx = r(varx)
		local typevarx = r(typevarx)		
	}
	if "`design'" == "comparative" {
		*local varx : word 1 of `regressors'
		gettoken varx catreg : catreg
		local typevarx = "i"
		/*		
		forvalues l=1(1)2 {	
			local label`l':label `varx' `l'	
		}
		replace `varx' = `varx' - 1
		label define `varx' 0 "`label1'" 1 "`label2'", replace
		label values `varx' `varx'
		
		local basecode = `=`basecode'-1'
		local baselab:label `varx' `basecode'
		*/
		if `basecode' == 1 {
			local indexcode "2"
		}
		else {
			local indexcode "1"
		}
		local indexlab:label `varx' `indexcode'
		
		if strpos("`outplot'", "r") != 0 {
			local varxlabs "`varx' `indexlab' `baselab'"
		}
	}
		 
	if "`design'" == "pcbnetwork" | "`design'" == "mcbnetwork" | "`design'" == "mpair" { 
		*local varx = "`ipair'"
		if "`inference'" == "bayesian" & "`design'" == "mcbnetwork"  {
			local varx = "`Comparator'"
		}
		else {
			local varx = "`ipair'"	
			*local varx = "`Comparator'"
		}
		local typevarx = "i"		
	}
	
	local pcont: word count `contreg'
	if "`typevarx'" != "" & "`typevarx'" == "c" {
		local ++pcont
	} 
	if `pcont' > 0 {
		local continuous = "continuous"
	}
	
	/*Model presenations*/
	if ("`design'" == "general" | "`design'" == "comparative" ) {
		local nu = "mu"
	}
	else if "`design'" == "pcbnetwork" | "`design'" == "mcbnetwork" {
		if "`interaction'" != "" {
			local nu = "Ipair*`Comparator' + `Index'"
		}
		else {
			local nu = "mu + Ipair + `Index'"
		}			
	}
	else if "`design'" == "mpair" {
		local nu = "mu + Ipair"
	}
	else if "`design'" == "abnetwork" {
		local nu = "mu.`first'"
	}
	if "`model'" == "cbbetabin" {
		local nu = "mu + b0"
	}
	
	local VarX: word 1 of `regressors'
	forvalues i=1/`p' {
		local c:word `i' of `regressors'
		local nu = "`nu' + `c'"		
		
		if "`interaction'" != "" & `i' > 1 {
				local nu = "`nu' + `c'*`VarX'"			
		}
	}
	
	if ("`catreg'" != " " | "`typevarx'" =="i" | ("`design'" == "comparative" | "`design'" == "mcbnetwork" | "`design'" == "pcbnetwork"))  {

		if "`design'" == "mcbnetwork" | "`design'" == "pcbnetwork" {
			local catregs = "`catreg' `Comparator' `Index'"
		}

		if "`design'" == "comparative" {
			local catregs = "`catreg' `varx'" 
		}
		if "`design'" == "abnetwork" {
			tokenize `catreg'
			macro shift
			local catregs "`*'"
		}
		if "`design'" == "general" | "`design'" == "mpair"  {
			local catregs "`catreg'"
		}
	}
	
	/*//Specify the grouping variable for proportions and/or ratios/difference
	if "`subgroup'" == "" & ("`catreg'" != "" | "`typevarx'" =="i" ) {
		if strpos("`outplot'", "abs") != 0 {
			if "`typevarx'" =="i" {
				//comparative
				local groupvar = "`varx'"
			}
			else {
				//other designs
				local groupvar : word 1 of `catreg'
			}
		}
		//Comparative
		if strpos("`outplot'", "r") != 0 & "`varx'" != ""  & "`catreg'" != "" {
			local groupvar : word 1 of `catreg'
		}
	}
	*/
	if "`by'" != "" {
		*local groupvar "`by'"
		local byvar "`by'"
		*How many times to loop
		*qui label list `by'
		*local nlevels = r(max)
		
		qui levelsof `by', local(codelevels)
		local nlevels = r(r)
					
	}
	
	/*if "`design'" == "abnetwork" {	
		//How many groups
		qui label list `first'
		local ngroups = r(max)
		
		if "`aliasdesign'" != "" {
			cap assert `ngroups' == 2
			cap assert "`aliasdesign'" == "comparative"
		}
		else {
			local groupvar "`first'"
			local overall "nooverall"
			
			if strpos("`outplot'", "r") != 0 {
				local itable "noitable"
			}
		}		
	} */
	
	*Stratify not allow in pcbnetwork, mcbnetwork or abnetwork analysis
	if "`stratify'" != "" {
		if ("`design'" == "pcbnetwork") | ("`design'" == "mcbnetwork") | ("`design'" == "abnetwork") {
			di as res "NOTE: The option stratify is ignored in `design' analysis"
			local stratify
		}
	}
	
	*Check by is active & that there are more than 1 levels
	if "`stratify'" != "" {
		if "`by'" == "" {
			di as error "The by() variable needs to be specified in stratified analysis"
			exit			
		}
		else {
			if `nlevels' < 2 {
				di as error "The by() variable should have atleast 2 categories in stratified analysis"	
				exit
			}
		}
	}	
	qui gen `use' = .
	
	//Replace population-averaged estimates with Conditional/exact estimates if model has issues e.g complete seperation etc
	if "`stratify'" != "" & `p' < 1  {
		local enhance "enhance"		
	}
	
	*Loop should begin here
	if "`stratify'" == "" {
		local nlevels = 0
	}
	local i = 1
	local byrownames 
	local bybirownames
	
	
	if "`design'" == "abnetwork"  {
		local hetdim 7
	}
	else {
		if (`p' == 0) & ("`model'" == "random") &  ("`design'" != "pcbnetwork" | "`design'" != "mcbnetwork" )  {
			local hetdim 5
		}
		else {
			local hetdim 7
		}
	}
	
	cap confirm string variable `studyid'
	if _rc == 0 {
		my_ncod `numsid', oldvar(`studyid')
		drop `studyid'
		rename `numsid' `studyid'
	}
	
	if strpos("`outplot'", "r") != 0 & "`design'" == "abnetwork" & "`aliasdesign'" == "" {
		local smooth
	}
	
	//Should run atleast once
	while `i' < `=`nlevels' + 2' {
		local modeli = "`model'"
		local modeloptsi = "`modelopts'"
		local smoothi = "`smooth'"
		local getmodel
		local optimizedi = 0
		local computewti = "computewt"
		local ilab
	
		//don't run last loop if stratify
		if (`i' > `nlevels') & ("`stratify'" != "") & ("`design'" == "comparative" | "`design'" == "mpair") {
			local overall "nooverall"
			continue, break
		}
		
		*Stratify except the last loop for the overall
		if (`i' < `=`nlevels' + 1') & ("`stratify'" != "") {
			local icode : word `i' of  `codelevels'
			local strataif `"if `by' == `icode'"'
			local ilab:label `by' `icode'
			local stratalab `":`by' = `ilab'"'
			local ilab = ustrregexra("`ilab'", " ", "_")
			local byrownames = "`byrownames' `by':`ilab'"
			local byrowname = "`by'|`ilab'"
			if ("`design'" == "comparative" | "`design'" == "mpair") & "`stratify'" != "" {
				local bybirownames = "`bybirownames' `ilab':`baselab' `ilab':`ilab'"
			}

			*Check if there is enough data in each strata
			//Number of obs in the analysis
			qui egen `obsid' = group(`rid') if `by' == `icode'
			qui summ `obsid'
			local Nobs= r(max)
			drop `obsid'	

			//Number of studies in the analysis
			qui egen `uniq' = group(`studyid') if `by' == `icode'
			qui summ `uniq'
			local Nuniq = r(max)
			drop `uniq'	
		}
		else {
			//Skip if overall not needed
			if "`overall'" != "" & (`i' > `=`nlevels'+1')  & ("`stratify'" != "" | "`design'" == "comparative" | "`design'" == "mpair" )   {
				continue, break
			}
			//Don't 1.smoothen after last loop if stratify 2.compute weights
			if (`i' > `nlevels') & ("`stratify'" != "") {
				local smoothi
				local computewti
			}
			
			//Nullify
			local strataif 
			local stratalab ": all studies"
			if "`stratify'" != "" {
				local byrownames = "`byrownames' Overall"
				local byrowname = "All_studies"				
			}
			
			//Number of obs in the analysis
			qui count
			local Nobs= r(N)
			if "`design'" == "mcbnetwork" | "`design'" == "pcbnetwork" {
				local Nobs = `Nobs'*0.5
			}
			qui egen `uniq' = group(`studyid')
			qui summ `uniq'
			local Nuniq = r(max)
			drop `uniq'
		}
		
		if "`model'" == "cbbetabin" {
			local Nobs = `Nobs'*0.5
		}
		
		if "`design'" == "comparative" {
			cap assert mod(`Nobs', 2) == 0 
			if _rc != 0 {
				di as error "Comparative analysis requires 2 observations per study"
				exit _rc
			}
		}
		if "`design'" == "abnetwork" {
			cap assert `Nobs'/`Nuniq' >= 2 
			if _rc != 0 {
				di as error "abnetwork design requires atleast 2 observations per study"
				exit _rc
			}
		}
		/*
		if ("`model'" == "cbbetabin" ) | ("`design'" == "general") {
			cap assert `Nobs'/`Nuniq' == 1 
			if _rc != 0 {
				di as error "Ensure no duplicate identifiers in studyid(`studyid')"
				exit _rc
			}
		}
		*/
		
		*if (`Nobs' < 3 & "`modeli'" != "hexact" & "`design'" != "comparative") | ((`Nobs' < 5 ) & ("`modeli'" == "random") & ("`design'" == "comparative")) {
		if `Nuniq' < 3 {
			if "`inference'" == "frequentist" {
				if  "`modeli'" == "random"  {
					local modeli fixed //If less than 3 studies, use exact model
					if "`modeloptsi'" != "" {
						local modeloptsi
						noi di as res _n  "Warning: `model'-effects model options ignored."
						noi di as res _n  "Warning: Homo-exact model fitted instead."
					}
				}
			}
		}
		
		//scimethod
		/*
		if "`cimethod'" != "" & (strpos("`cimethod'", ",") != 1) { 
			tokenize "`cimethod'", parse(",")
		
			if "`3'" != "" {
				local scimethod = strltrim("`3'")
			}
		}
		if "`cimethod'" != "" & (strpos("`cimethod'", ",") == 1) { 
			tokenize "`cimethod'", parse(",")		
			if "`2'" != "" {
				local scimethod = strltrim("`2'")
			}
		}
		*/
		if "`cimethod'" != "" {
			tokenize "`cimethod'", parse(",")
			if "`1'" == "," {
				local scimethod = strltrim("`2'")
			}
			
			if "`3'" != "" {
				local scimethod = strltrim("`3'")
			}
		}
		
		if "`model'" == "hexact" {
			local scimethod "exact"
		}
		else {
			//Summary
			if "`inference'" == "frequentist" {
				if "`scimethod'" != "" & ((strpos("`scimethod'", "t") != 1) &  (strpos("`scimethod'", "w") != 1) & (strpos("`scimethod'", "z") != 1) ){
					di as error "Option `scimethod' not allowed in cimethod(`cimethod')"
					exit	
				}
				else if "`scimethod'" == "z" | "`scimethod'" == "wald" {
					local scimethod "wald"
				}
				else {
					//Default
					local scimethod "t"
				}
			}
			else {
				//Bayesian credible intervals; default is equal-tailed interval
				if ("`scimethod'" != "") & ("`scimethod'" != "hpd") & ("`scimethod'" != "eti") {
					di as error "Option `scimethod' not allowed in cimethod(`cimethod')"
					exit
				}
				else if "`scimethod'" == "hpd" {
					local scimethod 
				}
				else {
					//Default
					local scimethod "eti"
				}
			}
		}
				
		if "`stratify'" != "" {
			di as res _n "*********************************** Model for `stratalab' ***************************************" 
		}
		else {
			di as res _n "**************************************************************************" 
		}
		
		*Run model if more than 1 study
		if (`Nobs' > 1) {
			if "`inference'" == "bayesian" {
				local bayesreps = "`bwd'\metapreg_bayesreps"
				local bayesest = "`bwd'\metapreg_bayesest"
			}
			
			preg `event' `nonevent' `total' `strataif', rid(`rid') sid(`studyid') studyid(`studyid') use(`use') regexpression(`regexpression') regexpression2(`regexpression2') nu(`nu')  ///
				regressors(`regressors')  catreg(`catreg') contreg(`contreg') level(`level') varx(`varx') typevarx(`typevarx')  /// 
				`progress' model(`modeli') modelopts(`modeloptsi') `mc' `interaction' `design' aliasdesign(`aliasdesign') by(`by') `stratify' baselevel(`basecode') ///
				comparator(`Comparator') scimethod(`scimethod') `gof' nsims(`nsims') link(`link') bayesrepsfilename(`metapregbayesreps') ///
				modeles(`modeles')  modellci(`modellci') modeluci(`modeluci') ///
				outplot(`outplot') `smoothi' cov(`cov') `computewti' ///
				inference(`inference') refsampling(`refsampling') stat(`stat') bayesreps(`bayesreps') bayesest(`bayesest')
				
			mat `rawesti' = r(rawest)
			if "`getmodeli'" != "crbetabin" {
				mat `popabsouti' = r(popabsout)
			}
			mat `exactabsouti' = r(exactabsout)

			local mdf = r(mdf) //mdf = 0 if saturated
			local getmodeli = r(model) //Returned model
			local rrsuccess = r(rrsuccess)
			
			if ("`catreg'" != " " | "`typevarx'" == "i") & `rrsuccess' {
				if "`getmodeli'" == "hexact" {
					mat `exactorouti' = r(exactorout)
					mat `exactlorouti' = r(exactlorout)	
				}
				else {
					mat `rrouti' = r(rrout)
					mat `rdouti' = r(rdout)
					mat `orouti' = r(orout)
					
					if "`getmodeli'" != "crbetabin" {
						mat `poprrouti' = r(poprrout)
						mat `poprdouti' = r(poprdout)
						mat `poplrrouti' = r(poplrrout)
						mat `poporouti' = r(poporout)
						mat `poplorouti' = r(poplorout)
					}
					
					local inltest = r(inltest)
					if "`inltest'" == "yes" & "`stratify'" == "" {
						mat `nltestRR' = r(nltestRR) 
						mat `nltestRD' = r(nltestRD)
						mat `nltestOR' = r(nltestOR) 
					}
				}
			}
			else {
				local ratios "noratios"
			}
						
			/*
			if (`p' > 0) & ("`mc'" =="") { 
				mat `mctest' = r(mctest) 
			}*/
			
			mat `absouti' = r(absout)
			mat `absoutpi' = r(absoutp)
			if strpos("`modeli'", "random") !=0 | strpos("`model'", "betabin") != 0 { 
				mat `hetouti' = r(hetout)
				mat `covmati' = r(covmat)
			}
			else {
				mat `hetouti' = J(1, `hetdim', .)
				mat `covmati' = J(1, 3, .)
			}			
		}
		*if 1 study or exact inference
		else {
			mat `rawesti' = J(1, 9, .)
			mat `popabsouti' = J(1, 6, .)
			mat `exactabsouti' = J(1, 11, .)
			mat `absouti' = J(1, 9, .)
			*mat `absoutpi' = J(1, 2, .)
			mat `hetouti' = J(1, `hetdim', .)
			mat `covmati' = J(1, 3, .)			
			mat `rrouti' = J(1, 6, 1)
			mat `rdouti' = J(1, 6, 0)
			mat `orouti' = J(1, 6, 1)
			
			mat `poprrouti' = J(1, 6, 1)
			mat `poprdouti' = J(1, 6, 0)
			mat `poplrrouti' = J(1, 6, 1)
			mat `poporouti' = J(1, 6, 1)
			mat `poplorouti' = J(1, 6, 1)

			mat `exactorouti' = J(1, 4, .)
			mat `exactlorouti' = J(1, 4, .)
			
			mat rownames `rawesti' = Overall
			mat rownames `popabsouti' = Overall
			mat rownames `exactabsouti' = Overall
			mat rownames `absouti' = Overall
			*mat rownames `absoutpi' = Overall
			mat rownames `hetouti' = Overall
			mat rownames `covmati' = Overall
			mat rownames `rrouti' = Overall
			mat rownames `poprrouti' = Overall
			mat rownames `rdouti' = Overall
			mat rownames `poprdouti' = Overall
			mat rownames `poplrrouti' = Overall
			mat rownames `orouti' = Overall
			mat rownames `poporouti' = Overall
			mat rownames `poplorouti' = Overall
			mat rownames `exactlorouti' = Overall
			mat rownames `exactorouti' = Overall
			
			qui replace `use' = 1 `strataif'
			local getmodeli = "none" //Returned model
		}
		
		if ("`stratify'" != "") {
			mat rownames `hetouti' = `byrowname'
			mat roweq `absouti' = `byrowname'
			mat roweq `popabsouti' = `byrowname'
			mat roweq `exactabsouti' = `byrowname'
			*mat roweq `absoutpi' = `byrowname'
			mat roweq `rawesti' = `byrowname'
			mat roweq `covmati' = `byrowname'
			
			if "`ratios'" == "" {
				if "`model'" == "hexact" {
					mat roweq `exactorouti' = `byrowname'
					mat roweq `exactlorouti' = `byrowname'
				}
				else {
					mat roweq `rrouti' = `byrowname'
					mat roweq `poprrouti' = `byrowname'
					mat roweq `rdouti' = `byrowname'
					mat roweq `poprdouti' = `byrowname'
					mat roweq `poplrrouti' = `byrowname'
					mat roweq `orouti' = `byrowname'
					mat roweq `poporouti' = `byrowname'
					mat roweq `poplorouti' = `byrowname'
				}
			}
		}
		*noi mat list `hetouti'
		*Stack the matrices
		if `i' == 1 {
			mat `absout' =	`absouti'
			if "`ratios'" == "" {
				if "`model'" != "hexact" {
					mat `rdout' =	`rdouti'
					mat `rrout' =	`rrouti'
					mat `orout' =	`orouti'
					
					if "`model'" != "crbetabin" {
						mat `poprdout' = `poprdouti'
						mat `poprrout' = `poprrouti'
						mat `poplrrout' = `poplrrouti'
						mat `poporout' = `poporouti'
						mat `poplorout' = `poplorouti'
					}
				}
				else {
					mat `exactorout' = `exactorouti'
					mat `exactlorout' = `exactlorouti'
				}
			}			
			mat `rawest' = `rawesti'
			mat `popabsout' = `popabsouti'
			mat `exactabsout' = `exactabsouti'
			mat `absoutp' = `absoutpi'
			mat `hetout' = `hetouti'
			mat `covmat' = `covmati'
		}
		else {
			mat `absout' = `absout' \ `absouti'
			mat `popabsout' = `popabsout' \ `popabsouti'
			mat `exactabsout' = `exactabsout' \ `exactabsouti'
			
			if "`ratios'" == "" {
				if "`model'" != "hexact" {
					mat `rdout' = `rdout' \ `rdouti'
					mat `rrout' = `rrout' \ `rrouti'
					mat `orout' = `orout' \ `orouti'
					
					if "`model'" != "crbetabin" {
						mat `poprdout' = `poprdout' \ `poprdouti'
						mat `poprrout' = `poprrout' \ `poprrouti'
						mat `poplrrout' = `poplrrout' \ `poplrrouti'
						mat `poporout' = `poporout' \ `poporouti'
						mat `poplorout' = `poplorout' \ `poplorouti'
					}
				}
				else {
					mat `exactorout' = `exactorout' \ `exactorouti'
					mat `exactlorout' = `exactlorout' \ `exactlorouti'
				}
			}
			mat `rawest' = `rawest' \ `rawesti'
			*mat `absoutp' = `absoutp' \ `absoutpi'
			mat `hetout' = `hetout' \ `hetouti'
			mat `covmat' = `covmat' \ `covmati'
		}
				
		tokenize `depvars'
		if "`design'" == "general" | "`design'" == "abnetwork" | "`design'" == "comparative" {
			if strpos("`model'", "betabin")!= 0 {
				di "{phang} `1' ~ beta-binomial(alpha, beta, `2') {p_end}"
				di "{phang}E(p) = alpha/(alpha + beta) {p_end}"
				di "{phang}phi = 1/(alpha*beta) {p_end}"				
			}
			else {
				di "{phang} `1' ~ binomial(p, `2'){p_end}"
			}
		}
		else if ( "`design'" == "mcbnetwork" ) {
			di "{phang} `1' + `2'  ~ binomial(p, `1' + `2' + `3' + `4'){p_end}"
			di "{phang} `1' + `3' ~ binomial(p, `1' + `2' + `3' + `4'){p_end}"
		}
		else if ( "`design'" == "mpair" ) {
			di "{phang} `1' + `2'  ~ binomial(p, `1' + `2' + `3' + `4')  << Response1{p_end}"
			di "{phang} `1' + `3' ~ binomial(p, `1' + `2' + `3' + `4') << Response2{p_end}"
		}
		else if "`design'" == "pcbnetwork" {
			di "{phang} `1' ~ binomial(p, `3'){p_end}"
			di "{phang} `2' ~ binomial(p, `3'){p_end}"
		}
		
		if strpos("`getmodeli'", "random") != 0 {
			if "`cov'" == "" {
				di "{phang} `link'(p) = `nu' + `studyid'{p_end}"	
			}
			else {
				if "`cov'" != "commonslope" {
					if ("`design'" == "mpair")  {
						di "{phang} `link'(p) = `nu' + Ipair.`studyid' + `studyid'{p_end}"
					}
					else {
						di "{phang} `link'(p) = `nu' + `first'.`studyid' + `studyid'{p_end}"
					}
				}
				else {
					di "{phang} `link'(p) = `nu' + `studyid'{p_end}"
				}	
			}
			
			if ("`cov'" == "commonslope" & ("`design'" == "comparative"  | "`design'" == "mpair" )) | ("`cov'" == "" & "`design'" != "comparative" )  {
				di "{phang}`studyid' ~ N(0, tau2){p_end}"
			}
			
			if (("`cov'" == "commonint" | "`cov'" == "freeint")& "`design'" == "comparative" ) {
				di "{phang}`first'.`studyid' ~ N(0, sigma2){p_end}"
			}
			
			if "`cov'" == "independent" {
				di "{phang}`studyid' ~ N(0, tau2){p_end}"
				if ("`design'" == "mpair") {
					di "{phang}Ipair.`studyid' ~ N(0, sigma2){p_end}"
				}
				else {
					di "{phang}`first'.`studyid' ~ N(0, sigma2){p_end}"
				}
			}
			
			if "`cov'" == "unstructured" {
				if ("`design'" == "mpair") {
					di "{phang}`studyid',Ipair.`studyid'  ~ biv.normal(0, Sigma){p_end}"
				}
				else {
					di "{phang}`studyid',`first'.`studyid'  ~ biv.normal(0, Sigma){p_end}"
				}
				
				di "{p 20}  Sigma = {c |}tau2, {space 1} rho{c |}{p_end}"
				di "{p 28}   {c |}rho, sigma2{c |}{p_end}"
			}
		}
		else if "`model'" == "crbetabin" {
			di "{phang} `link'(E(p)) = `nu'{p_end}"	
		}
		else if "`model'" == "cbbetabin" {
			di _n"{phang}Model fitted via conditional FE negative binomial regression where{p_end}"
			di "{phang}alpha = exp(`nu') {p_end}"
			di "{phang}beta = exp(b0) {p_end}"
		}
		else /*if "`getmodeli'" == "fixed" | "`getmodeli'" == "betabin"*/ {
			if "`getmodeli'" != "hexact" {
				di "{phang} `link'(p) = `nu'{p_end}"
			}			
		}
		if "`design'" == "pcbnetwork" | "`design'" == "mcbnetwork"  {
			di "{phang} Ipair = 0 if 1st pair{p_end}"
			di "{phang} Ipair = 1 if 2nd pair{p_end}"
		}
		if "`design'" == "mpair" {
			di "{phang} Ipair = 0 if Response1{p_end}"
			di "{phang} Ipair = 1 if Response2{p_end}"
		}
		if "`design'" == "abnetwork" {
			di "{phang}`first' ~ N(0, sigma2){p_end}"
			*qui label list `first'
			*local nfirst = r(max)
			
			qui levelsof `first', local(codelevels)
			local nfirst = r(r)
		}		
		if ("`catreg'" != " " | "`typevarx'" =="i" | ("`design'" == "comparative" | "`design'" == "mcbnetwork" | "`design'" == "pcbnetwork"))  {
			di _n "{phang}Base levels{p_end}"
			di _n as txt "{pmore} Variable  -- Base Level{p_end}"
		}
		foreach fv of local catregs  {
			/*if "`fv'" == "`varx'" & "`design'" == "comparative" {
				local lab:label `fv' 0
			}
			else {*/
				local lab:label `fv' 1  //look for the minimum
			*}
			if "`fv'" != "`studyid'" {
				di "{pmore} `fv'  -- `lab'{p_end}"
			}			
		}
		if "`design'" == "abnetwork" {
			local lab:label `first' `basecode'
			di "{pmore} `first'  -- `lab'{p_end}"
		}
		if "`design'" == "mpair" {
			local lab:label `varx' `basecode'
			di "{pmore} Ipair  -- `lab'{p_end}"
		}
			
		di _n
		di "{phang}" as txt "Number of observations = " as res "`Nobs'{p_end}"
		di "{phang}" as txt "Number of studies = " as res "`Nuniq'{p_end}"
		if "`design'" == "abnetwork" {
			di "{phang}" as txt "Number of `first's = " as res "`nfirst'{p_end}"
		}
				
		if `Nobs' > 1 & "`model'" != "hexact" {
		
			//Get GOF
			qui estimates restore metapreg_modest
			
			if "`inference'" == "frequentist" {
				qui estat ic
				mat `matgof' = r(S)
				local BIC =  `matgof'[1, 6]
				mat `matgof' = `matgof'[1..., 5..6]
				local widthc = 8

			}
			else {
				qui bayesstats ic
				mat `matgof' = r(ic)
				mat `matgof' = `matgof'[1..., 2..3]
				local DIC = `matgof'[1,1]
				mat rownames `matgof' = Value
				local widthc = 15
			}
			
			mat rownames `matgof' = Value
			//Display GOF	
			if ("`gof'" != "") {	
				di _n
				di as text "Goodness of Fit Criterion"
				#delimit ;
				noi matlist `matgof',  
							cspec(& %7s |   %8.`=`dp''f &  %`=`widthc''.`=`dp''f o2&) 
							rspec(&-&) underscore  nodotz
				;
				#delimit cr 
			}
			
			//Get the command, necessary for frequentist models 
			if "`inference'" != "bayesian" {
				qui capture estimates restore metapreg_modest
				if _rc == 0 {
					local command0 = e(cmd) 
				}
			}
			
			if "`stratify'" != "" {
				*qui estimates restore metapreg_modest
				di  _n	

				if (`i' < `=`nlevels' + 1') {					
					local ilab = ustrregexra("`ilab'", " ", "_")
					local ilab = ustrregexra("`ilab'", "-", "_")
					if strpos("`ilab'", "+") != 0 {
						local ilab = ustrregexra("`ilab'", "+", "")
					}
					if strpos("`ilab'", "/") != 0 {
						local ilab = ustrregexra("`ilab'", "/", "_")
					}
					local ilab = "`ilab'" + "$by_index_"
					
					//abbreviate if NECESSARY
					local ilablen = strlen("`ilab'")
					if `ilablen' > 20 {
						local ilab = abbrev("`ilab'", 15)
						if strpos("`ilab'", "~") != 0 {
							local ilab = ustrregexra("`ilab'", "~", "")
						}
					}
				}
				else {
					local ilab ="All_studies" + "$by_index_"
				}
				
				qui estimates store metapreg_`ilab'
				display `"{stata "estimates replay metapreg_`ilab'":Click to show the raw estimates}"'
			}
			else {
				if "$by_index_" !="" {
					di  _n	
					*qui estimates restore metapreg_modest
					local ilab = "$by_index_"
					qui estimates store metapreg_`ilab'
					display `"{stata "estimates replay metapreg_`ilab'":Click to show the raw estimates}"'
				}
				else {
					di  _n					
					display `"{stata "estimates replay metapreg_modest":Click to show the raw estimates}"'
				}
			}
			
			if ((`p' > 0 & "`abnetwork'" == "") | (`p' > 1 & "`abnetwork'" != "") | ("`interaction'" != "" & "`pcbnetwork'`mcbnetwork'" != "") ) & "`mc'" != "" {
				if "`inference'" == "bayesian" {
					local bayesnullest = "`bwd'\metapreg_bayesnullest"
				}
				
				capture noisily mcpreg `event' `nonevent' `total' `strataif', rid(`rid') sid(`studyid') studyid(`studyid') use(`use') regexpression(`regexpression') nu(`nu')  ///
					regressors(`regressors') level(`level') varx(`varx') typevarx(`typevarx')  `progress' /// 
					model(`modeli') modelopts(`modeloptsi') command0(`command0')  `interaction' `design' by(`by') `stratify' baselevel(`basecode') ///
					comparator(`Comparator')  link(`link')   ///
					 cov(`cov') inference(`inference') refsampling(`refsampling') bayesest(`bayesnullest')
										 
				if _rc == 0 {		
					mat `mctest' = r(mctest) 
				}
			}
		}
		local ++i
	}
	*End of loop
	cap drop  `numsid'
	
	qui keep if mu == 1
	
	//rownames for the matrix
	if "`stratify'" != "" & `i' > 1 {
		if "`abnetwork'`cov'" != "" {
			if "`cov'" == "unstructured" {
				if "`inference'" == "bayesian" {
					mat colnames `hetout' = Delta_ML log(BF) Post_prob tau2 sigma2 rho I2tau I2sigma 
				}
				else {
					mat colnames `hetout' = DF Chisq p tau2 sigma2 rho I2tau I2sigma 
				}
			}
			else {
				if "`inference'" == "bayesian" {
					mat colnames `hetout' = Delta_ML log(BF) Post_prob tau2 sigma2 I2tau I2sigma 
				}
				else {
					mat colnames `hetout' = DF Chisq p tau2 sigma2 I2tau I2sigma 
				}
			}
		}
		else {
			if (`p' == 0) & ("`model'" == "random") & "`pcbnetwork'`mcbnetwork'" == "" {
				if "`inference'" == "bayesian" {
					mat colnames `hetout' = Delta_ML log(BF) Post_prob tau2 I2tau
				}
				else {
					mat colnames `hetout' = DF Chisq p tau2 I2tau 
				}
			}
			else {
				if "`inference'" == "bayesian" {
					mat colnames `hetout' = Delta_ML log(BF) Post_prob tau2
				}
				else {
					mat colnames `hetout' = DF Chisq p tau2 
				}
			}
		}
	}
	//If stratify, no overall for comparative
	if "`stratify'" != "" & ("`design'" == "comparative" | "`design'" == "mpair") {
		local overall "nooverall"
	}
	
	//Subset matrices 
	if "`scimethod'" == "t" {
		mat `neorawest' = (`rawest'[1..., 1..3] , `rawest'[1..., 7..9])
		mat `neoabsout' = (`absout'[1..., 1..3] , `absout'[1..., 7..9])
		if "`ratios'" == "" {
			mat `neorrout' = (`rrout'[1..., 1..3] , `rrout'[1..., 7..9])
			mat `neordout' = (`rdout'[1..., 1..3] , `rdout'[1..., 7..9])
			mat `neoorout' = (`orout'[1..., 1..3] , `orout'[1..., 7..9])
		}
	}
	else if "`scimethod'" == "wald" {
		mat `neorawest' = `rawest'[1..., 1..6] 
		mat `neoabsout' = `absout'[1..., 1..6]
		if "`ratios'" == "" {
			mat `neorrout' 	= `rrout'[1..., 1..6]
			mat `neordout' 	= `rdout'[1..., 1..6]
			mat `neoorout' 	= `orout'[1..., 1..6]
		}
	}
	else if "`scimethod'" == "eti" {
		mat `neorawest' = `rawest'[1..., 1..6] 
		mat `neoabsout' = `absout'[1..., 1..6]
		if "`ratios'" == "" {
			mat `neorrout' 	= `rrout'[1..., 1..6]
			mat `neordout' 	= `rdout'[1..., 1..6]
			mat `neoorout' 	= `orout'[1..., 1..6]
		}
	}
	else if "`scimethod'" == "hpd" {
		mat `neorawest' = (`rawest'[1..., 1..4] , `rawest'[1..., 7..8])
		mat `neoabsout' = (`absout'[1..., 1..4] , `absout'[1..., 7..8])
		
		if "`ratios'" == "" {
			mat `neorrout' = (`rrout'[1..., 1..4] , `rrout'[1..., 7..8])
			mat `neordout' = (`rdout'[1..., 1..4] , `rdout'[1..., 7..8])
			mat `neoorout' = (`orout'[1..., 1..4] , `orout'[1..., 7..8])
		}
	}
	
	/*If cov(freeint) the conditional estimates are irrelevant. 
	Hence do not report them instead of using the average intercept in the calculations for logits, RR, RD 
	This is automatic for Bayesian estimation. Change applies to frequentist estimation where the margins
	command uses the average intercept*/
	
	
	if "`model'" == "hexact"  {
			//exactabs
			if  (("`sumtable'" == "all")  | (strpos("`sumtable'", "abs") != 0)) | ("`sumtable'" == "")  {
				printmat, matrixout(`exactabsout') type(exactabs) p(`p') dp(`dp') power(`power') `continuous'  model(`model') power(`power') inference(`inference')
			}
			
			//or
			if (("`sumtable'" == "all") | (strpos("`sumtable'", "or") != 0)) & (("`catreg'" != " ") | ("`typevarx'" == "i"))  | ("`sumtable'" == "")  {
				
				//rr
				cap confirm matrix `exactorout'
				if _rc == 0 {
					printmat, matrixout(`exactorout') type(exactor) p(`p') dp(`dp') power(`power')  model(`model')
				}			
			}
	}
	else {
		//Extra tables
		//het
		if strpos("`model'", "random") != 0 | strpos("`model'", "betabin") != 0  {			
			printmat, matrixout(`hetout') type(het) dp(`dp') `design' model(`model')
		}
		
		//Raw estimates
		if  ((("`sumtable'" == "all") | (strpos("`sumtable'", "logit") != 0)) | ("`sumtable'" == "")) & ("`cov'" != "freeint" ) {
			printmat, matrixout(`neorawest') type(raw) p(`p') dp(`dp') power(`power') `continuous' model(`model') `stratify' link(`link') inference(`inference')
		}
		//abs
		if  ((("`sumtable'" == "all")  |(strpos("`sumtable'", "abs") != 0)) | ("`sumtable'" == "")) & ("`cov'" != "freeint" )  {
			printmat, matrixout(`neoabsout') type(abs) p(`p') dp(`dp') power(`power') `continuous'  model(`model') inference(`inference')  power(`power')
		}
		//Pop p 
		if  (("`sumtable'" == "all")  | (strpos("`sumtable'", "abs") != 0) | ("`sumtable'" == "") ) & strpos("`model'", "betabin") == 0  {

			printmat, matrixout(`popabsout') type(popabs) dp(`dp') power(`power') nsims(`nsims') model(`model') 
		}
		
		//RD
		if (("`sumtable'" == "all")  |  ("`sumtable'" == "")  | (strpos("`sumtable'", "rd") != 0)) & (("`catreg'" != " ") | ("`typevarx'" == "i"))  {
			if "`cov'" != "freeint" { 
				//rd
				cap confirm matrix `neordout'
				if _rc == 0 {
					printmat, matrixout(`neordout') type(rd) p(`p') dp(`dp') power(`power')  model(`model') inference(`inference') 
				}
				
				//RD equal
				if "`inltest'" == "yes" {
					cap confirm matrix `nltestRD'
					if _rc == 0 {
						printmat, matrixout(`nltestRD') type(rde) dp(`dp') inference(`inference')
					}
				}
			}
			
			if strpos("`model'", "betabin") == 0 {
				cap confirm matrix `poprdout'
				if _rc == 0 {
					printmat, matrixout(`poprdout') type(poprd) p(`p') dp(`dp') power(`power')  model(`model') nsims(`nsims')
				}
			}
		}
		//rr
		if (("`sumtable'" == "all") | ("`sumtable'" == "")  |  (strpos("`sumtable'", "rr") != 0)) & (("`catreg'" != " ") | ("`typevarx'" == "i"))   {
			if "`cov'" != "freeint" {
				//rr
				cap confirm matrix `neorrout'
				if _rc == 0 {
					printmat, matrixout(`neorrout') type(rr) p(`p') dp(`dp') power(`power')  model(`model')  inference(`inference') power(`power')
				}
				
				//rr equal
				if "`inltest'" == "yes" {
					cap confirm matrix `nltestRR'
					if _rc == 0 {
						printmat, matrixout(`nltestRR') type(rre) dp(`dp') inference(`inference')
					}
				}
			}
			
			if strpos("`model'", "betabin") == 0 {
				cap confirm matrix `poprrout'
				if _rc == 0 {
					printmat, matrixout(`poprrout') type(poprr) p(`p') dp(`dp') power(`power')  model(`model') nsims(`nsims')
				}
			}
		}
		
		//or
		if (("`sumtable'" == "all") | ("`sumtable'" == "")   | (strpos("`sumtable'", "or") != 0)) & (("`catreg'" != " ") | ("`typevarx'" == "i"))   {
			
			//rr
			cap confirm matrix `neoorout'
			if _rc == 0 {
				printmat, matrixout(`neoorout') type(or) p(`p') dp(`dp') power(`power')  model(`model') inference(`inference') power(`power')
			}
			
			//or equal
			if "`inltest'" == "yes" {
				cap confirm matrix `nltestOR'
				if _rc == 0 {
					printmat, matrixout(`nltestOR') type(ore) dp(`dp') inference(`inference')
				}
			}
			if strpos("`model'", "betabin") == 0 {
				cap confirm matrix `poprrout'
				if _rc == 0 {
					printmat, matrixout(`poporout') type(popor) p(`p') dp(`dp') power(`power')  model(`model') nsims(`nsims')
				}
			}
		}	
	}
	
	tempfile master
	qui save "`master'"
	
	if "`itable'" == "" | "`graph'" == "" {
		local metricindex 0	
		
		foreach metric of local outplot {
		
			local ++metricindex
			local icimethod
			local scimethod
			local neorcols
			local neolcols
			local neoglcols
			local neogrcols
			local neoflcols
			local neofrcols
			local neoclcols
			local neocrcols	
			local groupvar
			local neologscale
			local compabs
			
			//Specify the grouping variable for proportions and/or ratios/difference
			if "`subgroup'" == "" & ("`catreg'" != "" | "`typevarx'" =="i" ) {
				if "`metric'"== "abs" {
					if "`typevarx'" =="i" {
						//comparative
						local groupvar = "`varx'"
					}
					else {
						//other designs
						local groupvar : word 1 of `catreg'
					}
				}
				//Comparative
				if strpos("`metric'", "r") != 0 & "`varx'" != ""  & "`catreg'" != "" {
					local groupvar : word 1 of `catreg'
				}
			}
			
			if "`by'" != "" {
				local groupvar "`by'"
				local byvar "`by'"
			}
			if "`design'" == "abnetwork" {	
				//How many groups
				*qui label list `first'
				*local ngroups = r(max)
				
				qui levelsof `first', local(codelevels)
				local ngroups = r(r)
				
				if "`aliasdesign'" != "" {
					cap assert `ngroups' == 2
					cap assert "`aliasdesign'" == "comparative"
				}
				else {
					local groupvar "`first'"
					local overall "nooverall"
					
					if strpos("`metric'", "r") != 0 {
						local itable "noitable"
					}
				}		
			} 
			
			//nullify groupvar if its the studyid
			if "`groupvar'" == "`studyid'" {
				local groupvar
			}	
			if "`groupvar'" == "" {
				local subgroup nosubgroup
			}
						
			//Logscale only in or/rr 
			if "`logscale'" != "" {
				*cap assert ("`metric'" != "abs" & "`metric'" != "lor" & "`metric'" != "lrr" & "`metric'" != "rd")
				if ("`metric'" != "abs" & "`metric'" != "lor" & "`metric'" != "lrr" & "`metric'" != "rd") {
					local neologscale "logscale"
				}
				/*if _rc != 0 {
					di as res "Option `logscale' not allowed"
					
				}*/
			}
			
			if "`sumstat'" == "" {
				if "`metric'" == "abs" {
					local neosumstat = "Proportion"
				}
				else if "`metric'" == "rd"  {
					local neosumstat = "Probability Difference"
				}
				else if "`metric'" == "rr"  {
					local neosumstat = "Proportion Ratio"
				}
				else if "`metric'" == "or"  {
					local neosumstat = "Odds Ratio"
				}
				else if "`metric'" == "lrr"  {
					local neosumstat = "Log Proportion Ratio"
				}
				else if "`metric'" == "lor"  {
					local neosumstat = "Log Odds Ratio"
				}
			}
			else {
				local neosumstat "`sumstat'"
			}
			
			//icimethod
			if "`cimethod'" != "" {
				tokenize "`cimethod'", parse(",")
				//inidividual CI method
				if "`1'" != "," {
					local icimethod = strltrim("`1'")
				}
				
				//Summary CI method
				if "`1'" == "," {
					local scimethod = strltrim("`2'")
				}
				
				//Indiviadual and summary CI method
				if "`3'" != "" {
					local scimethod = strltrim("`3'")
					local icimethod = strltrim("`1'")
				}
			}
			
			if "`metric'" != "abs"  {
				/*
				if "`cimethod'" != "" & (strpos("`cimethod'", ",") != 1) { 
					tokenize "`cimethod'", parse(",")
					*local icimethod "`1'"
					local icimethod : word `metricindex' of `1'
				}
				if "`cimethod'" != "" & (strpos("`cimethod'", ",") == 1) { 
					tokenize "`cimethod'", parse(",")		
					if "`2'" != "" {
						local scimethod = strltrim("`2'")
					}
				}
				*/
										
				//iOR/iLOR
				if ("`icimethod'" != "") & strpos("`metric'",  "or") != 0 {
					if (strpos("`icimethod'", "ex") != 1) & (strpos("`icimethod'", "wo") != 1) & (strpos("`icimethod'", "co") != 1)  {
						di as error "Option `icimethod' not allowed in cimethod(`cimethod')"
						exit	
					}
				}
				
				if ("`icimethod'" == "") & strpos("`metric'",  "or") != 0  {
					if ("`design'" == "mcbnetwork" | "`design'" == "mpair") {
						local icimethod "E.Fisher"
					}
					else {
						local icimethod "woolf"
					}
				}
				
				//iRR/iLRR
				if ("`design'" == "mcbnetwork" | "`design'" == "mpair") & strpos("`metric'", "rr") != 0  {
					local icimethod "CML"
				}
				
				if ("`design'" == "pcbnetwork" | "`design'" == "comparative" | "`aliasdesign'" == "comparative") & strpos("`metric'", "rr") != 0   {
					if ("`icimethod'" != "") {
						if (strpos("`icimethod'", "koo") != 1) & 	///
						(strpos("`icimethod'", "ka") != 1) & 		///
						(strpos("`icimethod'", "bail") != 1) & 		///
						(strpos("`icimethod'", "asin") != 1) & 		///	
						(strpos("`icimethod'", "noe") != 1) &		///				
						(strpos("`icimethod'", "adlo") != 1) {
							di as error "Option `icimethod' not allowed in cimethod(`cimethod')"
							exit	
						}
					}
					if ("`icimethod'" == "")  {
						local icimethod "koopman"
					}
				}
				//iRD
				if ("`icimethod'" == "") & strpos("`metric'",  "rd") != 0  {
					local icimethod "Newcombe"
				}
			}
			
			if "`metric'" == "abs"  {
				/*
				if "`cimethod'" != "" { 
					tokenize "`cimethod'", parse(",")
					local icimethod "`1'"
				}
				*/
				if "`icimethod'" != "" {
					if (strpos("`icimethod'", "ex") != 1) & (strpos("`icimethod'", "wi") != 1) &  (strpos("`icimethod'", "wa") != 1) & (strpos("`icimethod'", "e") != 1) & (strpos("`icimethod'", "ag") != 1) & (strpos("`icimethod'", "je") != 1)   {
						di as error "Option `icimethod' not allowed in cimethod(`cimethod')"
						exit	
					}
				}
				else {
					local icimethod "wilson"
				}
			}
			
				
			//CI
			if "`metric'" != "abs" {
				local se
			}
			
			//
			if "`metric'" != "abs" & "`design'" == "abnetwork" & "`aliasdesign'" == "" {
				local smooth
			}
			
			if "`metric'" == "abs" {
				local modelstats "`modelp' `modelplci' `modelpuci'"
			}
			else {
				local modelstats "`model`metric'' `model`metric'lci' `model`metric'uci'"
			}
			
			if "`metric'" == "abs" {
				local mes "`modelp'"
				local mlci "`modelplci'"
				local muci "`modelpuci'"
			}
			else {
				local mes "`model`metric''"
				local mlci "`model`metric'lci'"
				local muci "`model`metric'uci'"
			}
			
			if "`design'" == "mcbnetwork" | "`design'" == "pcbnetwork"  {
				*widesetup `event' `total', sid(`rid') idpair(`assignment')  jvar(`comparator')

				sort `rid'
				cap drop `assignment' `ipair' `nonevent'
				
				qui reshape wide `event' `total' _WT /// 
								`modelp' `modelplci' `modelpuci'  //// 
								`modelrr' `modelrrlci' `modelrruci'  ///
								`modelrd' `modelrdlci' `modelrduci'  ///
								`modellrr' `modellrrlci' `modellrruci'  ///
								`modellor' `modellorlci' `modelloruci'  ///
								`modelor' `modelorlci' `modeloruci'  ///
								, i(`rid') j(`idpair')

				//Add the weights
				qui gen _WT = _WT0 + _WT1
				qui drop _WT0  _WT1
				
				qui gen `mes' = `mes'1
				qui drop `mes'0 `mes'1
				
				qui gen `mlci' = `mlci'1
				qui drop `mlci'0 `mlci'1
				
				qui gen `muci' = `muci'1
				qui drop `muci'0 `muci'1
			}
			if "`design'" == "mpair" {
				sort `rid'
				cap drop `ipair' `nonevent'
				
				qui reshape wide `event' `total' _WT `modeles'   `modellci' `modeluci', i(`rid') j(`idpair')
				
				//Add the weights
				qui gen _WT = _WT0 + _WT1
				qui drop _WT0  _WT1
				
				qui gen `mes' = `mes'1
				qui drop `mes'0 `mes'1
				
				qui gen `mlci' = `mlci'1
				qui drop `mlci'0 `mlci'1
				
				qui gen `muci' = `muci'1
				qui drop `muci'0 `muci'1
			
			}
			 
						
			//Generate the CI		
			metapregci `depvars', studyid(`studyid') first(`first') es(`es') se(`se') uci(`uci') lci(`lci') `design' aliasdesign(`aliasdesign') ///
				id(`id') rid(`rid') regressors(`regressors') outplot(`metric') level(`level') ///
				icimethod(`icimethod') lcols(`lcols') rcols(`rcols')  sortby(`sortby') by(`by')  ///
				modeles(`mes') modellci(`mlci') modeluci(`muci') `smooth' 
			
			local neodepvars = r(depvars)
			local neorcols = r(rcols)
			local neolcols = r(lcols)
			local neosortby = r(sortby)
			if (`p' > 0) {
				local indvars = r(regressors)
			}
			
			/*
			if `mdf' == 0 {
				local smooth
			}
			*/
			
			//identifier as ordered  by ES
			qui gen `cid' = .	
						
			prep4show `id' `cid' `use' `neolabel' `es' `lci' `uci' `modelstats', `design' aliasdesign(`aliasdesign') ///
				sortby(`neosortby') groupvar(`groupvar') grptotal(`grptotal') se(`se') 	///
				outplot(`metric') rrout(`rrout') poprrout(`poprrout') rdout(`rdout') poprdout(`poprdout') poplrrout(`poplrrout') ///
				orout(`orout') poporout(`poporout') poplorout(`poplorout') exactorout(`exactorout') ///
				absout(`absout') popabsout(`popabsout') exactabsout(`exactabsout') absoutp(`absoutp') hetout(`hetout')	///
				`subgroup' `summaryonly' dp(`dp') pcont(`pcont') model(`model') `prediction' inference(`inference')	///
				`overall' download(`download') indvars(`indvars') depvars(`neodepvars') `stratify' level(`level') `enhance' stat(`stat')
				
			if "`itable'" == "" {
				sort `id'
				
				disptab `id'  `use' `neolabel' `es' `lci' `uci' `grptotal' `modelstats', ///
					`itable' dp(`dp') power(`power') design(`design') aliasdesign(`aliasdesign') `summaryonly' ///
					`subgroup' sumstat(`neosumstat') level(`level') `wt' `smooth' inference(`inference') ///
					scimethod(`scimethod') icimethod(`icimethod') model(`model') groupvar(`groupvar') outplot(`metric') stat(`stat')
			}
			
			if "`graph'" == "" {
				//Gather graph options
				if "`astext'" != "" {
					local neoastext "astext(`astext')"
				}
				if "`texts'" != "" {
					local neotexts "texts(`texts')"
				}
				if "`neolcols'" != "" {
					local neolcols "lcols(`neolcols')"
				}
				if "`neorcols'" != "" {
					local neorcols "rcols(`neorcols')"
				}
				if "`xlabel'" != "" {
					local neoxlabel "xlabel(`xlabel')"
				}
				else if strpos("`metric'", "abs") != 0 & "`pxlabel'" != ""  {
					local neoxlabel "xlabel(`pxlabel')"
				}
				else if strpos("`metric'", "r") != 0 & "`rxlabel'" != ""  {
					local neoxlabel "xlabel(`rxlabel')"
				}
				local goptions `" `fplot' `neolcols' `neorcols' `overall' `ovline' `stats' `box' `double' `neoastext' `ciopts' `diamopts' `olineopts' `pointopts' `boxopts' `predciopts' `prediction' `subline' `neotexts' `neoxlabel' `xline' `xtick'  `neologscale' `options'"'
				
				metapplotcheck,`summaryonly' `goptions'  //plots advance housekeeping
				local goptions = r(plotopts)
				local neoglcols = r(lcols)
				local neogrcols = r(rcols)
				if "`neoglcols'" == " " { //if empty
					local neoglcols
				}
				if "`neogrcols'" == " " { //if empty
					local neogrcols
				}
			}
						
			//Draw the forestplot
			if "`graph'`fplot'" == "" {
				if `"`foptions'"' != "" {
					metapplotcheck,`summaryonly' `goptions' `foptions'  //Forest plot advance housekeeping
					local neofoptions = r(plotopts)
					local neoflcols = r(lcols)
					local neofrcols = r(rcols)
					if "`neoflcols'" == " " { //if empty
						local neoflcols
					}
					if "`neofrcols'" == " " { //if empty
						local neofrcols
					}
				}
				else  {
					local neofoptions `"`goptions'"'
					if "`neoglcols'" != "" {
						local neoflcols = "`neoflcols' `neoglcols'"
					}
					if "`neogrcols'" != "" {
						local neofrcols = "`neofrcols' `neogrcols'"
					}
				}
				
				sort `id' 
										
				metapplot `es' `lci' `uci' `use' `neolabel' `grptotal' `id' `modelstats', model(`model') ///	
					studyid(`studyid') power(`power') dp(`dp') level(`level') ///
					groupvar(`groupvar')  type(fplot) sumstat(`neosumstat') ///
					outplot(`metric') lcols(`neoflcols') rcols(`neofrcols') ///
					`neofoptions'  design(`design') aliasdesign(`aliasdesign') `wt' `smooth' varxlabs(`varxlabs')
			}
			
			//Toggle for comparative
			if ("`design'" == "comparative" | "`aliasdesign'" == "comparative") & "`metric'" == "abs" {
				local compabs "compabs"
			}
			//Draw the catterpillar plot
			if "`graph'" == "" & "`catpplot'" != ""  & "`compabs'" == ""  {
				
				local neorcols
				local neolcols
			
				if "`coptions'" != "" {
					metapplotcheck,`summaryonly' `goptions' `coptions'  //Catterpillar plot advance housekeeping
					local neocoptions = r(plotopts)
					local neoclcols = r(lcols)
					local neocrcols = r(rcols)
					if "`neoclcols'" == " " { //if empty
						local neoclcols
					}
					if "`neocrcols'" == " " { //if empty
						local neocrcols
					}
				}
				else {
					local neocoptions = `" `goptions'"'
					if "`neoglcols'" != "" {
						local neoclcols = "`neoclcols' `neoglcols'"
					}
					if "`neogrcols'" != "" {
						local neocrcols = "`neocrcols' `neogrcols'"
					}
				}
				
				sort `cid'
				
				metapplot `es' `lci' `uci' `use' `neolabel' `grptotal' `cid' `modelstats', model(`model') ///	
					studyid(`studyid') power(`power') dp(`dp') level(`level') ///
					groupvar(`groupvar')  type(catpplot) sumstat(`neosumstat') ///
					outplot(`metric') lcols(`neoclcols') rcols(`neocrcols') ///
					`neocoptions' design(`design') aliasdesign(`aliasdesign') `wt' `smooth' varxlabs(`varxlabs')
				
			}
			
			use "`master'", clear		
		}
	}
		
	//model comparison
	if ((`p' > 0 & "`design'" != "abnetwork") | (`p' > 1 & "`design'" == "abnetwork")) & ("`mc'" != "") {
		cap confirm matrix `mctest'
		if _rc == 0 {
			printmat, matrixout(`mctest') type(mc) dp(`dp') 
			//Just initialize
			gettoken first confounders : regressors
			local p: word count `regressors'
			
			local redindex 0
		
			di as txt _n "Fitted reduced model(s) for comparison"
			if "`abnetwork'`interaction'" !="" {
				if "`mcbnetwork'`pcbnetwork'" != "" {
					local confariates "`comparator'"	
				}
				else {
					local confariates "`confounders'"
				}
			}
			else {
				local confariates "`regressors'"
			}
			local initial 1
			foreach c of local confariates {
				
				if ("`interaction'" != "" & "`pcbnetwork'`mcbnetwork'" != "")  {
						local omterm = "`c'*`ipair'"
						gettoken start end : regexpression
						local eqreduced = "Ipair + `end'"
						
				}
				else {						
					foreach term of local regexpression {
						if "`interaction'" != "" {
							if strpos("`term'", "`c'#") != 0 & strpos("`term'", "`first'") != 0 {
								local omterm = "`c'*`first'"
							}
						}
						else{
							if "`model'" == "cbbetabin" {
								if ("`term'" == "i.`c'#c.mu")|("`term'" == "c.`c'#c.mu") {
									local omterm = "`c'"
								}
							}
							else {
								if ("`term'" == "i.`c'")|("`term'" == "c.`c'")|("`term'" == "`c'") {
									local omterm = "`c'"
								} 
							}
						}
					}
					local eqreduced = subinstr("`nu'", "+ `omterm'", "", 1)
				}
				
				local ++redindex
				if "`model'" == "cbbetabin" {
					di as res _n "`redindex'. Ommitted `omterm' in alpha"
				}
				else {
					di as res _n "`redindex'. Ommitted `omterm' in `link'(p)"
				}
				if "`model'"  == "random" {
					di as res "{phang} `link'(p) = `eqreduced' + `studyid'{p_end}"
				}
				else if "`model'" == "cbbetabin" {
					di as res "{phang} alpha = exp(`eqreduced'){p_end}"
				}
				else {
					di as res "{phang} `link'(p) = `eqreduced'{p_end}"
				}
				//Ultimate null model
				if (`p' > 1 & "`abnetwork'" == "") | (`p' > 2 & "`abnetwork'" != "")  {
					local ++redindex 
					if "`model'" == "cbbetabin" {
						di as res _n "`redindex'. Ommitted all covariate effects in alpha"
					}
					else {
						di as res _n "`redindex'. Ommitted all covariate effects in `link'(p)"
					}
				}
			}
		}		
	}
	
	cap ereturn clear
	
	
	//Save the matrices if only relevent. For cov==freeint, the conditional matrices are irrelevant
	
	cap confirm matrix `nltestRR'
	if _rc == 0 {
		ereturn matrix rrtest = `nltestRR'
		ereturn matrix rdtest = `nltestRD'
		ereturn matrix ortest = `nltestOR'
	}
	
	cap confirm matrix `rawest'
	if _rc == 0 {
		ereturn matrix rawest = `rawest'
		ereturn matrix popabsout = `popabsout'
	}
	
	cap confirm matrix `absout'
	if _rc == 0 {
		ereturn matrix absout = `absout'
	}
			
	cap confirm matrix `rrout'
	if _rc == 0 {
		ereturn matrix rdout = `rdout'
		ereturn matrix rrout = `rrout'
		ereturn matrix orout = `orout'
	}
	
	cap confirm matrix `matgof'
	if _rc == 0 {
		ereturn matrix gof = `matgof'
	}

	cap confirm matrix `mctest'
	if _rc == 0 {
		ereturn matrix mctest = `mctest'
	}
	
	cap confirm matrix `hetout'
	if _rc == 0 {
		ereturn matrix hetout = `hetout'
	}
	
	cap confirm matrix `covmat'
	if _rc == 0 {
		ereturn matrix covmat = `covmat'
	}
	
	cap confirm matrix `exactabsout'
	if _rc == 0 {
		ereturn matrix exactabsout = `exactabsout'
	}
	
	cap confirm matrix `exactorout'
	if _rc == 0 {
		ereturn matrix exactorout = `exactorout'
	}	
	
	cap confirm matrix `poprdout'
	if _rc == 0 {
		ereturn matrix poprdout = `poprdout'
		ereturn matrix poprrout = `poprrout'
		ereturn matrix poplrrout = `poplrrout'
		ereturn matrix poporout = `poporout'
		ereturn matrix poplorout = `poplorout'
	}
		
	restore
	if "`inference'" == "bayesian" {
		cap erase `bayesreps'.dta
		cap erase `bayesreps'.ster 
		di _n
		di _n
		di _n
		di _n
		di as re "{pmore}| The bayesian estimation commands saved a dataset {p_end}" 
		di 		 "{pmore}| `bayesest'.dta {p_end}" 
		di as re "{pmore}| containing the MCMC samples of the parameters to the disk.{p_end}" 		
		di as re "{pmore}| It is your responsibility to erase the dataset {p_end}" 
		di as re "{pmore}| after it is no longer needed.{p_end}"
		di 		`"{pmore}{stata "erase `bayesest'.dta":Click to erase the dataset}"'
	}
end

/*+++++++++++++++++++++++++	SUPPORTING FUNCTIONS: ESTCOVAR +++++++++++++++++++++++++
							Compose the var-cov matrices
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
cap program drop estcovar
program define estcovar, rclass

	syntax, bmatrix(name)  sid(varname) [varx(varname) cov(string) predcmd(string) level(real 95) abnetwork]
	*matrix is colvector
	tempname matcoef rosevar coefvar covmat matvar bmat vmat
	mat `matcoef' = `bmatrix''
	local bnrows = rowsof(`matcoef')
	
	local critvalue -invnorm((100-`level')/200)

	if "`predcmd'" == "meqrlogit_p" {
		local scalefn1 "exp"
		local scalefn2 "tanh"
		local scalepow "2"
	}
	else {
		local scalepow "1"
	}				
				
	if strpos("`cov'", "uns") != 0 {
		if "`predcmd'" == "meqrlogit_p" {
			*local covarexpression  "(tausq:(exp(_b[lns1_1_1:_cons]))^2) (sigmasq:(exp(_b[lns1_1_2:_cons]))^2) (covar:exp(_b[lns1_1_1:_cons])*exp(_b[lns1_1_2:_cons])*tanh(_b[atr1_1_1_2:_cons]))" 
			local covarexpression  "(lntau:_b[lns1_1_2:_cons]) (lnsigma:_b[lns1_1_1:_cons]) (rho:exp(_b[lns1_1_1:_cons])*exp(_b[lns1_1_2:_cons])*tanh(_b[atr1_1_1_2:_cons]))" 
		}
		else {
			local covarexpression  "(lntausq:ln(_b[/var(mu[`sid'])])) (lnsigmasq:ln(_b[/var(2.`varx'[`sid'])])) (covar:_b[/cov(2.`varx'[`sid'],mu[`sid'])])"
		}		
				
		local k = 3	
	}
	
	else if strpos("`cov'", "ind") != 0 | "`abnetwork'" != "" {
		if "`predcmd'" == "meqrlogit_p" {
			local covarexpression  "(lntau:_b[lns1_1_2:_cons]) (lnsigma:_b[lns1_1_1:_cons])"
		}
		else {
			if "`abnetwork'" == "" { 
				local covarexpression  "(lntausq:ln(_b[/var(mu[`sid'])])) (lnsigmasq:ln(_b[/var(2.`varx'[`sid'])]))" 
			}
			else {
				local covarexpression  "(lntausq:ln(_b[/var(mu[`sid'])])) (lnsigmasq:ln(_b[/var(mu[`sid'>`varx'])]))" 
			}
		}					
		local k = 2
	}
	else if strpos("`cov'", "common") != 0 | strpos("`cov'", "free") != 0 | "`cov'" == "" {
		if "`predcmd'" == "meqrlogit_p" {
			local covarexpression  "(lnvar:_b[lns1_1_1:_cons])"
		}
		else {
			if "`cov'" == "commonslope" |  "`cov'" == "" {
				local covarexpression  "(lntausq:ln(_b[/var(mu[`sid'])]))"
			}
			if "`cov'" == "commonint" | "`cov'" == "freeint" {		
				local covarexpression  "(lnsigmasq:ln(_b[/var(2.`varx'[`sid'])]))" 
			}
		}
		local k = 1
	}

	mat `covmat' = J(`k', 3, .)
	
	capture nlcom `covarexpression'
	
	if _rc == 0 {
		mat `bmat' = r(b)
		mat `vmat' = r(V)
		mat `vmat' = vecdiag(`vmat')
				
		forvalues r=1(1)`k' {
			local scalepow 1
			
			if `r' != 3 {
				local scalefn "exp"
				
				if "`predcmd'" == "meqrlogit_p" {
					local scalepow 2
				}
			}
			else {
				local scalefn
			}
			
			mat `covmat'[`r', 1] = (`scalefn'(`bmat'[1,`r']))^`scalepow'
			mat `covmat'[`r', 2] = (`scalefn'(`bmat'[1, `r'] - `critvalue'*sqrt(`vmat'[1,`r'])))^`scalepow'
			mat `covmat'[`r', 3] = (`scalefn'(`bmat'[1, `r'] + `critvalue'*sqrt(`vmat'[1,`r'])))^`scalepow'
		}
	}
		
	if `k' == 3 {
		local rowids = "tausq sigmasq covar"
	}
	else {
		if "`cov'" != "" {
			local rowids = "tausq sigmasq"
		}
		else {
			local rowids = "tausq"
		}
	}
	
	if "`cov'" == "commonslope" {
			mat `covmat' = (`covmat' \ J(1, 3, 0))
	}
	if "`cov'" == "commonint"  {
		mat `covmat' = (J(1, 3, 0) \ `covmat')
	}
	if "`cov'" == "freeint"  {
		mat `covmat' = (J(1, 3, .) \ `covmat')
	}
	
	mat rownames `covmat' = `rowids'
	
	*return matrix rosevar = `rosevar' 
	return matrix covmat = `covmat' 
	return local k = `k' 
end

/**************************************************************************************************
							METAPREGCI - CONFIDENCE INTERVALS
**************************************************************************************************/
capture program drop metapregci
program define metapregci, rclass
	#delimit ;
	syntax varlist(min=2 max=4), studyid(varname) [first(varname) es(name) se(name) uci(name) lci(name)
		id(name) rid(varname) regressors(varlist) outplot(string) level(integer 95) by(varname)
		icimethod(string) lcols(varlist) rcols(varlist) mpair mcbnetwork pcbnetwork sortby(varlist) 
		comparative abnetwork general aliasdesign(string) modeles(varname) modellci(varname) modeluci(varname) smooth
		vlist(string asis) cc0(string asis) cc1(string asis)
		];
	#delimit cr
	tempvar uniq event event1 event2 total total1 total2 a b c d idpair
	*gettoken idpair confounders : regressors
	qui {
		tokenize `varlist'
		if "`mcbnetwork'`pcbnetwork'`mpair'" == "" {
			generate `event' = `1'
			generate `total' = `2'
			local depvars "`1' `2'"
		}
		else if "`mcbnetwork'`mpair'" != "" {
			gen `a' = `1'
			gen `b' = `2'
			gen `c' = `3'
			gen `d' = `4'
			local depvars "`1' `2' `3' `4'"
			gen `id' = _n
		}
		else {
			generate `event1' = `1'
			generate `event2' = `2'
			generate `total1' = `3'
			generate `total2' = `3'
			local depvars "`1' `2' `3'"
			gen `id' = _n
		}

		if "`outplot'" == "lor" | "`outplot'" == "lrr" {
			local transform "transform"
		}
		
		if "`outplot'" != "abs" {
			if "`abnetwork'" != "" & "`aliasdesign'" == "" {
				gen `id' = _n
				gen `es' = .
				gen `lci' = .
				gen `uci' = .
			}
			if "`mcbnetwork'`mpair'" != ""   { 
				if strpos("`outplot'", "rr") != 0 {
					//constrained maximum likelihood estimation
					cmlci `a' `b' `c' `d', r(`es') upperci(`uci') lowerci(`lci') alpha(`=1 - `level'*0.01') `transform'
				}
				if strpos("`outplot'", "or") != 0 {
					orccci `a' `b' `c' `d', r(`es') upperci(`uci') lowerci(`lci') level(`level') `mcbnetwork' `mpair' `transform'
				}
			}
			if "`pcbnetwork'" !="" {
				if strpos("`outplot'", "rr") != 0 {
					rrci `event1' `total1' `event2' `total2', r(`es') upperci(`uci') lowerci(`lci') alpha(`=1 - `level'*0.01') icimethod(`icimethod') `transform'
				}
				if strpos("`outplot'", "or") != 0 {
					orccci `event1' `total1' `event2' `total2', r(`es') upperci(`uci') lowerci(`lci') level(`level')  icimethod(`icimethod') `transform'
				}
			}
			if "`comparative'" != "" | "`aliasdesign'" == "comparative"  {
				egen `id' = group(`studyid' `by')
								
				sort `id' `rid'
				by `id': egen `idpair' = seq()

				count
				local Nobs = r(N) /*Check if the number of studies is half*/
				cap assert  mod(`Nobs', 2) == 0 
				if _rc != 0 {
					di as error "More than two observations per study for some studies"
					exit _rc, STATA
				}

				sort `id'  `idpair'
				
				if "`=`studyid'[1]'" != "`=`studyid'[2]'" {
					di as error "Data not properly sorted. `studyid' in row 1 and 2 should be the same. "
					exit _rc, STATA
				}
				
				widesetup `event' `total' `confounders' , sid(`id') idpair(`idpair') sortby(`sortby')
				local vlist = r(vlist)
				local cc0 = r(cc0)
				local cc1 = r(cc1)
				
				if "`comparative'" != "" | "`aliasdesign'" == "comparative" {
					if strpos("`outplot'", "rr") != 0 {
						rrci `event'1 `total'1 `event'0 `total'0, r(`es') upperci(`uci') lowerci(`lci') alpha(`=1 - `level'*0.01') icimethod(`icimethod') `transform'
					}
					if strpos("`outplot'", "rd") != 0 {
						rdci `event'1 `total'1 `event'0 `total'0, r(`es') upperci(`uci') lowerci(`lci') alpha(`=1 - `level'*0.01') icimethod(`icimethod') 
					}
					if strpos("`outplot'", "or") != 0 {
						orccci `event'1 `total'1 `event'0 `total'0, r(`es') upperci(`uci') lowerci(`lci') level(`level')  icimethod(`icimethod') `transform'
					}
				}
			}
			
			if "`comparative'" != "" | "`aliasdesign'" == "comparative" {	
					//Rename the varying columns
					local newcc0: label `first' `cc0'
					local newcc1: label `first' `cc1'
				
					foreach v of local vlist {
						rename `v'0 `v'_`cc0'
						label var `v'_`cc0' "`v'_`newcc0'"
						rename `v'1 `v'_`cc1'
						label var `v'_`cc1' "`v'_`newcc1'"
					}
					
					//Add the weights
					
					gen _WT = 0 
					replace _WT = _WT  +  _WT_1 if _WT_1 != .
					replace _WT = _WT  +  _WT_2 if _WT_2 != .
					
					*qui gen _WT = _WT_1 + _WT_2

					drop _WT_1  _WT_2
				
				//Remove unnecessary columns
				cap confirm variable `modeles'_1
				if _rc == 0 {
					drop `modeles'_1
					rename `modeles'_2 `modeles'
				}
				
				cap confirm variable `modellci'_1
				if _rc == 0 {
					drop `modellci'_1
					rename `modellci'_2 `modellci'
				}
				
				cap confirm variable `modeluci'_2
				if _rc == 0 {
					drop `modeluci'_1
					rename `modeluci'_2 `modeluci'
				}
								
				//make new lcols		
				foreach lcol of local lcols {
					local lenvar = strlen("`lcol'")
					
					foreach v of local vlist {
						local matchstr = substr("`v'", 1, `lenvar')
						
						if strmatch("`matchstr'", "`lcol'") == 1 {
							continue, break
						}
					}
					
					if strmatch("`matchstr'", "`lcol'") == 1 {
						local lcols_r "`lcols_r' `lcol'_`cc0' `lcol'_`cc1'"
					}
					else {
						local lcols_r "`lcols_r' `lcol'"
					}
				}
				local lcols "`lcols_r'"
				
				//make new rcols
				foreach rcol of local rcols {
					local lenvar = strlen("`rcol'")

					foreach v of local vlist {
						local matchstr = substr("`v'", 1, `lenvar')
						
						if strmatch("`matchstr'", "`rcol'") == 1 {
							continue, break
						}
					}
					
					if strmatch("`matchstr'", "`rcol'") == 1 {
						local rcols_r "`rcols_r' `rcol'_`cc0' `rcol'_`cc1'"
					}
					else {
						local rcols_r "`rcols_r' `rcol'"
					}
				}
				local rcols "`rcols_r'"
				
				//make new sortby
				foreach byv of local sortby {
					local lenvar = strlen("`byv'")

					foreach v of local vlist {
						local matchstr = substr("`v'", 1, `lenvar')
						
						if strmatch("`matchstr'", "`byv'") == 1 {
							continue, break
						}
					}
					
					if strmatch("`matchstr'", "`byv'") == 1 {
						local rcols_r "`sortby_r' `byv'_`cc0' `byv'_`cc1'"
					}
					else {
						local sortby_r "`sortby_r' `byv'"
					}
				}
				local sortyby "`sortby_r'"
				
				//make new depvars		
				foreach depvar of local depvars {
					local lenvar = strlen("`depvar'")

					foreach v of local vlist {
						local matchstr = substr("`v'", 1, `lenvar')
						
						if strmatch("`matchstr'", "`depvar'") == 1 {
							continue, break
						}
					}
					
					if strmatch("`matchstr'", "`depvar'") == 1 {
						local depvars_r "`depvars_r' `depvar'_`cc0' `depvar'_`cc1'"
					}
					else {
						local depvars_r "`depvars_r' `depvar'"
					}
				}
				
				local depvars "`depvars_r'"
				
				//make new indvars
				foreach indvar of local confounders {
					local lenvar = strlen("`indvar'")

					foreach v of local vlist {
						local matchstr = substr("`v'", 1, `lenvar')
						
						if strmatch("`matchstr'", "`indvar'") == 1 {
							continue, break
						}
					}
					
					if strmatch("`matchstr'", "`indvar'") == 1 {
						local indvars_r "`indvars_r' `indvar'_`cc0' `indvar'_`cc1'"
					}
					else {
						local indvars_r "`indvars_r' `indvar'"
					}
				}
				local regressors "`indvars_r'"
				local p: word count `confounders' 
				if `p' == 0 {
					local regressors = " "
				}
			}
		}
		else {
			metapreg_propci `total' `event', p(`es') se(`se') lowerci(`lci') upperci(`uci') icimethod(`icimethod') level(`level')
			gen `id' = _n
		}
	}
	if "`rcols'" =="" {
		local rcols = " "
	}
	if "`lcols'" =="" {
		local lcols = " "
	}
	if "`sortby'" == "" {
		local sortby = " "
	}
	return local regressors = "`regressors'"
	return local depvars = "`depvars'"
	return local rcols = "`rcols'"
	return local lcols = "`lcols'"
	return local sortby = "`sortby'"
end
/**************************************************************************************************
							PREG - MAIN REGRESSION 
**************************************************************************************************/
capture program drop preg
program define preg, rclass

	#delimit ;

	syntax varlist(min=3) [if] [in], sid(varname) studyid(varname) use(varname) [
		regexpression(string) regexpression2(string) nu(string) baselevel(passthru) rid(varname)
		regressors(varlist) varx(varname) typevarx(string) comparator(varname) 
		catreg(varlist) contreg(varlist)
		scimethod(string) cov(string)
		level(integer 95)
		DP(integer 2)
		progress
		model(string) modelopts(string asis) outplot(string)
		noMC noCONstant
		interaction	
		comparative mpair mcbnetwork pcbnetwork abnetwork
		aliasdesign(string)
		by(varname) stratify
		GOF nsims(string) link(string) bayesrepsfilename(string asis)
		modeles(varlist) modelse(varname) modellci(varlist) modeluci(varlist) 
		smooth computewt inference(string) refsampling(string) stat(string)
		bayesest(string asis) bayesreps(string asis)
			*];

	#delimit cr
	marksample touse, strok 
		
	tempvar event nonevent total invtotal predevent ill iw insample
	tempname coefmat coefvar testlr V logodds absout absoutp rrout rdout orout nltestRR nltestRD nltestOR  ///
			 hetout mctest absexact newobs matgof popabsout poprrout poplrrout poporout poplorout ///
			 rosevar covmat rawest rawestp coeflor coefor lorci exactabsout exactorout ///
			 exactlorout lnrho bayestest  bayesstats poprdout
	
	tokenize `varlist'
	local event = "`1'"
	local nonevent = "`2'"
	local total = "`3'"
	
	//fit the model
	if "`progress'" != "" {
		local echo noi
	}
	else {
		local echo qui
	}
	//Just initialize
	if "`mpair'" != "" {
		local first "`varx'"
	}
	else {
		gettoken first confounders : regressors
	}
	local p: word count `regressors'
	
	
	if "`mcbnetwork'`pcbnetwork'" != "" {		
		tokenize `regexpression'
		local one "`1'"
		local two "`2'"
		local three "`3'"
		
		if "`interaction'" != "" {
			tokenize `one', parse("#")
			tokenize `1', parse(".")
			local ipair "`3'"
			
			tokenize `two', parse(".")
			local index "`3'"
		}
		else {
			tokenize `two', parse(".")
			local ipair "`3'"
		
			tokenize `three', parse(".")
			local index "`3'"
		}
	}
	
	//Which outcome to model
	if "`link'" == "loglog" & "`model'" != "crbetabin" {
		local outcome "`nonevent'"
	}
	else {
		local outcome "`event'"
	}
	
	`echo' fitmodel `outcome' `total' if `touse', `modelopts' model(`model') regexpression(`regexpression') ///
		sid(`sid') studyid(`studyid') level(`level') nested(`first') `abnetwork' `comparative' `mpair' cov(`cov') ///
		link(`link') p(`p')  bayesest(`bayesest') inference(`inference') refsampling(`refsampling') `progress' `interaction'
		
	//Returned model	
	local getmodel = r(model)
	estimates store metapreg_modest
	local lnvar = r(lnvar)
	local sd = r(sd)
	
	qui {
		replace _ESAMPLE = e(sample) 
		replace `use' = 1 if (_ESAMPLE == 1)
	}
	
	local mdf = .

	if "`inference'" == "bayesian" {
		gen `insample' = e(sample) 
		
		gsort -`insample' `rid' `varx'
		//Generate the predictions
		qui bayespredict {_mu} if e(sample), saving("`bayesreps'", replace) rseed(1) 
	}
	else {
		local predcmd = e(predict)
		
		mat `coefmat' = e(b)
		mat `coefvar' = e(V)
		
		local DF = e(N) -  e(k)
		local mdf = e(df) //mdf = 0 if saturated model
	}
	
	if "`computewt'" !="" {
		qui {
			if "`model'" == "random" {
				//if random, needs atleast 7 studies to run predict command
				count
				local nobs = r(N)
				if ((`nobs' < 7) & ("`model'" == "random")) {
					local multipler = int(ceil(7/`nobs'))
					expand `multipler', gen(`newobs')
				}
			}
			
			if "`model'" == "cbbetabin" {
				predictnl `predevent' = invlogit(xb() - _b[_cons])*`total'
			}
			else if "`model'" == "crbetabin" {
				predict `predevent', n
			}
			else if strpos("`model'", "bayes") == 1 {
				bayespredict `predevent', mean rseed(1) 
			}
			else if "`predcmd'" == "poisso_p" {
				predict `predevent', n
			}
			else {
				//melogit | meqrlogit | binreg | mepoisson
				predict `predevent', mu
			}
	
			if "`model'" == "random" {
				//Revert to original data if filler data was generated
				if (`nobs' < 7)  {
					keep if !`newobs'
				}
			}
					
			//compute the weight
			gen `iw' = `total'*(`predevent'/`total')*(1 - `predevent'/`total') if  `predevent' !=`total'
			replace `iw' = `total' if  `predevent'==`total'  

		
			//compute the relative weight
			sum `iw' if (_ESAMPLE == 1 /*& mu == 1*/)
			local W = r(sum)
			
			//compute the weights
			replace _WT = (`iw'/`W')*100 if (_ESAMPLE == 1 /*& mu == 1*/) & (_WT == .) & (`iw' != .)
		}
	}
		
	//FE 
	local BHET = .
	local P_BHET = .
	local DF_BHET = .
	local MLdiff  = .
	local BF  = .
	local postprob = .

	if "`getmodel'" == "random" {
		local BHET = e(chi2_c)
		local P_BHET = e(p_c)
		if "`abnetwork'" == "" {
			local DF_BHET = 1
		}
		else {
			local DF_BHET = 2
		}
		
		capture estcovar, bmatrix(`coefmat') cov(`cov') predcmd(`predcmd') level(`level') sid(`sid') varx(`first') `abnetwork'
		local DF_BHET = r(k)
		mat `covmat' = r(covmat)  //var-cov matrix

		mat colnames `covmat' = Mean Lower Upper
	}
	else if "`getmodel'" == "crbetabin" {
		local BHET = e(chi2_c)
		local P_BHET = `=chi2tail(1, e(chi2_c))*0.5'
		local DF_BHET = 1
		local SIGMA = e(sigma)
		
		qui nlcom (rho:_b[/lnsigma]), level(`level')
		mat `covmat' = r(table)
		mat `covmat' = (exp(`covmat'[1,1]), exp(`covmat'[5,1]), exp(`covmat'[6,1]))
		mat rownames `covmat' = Sigma
		mat colnames `covmat' = Mean Lower Upper
	}
	else if "`getmodel'" == "cbbetabin" {
		//obtain ln rho  = ln (1/(1 + alpha + beta))
		qui {		
			margins , expression(-xb() - _b[_cons]) at(mu==1) 
			
			margins , expression(ln(1/(1 + exp(xb()) + exp(_b[_cons])))) at(mu==1) 

			mat `lnrho' = r(table)

			//fit the fixed effects model, to test if overdispersion
			fitmodel `outcome' `total' if `touse' & mu == 1, `modelopts' model(fixed) regexpression(`regexpression2') ///
				sid(`sid')  level(`level') nested(`first') `abnetwork' `comparative' `mpair' cov(`cov') link(`link') p(`p') ///
				inference(`inference') refsampling(`refsampling') `progress' 
				
			estimates store metapreg_Null
			
			//LR test the model
			capture lrtest metapreg_modest metapreg_Null, force
			if _rc == 0 {
				local BHET = r(chi2)
				local P_BHET = `=chi2tail(1, `BHET')*0.5'
				local SIGMA = exp(`lnrho'[1,1])
				local DF_BHET = 1
			}
			estimates drop metapreg_Null
		}		
	}
	else if "`getmodel'" == "bayesrandom" {	
		tempfile bayesfixedest
		local bayesfixedest = subinstr("`bayesfixedest'", ".tmp", ".dta", 1)
		
		if "`cov'" == "freeint" {
			tempvar holder
			my_ncod `holder', oldvar(`studyid')
			
			drop `studyid'
			gen `studyid' = `holder'

			local neoregexpression "`regexpression' i.`studyid'"
		}
		else {
			local neoregexpression "`regexpression'"
		}
		
		qui fitmodel `outcome' `total' if `touse', `modelopts' model(bayesfixed) regexpression(`neoregexpression') ///
		sid(`sid') studyid(`studyid') level(`level') nested(`first') `abnetwork' `comparative' `mpair' cov(`cov') link(`link') p(`p')  ///
		bayesest(`bayesfixedest') inference(`inference') refsampling(`refsampling') `progress' `interaction'
		
		qui estimates store metapreg_Null
		
		capture bayestest model metapreg_modest metapreg_Null
		
		if _rc == 0 {
			mat `bayestest' = r(test)
			local MLdiff = `bayestest'[2,2] - `bayestest'[1,2] //Difference in marginal likelihood
			local postprob = `bayestest'[1,4]  // posterior prob of the model
		}
		
		capture bayesstats ic  metapreg_Null metapreg_modest 
		if _rc == 0 {
			mat `bayesstats' = r(ic)
			local BF = `bayesstats'[2,4] //log Bayes factor
		}
	
		/*if _rc == 0 {
			mat `bayestest' = r(test)
			local BHET = `bayestest'[2,2] - `bayestest'[1,2] //Difference in marginal likelihood
			local P_BHET = `bayestest'[1,4]  // posterior prob of the model
			local DF_BHET = 1
		}*/

		estimates drop metapreg_Null
		rm "`bayesfixedest'" 
	}
	
	//FE
	local TAU21 = 0
	local TAU22 = 0
	local ISQ1 = .
	local ISQ2 = .
	local rho = .
	
	if "`getmodel'" == "random" {
	
		local npar = colsof(`coefmat')
		if "`predcmd'" == "meqrlogit_p" {
			local scalefn "exp"
			local scalepow "2"
		}
		else {
			local scalepow "1"
		}
		
		if "`abnetwork'`cov'" == "" {
			*local TAU21 = `scalefn'(`coefmat'[1, `npar'])^`scalepow' //Between study variance	1
			local TAU21 = `covmat'[1, 1] //Between study variance	1
			local TAU22 = 0
		}
		if "`abnetwork'`cov'" != ""  {
			local TAU21 = `covmat'[1, 1]
			local TAU22 = `covmat'[2, 1]
			/*
			local TAU21 = `scalefn'(`coefmat'[1, `=`npar'-1'])^`scalepow' //Between study variance	1
			local TAU22 = `scalefn'(`coefmat'[1, `npar'])^`scalepow' //Between study variance	2
			*/
		}
		if "`cov'" == "unstructured" {
			/*local TAU21 = `rosevar'[1, 1]
			local TAU22 = `rosevar'[2, 2]
			*/
			local rho = `covmat'[3, 1]/sqrt(`covmat'[2, 1]*`covmat'[1, 1])
			*local rho = `rosevar'[1, 2]/sqrt(`rosevar'[2, 2]*`rosevar'[1, 1])
			  
		}
	}
	
	if "`getmodel'" == "bayesrandom" {
		qui estimate restore metapreg_modest
		
		if `lnvar' {
			local parmsigma2 = `"(exp({lnsigma})^2)"'
			local parmtau2 = `"(exp({lntau})^2)"'
		}
		else if `sd' {
			local parmsigma2 = `"({sigma}^2)"'
			local parmtau2 = `"({tau}^2)"'
		}
		else  {
			local parmsigma2 = `"{sigmasq}"'
			local parmtau2 = `"{tausq}"'
		}
		
		if "`abnetwork'`cov'" == "" {
			qui bayesstats summary `parmtau2'
			mat `covmat' = r(summary)
		
			local TAU21 = `covmat'[1, 4] //Median Between study variance	1
			local TAU22 = 0
		}
		else if "`abnetwork'" != "" |  "`cov'" == "independent" {
			#delimit ;
			qui bayesstats summary `parmtau2' `parmsigma2'
				(isq1:`parmtau2'/(`parmtau2' + `parmsigma2'))
				(isq2:`parmsigma2'/(`parmtau2' + `parmsigma2'))
			;
			#delimit cr	
			mat `covmat' = r(summary) //Take medians
			local TAU21 = `covmat'[1, 4] 
			local TAU22 = `covmat'[2, 4] 
			local ISQ1 = `covmat'[3, 4]*100
			local ISQ2 = `covmat'[4, 4]*100 
		}
		else if (("`cov'" == "commonint") | ("`cov'" == "freeint")) {	
			qui bayesstats summary `parmsigma2'
			mat `covmat' = r(summary)
			local dim :colsof(`covmat')
			mat `covmat' = J(1, `dim', 0) \ `covmat'
			mat rownames `covmat' = `parmtau2' `parmsigma2'
			
			if "`cov'" == "commonint" {
				local TAU21 = 0
			}
			else {
				local TAU21 = .
			}
			local TAU22 = `covmat'[2, 4] //Median Between study variance	2
		}
		else if "`cov'" == "commonslope" {
			qui bayesstats summary `parmtau2'
			mat `covmat' = r(summary)
		
			local TAU21 = `covmat'[1, 4] //Median Between study variance	1
			local TAU22 = 0
			
			local dim :colsof(`covmat')
			mat `covmat' = `covmat' \ J(1, `dim', 0) 
			mat rownames `covmat' = `parmtau2' `parmsigma2'
		}
		else if "`cov'" == "unstructured" {			
			#delimit ;
			qui bayesstats summary 
				{Sigma_1_1} 
				{Sigma_2_2} 
				(rho:{Sigma_2_1}/sqrt({Sigma_2_2}*{Sigma_1_1}))
				(isq1:{Sigma_1_1}/({Sigma_1_1} + {Sigma_2_2}))
				(isq2:{Sigma_2_2}/({Sigma_1_1} + {Sigma_2_2}))
			;
			#delimit cr

			mat `covmat' = r(summary)	//Take medians			
			local TAU21 = `covmat'[1, 4]
			local TAU22 = `covmat'[2, 4]			
			local rho = `covmat'[3, 4]
			local ISQ1 = `covmat'[4, 4]*100
			local ISQ2 = `covmat'[5, 4]*100  
		}
	}

	if (`p' == 0) & (strpos("`getmodel'",  "random") != 0) & ("`pcbnetwork'`mcbnetwork'" == "") & "`link'" == "logit" {
		/*Compute I2*/				
		qui gen `invtotal' = 1/`total'
		qui summ `invtotal' if `touse'
		
		local avgN = r(mean)
		if (strpos("`getmodel'",  "bayes") != 0) {
			qui bayesstats summary {`event':mu}
			mat `coefmat' = r(summary)
			local Esigma = (exp(`TAU21'*0.5 + `coefmat'[1, 4]) + exp(`TAU21'*0.5 - `coefmat'[1, 4]) + 2) * `avgN'
		}
		else {
			local Esigma = (exp(`TAU21'*0.5 + `coefmat'[1, 1]) + exp(`TAU21'*0.5 - `coefmat'[1, 1]) + 2) * `avgN'
		}		
		local ISQ1 = `TAU21'/(`Esigma' + `TAU21')*100	
	}
	/* else */ if ("`abnetwork'`cov'" != "") & ("`cov'" != "commonslope") & (strpos("`getmodel'",  "random") != 0) {
		local ISQ1 = `TAU21'/(`TAU21' + `TAU22')*100
		local ISQ2 = `TAU22'/(`TAU21' + `TAU22')*100		
	}

	if "`inference'" == "frequentist" {
		//Raw estimates in logit/cloglog/loglog scale
		cap freqsummary, event(`event') total(`total') studyid(`studyid') estimates(metapreg_modest) ///
			`interaction' catreg(`catreg') contreg(`contreg') level(`level') model(`getmodel') scimethod(`scimethod')  ///
			varx(`varx') typevarx(`typevarx') by(`by') regexpression(`regexpression') ///
			`mcbnetwork' `comparative' `pcbnetwork' `abnetwork' `stratify' `mpair'  ///
			comparator(`comparator') link(`link')  total(`total') 
			
		mat `rawest' = r(outmatrix)
		mat `exactabsout' = r(exactabsout)
			
		//Conditional ABS
		estp, rawestmat(`rawest') link(`link') scimethod(`scimethod') model(`getmodel')
		mat `absout' = r(outmatrix)
	}
	else {
		//Conditional odds/abs
		cap bayessummary, event(`event') total(`total') studyid(`studyid') estimates(metapreg_modest) ///
			`interaction' catreg(`catreg') contreg(`contreg') level(`level') model(`getmodel') scimethod(`scimethod')  ///
			varx(`varx') typevarx(`typevarx') by(`by') regexpression(`regexpression') ///
			`mcbnetwork' `comparative' `pcbnetwork' `abnetwork' `stratify' `mpair'  ///
			comparator(`comparator') link(`link') cov(`cov')
		
		mat `rawest' = r(loddsout)
		mat `exactabsout' = r(exactabsout)
		mat `absout' = r(absout)
	}
		
	//Population abs
	if "`getmodel'" != "hexact" & "`getmodel'" != "crbetabin"  {
		//simulations ABS
		cap postsim , event(`event') total(`total') orderid(`rid') studyid(`studyid') todo(p) estimates(metapreg_modest) rawest(`rawest') ///
				level(`level')  model(`getmodel')  by(`by') `comparative' `interaction' `abnetwork' `mpair' catreg(`catreg') by(`by') ///
				modeles(`modeles') modellci(`modellci') modeluci(`modeluci') `stratify' ///
				`mcbnetwork' varx(`varx') cov(`cov') nsims(`nsims') link(`link') p(`p')  bayesreps(`bayesreps') 
				
		mat `popabsout' = r(outmatrix)
		
	}

	//R
	local rrsuccess 0
	if "`catreg'" != "" | "`typevarx'" == "i" {
		if "`inference'" == "frequentist" {
			cap freqestr, event(`event') total(`total') studyid(`studyid') estimates(metapreg_modest)  catreg(`catreg') ///
				level(`level') comparator(`comparator') `interaction' scimethod(`scimethod') ///
				varx(`varx') typevarx(`typevarx') by(`by') `mcbnetwork' `pcbnetwork' `mpair' ///
				`comparative' `abnetwork' aliasdesign(`aliasdesign') `stratify' model(`getmodel')  ///
				regexpression(`regexpression') `baselevel' link(`link') inference(`inference') total(`total') cov(`cov')
		}
		else {
			cap  bayesestr , event(`event') catreg(`catreg') ///
				level(`level') comparator(`comparator') `interaction' scimethod(`scimethod') ///
				varx(`varx') typevarx(`typevarx') by(`by') `mcbnetwork' `pcbnetwork'  `mpair' ///
				`comparative' `abnetwork' `stratify' model(`getmodel') cov(`cov') ///
				regexpression(`regexpression') `baselevel' link(`link') sid(`studyid')
		}	
			
		if _rc == 0 {
			if "`getmodel'" != "hexact" {
				mat `rrout' = r(rroutmatrix)
				mat `rdout' = r(rdoutmatrix)
				mat `orout' = r(oroutmatrix)
								
				if "`inference'" == "frequentist" {
					local inltest = r(inltest)
					if "`inltest'" == "yes" {
						mat `nltestRR' = r(nltestRR) //if RR by groups are equal
						mat `nltestRD' = r(nltestRD) //if RD by groups are equal
						mat `nltestOR' = r(nltestOR) //if RR by groups are equal
					}
				}
				if "`getmodel'" != "crbetabin" {
					//simulations
					cap postsim ,  event(`event') total(`total') orderid(`rid') studyid(`studyid') todo(r) estimates(metapreg_modest) rrout(`rrout') orout(`orout') ///
							 level(`level')  model(`getmodel')  by(`by') `comparative' `interaction' `abnetwork' `mpair' ///
							 modeles(`modeles') modellci(`modellci') modeluci(`modeluci') `stratify' ///
							 `baselevel' `mcbnetwork' aliasdesign(`aliasdesign') varx(`varx')  cov(`cov') nsims(`nsims') link(`link') p(`p') bayesreps(`bayesreps')
					
					mat `poprrout' = r(rroutmatrix)
					mat `poprdout' = r(rdoutmatrix)
					mat `poplrrout' = r(lrroutmatrix)
					mat `poporout' = r(oroutmatrix)
					mat `poplorout' = r(loroutmatrix)
					
				}
			}
			else {
				mat `exactorout' = r(exactorout)
				mat `exactlorout' = r(exactlorout)
			}
			
			local rrsuccess 1
		}
	}
	
	//Smooth estimates
	//simulations
	if "`smooth'" != "" {
		postsim, event(`event') total(`total') orderid(`rid') studyid(`studyid') todo(smooth) estimates(metapreg_modest)  ///
			level(`level')  model(`getmodel')  by(`by') `comparative' `abnetwork'  `mcbnetwork' `mpair' aliasdesign(`aliasdesign') ///
			modeles(`modeles') modellci(`modellci') modeluci(`modeluci') outplot(`outplot')	`stratify' ///
			`interaction'  varx(`varx')  cov(`cov') nsims(`nsims') link(`link')  bayesreps(`bayesreps') stat(`stat')
	}
	//===================================================================================
	//Return the matrices
	if "`abnetwork'`cov'" != "" {
		if "`inference'" == "bayesian" {
			if "`cov'" == "unstructured" {
				mat `hetout' = (`MLdiff', `BF', `postprob', `TAU21', `TAU22', `rho', `ISQ1', `ISQ2')
				mat colnames `hetout' = Delta_ML log(BF) Post_Prob tau2 sigma2 rho I2tau I2sigma 
			}
			else {
				mat `hetout' = (`MLdiff', `BF', `postprob', `TAU21', `TAU22', `ISQ1', `ISQ2')
				mat colnames `hetout' = Delta_ML log(BF) Post_Prob tau2 sigma2 I2tau I2sigma 
			}
		}
		else {
			if "`cov'" == "unstructured" {
				mat `hetout' = (`DF_BHET', `BHET' ,`P_BHET', `TAU21', `TAU22', `rho', `ISQ1', `ISQ2')
				mat colnames `hetout' = DF Chisq p tau2 sigma2 rho I2tau I2sigma 
			}
			else {
				mat `hetout' = (`DF_BHET', `BHET' ,`P_BHET', `TAU21', `TAU22', `ISQ1', `ISQ2')
				mat colnames `hetout' = DF Chisq p tau2 sigma2 I2tau I2sigma 
			}
		}
	}
	else {
		if (`p' == 0) & strpos("`model'", "random") !=0 & "`pcbnetwork'`mcbnetwork'" == "" {
			if "`inference'" == "bayesian" {
				mat `hetout' = (`MLdiff', `BF', `postprob', `TAU21', `ISQ1')
				mat colnames `hetout' = Delta_ML log(BF) Post_Prob tau2 I2tau 
			}
			else {
				mat `hetout' = (`DF_BHET', `BHET' ,`P_BHET', `TAU21', `ISQ1')
				mat colnames `hetout' = DF Chisq p tau2 I2tau 
			}
		}
		else {
			if "`inference'" == "bayesian" {
				mat `hetout' = (`MLdiff', `BF', `postprob', `TAU21')
				mat colnames `hetout' = Delta_ML log(BF) Post_Prob tau2 
			}
			else {
			
				mat `hetout' = (`DF_BHET', `BHET' ,`P_BHET', `TAU21')
				mat colnames `hetout' = DF Chisq p tau2 
			}
		}
	}
	
	if strpos("`getmodel'", "betabin") != 0 {
		mat `hetout' = (`DF_BHET', `BHET', `P_BHET', `SIGMA')
		if "`getmodel'" == "crbetabin" {
			mat colnames `hetout' = DF chibar2 p phi
		}
		else {
			mat colnames `hetout' = DF chibar2 p rho
		}
	}
	
	mat rownames `hetout' = Model
	return matrix hetout = `hetout'
	return local inltest = "`inltest'"
	
	cap confirm matrix `covmat'
	if _rc == 0 {
		return matrix covmat = `covmat'
	}
										
	cap confirm matrix `rawest'
	if _rc == 0 {
		return matrix rawest = `rawest'
	}
	
	cap confirm matrix `popabsout'
	if _rc == 0 {
		return matrix popabsout = `popabsout'
	}
	cap confirm matrix `absout'
	if _rc == 0 {
		return matrix absout = `absout'
	}
	cap confirm matrix `absoutp'
	if _rc == 0 {
		return matrix absoutp = `absoutp'
	}
	cap confirm matrix `exactabsout'
	if _rc == 0 {
		return matrix  exactabsout = `exactabsout'
	}
	cap confirm matrix `rrout'
	if _rc == 0 {
		return matrix rrout = `rrout'
	}
	cap confirm matrix `poprrout'
	if _rc == 0 {
		return matrix poprrout = `poprrout'
		return matrix poplrrout = `poplrrout'
	}
	cap confirm matrix `rdout'
	if _rc == 0 {
		return matrix rdout = `rdout'
	}
	cap confirm matrix `poprdout'
	if _rc == 0 {
		return matrix poprdout = `poprdout'
	}
	cap confirm matrix `orout'
	if _rc == 0 {
		return matrix orout = `orout'
	}
	cap confirm matrix `poporout'
	if _rc == 0 {
		return matrix poporout = `poporout'
		return matrix poplorout = `poplorout'
	}
	cap confirm matrix `exactorout'
	if _rc == 0 {
		return matrix exactorout = `exactorout'
		return matrix exactlorout = `exactlorout'
	}
	
	if "`inltest'" == "yes" {
		return matrix nltestRR = `nltestRR'
		return matrix nltestRD = `nltestRD'
		return matrix nltestOR = `nltestOR'
	}
	
	return scalar mdf = `mdf'
	return local model = "`getmodel'"
	return local rrsuccess = "`rrsuccess'"

end

/**************************************************************************************************
								MCPREG - BREGRESSIONS FOR MODEL COMPARISON
**************************************************************************************************/
capture program drop mcpreg
program define mcpreg, rclass
	#delimit ;

	syntax varlist(min=3) [if] [in], sid(varname)  [
		regexpression(string) nu(string) baselevel(passthru) 
		regressors(varlist) varx(varname)  typevarx(string) comparator(varname) catreg(varlist) contreg(varlist)
		scimethod(string) cov(string)
		level(integer 95)
		DP(integer 2)
		progress
		model(string) modelopts(string) command0(string)
		noMC noCONstant
		interaction	
		comparative mcbnetwork pcbnetwork abnetwork mpair
		link(string) inference(string) refsampling(string) 
		bayesest(string asis) 
			*];

	#delimit cr
	marksample touse, strok 
	
	tempvar event nonevent total
	
	tempname coefmat coefvar testlr V logodds absout absoutp rrout orout nltestRR nltestOR  ///
			 hetout mctest absexact newobs matgof popabsout poprrout poporout poplorout rosevar rawvar rawest bayesstats bayestest
			 
	tokenize `varlist'
	
	local event = "`1'"
	local nonevent = "`2'"
	local total = "`3'"
	//fit the model
	if "`progress'" != "" {
		local echo noi
	}
	else {
		local echo qui
	}
	//Just initialize
	gettoken first confounders : regressors
	local p: word count `regressors'
	
	if "`mcbnetwork'`pcbnetwork'" != "" {		
		tokenize `regexpression'
		local one "`1'"
		local two "`2'"
		local three "`3'"
		
		if "`interaction'" != "" {
			tokenize `one', parse("#")
			tokenize `1', parse(".")
			local ipair "`3'"
			
			tokenize `two', parse(".")
			local index "`3'"
		}
		else {
			tokenize `two', parse(".")
			local ipair "`3'"
		
			tokenize `three', parse(".")
			local index "`3'"
		}
	}
	//Which outcome to model
	if "`link'" == "loglog" & "`model'" != "crbetabin" {
		local outcome "`nonevent'"
	}
	else {
		local outcome "`event'"
	}
		
	if "`inference'" == "frequentist" {
		qui estimates restore metapreg_modest
		qui estat ic
		mat `matgof' = r(S)
		local BIC =  `matgof'[1, 6]
	}
	local redindex 0
	
	if "`abnetwork'`interaction'" !="" {
		if "`mcbnetwork'`pcbnetwork'" != "" {
			local confariates "`comparator'"	
		}
		else {
			local confariates "`confounders'"
		}
	}
	else {
		local confariates "`regressors'"
	}
	
	local initial 1
	foreach c of local confariates {
		local nureduced	
		if ("`interaction'" != "" & "`pcbnetwork'`mcbnetwork'" != "")  {
				local omterm = "`c'*`ipair'"
				
				gettoken start end : regexpression
				
				local nureduced "mu i.`ipair' `end'"
				
		}
		else {						
			foreach term of local regexpression {
				if "`interaction'" != "" {
					if strpos("`term'", "`c'#") != 0 & strpos("`term'", "`first'") != 0 {
						local omterm = "`c'*`first'"
					}
					else {
						local nureduced "`nureduced' `term'"
					}
				}
				else {
					if "`model'" == "cbbetabin" {
						if ("`term'" == "i.`c'#c.mu")|("`term'" == "c.`c'#c.mu") {
							local omterm = "`c'"
						} 
						else {
							local nureduced "`nureduced' `term'"
						}
					}
					else {
						if ("`term'" == "i.`c'")|("`term'" == "c.`c'")|("`term'" == "`c'") {
							local omterm = "`c'"
						} 
						else {
							local nureduced "`nureduced' `term'"
						}
					}
				}
			}
		}
		
		if "`omterm'" == "`first'" {
			local newcov 
		}
		else {
			local newcov "`cov'"
		}
		
		`echo' fitmodel `outcome' `total' if `touse', `modelopts' model(`model')  ///
			regexpression(`nureduced') sid(`sid') studyid(`studyid') level(`level')  nested(`first') `abnetwork'  `comparative' `mpair' ///
			cov(`newcov') link(`link') p(`p') inference(`inference') ///
			refsampling(`refsampling') bayesest(`bayesnullest') `progress'
		
		if "`inference'" == "frequentist" {
			qui estat ic
			mat `matgof' = r(S)
			local BICmc = `matgof'[1, 6]
		}
		estimates store metapreg_Null
		
		if "`inference'" == "frequentist" {
			//Do not compare if different commands
			
			//LR test the model
			capture lrtest metapreg_modest metapreg_Null
			if _rc == 0 {
				local lrp :di %10.`dp'f chi2tail(r(df), r(chi2))
				local lrchi2 = r(chi2)
				local lrdf = r(df)
			}
			else {
				local lrp = .
				local lrchi2 = .
				local lrdf = .
			}
		}
		else {
			capture bayestest model metapreg_modest metapreg_Null
		
			if _rc == 0 {
				mat `bayestest' = r(test)
				local MLdiff = `bayestest'[2,2] - `bayestest'[1,2] //Difference in marginal likelihood
				local postprob = `bayestest'[2,4]  // posterior prob of the model
				
			}
			
			capture bayesstats ic  metapreg_modest metapreg_Null 
			if _rc == 0 {
				mat `bayesstats' = r(ic)
				local DICdiff = `bayesstats'[1,2] - `bayesstats'[2,2] //Difference in DIC
				local BF = `bayesstats'[2,4]
			}
		}
		estimates drop metapreg_Null
		
		if `initial' == 1  {
			if "`inference'" == "frequentist" {
					mat `mctest' = [`lrchi2', `lrdf', `lrp', `=`BIC' -`BICmc'']
			}
			else {
				mat `mctest' = [`MLdiff', `BF', `postprob', `DICdiff']
			}
		}
		else {
			if "`inference'" == "frequentist" {
				mat `mctest' =  `mctest' \ [`lrchi2', `lrdf', `lrp', `=`BIC' -`BICmc'']
			}
			else {
				mat `mctest' =  `mctest' \ [`MLdiff', `BF', `postprob', `DICdiff']
			}
		}
		local rownameslr "`rownameslr' `omterm'"
		
		local initial 0
	}
	//Ultimate null model
	if (`p' > 1 & "`abnetwork'" == "") | (`p' > 2 & "`abnetwork'" != "")  {
		
		if "`abnetwork'" != ""  {
			local regexpression "ibn.`first'"
		}
		else if "`pcbnetwork'`mcbnetwork'" != "" {
			local regexpression "mu i.`ipair' i.`index'"
		}
		else {
			local regexpression "mu"
		}
		
		`echo' fitmodel `outcome' `total' if `touse', `modelopts' model(`model') regexpression(mu) ///
			sid(`sid') studyid(`studyid') level(`level')  nested(`first') `abnetwork' `comparative' `mpair' link(`link') p(`p')  ///
			inference(`inference') refsampling(`refsampling') bayesest(`bayesnullest') `progress'
		
		if "`inference'" == "frequentist" {
			qui estat ic
			mat `matgof' = r(S)
			local BICmc = `matgof'[1, 6]
		}
		
		estimates store metapreg_Null
		if "`inference'" == "frequentist" {
			capture lrtest metapreg_modest metapreg_Null
			
			if _rc == 0 {
				local lrchi2 = r(chi2)
				local lrdf = r(df)
				local lrp :di %10.`dp'f r(p)
			}
			else {
				local lrp = .
				local lrchi2 = .
				local lrdf = .
			}
			
			local lrchi2 = r(chi2)
			local lrdf = r(df)
			local lrp :di %10.`dp'f r(p)
		}
		else {
			capture bayestest model metapreg_modest metapreg_Null
		
			if _rc == 0 {
				mat `bayestest' = r(test)
				local MLdiff = `bayestest'[2,2] - `bayestest'[1,2] //Difference in marginal likelihood
				local postprob = `bayestest'[2,4]  // posterior prob of the model
			}
			
			capture bayesstats ic metapreg_modest metapreg_Null
			if _rc == 0 {
				mat `bayesstats' = r(ic)
				local DICdiff = `bayesstats'[1,2] - `bayesstats'[2,2] //Difference in DIC
				local BF = `bayesstats'[2,4]
			}
		}
		
		estimates drop metapreg_Null
		if "`inference'" == "frequentist" {
			mat `mctest' = `mctest' \ [`lrchi2', `lrdf', `lrp', `=`BIC' -`BICmc'']
		}
		else {
			mat `mctest' =  `mctest' \ [`MLdiff', `BF', `postprob', `DICdiff']
		}
		local rownameslr "`rownameslr' All"
	}
	
	mat rownames `mctest' = `rownameslr'
	if "`inference'" == "frequentist" {
		mat colnames `mctest' =  chi2 df p Delta_BIC
	}
	else {
		mat colnames `mctest' =  Delta_ML log(BF) Post_Prob Delta_BIC
	}
	
	cap confirm matrix `mctest'
	if _rc == 0 {
		return matrix mctest = `mctest'
	}
end
/*+++++++++++++++++++++++++	SUPPORTING FUNCTIONS: myncod +++++++++++++++++++++++++
								Decode by order of data
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/	
cap program drop my_ncod
program define my_ncod
	syntax newvarname(gen), oldvar(varname)
	
	qui {
		cap confirm numeric var `oldvar'
		tempvar by_num 
		
		if _rc == 0 {
			decode `oldvar', gen(`by_num')
			drop `oldvar'
			rename `by_num' `oldvar'
		}

		//Remove spaces if present
		replace `oldvar' =  strtrim(`oldvar') 
		
		* The _by variable is generated according to the original
		* sort order of the data, and not done alpha-numerically

		qui count
		local N = r(N)
		cap drop `varlist'
		gen `varlist' = 1 in 1
		local lab = `oldvar'[1]
		cap label drop `oldvar'
		if "`lab'" != ""{
			label define `oldvar' 1 "`lab'"
		}
		local found1 "`lab'"
		local max = 1
		forvalues i = 2/`N'{
			local thisval = `oldvar'[`i']
			local already = 0
			forvalues j = 1/`max'{
				if "`thisval'" == "`found`j''"{
					local already = `j'
				}
			}
			if `already' > 0{
				replace `varlist' = `already' in `i'
			}
			else{
				local max = `max' + 1
				replace `varlist' = `max' in `i'
				local lab = `oldvar'[`i']
				if "`lab'" != ""{
					label define `oldvar' `max' "`lab'", modify
				}
				local found`max' "`lab'"
			}
		}

		label values `varlist' `oldvar'
		label copy `oldvar' `varlist', replace
	}
end

	
/*+++++++++++++++++++++++++	SUPPORTING FUNCTIONS: FITMODEL +++++++++++++++++++++++++
							Fit the regression model
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
	 
cap program drop fitmodel
program define fitmodel, rclass
	#delimit ;
	syntax varlist [if] [in], 
		[model(string) progress regexpression(string) sid(varname) studyid(varname) p(string) interaction
		level(integer 95) mpair mcbnetwork pcbnetwork abnetwork general comparative nested(string) cov(string) link(string) 
		bayesest(string asis) inference(string) 
		refsampling(string) 
		feprior(string asis)
		varprior(string asis)
		nchains(integer 3)
		thinning(integer 5) /*5*/
		burnin(integer 5000) /*5000*/
		mcmcsize(integer 3000) /*3000*/
		rseed(integer 1)		
		*]
	;
	#delimit cr
	marksample touse, strok 
	
	local modelopts `"`options'"'
	
	if "`progress'" != "" {
		local echo noisily
	}
	
	tokenize `varlist'
	local events = "`1'"
	local total = "`2'"
	
	if "`cov'" != "" {	
		if ("`comparative'`mpair'" != "" ) & ("`cov'" != "commonslope")  {
			local slope "2.`nested'"
		}

		local varx "`nested'"
	}
	else {
		local varx
	}
	
	if (strpos("`regexpression'", "mu") != 0) | ("`abnetwork'" != "" & "`model'" == "random"){
		local intercept "mu"
	}
	
	if "`comparative'`mpair'" != "" & ("`cov'" == "commonint" | "`cov'" == "freeint")  { 
		local intercept
	}

	if "`inference'" == "frequentist" {
		if "`cov'" != "" {
			if "`cov'" == "unstructured" {
				local cov "cov(`cov')"
			}
			else {				
				//if free intercepts, then remove mu and base for study
				if "`cov'" == "freeint" {			
					fvset base none `sid'	
					gettoken mu regexpression: regexpression
					local regexpression = "i.`sid' `regexpression'"
				}
				local cov
			}
		}
		else {
			local cov
		}
	
		if "`abnetwork'" != ""  {
			fvset base none `nested'
			
			if "`model'" == "random" {
				local nested = `"|| (`nested': mu, noconstant)"'
			}
		}
		else {
			local nested
		}
		//Specify the engine
		if "`link'" == "logit" {	
			if _caller() >= 16 {
				local fitcommand "melogit"
			}
			else {
				local fitcommand "meqrlogit"
			}
		}
		else if "`link'" == "log" {
			local fitcommand "mepoisson"
		} 
		else {
			local fitcommand "mecloglog"
		}
	}
	
	if "`inference'" == "bayesian" {
		if  "`bayesest'" != "" {
			local saving = `"saving(`bayesest', replace)"'
		}
				
		if "`feprior'" == "" {
			local feprior = "normal(0, 10)"
			*local feprior = "uniform(-10, 10)"
		}
		
		if "`link'" == "log" {
			local lkll "likelihood(poisson, exposure(`total'))"
		}
		
		if "`link'" == "logit" {
			local lkll "likelihood(binomial(`total'))"	
		}
	}
	
	//fit bayesian fixed
	if "`model'" == "bayesfixed" {
		if "`link'" == "log" {
			local inits "init1({`events':}  runiform(0, 0.1)) init2({`events':}  runiform(-0.1, 0)) init3({`events':}  runiform(-0.01, 0.01))"
		}
		
		//if saturated, then remove mu and base for study
		if strpos("`regexpression'" , "`studyid'") != 0 {			
			fvset base none `studyid'	
			gettoken mu regexpression: regexpression
		}
		
		local priorfe = `" prior({`events':`regexpression'}, `feprior')"'
		local blockfe = `"block({`events':`regexpression'})"'
			
		#delim ;
		capture `echo' bayesmh `events' `regexpression' if `touse', noconstant  `lkll'  /*likelihood(binomial(`total'))*/ 
			`priormu' `priorfe' 
			`blockfe' `blockmu' `inits'
			nchains(`nchains') thinning(`thinning') burnin(`burnin') mcmcsize(`mcmcsize') rseed(`rseed')
			`saving'
			;
		#delimit cr 
		local success = _rc			
	}
	
	//fit bayesian re
	local lnvar = 0
	local sd = 0
	if "`model'" == "bayesrandom" {
			
		if "`varprior'" == "" {
			if "`cov'" == "unstructured" {
				local varprior = "iwishart(2,3,I(2))"
				*local varprior = "wishart(2,3,I(2))"
			}
			else {
				local varprior "igamma(0.01, 0.01)"
				*local varprior "igamma(0.1, 0.1)"
				
				if "`cov'" != "commonslope" &  "`comparative'" != "" {
					local expsigma = `"{sigmasq}"'
					local parmsigma = `"{sigmasq}"'
					local priorsigma =  "`varprior'"
				}
				if "`cov'" != "commonint" & "`cov'" != "freeint" {
					local exptau = `"{tausq}"'
					local parmtau = `"{tausq}"'
					local priortau =  "`varprior'"
				}
			}
		}
		else {
			if ("`cov'" != "unstructured") & ("`cov'" != "") {
				if "`cov'" != "commonslope" {
					local expsigma = `"{sigmasq}"'
					local parmsigma = `"{sigmasq}"'
					local priorsigma =  "`varprior'"
				}
				if "`cov'" != "commonint" & "`cov'" != "freeint" {				
					local exptau = `"{tausq}"'
					local parmtau = `"{tausq}"'
					local priortau =  "`varprior'"
				}
			}
		}

		fvset base none `sid'
			
		gettoken mu variableterms: regexpression
						
		if "`cov'" != "" {
			tokenize `variableterms'
			macro shift
			local variableterms `*'
			
			if "`interaction'" != "" {
				//Rewrite the regexpression if interaction terms are present; discard the main terms and add hash to the interaction term
				local neoterms
				foreach term of local variableterms {
					if strpos("`term'", "#") != 0 {					
						local term = subinstr("`term'", "#", "##", 1)					
						local neoterms "`neoterms' `term'"
					}
				}
				local variableterms "`neoterms'"
			}
			
			if 	"`cov'" == "commonint" {
				local neoregexpression = `"mu i.`sid'#2.`varx' `variableterms'"'
			}
			else if "`cov'" == "commonslope" {
				if "`interaction'" == "" { 
					local neoregexpression = `"i.`sid' 2.`varx' `variableterms'"'
				}
				else {
					local neoregexpression = `"i.`sid' `variableterms'"'
				}
			}
			else {
				//freeint, independent or unstructured
				local neoregexpression = `"i.`sid' i.`sid'#2.`varx' `variableterms'"'
			}
			
			//Assign re priors
			if "`cov'" == "commonint"  |  "`cov'" == "freeint"{
				local priorslopes = `"prior({`events':i.`sid'#2.`varx'}, normal({`events':2.`varx'}, `expsigma'))"'
				local priorsigma = `"prior(`parmsigma', `priorsigma')"'
			}
			else if "`cov'" == "commonslope" {
				local priorsid = `"prior({`events':i.`sid'}, normal({`events':mu}, `exptau'))"'			
				local priortau = `"prior(`parmtau', `priortau')"'
			}
			else if "`cov'" == "independent" {
				local priorslopes = `"prior({`events':i.`sid'#2.`varx'}, normal({`events':2.`varx'}, `expsigma'))"'
				local priorsid = `"prior({`events':i.`sid'}, normal({`events':mu}, `exptau'))"'	
				
				local priorsigma = `"prior(`parmsigma', `priorsigma')"'
				local priortau = `"prior(`parmtau', `priortau')"'
			}
			else if "`cov'" == "unstructured" {
				local priorre = `"prior({`events':i.`sid' i.`sid'#2.`varx'}, mvnormal(2, {`events':mu}, {`events':2.`varx'}, {Sigma, matrix}))"'
				local priorvarcov = `"prior({Sigma, matrix}, `varprior')"'
			}
						
			if "`interaction'" == "" {
				//prior
				local priorvarx = `"prior({`events':2.`varx'}, `feprior')"'
				
				//block
				local blockvarx = `"block({`events':2.`varx'})"'
			}
			
			if "`cov'" == "commonint" | "`cov'" == "freeint"  {
				if	`refsampling' == 1 {
					local blockslopes = `"block({`events':i.`sid'#2.`varx'}, split)"'
				}
				else if `refsampling' == 2  {
					local blockslopes = `"block({`events':i.`sid'#2.`varx'}, reffects)"'
				}
				if strpos("`varprior'", "gamma") != 0 { 
					local blocksigma = `"block(`parmsigma', gibbs)"'
				}
				else {
					local blocksigma = `"block(`parmsigma')"'
				}
			}			
			else if "`cov'" == "commonslope" {
				if	`refsampling' == 1 {
					local blocksid = `"block({`events':i.`sid'}, split)"'
				}
				else if `refsampling' == 2 {
					local blocksid = `"block({`events':i.`sid'}, reffects)"'
				}
				else if `refsampling' == 3 {
					local blocksid = `"reffects(`sid')"'
				}
				if strpos("`varprior'", "gamma") != 0 { 
					local blocktau = `"block(`parmtau', gibbs)"'
				}
				else {
					local blocktau = `"block(`parmtau')"'
				}
			}
			else if "`cov'" == "independent" {
				if	`refsampling' == 1 {
					local blockslopes = `"block({`events':i.`sid'#2.`varx'}, split)"'
				}
				else if `refsampling' == 2  {
					local blockslopes = `"block({`events':i.`sid'#2.`varx'}, reffects)"'
				}
				if strpos("`varprior'", "gamma") != 0 { 
					local blocksigma = `"block(`parmsigma', gibbs)"'
					local blocktau = `"block(`parmtau', gibbs)"'
				}
				else {
					local blocksigma = `"block(`parmsigma')"'
					local blocktau = `"block(`parmtau')"'
				}
			}
			else if "`cov'" == "unstructured" {
				local blockvarcov = `"block({Sigma, matrix}, gibbs)"'
			}
		}
		//general design
		if "`cov'" == "" & strpos("`model'", "random") != 0 {
			*tokenize `variableterms'
			*macro shift
			*local variableterms `*'
			
			if "`interaction'" != "" {
				//Rewrite the regexpression if interaction terms are present; discard the main terms and add hash to the interaction term
				local neoterms
				foreach term of local variableterms {
					if strpos("`term'", "#") != 0 {					
						local term = subinstr("`term'", "#", "##", 1)					
						local neoterms "`neoterms' `term'"
					}
				}
				local variableterms "`neoterms'"
			}
			
			local neoregexpression = `"i.`sid' `variableterms'"'
			
			//Assign re priors
			local priorsid = `"prior({`events':i.`sid'}, normal({`events':mu}, `exptau'))"'			
			local priortau = `"prior(`parmtau', `priortau')"'
													
			if	`refsampling' == 1 {
				local blocksid = `"block({`events':i.`sid'}, split)"'
			}
			else if `refsampling' == 2 {
				local blocksid = `"block({`events':i.`sid'}, reffects)"'
			}
			else if `refsampling' == 3 {
				local blocksid = `"reffects(`sid')"'
			}
			if strpos("`varprior'", "gamma") != 0 { 
				local blocktau = `"block(`parmtau', gibbs)"'
			}
			else {
				local blocktau = `"block(`parmtau')"'
			}
		}

		//assign fe priors			
		if "`variableterms'" != "" {
			local priorfe = `"prior({`events':`variableterms'}, `feprior')"'
		}
		
		if "`cov'" != "freeint" {
			local priormu = `"prior({`events':mu}, `feprior')"'
		}
		if "`cov'" == "freeint" {
			local priormu = `"prior({`events':i.`sid'}, `feprior')"'
		}
		
		//block
		if "`cov'" == "unstructured" {
			local blockmu = `"block({`events':mu})"'
		}
		else if "`cov'" == "freeint" {
			local blockmu = `"block({`events':i.`sid'})"'
		}
		else  { 
			local blockmu = `"block({`events':mu})"'
		}
		
		if "`variableterms'" != "" {
			local blockfe = `"block({`events':`variableterms'})"'
		}
		
		if strpos("`modelopts'", "adapt") == 0 {
			local adapt "adaptation(maxiter(50))"
		}
		if "`link'" == "log" {
			local inits "init1({`events':}  runiform(0, 0.1)) init2({`events':}  runiform(-0.1, 0)) init3({`events':}  runiform(-0.01, 0.01))"
		}
						
		#delim ;
		capture `echo' bayesmh `events' `neoregexpression' if `touse', noconstant `lkll' /*likelihood(binomial(`total')) */ 
			`priorre' `priorvarcov' `priorslopes' `priorsid'  `priormu' `priorvarx' `priortau' `priorsigma' `priorfe'
			`blockfe' `blockvarx' `blockmu' `blockslopes' `blocksid' `blocksigma' `blocktau' `inits'
			nchains(`nchains') thinning(`thinning') burnin(`burnin') mcmcsize(`mcmcsize') rseed(`rseed')
			`saving' `adapt' `modelopts'
			;
		#delimit cr 
		local success = _rc	
	}

	/*===========================Frequentist models===========================================*/
	//Fit cbbetabin - common beta beta-binomial
	if "`model'" == "cbbetabin" {
		//Default iterations
		if strpos(`"`modelopts'"', "iterate") == 0  {
			local iterate = `"iterate(100)"'
		}
			
		qui xtset `sid'
		
		capture `echo' xtnbreg `events' `regexpression' if `touse', fe `modelopts' `iterate' l(`level')	
		local success = _rc	
	}
	
	//Fit crbetabin - common rho beta-binomial
	if "`model'" == "crbetabin" {
		//Default iterations
		if strpos(`"`modelopts'"', "iterate") == 0  {
			local iterate = `"iterate(100)"'
		}
				
		if "`comparative'" != "" {
			local nterms = wordcount("`regexpression'")
			if `nterms' > 1 {
				local term2 : word 3 of `regexpression'
				if strpos("`term2'", ".`studyid'") == 2 {
					qui xtset `sid'				
					local regexpression = subinword("`regexpression'", "`term2'", "`sid'", 1)
				}
			}
		}
		
		//if saturated, then remove mu and base for study
		if strpos("`regexpression'" , "`studyid'") != 0 & "`comparative'" == "" {			
			fvset base none `studyid'	
			gettoken mu regexpression: regexpression
		}
				
		capture `echo' betabin `events' `regexpression' if `touse', noconstant n(`total') link(`link') `modelopts' `iterate' l(`level')	
		local success = _rc
	}
			
	//Fit the FE model
	if ("`model'" == "fixed") |("`model'" == "hexact"){
		
		//if saturated, then remove mu and base for study : Why?
		/*if (strpos("`regexpression'" , "`studyid'") != 0) &  ("`studyid'" != "") {			
			fvset base none `studyid'	
			gettoken mu regexpression: regexpression
		}*/
		
		if "`link'" == "logit" { 
			capture `echo' binreg `events' `regexpression' if `touse', noconstant n(`total') ml `modelopts' l(`level')
		}
		else if "`link'" == "log" { 
			*capture `echo' glm `events' `regexpression' if `touse', noconstant family(binomial `total') link(log) ml `modelopts' l(`level')
			*capture `echo' glm `events' `regexpression' if `touse', noconstant family(poisson) exposure(`total') link(log) ml `modelopts' l(`level')
			capture `echo' poisson `events' `regexpression' if `touse', noconstant  exposure(`total') `modelopts' l(`level')	
		}
		else {
			capture `echo' glm `events' `regexpression' if `touse', noconstant family(binomial `total') link(cloglog) ml `modelopts' l(`level')	
		}
		local success = _rc
	}
	
	//Fit the ME model
	if ("`model'" == "random") {
		if (strpos(`"`modelopts'"', "intpoi") == 0) & (strpos(`"`modelopts'"', "lapl") == 0)  {
			qui count if `touse'
			if `=r(N)' < 7 {
				local ipoints = `"intpoints(`=r(N)')"'
			}
		}
		else {
			qui count if `touse'
			local nobs =r(N)
			
			local oldopts `"`modelopts'"'
			local modelopts
			while `"`oldopts'"' != "" {
				gettoken first oldopts : oldopts
				if strpos(`"`first'"', "intpoi")!= 0  {
					local b1 = strpos(`"`first'"', "(") + 1
					local b2 = strpos(`"`first'"', ")") 
					local oldpoints = substr(`"`first'"', `b1', `=`b2'-`b1'')
					if `oldpoints' > `nobs' {
						local ipoints = `"intpoints(`nobs')"'
					}
					local first
				}
				local modelopts "`modelopts' `first'"
				if "`first'" == "" {
					local local modelopts "`modelopts' `oldopts'"
					continue, break
				}
			}
		}
		//Default iterations
		if strpos(`"`modelopts'"', "iterate") == 0  {
			if "`fitcommand'" == "meqrlogit" {
				local iterate = `"iterate(30)"'
			}
			else {
				local iterate = `"iterate(100)"'
			}
		}
		
		//First trial
		local try = 1
		if "`link'" == "log" {
			#delim ;
			capture `echo' `fitcommand' (`events' `regexpression' if `touse', noconstant exposure(`total'))|| 
			  (`sid': `slope'  `intercept', `cov' noconstant) `nested',
			  `ipoints' `modelopts' l(`level') `iterate';
			#delimit cr 

		}
		else {
			#delim ;
			capture `echo' `fitcommand' (`events' `regexpression' if `touse', noconstant)|| 
			  (`sid': `slope'  `intercept', `cov' noconstant) `nested' ,
			  binomial(`total') `ipoints' `modelopts' l(`level') `iterate';
			#delimit cr 
		}
		
		local success = _rc
		local converged = e(converged)
		
		//Try dnumerical and intmethod(gh) second time
		if  ("`fitcommand'" == "melogit") & (`success' != 0) & strpos(`"`modelopts'"', "dnumerical") == 0 & strpos(`"`modelopts'"', "intme") == 0 {
			
			local ++try
			#delim ;
			capture `echo' `fitcommand' (`events' `regexpression' if `touse', noconstant)|| 
			  (`sid': `slope' `intercept', `cov' noconstant) `nested' ,
			  binomial(`total') `ipoints' dnumerical intmethod(gh) `modelopts' l(`level') `iterate';
			#delimit cr 
			
			local success = _rc
			local converged = e(converged)
		}
		
		//Got to meqrlogit if melogit fails
		if  ("`fitcommand'" == "melogit") & (`success' != 0)  {
			local fitcommand = "meqrlogit"
			local iterate = `"iterate(30)"'
			
			local ++try
			#delim ;
			capture `echo' `fitcommand' (`events' `regexpression' if `touse', noconstant)|| 
			  (`sid': `slope' `intercept', `cov' noconstant) `nested' ,
			  binomial(`total') `ipoints' `modelopts' l(`level') `iterate';
			#delimit cr 
			
			local success = _rc
			local converged = e(converged)
		}
		
		if (`success' != 0) & ("`fitcommand'" == "meqrlogit") & (strpos(`"`modelopts'"', "from") == 0) {
			//First fit laplace to get better starting values
			noi di _n"*********************************** ************* ***************************************" 
			noi di as txt _n "Just a moment - Obtaining better initial values "
			noi di   "*********************************** ************* ***************************************" 
			local lapsuccess 1
			
			local ++try	
			#delim ;
			capture `echo'  `fitcommand' (`events' `regexpression' if `touse', noconstant)|| 
				(`sid': `slope' `intercept', `cov' noconstant) `nested' ,
				binomial(`total') laplace l(`level') `iterate';
			#delimit cr 
			
			local lapsuccess = _rc //0 is success
			local converged = e(converged)
			
			if `lapsuccess' == 0 {
				qui estimates table
				tempname initmat
				mat `initmat' = r(coef)

				local ninits = rowsof(`initmat')
				forvalues e = 1(1)`ninits' {
					local init = `initmat'[`e', 1]
					if `init' != .b {
						if `e' == 1 {
							local inits = `"`init'"'
						}
						else {
							local inits = `"`inits', `init'"'
						}
					}
				}
				mat `initmat' = (`inits')
			
				local inits = `"from(`initmat', copy)"'
				
					
				//second trial with initial values
				local ++try
				#delim ;
				capture `echo'  `fitcommand' (`events' `regexpression' if `touse', noconstant)|| 
				  (`sid': `slope' `intercept', `cov' noconstant) `nested' ,
				  binomial(`total') `ipoints'  `inits' l(`level') `iterate';
				#delimit cr 
				
				local success = _rc
				local converged = e(converged)
			}
		}
		
		/*//Try to refineopts 3 times
		if strpos(`"`modelopts'"', "refineopts") == 0 & ("`fitcommand'" == "meqrlogit") {
			local try = 1
			while `try' < 3 & `converged' == 0 {
			
				#delim ;					
				capture noisily  `fitcommand' (`events' `regexpression' if `touse', noconstant)|| 
					(`sid': `varx' , `cov') `nested' ,
					binomial(`total') `ipoints'  l(`level') refineopts(iterate(`=10 * `try'')) `iterate';
				#delimit cr 
				
				local success = _rc
				local converged = e(converged)
				local try = `try' + 1
			}
		}*/
		
		*Try matlog + refineopts if still difficult
		if (strpos(`"`modelopts'"', "matlog") == 0) & ("`fitcommand'" == "meqrlogit") & ((`converged' == 0) | (`success' != 0)) {
			if strpos(`"`modelopts'"', "refineopts") == 0 {
				local refineopts = "refineopts(iterate(50))"
			}
			local ++try
			#delim ;
			capture `echo'  `fitcommand' (`events' `regexpression' if `touse', noconstant )|| 
				(`sid': `slope' `intercept', `cov' noconstant) `nested' ,
				binomial(`total') `ipoints'  l(`level') `refineopts' matlog `iterate';
			#delimit cr
			
			local success = _rc 
			
			local converged = e(converged)
		}
		
		*Try laplace if not for other commands
		if (`success' != 0) & ("`fitcommand'" != "meqrlogit") & (strpos(`"`modelopts'"', "laplac") == 0) {
			#delim ;
			capture `echo'  `fitcommand' (`events' `regexpression' if `touse', noconstant )|| 
				(`sid': `slope'  `intercept', `cov' noconstant) `nested' ,
				binomial(`total') `modelopts' l(`level') intmethod(laplace) `iterate';
			#delimit cr
			
			local success = _rc 
			local converged = e(converged)		
		}
	}
	//Revert to FE if ME fails
	if (`success' != 0) & ("`model'" == "random") {	
		if "`link'" == "logit" { 
			capture `echo' binreg `events' `regexpression' if `touse', noconstant n(`total') ml `modelopts' l(`level')
		}
		else if "`link'" == "log"   {
			capture `echo' glm `events' `regexpression' if `touse', noconstant family(poisson) exposure(`total') link(log) ml `modelopts' l(`level')	
		}
		else {
			capture `echo' glm `events' `regexpression' if `touse', noconstant family(binomial `total') link(cloglog) ml `modelopts' l(`level')	
		}
		local success = _rc
		local model "fixed"
	}
	*If not converged, exit and offer possible solutions
	if `success' != 0 {
		di as error "Model fitting failed"
		di as error "Try fitting a simpler model or better model option specifications"
		exit `success'
	}

	return local model "`model'"
	return local lnvar = "`lnvar'"
	return local sd = "`sd'"
	/*if "`model'" == "hexact" {
		return matrix absexact = `absexact'
	}*/
end

	/*+++++++++++++++++++++++++	SUPPORTING FUNCTIONS: metadta_PROPCI +++++++++++++++++++++++++
								CI for proportions
	++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
	cap program drop metapreg_propci
	program define metapreg_propci

		syntax varlist [if] [in], p(name) se(name)lowerci(name) upperci(name) [icimethod(string) level(real 95)]
		
		qui {	
			tokenize `varlist'
			gen `p' = .
			gen `lowerci' = .
			gen `upperci' = .
			gen `se' = .
			
			count `if' `in'
			forvalues i = 1/`r(N)' {
				local N = `1'[`i']
				local n = `2'[`i']

				cii proportions `N' `n', `icimethod' level(`level')
				
				replace `p' = r(proportion) in `i'
				replace `lowerci' = r(lb) in `i'
				replace `upperci' = r(ub) in `i'
				replace `se' =  r(se) in `i'
			}
		}
	end
/*+++++++++++++++++++++++++	SUPPORTING FUNCTIONS: LONGSETUP +++++++++++++++++++++++++
							Transform data to long format
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
cap program drop longsetup
program define longsetup

syntax varlist, rid(name) [assignment(name) event(name) total(name) idpair(name) mpair  mcbnetwork pcbnetwork abnetwork general comparative panelize]

	qui {
	
		tokenize `varlist'
		
		if "`mcbnetwork'" != "" {		
			/*4 variables per study : a b c d*/
			gen `event'1 = `1' + `2'  /* a + b */
			gen `event'0 = `1' + `3'  /* a + c */
			gen `total'1 = `1' + `2' + `3' + `4'  /* n */
			gen `total'0 = `1' + `2' + `3' + `4'  /* n */
			gen `assignment'1 = `5'
			gen `assignment'0 = `6'
		}
		if "`mpair'" != "" {		
			/*4 variables per study : a b c d*/
			gen `event'1 = `1' + `2'  /* a + b */
			gen `event'0 = `1' + `3'  /* a + c */
			gen `total'1 = `1' + `2' + `3' + `4'  /* n */
			gen `total'0 = `1' + `2' + `3' + `4'  /* n */
		}
		else if "`pcbnetwork'" != ""  {
			/*3 variables per study : n1 n2 N*/
			gen `event'1 = `1'  /* n1 */
			gen `event'0 = `2'  /* n2 */
			gen `total'1 = `3'  /* N */
			gen `total'0 = `3'  /* N */
			gen `assignment'1 = `4'
			gen `assignment'0 = `5'
		}
		else if "`panelize'" != "" {
			gen `event'1 = `1'  /* successes */
			gen `event'0 = `2'  /* failures */
		}
		
		gen `rid' = _n	
		if "`panelize'" != ""  {
			reshape long `event', i(`rid') j(`idpair')
		}
		else {		
			reshape long `event' `total' `assignment', i(`rid') j(`idpair')
		}
	}	
end	


/*+++++++++++++++++++++++++	SUPPORTING FUNCTIONS: WIDESETUP +++++++++++++++++++++++++
							Transform data to wide format
	+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
	cap program drop widesetup
	program define widesetup, rclass

	syntax varlist, sid(varlist) idpair(varname) [sortby(varlist) jvar(varname) mpair mcbnetwork pcbnetwork abnetwork general comparative]

		qui {
			tokenize `varlist'

			tempvar modey diffy
			*if "`mcbnetwork'" == "" {
				tempvar jvar
				gen `jvar' = `idpair' - 1
			*}
			
			/*Check for varying variable and store them*/
			ds
			local vnames = r(varlist)
			local vlist
			foreach v of local vnames {	
				cap drop `modey' `diffy'
				bysort `sid': egen `modey' = mode(`v'), minmode
				egen `diffy' = diff(`v' `modey')
				sum `diffy'
				local sumy = r(sum)

				local in 0
				foreach str of local varlist {
					if `in' {
					continue, break	
					}
					if "`v'" =="`str'" {
						local in 1
					}
				}
				
				if (!`in') & (`sumy' > 0) & "`v'" != "`jvar'" & "`v'" != "`idpair'" {
					local vlist "`vlist' `v'"
				}
			}
			cap drop `modey' `diffy'
			
			sort `sid' `jvar' `sortby'
			
			/*2 variables per study : n N*/			
			reshape wide `1' `2'  `idpair' `vlist', i(`sid') j(`jvar')
			local cc0 = `idpair'0[1]
			local cc1 = `idpair'1[1]
			local idpair0 : lab `idpair' `cc0'
			local idpair1 : lab `idpair' `cc1'
			
			return local vlist = "`vlist'"
			return local cc0 = "`idpair0'"
			return local cc1 = "`idpair1'"
		}
	end	
/*+++++++++++++++++++++++++	SUPPORTING FUNCTIONS: PREP4SHOW +++++++++++++++++++++++++
							Prepare data for display table and graph
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
cap program drop prep4show
program define prep4show

	#delimit ;
	syntax varlist, [exactorout(name) poprrout(name) poprdout(name) poplrrout(name) rrout(name) rdout(name)  poporout(name) poplorout(name) orout(name) 
		popabsout(name)  exactabsout(name) absout(name) absoutp(name) sortby(varlist)  by(varname) hetout(name) model(string) 
		groupvar(varname) se(varname) summaryonly nooverall nosubgroup outplot(string) grptotal(name) download(string asis) 
		indvars(varlist) depvars(varlist) dp(integer 2) stratify pcont(integer 0) level(integer 95) prediction inference(string)
		comparative abnetwork general pcbnetwork mpair mcbnetwork aliasdesign(string) enhance stat(string)
		]
	;
	#delimit cr
	tempvar  expand serror
	tokenize `varlist'
	 
	local id = "`1'"
	local cid = "`2'"
	local use = "`3'"
	local label = "`4'"
	local es = "`5'"
	local lci = "`6'"
	local uci = "`7'"
	local modeles = "`8'"
	local modellci = "`9'"
	local modeluci = "`10'"
	
	if "`se'" !="" {
		gen `serror' = `se'
	}
	else{
		gen `serror' = 0
	}
	
	if "`outplot'" != "abs" & "`design'" == "abnetwork" & "`aliasdesign'" == "" {
		local summaryonly "summaryonly"
	}
	
	if "`stat'" == "Median" {
		local statscol = 3
	}
	else {
		local statscol = 1
	}
	
	if "`outplot'" == "lor" | "`outplot'" == "lrr"  {
		local transform "ln"
	}
	
	if "`inference'" == "bayesian" {
		local enhance 
	}
	
	qui {	
		gen `expand' = 1
	
		//Groups
		if "`groupvar'" != "" {	
			
			bys `groupvar' : egen `grptotal' = count(`id') //# studies in each group
			gsort `groupvar' `sortby' `id'
			bys `groupvar' : replace `expand' = 1 + 1*(_n==1) + 3*(_n==_N) 
			expand `expand'
			
			gsort `groupvar' `sortby' `id' `expand'
			bys `groupvar' : replace `use' = -2 if _n==1  //group label
			bys `groupvar' : replace `use' = 2 if _n==_N-2  //summary
			bys `groupvar' : replace `use' = 4 if _n==_N-1  //prediction
			bys `groupvar' : replace `use' = 3 if _n==_N //blank */
			replace `id' = `id' + 1 if `use' == 1
			replace `id' = `id' + 2 if `use' == 2  //summary 
			replace `id' = `id' + 3 if `use' == 4  //Prediction
			replace `id' = `id' + 4 if `use' == 3 //blank
			
			*replace `label' = "Summary" if `use' == 2 
			replace `label' = "Group Summary" if `use' == 2 
			replace _WT = . if `use' == 2 
						
			*qui label list `groupvar'
			*local nlevels = r(max)
			
			qui levelsof `groupvar', local(codelevels)
			local nlevels = r(r)
			
			*foreach l of local codelevels {
			forvalues l = 1/`nlevels' {
				local groupcode : word `l' of `codelevels'
				local lab:label `groupvar' `groupcode'
				
				if "`outplot'" == "abs" {
					if "`model'" == "hexact" {
						local S_1 = `exactabsout'[`l', 1]
						local S_3 = `exactabsout'[`l', 5]
						local S_4 = `exactabsout'[`l', 6]
					}
					else {
						if "`model'" != "crbetabin" {
							local S_1 = `popabsout'[`l', `statscol']
							local S_3 = `popabsout'[`l', 4]
							local S_4 = `popabsout'[`l', 5]
							
							
							//Get Conditional CI
							local C_1 = `absout'[`=`pcont' +`l'', 1]
							local C_3 = `absout'[`=`pcont' +`l'', 5]
							local C_4 = `absout'[`=`pcont' +`l'', 6]
						}
						else {
							//Get Conditional CI
							local S_1 = `absout'[`=`pcont' +`l'', 1]
							local S_3 = `absout'[`=`pcont' +`l'', 5]
							local S_4 = `absout'[`=`pcont' +`l'', 6]
						}
						
						//Get exact estimates
						local E_1 = `exactabsout'[`l', 1]
						local E_3 = `exactabsout'[`l', 5]
						local E_4 = `exactabsout'[`l', 6]
						local nstudies = `exactabsout'[`l', 9]   //number of studies
						local np0 = `exactabsout'[`l', 10]    //number of zero p's
						local np1 = `exactabsout'[`l', 11]   //number of 1 p's
						
						if "`enhance'" != "" {
							if "`model'" != "crbetabin" {
								//if simulated RE more than 5 times larger, replace with conditional stats 
								if (`=(`S_4' - `S_3')/(`C_4' - `C_3')' > 5) & (`C_4' == .) & (`C_3' == .)  {
									local S_3 = `C_3'
									local S_4 = `C_4'
								}
							}
							
							//if simulated FE variance 5 times larger replace with exact
							if /*(`=(`S_4' - `S_3')/(`E_4' - `E_3')' > 5) &*/ (`nstudies' == `np0') | (`nstudies' == `np1')  {
								local S_3 = `E_3'
								local S_4 = `E_4'
							}
						}
					}
					
					if "`prediction'" != "" {
						local S_5 = `absoutp'[`l', 1]
						local S_6 = `absoutp'[`l', 2]
					}
					if "`model'" == "random" & "`indvars'" == "" & "`stratify'" !="" {
						local isq = `hetout'[`l', 5]
						local phet = `hetout'[`l', 3]
						replace `label' = "Group Summary" + " (Isq = " + string(`isq', "%10.`=`dp''f") + "%, p = " + string(`phet', "%10.`=`dp''f") + ")" if (`use' == 2) & (`groupvar' == `groupcode') & (`grptotal' > 2)  & (`isq' != .)	
					}	 
				}
				else {
					if "`model'" == "hexact" {
						local S_1 = `transform'(`exactorout'[`l', 1])
						local S_3 = `transform'(`exactorout'[`l', 3])
						local S_4 = `transform'(`exactorout'[`l', 4])
					}
					else {
						if strpos("`outplot'", "rd") != 0 {
							if "`model'" != "crbetabin" {
								local S_1 = `poprdout'[`l', `statscol']
								local S_3 = `poprdout'[`l', 4]
								local S_4 = `poprdout'[`l', 5]
																						
								//Conditional stats
								local C_1 = (`rdout'[`l', 1])
								local C_3 = (`rdout'[`l', 5])
								local C_4 = (`rdout'[`l', 6])
							}
							else {
								//Conditional stats
								local S_1 = (`rdout'[`l', 1])
								local S_3 = (`rdout'[`l', 5])
								local S_4 = (`rdout'[`l', 6])
							}
						}
						if strpos("`outplot'", "rr") != 0 {
							if "`model'" != "crbetabin" {
								if "`outplot'" == "rr" {
									local S_1 = `poprrout'[`l', `statscol']
									local S_3 = `poprrout'[`l', 4]
									local S_4 = `poprrout'[`l', 5]
								}
								else {
									local S_1 = `poplrrout'[`l', `statscol']
									local S_3 = `poplrrout'[`l', 4]
									local S_4 = `poplrrout'[`l', 5]
								}
								
								//Conditional stats
								local C_1 = `transform'(`rrout'[`l', 1])
								local C_3 = `transform'(`rrout'[`l', 5])
								local C_4 = `transform'(`rrout'[`l', 6])
							}
							else {
								//Conditional stats
								local S_1 = `transform'(`rrout'[`l', 1])
								local S_3 = `transform'(`rrout'[`l', 5])
								local S_4 = `transform'(`rrout'[`l', 6])
							}
						}
						if strpos("`outplot'", "or") != 0 {
							if "`model'" != "crbetabin" {
								if "`outplot'" == "or" {
									local S_1 = `poporout'[`l', `statscol']
									local S_3 = `poporout'[`l', 4]
									local S_4 = `poporout'[`l', 5]
								}
								else {
									local S_1 = `poplorout'[`l', `statscol']
									local S_3 = `poplorout'[`l', 4]
									local S_4 = `poplorout'[`l', 5]
								}
								
								//Conditional stats
								local C_1 = `transform'(`orout'[`l', 1])
								local C_3 = `transform'(`orout'[`l', 5])
								local C_4 = `transform'(`orout'[`l', 6])
							}
							else {
								//Conditional stats
								local S_1 = `transform'(`orout'[`l', 1])
								local S_3 = `transform'(`orout'[`l', 5])
								local S_4 = `transform'(`orout'[`l', 6])
							}
						}
					}
					
					if "`enhance'" != "" & "`model'" != "crbetabin" {
					//if simulated more than 5 times larger, replace with conditional stats 
						if `=(`S_4' - `S_3')/(`C_4' - `C_3')' > 5 {
							local S_3 = `C_3'
							local S_4 = `C_4'
						}
					}
				}
				
				replace `label' = "`groupvar' = `lab'" if `use' == -2 & `groupvar' == `groupcode' & (("`abnetwork'" == "") |("`outplot'" == "abs" & "`abnetwork'" != ""))	
				replace `label' = "`lab'" if `use' == 2 & `groupvar' == `groupcode' & "`outplot'" != "abs" & "`abnetwork'" != ""		
				replace `es'  = `S_1' if `use' == 2 & `groupvar' == `groupcode'	
				replace `lci' = `S_3' if `use' == 2 & `groupvar' == `groupcode'	
				replace `uci' = `S_4' if `use' == 2 & `groupvar' == `groupcode'	
				//Predictions
				/*
				if "`outplot'" == "abs" & "`prediction'" != "" {
					replace `lci' = `S_5' if `use' == 4 & `groupvar' == `l'	
					replace `uci' = `S_6' if `use' == 4 & `groupvar' == `l'	
				}
				*/
				//Weights
				sum _WT if `use' == 1 & `groupvar' == `groupcode'
				local groupwt = r(sum)
				replace _WT = `groupwt' if `use' == 2 & `groupvar' == `groupcode'	
			}	
		}
		else {
			egen `grptotal' = count(`id') //# studies total
		}
		//Overall
		if "`overall'" == "" {		
			gsort  `groupvar' `sortby' `id' 
			replace `expand' = 1 + 3*(_n==_N)
			expand `expand'
			gsort  `groupvar' `sortby' `id' `expand'
			replace `use' = 4 if _n==_N  //Prediction
			replace `use' = 5 if _n==_N-1  //Overall
			replace `use' = 3 if _n==_N-2 //blank
			replace `id' = `id' + 3 if _n==_N  //Prediction
			replace `id' = `id' + 2 if _n==_N-1  //Overall
			replace `id' = `id' + 1 if _n==_N-2 //blank

			//Fill in the right info
			if "`outplot'" == "abs" {				
				if "`model'" == "hexact" {
					local nrows = rowsof(`exactabsout')
					local S_1 = `exactabsout'[`nrows', 1]
					local S_3 = `exactabsout'[`nrows', 5]
					local S_4 = `exactabsout'[`nrows', 6]
				}
				else if "`model'" == "crbetabin" {
					local nrows = rowsof(`absout')
					local S_1 = `absout'[`nrows', 1]
					local S_3 = `absout'[`nrows', 5]
					local S_4 = `absout'[`nrows', 6]
				}
				else {
					local nrows = rowsof(`popabsout')
					local S_1 = `popabsout'[`nrows', `statscol']
					local S_3 = `popabsout'[`nrows', 4]
					local S_4 = `popabsout'[`nrows', 5]
				}
				/*
				//predictions
				if "`prediction'" != "" {
					local nrows = rowsof(`absoutp')
					local S_5 = `absoutp'[`nrows', 1]
					local S_6 = `absoutp'[`nrows', 2]
				}*/			
			}
			else {
				if "`model'" == "hexact" {
					local nrows = rowsof(`exactorout')
					local S_1 = `transform'(`exactorout'[`nrows', 1])
					local S_3 = `transform'(`exactorout'[`nrows', 3])
					local S_4 = `transform'(`exactorout'[`nrows', 4])
				}
				else if "`model'" == "crbetabin" {
					if strpos("`outplot'", "rd") != 0 {
						local nrows = rowsof(`rdout')
						
						local S_1 = `rdout'[`nrows', 1]
						local S_3 = `rdout'[`nrows', 4]
						local S_4 = `rdout'[`nrows', 5]
					}
					
					if strpos("`outplot'", "rr") != 0 {
						local nrows = rowsof(`rrout')
						local S_1 = `transform'(`rrout'[`nrows', 1])
						local S_3 = `transform'(`rrout'[`nrows', 5])
						local S_4 = `transform'(`rrout'[`nrows', 6])
					}
					
					if strpos("`outplot'", "or") != 0 {
						local nrows = rowsof(`orout')
						
						local S_1 = `transform'(`orout'[`nrows', 1])
						local S_3 = `transform'(`orout'[`nrows', 5])
						local S_4 = `transform'(`orout'[`nrows', 6])
					}
				}
				else {
					if strpos("`outplot'", "rd") != 0 {
						local nrows = rowsof(`poprdout')
						
						local S_1 = `poprdout'[`nrows', `statscol']
						local S_3 = `poprdout'[`nrows', 4]
						local S_4 = `poprdout'[`nrows', 5]
					}
					
					if strpos("`outplot'", "rr") != 0 {
						local nrows = rowsof(`poprrout')
						if "`outplot'" == "rr" {
							local S_1 = `poprrout'[`nrows', `statscol']
							local S_3 = `poprrout'[`nrows', 4]
							local S_4 = `poprrout'[`nrows', 5]
						}
						else {
							local S_1 = `poplrrout'[`nrows', `statscol']
							local S_3 = `poplrrout'[`nrows', 4]
							local S_4 = `poplrrout'[`nrows', 5]
						}
					}
					
					if strpos("`outplot'", "or") != 0 {
						local nrows = rowsof(`poporout')
						if "`outplot'" == "or" {
							local S_1 = `poporout'[`nrows', `statscol']
							local S_3 = `poporout'[`nrows', 4]
							local S_4 = `poporout'[`nrows', 5]
						}
						else {
							local S_1 = `poplorout'[`nrows', `statscol']
							local S_3 = `poplorout'[`nrows', 4]
							local S_4 = `poplorout'[`nrows', 5]
						}
					}
				}
			}
			replace `label' = "Population Summary" if `use' == 5
			*replace `label' = "Population `stat'" if `use' == 5
			
			if "`model'" == "random" & "`indvars'" == ""  & "`outplot'" == "abs" {
				local nrows = rowsof(`hetout')
				local isq = `hetout'[`nrows', 5]
				local phet = `hetout'[`nrows', 3]
				*replace `label' = "Overall (Isq = " + string(`isq', "%10.`=`dp''f") + "%, p = " + string(`phet', "%10.`=`dp''f") + ")" if `use' == 3
				*replace `label' = "Population `stat' (Isq = " + string(`isq', "%10.`=`dp''f") + "%, p = " + string(`phet', "%10.`=`dp''f") + ")" if `use' == 5 & (`isq' != .)	
				replace `label' = "Population Summary (Isq = " + string(`isq', "%10.`=`dp''f") + "%, p = " + string(`phet', "%10.`=`dp''f") + ")" if `use' == 5 & (`isq' != .)	
			}
					
			replace `es' = `S_1' if `use' == 5	
			replace `lci' = `S_3' if `use' == 5
			replace `uci' = `S_4' if `use' == 5
			replace _WT = . if (`use' == 5) & ("`stratify'" != "")
			replace _WT = 100 if (`use' == 5) & ("`stratify'" == "")
			//Predictions
			if "`outplot'" == "abs" & "`prediction'" != "" {
				replace `lci' = `S_5' if _n==_N
				replace `uci' = `S_6' if _n==_N
			}
		}
		count if `use'==1 
		replace `grptotal' = `=r(N)' if `use'==5
		replace `grptotal' = `=r(N)' if _n==_N
		
		replace `label' = "" if `use' == 3 | `use' == 4
		replace `es' = . if `use' == 3 | `use' == -2 | `use' == 4  //4 is prediction 
		replace `lci' = . if `use' == 3 | `use' == -2
		replace `uci' = . if `use' == 3 | `use' == -2
				
		gsort `groupvar' `sortby'  `id' 
				
		*replace `label' = "Predictive t Interval" if `use' == 4 & "`model'" == "random"
		*replace `label' = "t Interval" if `use' == 4 & "`model'" != "random"
	}
	qui {
		replace `modeles' = . if `use' != 1
		replace `modellci' = . if `use' != 1
		replace `modeluci' = . if `use' != 1
		replace _WT = . if `use'==3 | `use'==-2 | `use'==4
	}	
	if `"`download'"' != "" {
		local ZOVE -invnorm((100-`level')/200)
		preserve
		qui {
			cap drop _ES  _SE _LCI _UCI _USE _LABEL _MODELES _MODELLCI _MODELUCI
			gen _ES = `es'
			gen _SE = `serror'
			gen _LCI = `lci'
			gen _UCI = `uci'
			gen _USE = `use'
			gen _LABEL = `label'
			gen _ID = `id'
			gen _MODELES = `modeles'
			gen _MODELLCI = `modellci'
			gen _MODELUCI = `modeluci'
			replace _ID = _n
			replace _SE = ( `uci' - `lci')/(2*`ZOVE') if _SE == 0
			
			*keep if _USE == 1
			keep `depvars' `indvars' `groupvar' _ES _SE _LCI _UCI _USE _ESAMPLE _WT _LABEL _ID _MODELES _MODELLCI _MODELUCI 
		}
		di "*********************************************************************"
		di _n "Saving data....."
		di "Note: For n=N or n=0, _SE=0"
		di "and approximated with _SE = (_UCI – _LCI)/(2*Z(`level'))"
		noi save `download', replace
		
		restore
	}
	qui {
				
		if "`abnetwork'" == "" | ("`abnetwork'" != "" & "`outplot'" == "abs") {
			drop if (`use' == 2 | `use' == 5) & (`grptotal' == 1)  //drop summary if 1 study
		}
		drop if (`use' == 1 & "`summaryonly'" != "" & `grptotal' > 1) | (`use' == 2 & "`subgroup'" != "") | (`use' == 5 & "`overall'" != "") | (`use' == 4 & "`prediction'" == "") //Drop unnecessary rows
		
		if "`abnetwork'" != "" & "`outplot'" != "abs" & "`aliasdesign'" != "comparative" {
			drop if `use' == 1 | `use' == -2
			replace `use' = 1 if `use' == 2
		}

		gsort `groupvar' `sortby' `id'
				
		replace `id' = _n
		
		gsort `groupvar' `use' -`es'
		
		replace `cid' = _n		
	}
end	
/*+++++++++++++++++++++++++	SUPPORTING FUNCTIONS: DISPTAB +++++++++++++++++++++++++
							Prepare data for display table and graph
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
cap program drop disptab
program define disptab
	#delimit ;
	syntax varlist, [nosubgroup nooverall level(integer 95) sumstat(string asis) model(string)
	dp(integer 2) power(integer 0) nowt smooth icimethod(string) scimethod(string) groupvar(varname) 
	design(string) aliasdesign(string) outplot(string) summaryonly stat(string) inference(string)]
	;
	#delimit cr
	
	tempvar rid id use label es lci uci grptotal modeles modellci modeluci
	tokenize `varlist'
	qui {
		gen `id' = `1'
		gen `use' = `2'
		gen `label' = `3'
		gen `es' = `4'
		gen `lci' = `5'
		gen `uci' = `6'
		gen `grptotal' = `7'
		
		if "`smooth'" !="" {
			gen `modeles' = `8'
			gen `modellci' = `9'
			gen `modeluci' = `10'
		}
	}
	
	if "`outplot'" != "abs" & "`design'" == "abnetwork" & "`aliasdesign'" == "" {
		local summaryonly "summaryonly"
	}
		
	preserve
		tempvar tlabellen 
		//study label
		local studylb: variable label `label'
		if "`studylb'" == "" {
			local studylb "Study"
		}		
		local studylen = strlen("`studylb'")
	 
		qui gen `tlabellen' = strlen(`label')
		qui summ `tlabellen' if `use' == 1 
			
		local nlen = `=max(r(max), 15) + 2' 
		local nlens = strlen("`sumstat'")
		
		local level: displ %2.0f `level'
		local start: displ %2.0f `=`nlen'/2 - `studylen'/2 + 2'
		
		di as res _n "***********************************************************************"
		if "`smooth'" != "" {
			di as res "{pmore2} Study specific `sumstat' :  Observed (Smoothed) {p_end}"
		}
		else {
			di as res "{pmore2} Study specific `sumstat'  {p_end}"
		}
		di as res    "***********************************************************************" 
		
		local colstat = int(`=`nlen' + `nlens'*0.5')
		
		
		if "`wt'" =="" {
			local dispwt "% Weight"
		}
		
		//Find the length of the estimates
		qui {
			tempvar hold holdstr slimest
			
			gen `hold' = `uci'*(10^`power')
			replace `hold' = `lci'*(10^`power') if `hold' == .
			replace `hold' = `es'*(10^`power') if `hold' == .
			
			tostring `hold', gen(`holdstr') format(%10.`dp'f) force
			gen `slimest' = strlen(strltrim(`holdstr'))
			sum `slimest'
			local est_i_len = r(max)
		}
		
		if "`smooth'" !="" {
			local open " ("
			local close ")"
			local aesclose "%1s"
			local aes "%`=`est_i_len''.`=`dp''f"
		}
		if "`summaryonly'"  == "" {
			if "`smooth'" == "" {
				local citext "- `icimethod' CI -"
			}
			else  {
				local citext "- `icimethod' (Wald) CI - "
			}
		}
		else {
			local citext "- Centile CI - "
		}
		if "`outplot'" == "abs" {
			if "`summaryonly'"  == "" {
				if "`smooth'" == "" {
					local citext "- `icimethod' CI -"
				}
				else  {
					if "`inference'" == "bayesian" {
						local citext "- `icimethod' (Centile) CI - "
					}
					else {
						local citext "- `icimethod' (Wald) CI - "
					}
				}
			}
			else {
				local citext "- Centile CI - "
			}
		}
		else {
			if "`smooth'" == "" {
				local citext "- `icimethod' CI -"
			}
			else  {
				local citext "- `icimethod' (Centile) CI - "
			}
		}	
		
		if "`outplot'" =="abs" & ("`design'" == "mpair" | "`design'" == "comparative" | "`aliasdesign'" == "comparative")  & "`groupvar'" != "" {
			qui {
				if "`smooth'" == "" {
					keep `id'  `use' `label' `es' `lci' `uci' `grptotal'  `groupvar' _WT
				}
				else {
					keep `id'  `use' `label' `es' `lci' `uci' `grptotal' `modeles' `modellci' `modeluci' `groupvar' _WT
				}
				bys `groupvar': egen `rid' = seq()
				replace `groupvar' = `groupvar' - 1
				duplicates drop `groupvar' if `use'==3, force
				if "`smooth'" == "" {
					reshape wide `id' `label' `es' `lci' `uci' `grptotal' _WT, i(`rid') j(`groupvar')
				}
				else {
					reshape wide `id' `label' `es' `lci' `uci' `grptotal' `modeles' `modellci' `modeluci' _WT, i(`rid') j(`groupvar')
				}
				local label0 = `label'0[1]
				local label1 = `label'1[1]

				local nlen0 = strlen("`label0'")
				local nlen1 = strlen("`label1'")
				
				replace _WT0 = _WT0 + _WT1
				rename _WT0 _WT				
			}
			
			di _n as txt _col(`nlen') "| "   _skip(`=21 - round(`nlen0'/2)') "`label0'" ///
					  _skip(`=47 - (21 - round(`nlen0'/2)) - `nlen0' - 1')	"| " _skip(`=21 - round(`nlen1'/2)') "`label1'" _cont
			
			di  _n  as txt _col(`start') "`studylb'" _col(`nlen') "| "   _skip(5) "Estimate" ///
					  _skip(5) "[`level'% Conf. Interval]"  ///
					  _skip(9)	"| " _skip(5) "Estimate" ///
					  _skip(5) "[`level'% Conf. Interval]" ///
					  _skip(12) "| " _skip(3) "`dispwt'"
					  
			di  _dup(`=`nlen'-1') "-" "+" _dup(48) "-" "+" _dup(51) "-" "+" _dup(10) "-"
			qui count
			local N = r(N)
			

			
			forvalues i = 1(1)`N' {
				//Weight
				if "`wt'" =="" {
					local ww = _WT[`i']
				}
				
				//Studies -- Control
				if ((`use'[`i'] ==1)) {
					//Smooth estimates
					if "`smooth'" !="" {					
						local mes0 "`modeles'0[`i']*(10^`power')"
						local mlci0 "`modellci'0[`i']*(10^`power')"
						local muci0 "`modeluci'0[`i']*(10^`power')"
					}
					
					di _col(2) as txt `label'0[`i'] _col(`nlen') "|  "  ///
					_skip(2) as res  %5.`=`dp''f  `es'0[`i']*(10^`power') "`open'" `aes' `mes0' "`close'" /// 
					_col(`=`nlen' + 20') %5.`=`dp''f `lci'0[`i']*(10^`power') "`open'" `aes' `mlci0'  `aesclose' "`close'" ///
					_skip(5) %5.`=`dp''f `uci'0[`i']*(10^`power') "`open'" `aes' `muci0'  `aesclose' "`close'" _cont
				}
				//studies - Treatment
				if (`use'[`i'] ==1 )   { 
					//Smooth estimates
					if "`smooth'" !="" {					
						local mes1 "`modeles'1[`i']*(10^`power')"
						local mlci1 "`modellci'1[`i']*(10^`power')"
						local muci1 "`modeluci'1[`i']*(10^`power')"
					}
					
					di as txt _col(`=`nlen' + 45') "|  "  ///
					_skip(2) as res  %5.`=`dp''f  `es'1[`i']*(10^`power') "`open'" `aes' `mes1' "`close'" /// 
					_col(`=`nlen' + 72') %5.`=`dp''f `lci'1[`i']*(10^`power') "`open'" `aes' `mlci1'  `aesclose' "`close'" ///
					_skip(5) %5.`=`dp''f `uci'1[`i']*(10^`power') "`open'" `aes' `muci1'  `aesclose' "`close'"   _col(`=`nlen' + 90') as txt "|  " _skip(3) as res %5.`=`dp''f `ww'
				}
				//Summaries
				if (`use'[`i']== 2) {
					di as res _dup(`=`nlen'-1') "-" "+" _dup(48) "-" "+" _dup(51) "-" "+" _dup(10) "-"
					
					di _col(2) as txt `label'0[`i'] _col(`nlen') "|  "  ///
					_skip(`=3 + 5') as res  %5.`=`dp''f  `es'0[`i']*(10^`power') /// 
					_col(`=`nlen' + 20 + 6') %5.`=`dp''f `lci'0[`i']*(10^`power') ///
					_skip(`=5 + 7') %5.`=`dp''f `uci'0[`i']*(10^`power') ///
					as txt _col(`=`nlen' + 45 + 4') "|  " ///
					_skip(`=3 + 5') as res  %5.`=`dp''f  `es'1[`i']*(10^`power') /// 
					_col(`=`nlen' + 66 + 12') %5.`=`dp''f `lci'1[`i']*(10^`power') ///
					_skip(`=5 + 7') %5.`=`dp''f `uci'1[`i']*(10^`power')  ///
					_col(`=`nlen' + 90 + 8') as txt " |  " _skip(2) as res %5.`=`dp''f `ww'
				}
				//Blanks
				if (`use'[`i'] == 0 ){
						di as res _dup(`=`nlen'-1') "-" "+" _dup(48) "-" "+" _dup(51) "-" "+" _dup(10) "-"
				}
			}
		}
		else {		
		
			if "`smooth'" !=""  {			  
				di  _n  as txt _col(`start') "`studylb'" _col(`nlen') "|  "   _skip(5) "Estimate" ///
				  _col(`=`nlen' + `nlens' + 20') "`=(100-`level')/2'% "`"`citext'"'" `=100 - (100-`level')/2'%" _col(`=`nlen' + `nlens' + 60') "`dispwt'"  
				
				local colwt = int(`=`nlen' + `nlens' + 55')
			}
			else{
				di  _n  as txt _col(`start') "`studylb'" _col(`nlen') "|  "   _skip(5) "Estimate" ///
				  _col(`=`nlen' + `nlens' + 10') "`=(100-`level')/2'% "`"`citext'"'" `=100 - (100-`level')/2'%" _col(`=`nlen' + `nlens' + 40') "`dispwt'"
			
				local colwt = int(`=`nlen' + `nlens' + 35')
			}
			di  _dup(`=`nlen'-1') "-" "+" _dup(57) "-" 
			
			qui count
			local N = r(N)
					
			forvalues i = 1(1)`N' {
				//Weight
				if "`wt'" =="" {
					local ww = _WT[`i']
				}
				//Group labels
				if ((`use'[`i']== -2)){ 
					di _col(2) as txt `label'[`i'] _col(`nlen') /*"|  "*/
				}
				
				//Studies 
				if ((`use'[`i'] ==1)) {
							
					//Smooth estimates
					if "`smooth'" !="" {						
						local mes = "`=`modeles'[`i']*(10^`power')'"
						local mlci = "`=`modellci'[`i']*(10^`power')'"
						local muci = "`=`modeluci'[`i']*(10^`power')'"
					}
					
					di _col(2) as txt `label'[`i'] _col(`nlen') "|  "  ///
					_col(`colstat')  as res  %10.`=`dp''f  `es'[`i']*(10^`power')  "`open'" `aes' `mes' "`close'"  /// 
					_col(`=`nlen' + `nlens' + 5') %10.`=`dp''f `lci'[`i']*(10^`power') "`open'" `aes' `mlci'  `aesclose' "`close'"  ///
					_skip(5) %10.`=`dp''f `uci'[`i']*(10^`power') "`open'" `aes' `muci' "`close'"   _col(`colwt') %10.`=`dp''f `ww'
				}
				//Summaries
				if ( (`use'[`i']== 5) | ((`use'[`i']== 2) & (`grptotal'[`i'] > 1))){
					if ((`use'[`i']== 2) & (`grptotal'[`i'] > 1)) {
						di _col(2) as txt _col(`nlen') "|  " 
					}
					if (`use'[`i']== 2)	{
						local sumtext = "Group Summary"
					}
					else {
						local sumtext = "Population Summary"			
					}
					if "`smooth'" != "" {
						di _col(2) as txt "`sumtext'" _col(`nlen') "|  "  ///
							_col(`=`colstat'+8') as res  %`=`est_i_len''.`=`dp''f  `es'[`i']*(10^`power') /// 
							_col(`=`nlen' + `nlens' + 22') %`=`est_i_len''.`=`dp''f `lci'[`i']*(10^`power') ///
							_skip(14) %`=`est_i_len''.`=`dp''f `uci'[`i']*(10^`power') _col(`colwt')  %10.`=`dp''f `ww'
					}
					else {
						di _col(2) as txt "`sumtext'" _col(`nlen') "|  "  ///
						_col(`colstat') as res  %10.`=`dp''f  `es'[`i']*(10^`power') /// 
						_col(`=`nlen' + `nlens' + 5') %10.`=`dp''f `lci'[`i']*(10^`power') ///
						_skip(5) %10.`=`dp''f `uci'[`i']*(10^`power') _col(`colwt')  %10.`=`dp''f `ww'
					}
				}
				//Blanks
				if (`use'[`i'] == 0 | `use'[`i'] == 3  ){
					di as txt _dup(`=`nlen'-1') "-" "+" _dup(57) "-"		
				}
			}
		}		
	restore
end

	/*++++++++++++++++	SUPPORTING FUNCTIONS: BUILDEXPRESSIONS +++++++++++++++++++++
				buildexpressions the regression and estimation expressions
	+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
	cap program drop buildregexpr
	program define buildregexpr, rclass
		
		syntax varlist, [interaction alphasort mcbnetwork pcbnetwork abnetwork general mpair comparative ipair(varname) baselevel(string) studyid(varname) model(string)]
		
		tempvar holder
		tokenize `varlist'

		if "`mcbnetwork'`pcbnetwork'`mpair'" == "" {
			macro shift 2
			local regressors "`*'"
		}
		else if "`mpair'"  != "" {
			macro shift 4
			local regressors "`*'"
			
			my_ncod `holder', oldvar(`ipair')
			drop `ipair'
			rename `holder' `ipair'
		}
		else {
			if "`mcbnetwork'" != "" {
				local index = "`5'"
				local comparator = "`6'"
				macro shift 6
				}
			else {
				local index = "`4'"
				local comparator = "`5'"
				macro shift 5
			}			
			local regressors "`*'"
			
			my_ncod `holder', oldvar(`index')
			drop `index'
			rename `holder' `index'

			my_ncod `holder', oldvar(`ipair')
			drop `ipair'
			rename `holder' `ipair'
			
			my_ncod `holder', oldvar(`comparator')
			drop `comparator'
			rename `holder' `comparator'
		}
		
		local p: word count `regressors'
		
		local catreg " "
		local contreg " "
		
		if ("`general'`comparative'" != "") {
			local regexpression = "mu"
		}
		else if "`mpair'" != "" {
			local regexpression = "mu i.`ipair'"	
		}
		else if "`mcbnetwork'`pcbnetwork'" != "" {
			if "`interaction'" != "" {				
				local regexpression = "ibn.`ipair'#ibn.`comparator' i.`index'"
				//nulllify
				local interaction
			}
			else {				
				if "`mcbnetwork'" != "" {
					*local regexpression = "mu i.`comparator' i.`index'"
					local regexpression = "mu i.`ipair' i.`index'"						
				}
				else {
					local regexpression = "mu i.`ipair' i.`index'"	
				}
			}
		}
		else { 
			*abnetwork 
			local regexpression 
		}
		if ("`model'" == "cbbetabin") {
			local regexpression2 = "mu"
			
			if "`abnetwork'" != "" {
				local regexpression2
			}
		}
		
		local basecode 1
		tokenize `regressors'
		forvalues i = 1(1)`p' {	
		
			cap label list ``i'' //see if labelled
			
			if _rc != 0 {
				cap confirm string variable ``i''  //check if string
				if _rc == 0 {
					if "`alphasort'" != "" {
						sort ``i''
					}
					my_ncod `holder', oldvar(``i'')
					drop ``i''
					rename `holder' ``i''
					local prefix_`i' "i"
				}
				else {
					capture confirm numeric var ``i''
					if _rc == 0 {
						local prefix_`i' "c"
					}
				}
			}
			else {
				local prefix_`i' "i"
			}
			
			if "`baselevel'" != "" & `i'==1 {
				//Find the base level
				*qui label list ``i''
				*local nlevels = r(max)
				
				qui levelsof ``i'', local(codelevels)
				local nlevels = r(r)
				
				local found = 0
				
				foreach l of local codelevels {
					local lab:label ``i'' `l'
					if "`lab'" == "`baselevel'" {
						local found = 1
						local basecode `l'
					}
					if `found' {
						continue, break
					}
				}
				/*
				local level 1
				while !`found' & `level' < `=`nlevels'+1' {
					local lab:label ``i'' `level'
					if "`lab'" == "`baselevel'" {
						local found = 1
						local basecode `level'
					}
					local ++level
				}
				*/
				
				/*if "`abnetwork'" != "" {
					local prefix_`i' "ibn"
				}*/
				if "`general'`comparative'" != "" {
					local prefix_`i' "ib`basecode'"
				}
			}
			/*Add the proper expression for regression*/
			local regexpression2 = "`regexpression2' `prefix_`i''.``i''#c.mu"     //for cbbetabin
			local regexpression = "`regexpression' `prefix_`i''.``i''"   //for other models
				
			if `i' > 1 & "`interaction'" != "" {
				local regexpression = "`regexpression' `prefix_`i''.``i''#`prefix_1'.`1'"   //for cbbetabin
				local regexpression2 = "`regexpression2' `prefix_`i''.``i''#`prefix_1'.`1'#c.mu"	 //for other models			
			}
						
			if "``i''" == "`studyid'" {
				continue
			}
			//Pick out the interactor variable
			if `i' == 1 & "`interaction'" != "" {
				local varx = "``i''"
				local typevarx = "`prefix_`i''"
			}
			if strpos("`prefix_`i''","i")  != 0 {
				local catreg "`catreg' ``i''"
			}
			else {
				local contreg "`contreg' ``i''"
			}
		}
		
		if "`interaction'" != "" {
			return local varx = "`varx'"
			return local typevarx  = "`typevarx'"
		}				
		return local  regexpression = "`regexpression'"
		return local  regexpression2 = "`regexpression2'"
		return local  catreg = "`catreg'"
		return local  contreg = "`contreg'"
		return local basecode = "`basecode'"
	end

/*+++++++++++++++++++++++++	SUPPORTING FUNCTIONS:  estp +++++++++++++++++++++++++
							Proportions after modelling
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
cap program drop estp
	program define estp, rclass
	
	syntax, rawestmat(name)  [link(string) scimethod(string) model(string)]
	
	if "`link'" == "cloglog" {
		local invfn "invcloglog"
	}
	else if "`link'" == "loglog" {
		if "`model'" == "crbetabin" {
			local invfn "exp(-exp(-"
			local closebracket "))"
		}
		else {
			local invfn "1 - invcloglog"
		}
	}
	else {
		local invfn "invlogit"
	}
	
	if "`scimethod'"== "t" {
		local statistic "t"
	}
	else {
		local statistic "z"
	}
	
	tempname matrixout
	mat `matrixout' = `rawestmat'
	
	local nrows = rowsof(`matrixout')
	local ncols = colsof(`matrixout')
			
	if `ncols' > 2 {	
		forvalues r = 1(1)`nrows' {
			mat `matrixout'[`r', 1] = `invfn'(`matrixout'[`r', 1])`closebracket' //p
			
			if "`link'" == "loglog" & "`model'" != "crbetabin" {
				mat `matrixout'[`r', 5] = `invfn'(`rawestmat'[`r', 6])`closebracket' //lower
				mat `matrixout'[`r', 6] = `invfn'(`rawestmat'[`r', 5])`closebracket' //upper
				mat `matrixout'[`r', 8] = `invfn'(`rawestmat'[`r', 9])`closebracket' //lower
				mat `matrixout'[`r', 9] = `invfn'(`rawestmat'[`r', 8])`closebracket' //upper
			}
			else {
				mat `matrixout'[`r', 5] = `invfn'(`matrixout'[`r', 5])`closebracket' //lower
				mat `matrixout'[`r', 6] = `invfn'(`matrixout'[`r', 6])`closebracket' //upper
				mat `matrixout'[`r', 8] = `invfn'(`matrixout'[`r', 8])`closebracket' //lower
				mat `matrixout'[`r', 9] = `invfn'(`matrixout'[`r', 9])`closebracket' //upper
			}
		}
		mat colnames `matrixout' = Mean SE(`link') `statistic'(`link') P>|z| z_Lower z_Upper P>|t| t_Lower t_Upper
	}
	else {
		forvalues r = 1(1)`nrows' {
			mat `matrixout'[`r', 1] = `invfn'(`matrixout'[`r', 1])`closebracket'  //lower
			mat `matrixout'[`r', 2] = `invfn'(`matrixout'[`r', 2])`closebracket'  //upper
		}
	}
	
	return matrix outmatrix = `matrixout' 
end	
	
	
/*+++++++++++++++++++++++++	SUPPORTING FUNCTIONS:  BAYEssummary +++++++++++++++++++++++++
							estimate log odds or proportions after modelling
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/	
	cap program drop bayessummary
	program define bayessummary, rclass
		syntax, estimates(string) studyid(varname) [event(varname) total(varname) DP(integer 2) model(string) varx(varname) typevarx(string) regexpression(string) ///
			comparator(varname) scimethod(string) mpair mcbnetwork pcbnetwork abnetwork general comparative stratify interaction ///
			catreg(varlist) contreg(varlist) power(integer 0) level(integer 95) by(varname) link(string) cov(string) baselevel(integer 1)]
		
		tempname absout loddsout exactabsout exactabsouti absexact etimat hpdmat
		tempvar subset insample hold holdleft holdright
		
		tokenize `regexpression'
		if "`mcbnetwork'`pcbnetwork'" != "" {
			 if "`interaction'" != "" {
				tokenize `2', parse(".")
			 }
			 else {
				tokenize `3', parse(".")
			 }
			local index "`3'"
			local catreg = "`comparator' `3' `catreg'"
			local varx //nullify
			if "`by'" != "`comparator'" {
				*local catreg = "`comparator' `catreg'"
			}
		}
		if "`abnetwork'" != "" {
			tokenize `2', parse(".")
			local assignment "`3'"
			local catreg = "`3' `catreg'"
		}
		
		if "`interaction'" != "" & "`typevarx'" == "i" {
			local idpairconcat "#`varx'"
		}
		
		if "`typevarx'" == "i"  {
			if "`catreg'" == "" {
				local catreg = "`varx'"
			}
			else {
				if "`idpairconcat'" == "" & "`by'" != "`varx'" {
					local catreg "`varx' `catreg'"
				}
			}
		}
		else {
			if "`contreg'" == "" {
				local contreg = "`varx'"
			}
		}
		
		if "`mcbnetwork'`pcbnetwork'" != "" {
			local rownamesmaxlen : strlen local Index
			local rownamesmaxlen = max(`rownamesmaxlen', 10)
		}
		else {
			local rownamesmaxlen = 10 /*Default*/
		}
		
		local invfn "invlogit"
						
		local ncatreg 0
		local parmlodds
		local parmp
		
		local catvars = "`catreg'"
		if strpos("`regexpression'", "`studyid'") != 0 {
			local mu "muoff"
		}  
		
		local ncatreg 0
		if  "`cov'" == "freeint" | "`mu'" == "muoff" {
			foreach c of local catvars {
				
				*qui label list `c'
				*local nlevels = r(max)
				
				qui levelsof `c', local(codelevels)
				local nlevels = r(r)
				
				foreach l of local 	codelevels {
				*forvalues l = 1/`nlevels' {
					local lab:label `c' `l'
					local lab = ustrregexra("`lab'", " ", "_")
					local nlen : strlen local lab
					local rownamesmaxlen = max(`rownamesmaxlen', min(`nlen', 32)) //Check if there is a longer name
					local rownames = "`c':`lab' `rownames'" 
				}
				local ncatreg = `nlevels' + `ncatreg'
			}
			
			mat `loddsout' = J(`ncatreg', 8, .)
			mat `absout' = J(`ncatreg', 8, .)
			mat colnames `loddsout' = Mean SD MCSE Median eti_Lower eti_Upper hpd_Lower hpd_Upper
			mat colnames `absout' = Mean SD MCSE Median eti_Lower eti_Upper hpd_Lower hpd_Upper
		}
		else {
			if "`catvars'" != "" {
				foreach c of local catvars {
					*qui label list `c'
					*local nlevels = r(max)
					
					qui levelsof `c', local(codelevels)
					local nlevels = r(r)
					
					foreach l of local 	codelevels {
						if "`interaction'" != "" {								
							local xterm = " + {`event':`l'.`c'#2.`varx'}"
													
							if `l' == `baselevel' {
								local parmp = "`parmp' (`c'_`l'_`varx'_1:`invfn'({`event':mu})) (`c'_`l'_`varx'_2:`invfn'({`event':mu} + {`event':2.`varx'}))  "
								local parmlodds = "`parmlodds' (`c'_`l'_`varx'_1:({`event':mu})) (`c'_`l'_`varx'_2:({`event':mu} + {`event':2.`varx'}))"
							} 
							else {
								local parmp = "`parmp' (`c'_`l'_`varx'_1:`invfn'({`event':`l'.`c'} + {`event':mu})) (`c'_`l'_`varx'_2:`invfn'({`event':`l'.`c'} + {`event':mu} + {`event':2.`varx'} `xterm'))"
								local parmlodds = "`parmlodds' (`c'_`l'_`varx'_1:({`event':`l'.`c'} + {`event':mu})) (`c'_`l'_`varx'_2:({`event':`l'.`c'} + {`event':mu} + {`event':2.`varx'} `xterm'))"
							}
						}
						else {
							if `l' == `baselevel' {
								local parmp = "`parmp' (`c'_`l':`invfn'({`event':mu}))"
								local parmlodds = "`parmlodds' (`c'_`l':({`event':mu}))"
							} 
							else {
								local parmp = "`parmp' (`c'_`l':`invfn'({`event':`l'.`c'} + {`event':mu}))"
								local parmlodds = "`parmlodds' (`c'_`l':({`event':`l'.`c'} + {`event':mu}))"
							}
						}
					}
				}
			}
			else {
				//mu
				local parmlodds = "(Overall:{`event':mu})"
				local parmp = "(Overall:`invfn'({`event':mu}))"
			}
			
			if "`parmlodds'" != "" | (  "`parmlodds'" == "" & ("`abnetwork'`mcbnetwork'`pcbnetwork'" == "" )) {
				//eti
				bayesstats summary `parmlodds', clevel(`level') 
				mat `etimat' = r(summary)	
				
				//hpd
				bayesstats summary `parmlodds', clevel(`level') hpd
				mat `hpdmat' = r(summary)
				mat `hpdmat' = `hpdmat'[1..., 5..6]
				
				mat `loddsout' = (`etimat' , `hpdmat') 
				local rnames :rownames `loddsout'
				local ncatreg = rowsof(`loddsout')
				
				//eti
				bayesstats summary `parmp', clevel(`level') 
				mat `etimat'  = r(summary)
				
				//hpd
				bayesstats summary `parmp', clevel(`level') hpd
				mat `hpdmat'  = r(summary)
				mat `hpdmat' = `hpdmat'[1..., 5..6]

				mat `absout' = (`etimat' , `hpdmat')			
			}
																		
			//Nice labels
			forvalues r = 1(1)`ncatreg' {
				local rname`r':word `r' of `rnames'
				
				tokenize `rname`r'', parse("_")	
								
				if "`mpair'" == "" {
					if "`7'" != "" {
						local leftvar = "`1'"
						local leftlab :label `1' `3' 
						local leftlab = ustrregexra("`leftlab'", " ", "_")
						local rightvar = "`5'"
						local rightlab : label `5' `7'
						local rightlab = ustrregexra("`rightlab'", " ", "_")
						
						local nlen1l:strlen local leftlab
						local nlenrl:strlen local rightlab
						local nlen1v:strlen local leftvar
						local nlenrv:strlen local rightvar
						
						local nlen = `nlen1l' + `nlenrl' + `nlen1v' + `nlenrv'
						local rownamesmaxlen = max(`rownamesmaxlen', min(`nlen', 32)) //Check if there is a longer name
						local rownames = "`rownames' `leftvar'*`rightvar':`leftlab'|`rightlab'" 
					}
					else if "`3'" != "" {
						local left = "`1'"
						local right = "`3'"
						local lab:label `left' `right'
						local lab = ustrregexra("`lab'", " ", "_")
						local nlen : strlen local lab
						local rownamesmaxlen = max(`rownamesmaxlen', min(`nlen', 32)) //Check if there is a longer name
						local rownames = "`rownames' `left':`lab'" 
					}
					else if "`3'" == "" {
						local rownames = "`rownames' `rname`r''" 
					}
				}
				else {
					local left = "`1'`2'`3'"
					local right = "`5'"
					local lab:label `left' `right'
					local lab = ustrregexra("`lab'", " ", "_")
					local nlen : strlen local lab
					local rownamesmaxlen = max(`rownamesmaxlen', min(`nlen', 32)) //Check if there is a longer name
					local rownames = "`rownames' `left':`lab'" 
				}
			}
		}

		mat rownames `loddsout' = `rownames'
		mat rownames `absout' = `rownames'
		
		mat colnames `loddsout' = Mean SD MCSE Median eti_Lower eti_Upper hpd_Lower hpd_Upper
		mat colnames `absout' = Mean SD MCSE Median eti_Lower eti_Upper hpd_Lower hpd_Upper
										
		//Get exact stats
		qui {
			gen `insample' = e(sample)
			local nrows = rowsof(`absout') //length of the vector
			local rnames :rownames `absout'
			local eqnames :roweq `absout'
			local newnrows = 0
			local mindex = 0
				
			foreach vari of local eqnames {		
				local ++mindex
				local group : word `mindex' of `rnames'
				
				//Skip if continous variable
				if (strpos("`vari'", "_") == 1) & ("`group'" != "Overall") & "`mpair'" == ""{
					continue
				}
				
				cap drop `subset' 
				
				if "`group'" != "Overall" {
					if strpos("`vari'", "*") == 0 {
					*if "`interaction'" == "" {
						cap drop `hold'
						decode `vari', gen(`hold')
						cap drop `subset'
						local latentgroup = ustrregexra("`group'", "_", " ")
						gen `subset' = 1 if `hold' == "`latentgroup'" & `insample' == 1 
					}
					else {
						tokenize `vari', parse("*")
						local leftvar = "`1'"
						local rightvar = "`3'"
						
						tokenize `group', parse("|")
						local leftgroup = "`1'"
						local rightgroup = "`3'"
						
						cap drop `holdleft' `holdright'
						decode `leftvar', gen(`holdleft')
						decode `rightvar', gen(`holdright')
						cap drop `subset'
						local latentleftgroup = ustrregexra("`leftgroup'", "_", " ")
						local latentrightgroup = ustrregexra("`rightgroup'", "_", " ")
						gen `subset' = 1 if (`holdleft' == "`latentleftgroup'") & (`holdright' == "`latentrightgroup'") & (`insample' == 1)
					}					
				}
				else {
					//All
					gen `subset' = 1 if `insample' == 1 
				}
				
				count if `subset' == 1 
				local nsubset = r(N)
				
				//Get the strata total
				sum `total' if `subset' == 1
				local stratatotal = r(sum)
				
				//Get the strata events
				sum `event' if `subset' == 1
				local strataevents = r(sum)
				
				//Get the strata size
				count if `subset' == 1
				local stratasize = r(N)
				
				tempvar ones zeros
				gen `ones' = `event'==`total' 
				sum `ones' if `subset' == 1 
				local sumones = r(sum)
				
				gen `zeros' = `event'==0 
				sum `zeros' if `subset' == 1  
				local sumzeros = r(sum)
				
				//Obtain the exact stats
				absexactci `stratatotal' `strataevents',  level(`level') //exact ci
				mat `absexact' = r(absexact)
				local modelp = `absexact'[1, 1]
				local postse = `absexact'[1, 2]
				local lowerp = `absexact'[1, 5]
				local upperp = `absexact'[1, 6]
				
				mat `exactabsouti' = (r(absexact), `stratatotal', `strataevents', `stratasize', `sumones', `sumzeros')
				mat rownames `exactabsouti' = `vari':`group'
				
				//Stack the matrices
				local ++newnrows
				if `newnrows' == 1 {
					mat `exactabsout' = `exactabsouti'	
				}
				else {
					mat `exactabsout' = `exactabsout'	\  `exactabsouti'
				}
			}
			mat colnames `exactabsout' = Mean SE z P>|z| Lower Upper Total Events Studies Ones Zeros
		}
		
		return matrix exactabsout = `exactabsout'		
		return matrix loddsout = `loddsout'
		return matrix absout = `absout'
	end	
	
/*+++++++++++++++++++++++++	SUPPORTING FUNCTIONS:  freqsummary +++++++++++++++++++++++++
							estimate raw estimates after modelling
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/	
	cap program drop freqsummary
	program define freqsummary, rclass

		syntax, estimates(string) studyid(varname) [event(varname) total(varname) abs DP(integer 2) model(string) varx(varname) typevarx(string) regexpression(string) ///
			comparator(varname) scimethod(string) mpair mcbnetwork pcbnetwork abnetwork general comparative  stratify interaction ///
			catreg(varlist) contreg(varlist) power(integer 0) level(integer 95) by(varname) link(string) total(varname)]
		
		tempname coefmat outmatrix outmatrixp matrixout bycatregmatrixout catregmatrixout contregmatrixout row ///
		outmatrixr overall Vmatrix byVmatrix exactabsouti exactabsout absexact tstats
		tempvar subset insample hold holdleft holdright
		
		tokenize `regexpression'
		if "`mcbnetwork'`pcbnetwork'" != "" {
			if "`interaction'" != "" {
				tokenize `2', parse(".")
			}
			 else {
				tokenize `3', parse(".")
			}
			local index "`3'"
			local catreg = "`3' `catreg'"
			local varx //nullify
			if "`by'" != "`comparator'" {
				*local catreg = "`comparator' `catreg'"
			}
		}
		if "`abnetwork'" != "" {
			tokenize `2', parse(".")
			local assignment "`3'"
			local catreg = "`3' `catreg'"
		}
		
		if "`interaction'" != "" & "`typevarx'" == "i" {
			local idpairconcat "#`varx'"
		}
		
		if "`typevarx'" == "i"  {
			if "`catreg'" == "" {
				local catreg = "`varx'"
			}
			else {
				if "`idpairconcat'" == "" & "`by'" !=  "`varx'" {
					local catreg "`varx' `catreg'"
				}
			}
		}
		else {
			if "`contreg'" == "" {
				local contreg = "`varx'"
			}
		}
		
		if "`model'" == "cbbetabin" {
			local at "at(mu==1)"
			local atexp "mu==1"
			local expression "expression(xb() - _b[_cons])"
		}
		else if "`link'" == "log" {
			local expression "expression(logit(predict(ir)))"
		}
		else {
			local expression "expression(predict(xb))"
		}
		
		if "`scimethod'"== "t" {
			local statistic "t"
		}
		else {
			local statistic "z"
		}
		/*
		if "`link'" == "cloglog" {
			local invfn "invcloglog"
		}
		else if "`link'" == "loglog" {
			local invfn "1 - invcloglog"
		}
		else if "`link'" == "log" {
			local invfn "exp"
		}
		else {
			local invfn "invlogit"
		}
		*/
		
		if "`idpairconcat'" != "" {
			local marginlist = `"`varx'"'
		}
		else {
			local marginlist
		}
		
		while "`catreg'" != ""  {
			tokenize `catreg'
			if ("`1'" != "`by'" & "`by'" != "") | "`by'" =="" {
				if ("`1'" != "`studyid'") {
				
					if "`idpairconcat'" != ""{
						local marginlist = `"`marginlist' `1'"'
					}
					local marginlist = `"`marginlist' `1'`idpairconcat'"'
				}
			}
			macro shift 
			local catreg `*'
		}
		
		estimates restore `estimates'
		if "`model'" == "random" {
			local df = e(N) -  e(k_f) - e(k_r)
		}
		else {
			local df = e(N) -  e(k)
		}

		mat `coefmat' = e(b)
		local predcmd = e(predict)
		
		if "`predcmd'" == "mepoisson_p" {
			local expression "expression(logit(exp(predict(xb))/`total'))"
		}
		
		local byncatreg 0
		if ("`by'" != "" & "`stratify'"  == "")  {
			margin , `expression' `at' over(`by') level(`level')
			
			
			mat `bycatregmatrixout' = r(table)'
			mat `byVmatrix' = r(V)
			mat `bycatregmatrixout' = `bycatregmatrixout'[1..., 1..6]
			
			local byrnames :rownames `bycatregmatrixout'
			local byncatreg = rowsof(`bycatregmatrixout')
		}
		
		if "`abnetwork'`mcbnetwork'`pcbnetwork'" == ""  {
			local grand "grand"
			local Overall "Overall"			
		}
		
		if "`comparative'" != "" & "`stratify'" != "" {
			local grand
			local Overall
		}
		local ncatreg 0
		
		//Overall
		if "`marginlist'" != "" | (  "`marginlist'" == "" & ("`mpair'`abnetwork'`mcbnetwork'`pcbnetwork'" == "" )) {
			margin `marginlist', `expression' `at' `grand' level(`level')
						
			mat `catregmatrixout' = r(table)'
			mat `Vmatrix' = r(V)
			mat `catregmatrixout' = `catregmatrixout'[1..., 1..6]
			
			local rnames :rownames `catregmatrixout'	
			local ncatreg = rowsof(`catregmatrixout')
		}
				
		local init 1
		local ncontreg 0
		local contrownames = ""
		if "`contreg'" != "" {
			foreach v of local contreg {
				summ `v', meanonly
				local vmean = r(mean)
				qui margin, `expression' at(`v'=`vmean' `atexp') level(`level')
				mat `matrixout' = r(table)'
				mat `matrixout' = `matrixout'[1..., 1..6]
				if `init' {
					local init 0
					mat `contregmatrixout' = `matrixout' 
				}
				else {
					mat `contregmatrixout' =  `contregmatrixout' \ `matrixout'
				}
				local contrownames = "`contrownames' `v'"
				local ++ncontreg
			}
		}
						
		//Stack the matrices
		if (`ncatreg' > 0 & `byncatreg' > 0) {
			mat `matrixout' =  `bycatregmatrixout' \ `catregmatrixout'
		}
		if (`ncatreg' > 0 & `byncatreg' == 0) {
			mat `matrixout' =  `catregmatrixout' 
		}
		
		if (`ncatreg' == 0 & `byncatreg' > 0) {
			mat `matrixout' =  `bycatregmatrixout' 
		}
		
		if (`ncontreg' > 0) {
			mat `matrixout' =  `contregmatrixout' \ `matrixout'
		}
		
		mat `tstats' = J(`=`byncatreg' + `ncatreg' + `ncontreg'', 3,.)
		
		mat `matrixout' = (`matrixout', `tstats')
		
		forvalues r = 1(1)`=`byncatreg' + `ncatreg' + `ncontreg''  {
				local tstat = `matrixout'[`r', 3]
				mat `matrixout'[`r', 7] = ttail(`df', abs(`tstat'))*2
				mat `matrixout'[`r', 8] = `matrixout'[`r', 1] - invttail((`df'), 0.5-`level'/200) * `matrixout'[`r', 2]
				mat `matrixout'[`r', 9] = `matrixout'[`r', 1] + invttail((`df'), 0.5-`level'/200) * `matrixout'[`r', 2]
		}
		
		local catrownames = ""
		if "`mcbnetwork'`pcbnetwork'" != "" {
			local rownamesmaxlen : strlen local Index
			local rownamesmaxlen = max(`rownamesmaxlen', 10)
		}
		else {
			local rownamesmaxlen = 10 /*Default*/
		}
		
		//# equations
		local init 0

		local rnames = "`byrnames' `rnames'" //attach the bynames
		
		//Except the grand rows	
		forvalues r = 1(1)`=`byncatreg' + `ncatreg'  - `="`grand'"!=""'' {
			//Labels
			local rname`r':word `r' of `rnames'
			tokenize `rname`r'', parse("#")					
			local left = "`1'"
			local right = "`3'"
			
			tokenize `left', parse(.)
			local leftv = "`3'"
			local leftlabel = "`1'"
			
			//no Interaction
			if "`right'" == "" {
				if "`leftv'" != "" {
					if strpos("`leftlabel'", "o") != 0 {
						local indexo = strlen("`leftlabel'") - 1
						local leftlabel = substr("`leftlabel'", 1, `indexo')
					}
					if strpos("`rname`r''", "b") != 0 {
						local leftlabel = ustrregexra("`leftlabel'", "bn", "")	
					}
					local lab:label `leftv' `leftlabel'
					local eqlab "`leftv'"
				}
				else {
					local lab "`leftlabel'"
					local eqlab ""
				}
				local nlencovl : strlen local llab
				local nlencov = `nlencovl' + 1					
			}
			else {
				//Interaction
				tokenize `right', parse(.)
				local rightv = "`3'"
				local rightlabel = "`1'"
				
				if strpos("`leftlabel'", "c") == 0 {
					if strpos("`leftlabel'", "o") != 0 {
						local indexo = strlen("`leftlabel'") - 1
						local leftlabel = substr("`leftlabel'", 1, `indexo')
					}
					if strpos("`leftlabel'", "b") != 0 {
						local leftlabel = ustrregexra("`leftlabel'", "bn", "")
					}
					local llab:label `leftv' `leftlabel'
				} 
				else {
					local llab
				}
				
				if strpos("`rightlabel'", "c") == 0 {
					if strpos("`rightlabel'", "o") != 0 {
						local indexo = strlen("`rightlabel'") - 1
						local rightlabel = substr("`rightlabel'", 1, `indexo')
					}
					if strpos("`rightlabel'", "b") != 0 {
						local rightlabel = ustrregexra("`rightlabel'", "bn", "")
					}
					local rlab:label `rightv' `rightlabel'
				} 
				else {
					local rlab
				}
				
				if (("`rlab'" != "") + ("`llab'" != "")) ==  0 {
					local lab = "`leftv'#`rightv'"
					local eqlab = ""
				}
				if (("`rlab'" != "") + ("`llab'" != "")) ==  1 {
					local lab = "`llab'`rlab'" 
					local eqlab = "`leftv'*`rightv'"
				}
				if (("`rlab'" != "") + ("`llab'" != "")) ==  2 {
					local lab = "`llab'|`rlab'" 
					local eqlab = "`leftv'*`rightv'"
				}
				local nlencovl : strlen local leftv
				local nlencovr : strlen local rightv
				local nlencov = `nlencovl' + `nlencovr' + 1
			}
			//check no underscore in the group names, replace with -
			if strpos("`lab'", "_") != 0 {
				local lab = ustrregexra("`lab'", "_", "-")
			}
			local lab = ustrregexra("`lab'", " ", "_")
			
			local nlenlab : strlen local lab
			if "`eqlab'" != "" {
				local nlencov = `nlencov'
			}
			else {
				local nlencov = 0
			}
			local rownamesmaxlen = max(`rownamesmaxlen', min(`=`nlenlab' + `nlencov' + 1', 32)) /*Check if there is a longer name*/
			local catrownames = "`catrownames' `eqlab':`lab'"
		}
		
		local rownames = "`contrownames' `catrownames' `Overall'"
		mat rownames `matrixout' = `rownames'
		mat colnames `matrixout' = Mean SE `statistic' P>|z| z_Lower z_Upper P>|t| t_Lower t_Upper

		//Get exact stats
		qui {
			gen `insample' = e(sample)
			local nrows = rowsof(`matrixout') //length of the vector
			local rnames :rownames `matrixout'
			local eqnames :roweq `matrixout'
			local newnrows = 0
			local mindex = 0
				
			foreach vari of local eqnames {		
				local ++mindex
				local group : word `mindex' of `rnames'
				
				//Skip if continous variable
				if (strpos("`vari'", "_") == 1) & ("`group'" != "Overall"){
					continue
				}
				
				cap drop `subset' 
				
				if "`group'" != "Overall" {
					if strpos("`vari'", "*") == 0 {
						cap drop `hold'
						noi decode `vari', gen(`hold')
						cap drop `subset'
						local latentgroup = ustrregexra("`group'", "_", " ")
						gen `subset' = 1 if `hold' == "`latentgroup'" & `insample' == 1 
					}
					else {
						tokenize `vari', parse("*")
						local leftvar = "`1'"
						local rightvar = "`3'"
						
						tokenize `group', parse("|")
						local leftgroup = "`1'"
						local rightgroup = "`3'"
						
						cap drop `holdleft' `holdright'
						decode `leftvar', gen(`holdleft')
						decode `rightvar', gen(`holdright')
						cap drop `subset'
						local latentleftgroup = ustrregexra("`leftgroup'", "_", " ")
						local latentrightgroup = ustrregexra("`rightgroup'", "_", " ")
						gen `subset' = 1 if (`holdleft' == "`latentleftgroup'") & (`holdright' == "`latentrightgroup'") & (`insample' == 1)
					}					
				}
				else {
					//All
					gen `subset' = 1 if `insample' == 1 
				}
				
				count if `subset' == 1 
				local nsubset = r(N)
				
				//Get the strata total
				sum `total' if `subset' == 1
				local stratatotal = r(sum)
				
				//Get the strata events
				sum `event' if `subset' == 1
				local strataevents = r(sum)
				
				//Get the strata size
				count if `subset' == 1
				local stratasize = r(N)
				
				tempvar ones zeros
				gen `ones' = `event'==`total' 
				sum `ones' if `subset' == 1 
				local sumones = r(sum)
				
				gen `zeros' = `event'==0 
				sum `zeros' if `subset' == 1  
				local sumzeros = r(sum)
				
				//Obtain the exact stats
				absexactci `stratatotal' `strataevents',  level(`level') //exact ci
				mat `absexact' = r(absexact)
				local modelp = `absexact'[1, 1]
				local postse = `absexact'[1, 2]
				local lowerp = `absexact'[1, 5]
				local upperp = `absexact'[1, 6]
				
				mat `exactabsouti' = (r(absexact), `stratatotal', `strataevents', `stratasize', `sumones', `sumzeros')
				mat rownames `exactabsouti' = `vari':`group'
				
				//Stack the matrices
				local ++newnrows
				if `newnrows' == 1 {
					mat `exactabsout' = `exactabsouti'	
				}
				else {
					mat `exactabsout' = `exactabsout'	\  `exactabsouti'
				}
			}
			mat colnames `exactabsout' = Mean SE z P>|z| Lower Upper Total Events Studies Ones Zeros
		}
		
		return matrix exactabsout = `exactabsout'		
		*return matrix outmatrixp = `outmatrixp'	
		return matrix outmatrix = `matrixout'
	end	

/*+++++++++++++++++++++++++	SUPPORTING FUNCTIONS:  postsim +++++++++++++++++++++++++
							Simulate and/or summarize posterior distribution
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/	
cap program drop postsim
program define postsim, rclass
	#delimit ;
	syntax  [if] [in], todo(string) orderid(varname) studyid(varname) estimates(name) 
	[bayesreps(string asis) event(varname) total(varname) rawest(name) rrout(name) orout(name) link(string) 
	modeles(varlist) modellci(varlist) modeluci(varlist) outplot(string) baselevel(integer 1) cov(string)
	model(string) comparative aliasdesign(string) by(varname) level(real 95) interaction abnetwork mcbnetwork mpair varx(varname) 
	stat(string) nsims(string) p(integer 0) catreg(varlist) by(varname) stratify]
	;
	#delimit cr 
	
	marksample touse, strok
		
	tempname betacoef rawcoef varrawcoef ///
				fullrawcoef fullvarrawcoef X beta sims ///
				popabsout popabsouti poprrout poprdout poprrouti poprdouti  poplrrout poplrrouti ///
				poporout poporouti poplorout poplorouti simvar absexact
	
	tempvar feff sfeff reff sreff reff1 sreff1 reff2 sreff2 eta insample ///
			newobs idpair gid rid hold holdleft holdright ///
			simmu sumphat meanphat subset subsetid subsetid1 sumphat1 ///
			meanphat1 gid1 modelse sumrrhat ///
			meanrrhat meanrdhat meanlrrhat sumorhat meanorhat meanlorhat sumlorhat varint varslope fisherrho ///
			simvarint simvarslope simfisherrho simrho covar simcovar lnsigma simlnsigma
	
	tokenize `modeles'
	local modelp "`1'"
	local modelrr "`2'"
	local  modelrd "`3'"
	local  modellrr  "`4'"
	local  modelor "`5'"
	local  modellor	"`6'"
	
	tokenize `modellci'
	local modelplci "`1'"
	local modelrrlci "`2'"
	local  modelrdlci "`3'"
	local  modellrrlci  "`4'"
	local  modelorlci "`5'"
	local  modellorlci	"`6'"
	
	tokenize `modeluci'
	local modelpuci "`1'"
	local modelrruci "`2'"
	local  modelrduci "`3'"
	local  modellrruci  "`4'"
	local  modeloruci "`5'"
	local  modelloruci	"`6'"
		
	if "`link'" == "cloglog" {
		local invfn "invcloglog"
	}
	else if "`link'" == "loglog" {
		if "`model'" == "crbetabin" {
			local invfn "exp(-exp(-"
			local closebracket "))"
		}
		else {
			local invfn "1-invcloglog"
		}
		local sign -
	}
	else if "`link'" == "log" {
		local invfn "exp"
	}
	else {
		local invfn "invlogit"
	}
	
	//if fixed, nullify covariances
	if "`model'" != "random" & "`cov'" != "" {
		local cov
	}
	
	//Restore 
	qui {
		estimates restore `estimates'
		gen `insample' = e(sample) /** mu*/
		
		if strpos("`model'", "bayes") == 0 {
			local predcmd = e(predict)
			
			//Coefficients estimates and varcov
			mat `fullrawcoef' = e(b)
			mat `fullvarrawcoef' = e(V)		
			
			local ncoef = colsof(`fullrawcoef')
			local rho = 0
			if "`model'" == "random" {
				if "`abnetwork'`cov'" == "" | "`cov'" =="commonslope"  {
					local nfeff = `=`ncoef' - 1'
					local varnames "`varint'"
					local simvarnames "`simvarint'"
				}
				if "`cov'" =="commonint" | "`cov'" =="freeint"  {
					local nfeff = `=`ncoef' - 1'
					local varnames "`varslope'"
					local simvarnames "`simvarslope'"
				}
				else if ("`abnetwork'" != "") | ("`cov'" =="independent") {
					local nfeff = `=`ncoef' - 2'
					local varnames "`varslope' `varint'"
					local simvarnames "`simvarslope' `simvarint'"
				}
				else if "`cov'" =="unstructured" {
					local nfeff = `=`ncoef' - 3'
					
					if "`predcmd'" == "meqrlogit_p" {
						local rho = tanh(`fullrawcoef'[1, `ncoef'])
						local varnames "`varslope' `varint' `fisherrho'"
						local simvarnames "`simvarslope'  `simvarint' `simfisherrho'"
					}
					else {
						local rho = `fullrawcoef'[1, `ncoef']/(sqrt(`fullrawcoef'[1, `=`ncoef'-1']*`fullrawcoef'[1, `=`ncoef'-2']))
						local varnames "`varslope' `varint' `covar'"
						local simvarnames "`simvarslope'  `simvarint' `simcovar'"
					}
				}
			}
			else if strpos("`model'", "betabin") == 1 {
				local nfeff = `=`ncoef' - 1'
				local varnames "`lnsigma'"
				local simvarnames "`simlnsigma'"
			}
			else {
				local nfeff = `ncoef'
			}
			//Get the FE parameters and their covariances
			mat `betacoef' = `fullrawcoef'[1, 1..`nfeff']
			//Predict		
			//Fill data if less than 7
			count
			local nobs = r(N)
			if ((`nobs' < 7) & ("`model'" == "random")) {
				local multipler = int(ceil(7/`nobs'))
				qui expand `multipler', gen(`newobs')
			}
			if "`model'" == "cbbetabin" {
				predictnl `feff' = xb() - _b[_cons], se(`sfeff')
			}
			else {
				if "`link'" == "log" {
					local offset "nooffset"
				}
				if "`predcmd'" == "mepoisson_p" {
					predictnl `feff' = log(exp(predict(xb))/`total'), se(`sfeff')  //logscale
				}
				else {
					//depends on the link used
					predict `feff' if `insample'==1, xb `offset' //FE
					predict `sfeff' if `insample'==1, stdp `offset' //se of FE	
				}
			}
			
			if "`model'" == "random" {
				if "`abnetwork'`cov'" == "" {
					predict `reff' if `insample'==1, reffects reses(`sreff')
					gen `modelse' = sqrt(`sreff'^2 + `sfeff'^2) if `insample'==1
				}
				else if ("`cov'" == "commonslope") {
					predict `reff2' if `insample'==1, reffects reses(`sreff2')
					gen `reff' = `reff2'
					gen `sreff' = `sreff2'	
					gen `modelse' = sqrt(`sreff'^2 + `sfeff'^2) if `insample'==1
				}
				else if ("`cov'" == "commonint") | ("`cov'" == "freeint") {
					predict `reff1' if `insample'==1, reffects reses(`sreff1')
					gen `reff' = `reff1'*2.`varx'
					gen `sreff' = 2.`varx'*`sreff1'	
					gen `modelse' = sqrt(`sreff'^2 + `sfeff'^2) if `insample'==1
				}
				else if "`abnetwork'" !="" | ("`cov'" == "unstructured") | ("`cov'" == "independent") {
					predict `reff1' `reff2' if `insample'==1, reffects reses(`sreff1' `sreff2')  //slope=1  int=2 
					
					if "`abnetwork'"  == "" {
						
						gen `reff' = `reff1'*2.`varx' + `reff2' 
						gen `sreff' = sqrt(`sreff1'^2 + `sreff2'^2)	
					}
					else {
						gen `reff' = `reff1' + `reff2' 
						gen `sreff' = sqrt(`sreff1'^2 + `sreff2'^2)	
					}
															
					gen `modelse' = sqrt(`sreff1'^2 + `sreff2'^2 + `sfeff'^2) if `insample'==1
				}
				
				gen `eta' = `feff' + `reff' //linear predictor
			}
			else {
				gen `eta' = `feff' if `insample'==1
				gen `modelse' = `sfeff' if `insample'==1
			}
					
			//Revert to original data if filler data was generated
			if (("`model'" == "random") & (`nobs' < 7))  {
				keep if !`newobs'
			}
			
			//Smooth p estimates
			replace `modelp' = `invfn'(`eta')`closebracket' if `insample'==1
			
			//identifiers
			sort `insample' `orderid'
			*egen `rid' = seq() if `insample'==1  //rowid
			
			if "`comparative'`mcbnetwork'`abnetwork'`mpair'" != ""  {
				egen `gid' = group(`studyid' `by') if `insample'==1  
				sort `gid' `orderid' `varx'
				by `gid': egen `idpair' = seq()
				egen `rid' = seq() if `insample'==1  //rowid
				
				if "`comparative'`mcbnetwork'`aliasdesign'`mpair'" != ""  {
					replace `modelrr' = `modelp'[_n] / `modelp'[_n-1] if (`gid'[_n]==`gid'[_n-1]) & (`idpair' == 2)
					replace `modelrd' = -`modelp'[_n] + `modelp'[_n-1] if (`gid'[_n]==`gid'[_n-1]) & (`idpair' == 2)
					replace `modellrr' = ln(`modelrr')
					replace `modelor' = (`modelp'[_n]/(1 - `modelp'[_n])) / (`modelp'[_n-1]/(1 - `modelp'[_n-1])) if (`gid'[_n]==`gid'[_n-1]) & (`idpair' == 2)
					replace `modellor' = ln(`modelor')
				}
			}
			else {
				egen `rid' = seq() if `insample'==1  //rowid
				gen `gid' = `rid'
			}
					
			//Generate designmatrix
			local colnames :colnames `betacoef'
			local nvars: word count `colnames'
			forvalues i=1(1)`nvars' {
				tempvar v`i' beta`i'
				
				local var`i' : word `i' of `colnames'

				local left
				local right
				local outcome
				local rightleft
				local rightright
				local leftleft
				local leftright
				
				if "`model'" == "cbbetabin" {
					//Split the term
					if strpos("`var`i''", "#") != 0 {
						tokenize `var`i'', parse("#")
						if "`5'" != "" {
							local left = "`1'"
							local right = "`3'"
							local outcome = "`5'"
						}
						else {
							local right = "`1'"
							local outcome = "`3'"
						}
						
						tokenize `outcome', parse(.)
						local outcome  = "`3'"
						
						tokenize `right', parse(.)
						local rightleft = "`1'"
						local rightright = "`3'"
						
						if "`left'" != "" {
							tokenize `left', parse(.)
							local leftleft = "`1'"
							local leftright = "`3'"
						}
					}
					else {
						local outcome "`var`i''"
					}
					
					//Constants
					if "`right'" == "" {
						cap confirm var `outcome'
						if _rc!=0  {	
							gen `v`i'' = 0
						}
						else {				
							gen `v`i'' = `outcome'
						}
					}
					
					//Main effects
					if "`right'" != "" & "`left'" == ""  {
						//Continous  
						if strpos("`rightleft'", "c") != 0 {
							gen `v`i'' = `outcome'*`rightright'
						}
						else {
							//Categorical
							if strpos("`rightleft'", "bn") != 0 {
								local rightleft = ustrregexra("`rightleft'", "bn", "")
							}
							if strpos("`rightleft'", "b") != 0 {
								local rightleft = ustrregexra("`rightleft'", "b", "")
							}		
							gen `v`i'' = 0 +  1*`outcome'*(`rightright' == `rightleft')
						}
					}
					
					//Interactions
					if "`left'" != "" {	
						//continous left
						if strpos("`leftleft'", "c") == 1 {
							local factorleft 0
						}
						else {
							local factorleft 1
						}
						
						//continous right
						if strpos("`rightleft'", "c") == 1 {
							local factorright 0
						}
						else {
							local factorright 1
						}
								
						local part 1
						local prefices "`leftleft' `rightleft'"
						foreach prefix of local prefices {
							//Continous  
							if strpos("`prefix'", "c") != 0 {
								local level`part' = `part'
							}
							else {
								//Categorical
								if strpos("`prefix'", "bn") != 0 {
									local level`part' = ustrregexra("`prefix'", "bn", "")
								}
								else if strpos("`prefix'", "b") != 0 {
									local level`part' = ustrregexra("`prefix'", "b", "")
								}
								else if strpos("`prefix'", "o") != 0 {
									local level`part' = ustrregexra("`prefix'", "o", "")
								}
								else {
									local level`part' = `prefix'
								}
							}
							local ++part
							
						}
						gen `v`i'' = ((`factorleft'*(`leftright'==`level1') + !`factorleft'*`leftright') * (`factorright'*(`rightright'==`level2') + !`factorright'*`rightright'))*`outcome'
					}
				}
				else {
					//Interaction
					tokenize `var`i'', parse("#")
					local left = "`1'"
					local right = "`3'"
					
					tokenize `left', parse(.)
					local leftleft = "`1'"
					local leftright = "`3'"
					
					//Constant or continous
					if "`right'" == ""  & "`leftright'" == "" {
						cap confirm var `leftleft'
						if _rc!=0  {	
							gen `v`i'' = 0
						}
						else {				
							gen `v`i'' = `leftleft'
						}
					}
					
					//Main categorical effects
					if "`right'" == "" & "`leftright'" != "" {
						if strpos("`leftleft'", "bn") != 0 {
							local leftleft = ustrregexra("`leftleft'", "bn", "")
						}
						if strpos("`leftleft'", "b") != 0 {
							local leftleft = ustrregexra("`leftleft'", "b", "")
						}		
						
						gen `v`i'' = 0 
						replace `v`i'' = 1 if `leftright' == `leftleft'
					}
					
					//Interactions
					if "`right'" != "" {
						tokenize `right', parse(.)
						local rightleft = "`1'"
						local rightright = "`3'"
					
						//continous left
						if strpos("`leftleft'", "c") == 1 {
							local factorleft 0
						}
						else {
							local factorleft 1
						}
						
						//continous right
						if strpos("`rightleft'", "c") == 1 {
							local factorright 0
						}
						else {
							local factorright 1
						}
						
						if `factorleft' == 1  {		
							//Categorical
							if strpos("`leftleft'", "bn") != 0 {
								local leftleft = ustrregexra("`leftleft'", "bn", "")
							}
							if strpos("`leftleft'", "b") != 0 {
								local leftleft = ustrregexra("`leftleft'", "b", "")
							}
							if strpos("`leftleft'", "o") != 0 {
								local leftleft = ustrregexra("`leftleft'", "o", "")
							}
							
							if `factorright' == 1 {
								if strpos("`rightleft'", "bn") != 0 {
									local rightleft = ustrregexra("`rightleft'", "bn", "")
								}
								if strpos("`rightleft'", "b") != 0 {
									local rightleft = ustrregexra("`rightleft'", "b", "")
								}
								if strpos("`rightleft'", "o") != 0 {
									local rightleft = ustrregexra("`rightleft'", "o", "")
								}
								
								gen `v`i'' = 0
								replace `v`i'' = 1 if (`leftright' == `leftleft') & (`rightright' == `rightleft')
						
							}
							else {
								gen `v`i'' = 0
								replace `v`i'' = 1*`rightright' if (`leftright' == `leftleft') 
							}
						}
						else {
							//Continous
							if `factorright' == 1 {
								if strpos("`rightleft'", "bn") != 0 {
									local rightleft = ustrregexra("`rightleft'", "bn", "")
								}
								if strpos("`rightleft'", "b") != 0 {
									local rightleft = ustrregexra("`rightleft'", "b", "")
								}
								gen `v`i'' = 0
								replace `v`i'' = 1*`leftright' if (`rightright' == `rightleft')
								
							}
							else {
								gen `v`i'' = .
								replace `v`i'' = `leftright'*`rightright'
							}
						}	
					}
				}
				local vnamelist "`vnamelist' `v`i''"
				local bnamelist "`bnamelist' `beta`i''"
			}
			
			if "`model'" == "random" | strpos("`model'", "betabin") == 1 {
				//Add varnames
				local bnamelist "`bnamelist' `varnames'"
			}
			set matsize `nsims'

			//make matrices from the dataset
			//roweq(`idpair')
			mkmat `vnamelist' if `insample'==1, matrix(`X')  rownames(`rid')
						
			tempvar present
			gen `present' = 1	
			
			//Simulate the parameters
			if `nobs' < `nsims' {
				set obs `nsims'
			}
			
			drawnorm `bnamelist', n(`nsims') cov(`fullvarrawcoef') means(`fullrawcoef') seed(1)

			mkmat `bnamelist', matrix(`beta')
			
			//Subset the matrix
			if `ncoef' > `nfeff' {
				mat `simvar' = `beta'[1..`nsims', `=`nfeff'+1'..`ncoef']
				mat `beta' = `beta'[1..`nsims', 1..`nfeff']
			}
			
			mat `sims' = `beta'*`X''
			
			//Construct the names
			local ncols = colsof(`sims') //length of the vector
			local cnames :colnames `sims'
			
			local matcolnames
			
			forvalues c=1(1)`ncols' {
				local matrid : word `c' of `cnames'

				tempvar festudy`matrid'
				local matcolnames = "`matcolnames' `festudy`matrid''"
			}
			
			if "`model'" == "random" {
				//Append the var matrix
				mat `sims' = (`sims', `simvar')
				
				//Add varnames
				local matcolnames "`matcolnames' `simvarnames'"
			}
			
			//pass the names
			matname `sims' `matcolnames', col(.) explicit

			//Bring the matrix to the dataset
			svmat `sims', names(col)
			
			if "`model'" == "random" {
				if ("`abnetwork'" !="" | "`cov'" !="")  {
					if "`cov'" !="unstructured" {
						gen `simrho' = 0
					}
					else {
						if "`predcmd'" == "meqrlogit_p" {
							gen `simrho' = tanh(`simfisherrho')
							replace `sreff1' = exp(`simvarslope') //Marginal se
							replace `sreff2' = sqrt((1 - (`simrho')^2)*(exp(`simvarint')^2)) //Conditional se
						}
						else {
							//Truncate the values to zero
							replace `simvarslope' = 0 if `simvarslope' < 0
							replace `simvarint' = 0 if `simvarint' < 0
							replace `simcovar' = 0 if `simvarslope' == 0 & `simvarint' == 0
							
							gen `simrho' = `simcovar'/sqrt(`simvarint'*`simvarslope')
							replace `sreff1' = sqrt(`simvarslope') //Marginal se
							replace `sreff2' = sqrt((1 - (`simrho')^2)*(`simvarint')) //Conditional se
						}
					}
					
					if "`predcmd'" == "meqrlogit_p" {
						if "`cov'" == "commonslope" | "`cov'" == "independent" | "`abnetwork'" !=""  {	
							replace `sreff2' = exp(`simvarint') //Marginal se
							
							if "`cov'" == "commonslope"  {
								gen `sreff1' = 0
								gen `reff1' = 0
							}
						}
						
						if "`cov'" == "commonint"| "`cov'" == "freeint"  | "`cov'" == "independent" | "`abnetwork'" !="" {	
							replace `sreff1' = exp(`simvarslope') //Marginal se
							
							if "`cov'" == "commonint" | "`cov'" == "freeint"  {
								gen `sreff2' = 0
								gen `reff2' = 0
							}
						}						
					}
					else {					
						if "`cov'" == "commonslope" | "`cov'" == "independent" | "`abnetwork'" !=""  {	
							replace `simvarint' = 0 if `simvarint' < 0
							replace `sreff2' = sqrt(`simvarint') //Marginal se
							
							if "`cov'" == "commonslope" {
								gen `sreff1' = 0
								gen `reff1' = 0
							}
						}
						
						if "`cov'" == "commonint" | "`cov'" == "freeint"  | "`cov'" == "independent" | "`abnetwork'" !="" {
							replace `simvarslope' = 0 if `simvarslope' < 0						
							replace `sreff1' = sqrt(`simvarslope') //Marginal se
							
							if "`cov'" == "commonint" | "`cov'" == "freeint" {
								gen `sreff2' = 0
								gen `reff2' = 0
							}
						}
					}
				}
				else {
					//Truncate the values to zero
					if "`predcmd'" == "melogit_p" {
						replace `simvarint' = 0 if `simvarint' < 0
					}
				}
			}
			
			//# of obs
			count if `insample' == 1
			local nobs = r(N)
			
			//Generate the p's and r's
			forvalues j=1(1)`nobs' { 
				tempvar phat`j' 
					
				if "`comparative'`mcbnetwork'`abnetwork'`mpair'" == "" {
					*tempvar phat`j' 
				
					if "`model'" == "random" {
						tempvar restudy`j'
						//EB re 
						sum `reff' if `rid' == `j' 
						local reff_`j' = r(mean)
						
						if "`predcmd'" == "meqrlogit_p" {
							gen `restudy`j'' = rnormal(0, exp(`simvarint'))
						}
						else {
							gen `restudy`j'' = rnormal(0, sqrt(`simvarint'))
						}
						
						gen `phat`j'' = `invfn'(`reff_`j'' + `restudy`j'' + `festudy`j'')`closebracket'
					}
					else {
						gen `phat`j'' = `invfn'(`festudy`j'')`closebracket'
					}
				}
				if "`comparative'`mcbnetwork'`abnetwork'`mpair'" != "" {
					sum `gid' if `rid' == `j'
					local index = r(min)
					
					sum `idpair' if `rid' == `j'
					local pair = r(min)
					
					//Generate the variables
					tempvar phat_`pair'`index'
					
					if "`model'" == "random" {
						//EB re 
						sum `reff' if `rid' == `j' 
						local reff_`j' = r(mean)
							
						if `pair' == 1 {
							//re - same per study				
							tempvar restudy`index'
							
							if "`abnetwork'`cov'" == "" {
								if "`predcmd'" == "meqrlogit_p" {
									gen `restudy`index'' = rnormal(0, exp(`simvarint'))
								}
								else {
									gen `restudy`index'' = rnormal(0, sqrt(`simvarint'))
								}
							}
							else if "`abnetwork'" !="" | "`cov'" !="" {
								replace `reff1' = rnormal(0, `sreff1')
								
								if "`cov'" =="unstructured" { 
									replace `reff2' = rnormal(`simrho'*`sreff2'*(`reff1'/`sreff1'), `sreff2')
								}
								else {
									replace `reff2' = rnormal(0, `sreff2')
								}
								
								gen `restudy`index'' = `reff1' + `reff2'
							}
						}					
						gen `phat`j'' = `invfn'(`reff_`j'' + `restudy`index'' +  `festudy`j'')`closebracket'
					}
					else {
						gen `phat`j'' = `invfn'(`festudy`j'')`closebracket'
					}
										
					if "`comparative'`mcbnetwork'`aliasdesign'`mpair'" != ""  {
						//Create the pairs
						gen `phat_`pair'`index'' = `phat`j''
						
						if `pair' == 2 {
							tempvar rrhat`index' rdhat`index' lrrhat`index' orhat`index' lorhat`index'
							
							gen `rrhat`index'' = `phat_2`index'' / `phat_1`index''
							gen `rdhat`index'' = -`sign'`phat_2`index'' + `sign'`phat_1`index''
							gen `lrrhat`index''  = ln(`rrhat`index'')
							gen `orhat`index'' = (`phat_2`index''/(1 - `phat_2`index'')) / (`phat_1`index'' /(1 - `phat_1`index''))
							gen `lorhat`index'' = ln(`orhat`index'')
						}
					}
				}
			}
		}	
		else {
			//identifiers
			gsort -`insample' `orderid' `varx'
			egen `rid' = seq() if `insample'==1  //rowid
			
			if "`comparative'`mcbnetwork'`abnetwork'`mpair'" != ""  {
				egen `gid' = group(`studyid' `by') if `insample'==1  
				sort `gid' `orderid' `varx'
				by `gid': egen `idpair' = seq()
			}
			else {
				gen `gid' = `rid'
			}
			
			//Generate the bayesian parameters
			tempvar present
			gen `present' = 1
			
			//merge with simulations
			sort `rid'
			merge 1:1 _n  using `bayesreps',  nogenerate 
						
			//# of obs
			count if `insample' == 1
			local nobs = r(N)
			
			//Generate the p's and r's
		
			forvalues j=1(1)`nobs' { 
				tempvar phat`j'  
			
				//total 
				sum `total' if `rid' == `j' 
				local total_`j' = r(mean)
											
				gen `phat`j'' = _mu1_`j'/`total_`j'' //most important value
	
				if "`comparative'`mcbnetwork'`abnetwork'`mpair'" != "" {
					sum `gid' if `rid' == `j'
					local index = r(min)
					
					sum `idpair' if `rid' == `j'
					local pair = r(min)
					
					//Generate the variables
					tempvar phat_`pair'`index'
					
					if "`comparative'`mcbnetwork'`aliasdesign'`mpair'" != ""  {
						gen `phat_`pair'`index'' = `phat`j''
					}
				}	
			}
			
			forvalues j=1(1)`nobs' { 
				if "`comparative'`mcbnetwork'`abnetwork'`mpair'" != "" {
					sum `gid' if `rid' == `j'
					local index = r(min)
					
					sum `idpair' if `rid' == `j'
					local pair = r(min)
										
					if "`comparative'`mcbnetwork'`aliasdesign'`mpair'" != ""  {						
						if `pair' == 2 {
							tempvar rrhat`index' rdhat`index' lrrhat`index' orhat`index' lorhat`index'
							
							gen `rrhat`index'' = `phat_2`index'' / `phat_1`index''
							gen `rdhat`index'' = -`sign'`phat_2`index'' + `sign'`phat_1`index''
							gen `lrrhat`index'' = ln(`rrhat`index'')
							gen `orhat`index'' = (`phat_2`index''/(1 - `phat_2`index'')) / (`phat_1`index'' /(1 - `phat_1`index''))
							gen `lorhat`index'' = ln(`orhat`index'')
						}
					}
				}	
			}
			
		}
		
		//Summarize	
		if "`todo'" == "p" {
			//Summarize p
			local nrows = rowsof(`rawest') //length of the vector
			local rnames :rownames `rawest'
			local eqnames :roweq `rawest'
			local newnrows = 0
			local mindex = 0
			
			
			if strpos("`model'", "bayes") == 1 {
				//Add overall
				if `nrows' > 1 {
					local eqnames = "`eqnames' _"
					local rnames = "`rnames' Overall"
				}
			
				//Add main effects if absent
				if "`interaction'" != "" {
					local catvars "`catreg' `varx'"
					foreach c of local catvars {
						*qui label list `c'
						*local nlevels = r(max)
						
						qui levelsof `c', local(codelevels)
						local nlevels = r(r)
						
						foreach l of local codelevels {
						*forvalues l = 1/`nlevels' {
							local lab:label `c' `l'
							local lab = ustrregexra("`lab'", " ", "_")
							local eqnames = "`c' `eqnames'"
							local rnames = "`lab' `rnames'"						
						}
					}
				}
				
				//Add by
				if ("`by'" != "" & "`stratify'"  == "")  {
					*qui label list `by'
					*local nlevels = r(max)
					
					qui levelsof `by', local(codelevels)
					local nlevels = r(r)
						
					foreach l of local codelevels {
						local lab:label `by' `l'
						local lab = ustrregexra("`lab'", " ", "_")
						local eqnames = "`by' `eqnames'"
						local rnames = "`lab' `rnames'"						
					}
				}
			}
	
			foreach vari of local eqnames {		
				local ++mindex
				local group : word `mindex' of `rnames'
				
				//Skip if continous variable
				if (strpos("`vari'", "_") == 1) & ("`group'" != "Overall") & "`mpair'" == "" {
					continue
				}
				
				cap drop `subset' `subsetid'
				
				if "`group'" != "Overall" {
					if strpos("`vari'", "*") == 0 {
						cap drop `hold'
						decode `vari', gen(`hold')
						cap drop `subset'
						local latentgroup = ustrregexra("`group'", "_", " ")
						gen `subset' = 1 if `hold' == "`latentgroup'" & `insample' == 1 
					}
					else {
						tokenize `vari', parse("*")
						local leftvar = "`1'"
						local rightvar = "`3'"
						
						tokenize `group', parse("|")
						local leftgroup = "`1'"
						local rightgroup = "`3'"
						
						cap drop `holdleft' `holdright'
						decode `leftvar', gen(`holdleft')
						decode `rightvar', gen(`holdright')
						cap drop `subset'
						local latentleftgroup = ustrregexra("`leftgroup'", "_", " ")
						local latentrightgroup = ustrregexra("`rightgroup'", "_", " ")
						gen `subset' = 1 if (`holdleft' == "`latentleftgroup'") & (`holdright' == "`latentrightgroup'") & (`insample' == 1)
					}					
					egen `subsetid' = seq() if `subset' == 1
				}
				else {
					//All
					gen `subset' = 1 if `insample' == 1 
					gen `subsetid' = `rid'
				}
				
				count if `subset' == 1 
				local nsubset = r(N)
				
				//Get the strata total
				sum `total' if `subset' == 1
				local stratatotal = r(sum)
				
				//Get the strata events
				sum `event' if `subset' == 1
				local strataevents = r(sum)
								
				tempvar ones zeros
				gen `ones' = `event'==`total' 
				sum `ones' if `subset' == 1 
				local sumones = r(sum)
				
				gen `zeros' = `event'==0 
				sum `zeros' if `subset' == 1  
				local sumzeros = r(sum)
				
				cap drop `sumphat' `meanphat'
				
				local plistvar
				forvalues j=1(1)`nsubset' { 
				
					sum `rid' if `subsetid' == `j'
					local index = r(min)
					
					local plistvar = "`plistvar' `phat`index''"

					//Replace ones/zeros if seperated				
					if ((`sumones' == `nsubset') | (`sumzeros' == `nsubset')) & `p'== 1 & strpos("`model'", "bayes") == 0 {
						replace `phat`index'' = `strataevents'/`stratatotal'	
					}
					
					if `j'== 1 {
						gen `sumphat' = `phat`index''	
					}
					else {					
						replace `sumphat' = `sumphat' + `phat`index''
					}
				}
				
				if strpos("`model'", "bayes") == 0 {
					//Obtain mean of modelled estimates
					sum `modelp' if `subset' == 1
					local meanmodelp = r(mean)
				}
								
				*gen `meanphat' = `sumphat'/`nsubset'
				egen `meanphat' = rowmean(`plistvar')
							
				//Standard error
				sum `meanphat'	
				local postse = r(sd)
				local postmean =  r(mean)
				
				//Obtain the quantiles
				centile `meanphat', centile(50 `=(100-`level')/2' `=100 - (100-`level')/2')
				local median = r(c_1) //Median
				local lowerp = r(c_2) //Lower centile
				local upperp = r(c_3) //Upper centile
				local nreps = r(N)
				
				if strpos("`model'", "bayes") == 0 {
					mat `popabsouti' = (`meanmodelp', `postse', `median', `lowerp', `upperp', `nreps')
				}
				else {
					mat `popabsouti' = (`postmean', `postse', `median', `lowerp', `upperp', `nreps')
				}
				mat rownames `popabsouti' = `vari':`group'
				
				//Stack the matrices
				local ++newnrows
				if `newnrows' == 1 {
					mat `popabsout' = `popabsouti'	
				}
				else {
					mat `popabsout' = `popabsout'	\  `popabsouti'
				}
			}
		}
		
		if "`todo'" == "r" {
			//Summarize RR
			local nrows = rowsof(`rrout') //length of the vector
			local rnames :rownames `rrout'
			local eqnames :roweq  `rrout'
			local newnrows 0
						
			if "`comparative'`mcbnetwork'`aliasdesign'`mpair'" == "" {
				local catvars : list uniq eqnames	
				foreach vari of local catvars {
					
					cap drop `hold'	
					decode `vari' if `insample' == 1, gen(`hold')
					
					*label list `vari'
					*local ngroups = r(max)
					levelsof `vari' if `insample' == 1, local(groupcodes)
					local ngroups = r(r)
			
					local baselab:label `vari' `baselevel'
					
					//count in basegroup
					tempvar meanphat`baselevel' meanrrhat`baselevel' meanrdhat`baselevel' meanlrrhat`baselevel' meanorhat`baselevel' ///
							meanlorhat`baselevel' gid`baselevel' sumphat`baselevel' subsetid`baselevel'
							
					tempname poprrouti`baselevel' poprdouti`baselevel' poplrrouti`baselevel' poporouti`baselevel'  poplorouti`baselevel'
					
					count if `vari' == `baselevel' & `insample' == 1
					local ngroup`baselevel' = r(N)
					
					egen `subsetid`baselevel'' = group(`rid') if `vari' == `baselevel' & `insample' == 1
					
					cap drop `sumphat`baselevel'' `meanphat`baselevel''
					
					//Get the strata total
					sum `total' if `vari' == `baselevel' & `insample' == 1
					local stratatotal = r(sum)
					
					//Get the strata events
					sum `event' if `vari' == `baselevel' & `insample' == 1
					local strataevents = r(sum)
					
					//Get the strata size
					count if `vari' == `baselevel' & `insample' == 1
					local stratasize = r(N)
					
					tempvar ones zeros
					gen `ones' = `event'==`total' 
					sum `ones' if `vari' == `baselevel' & `insample' == 1
					local sumones = r(sum)
					
					gen `zeros' = `event'==0 
					sum `zeros' if `vari' == `baselevel' & `insample' == 1 
					local sumzeros = r(sum)
					
					//basegroup				
					forvalues j=1(1)`ngroup`baselevel'' {
						sum `rid' if `subsetid`baselevel'' == `j'
						local index = r(min)
						
						//Replace ones/zeros if seperated				
						if ((`sumones' == `stratasize') | (`sumzeros' == `stratasize')) & `p'== 1  {
							replace `phat`index'' = `strataevents'/`stratatotal'	
						}
						
						if 	`j' == 1 {
							gen `sumphat`baselevel''  = `phat`index''
						}
						else {						
							replace `sumphat`baselevel'' = `sumphat`baselevel'' + `phat`index''
						}
						
					}
					gen `meanphat`baselevel'' = `sumphat`baselevel''/`ngroup`baselevel''
					
					if strpos("`model'", "bayes") == 0 {
						sum `modelp' if `vari' == `baselevel' & `insample' == 1
						local meanmodelp`baselevel' = r(mean)
					}
					else {
						sum `meanphat`baselevel''
						local meanmodelp`baselevel' =  r(mean)
					}
					
					mat `poprrouti`baselevel'' = (1, 0, 1, 1, 1, .)
					mat `poprdouti`baselevel'' = (0, 0, 0, 0, 0, .)
					mat `poplrrouti`baselevel'' = (1, 0, 1, 1, 1, .)
					mat `poporouti`baselevel'' = (1, 0, 1, 1, 1, .)
					mat `poplorouti`baselevel'' = (1, 0, 1, 1, 1, .)
					
					local baselab = ustrregexra("`baselab'", " ", "_")
					mat rownames `poprrouti`baselevel'' = `vari':`baselab'
					mat rownames `poprdouti`baselevel'' = `vari':`baselab'
					mat rownames `poplrrouti`baselevel'' = `vari':`baselab'
					mat rownames `poporouti`baselevel'' = `vari':`baselab'
					mat rownames `poplorouti`baselevel'' = `vari':`baselab'
					
					//Other groups
					foreach g of local groupcodes {
					*forvalues g=1(1)`ngroups' {
						*local g: word `l' of `groupcodes'
						if `g' != `baselevel' {
							tempvar meanphat`g' meanrrhat`g' meanrdhat`g' meanlrrhat`g' meanorhat`g' meanlorhat`g' gid`g' sumphat`g' subsetid`g'
							tempname poprrouti`g' poprdouti`g' poplrrouti`g' poporouti`g' poplorouti`g'
							
							local glab:label `vari' `g'
							count if `vari' == `g' & `insample' == 1
							local ngroup`g' = r(N)	
							egen `subsetid`g'' = group(`rid') if `vari' == `g' & `insample' == 1
							
							//Get the strata total
							sum `total' if `vari' == `g' & `insample' == 1
							local stratatotal = r(sum)
							
							//Get the strata events
							sum `event' if `vari' == `g' & `insample' == 1
							local strataevents = r(sum)
							
							//Get the strata size
							count if `vari' == `g' & `insample' == 1
							local stratasize = r(N)
							
							tempvar ones zeros
							gen `ones' = `event'==`total' 
							sum `ones' if `vari' == `g' & `insample' == 1
							local sumones = r(sum)
							
							gen `zeros' = `event'==0 
							sum `zeros' if `vari' == `g' & `insample' == 1 
							local sumzeros = r(sum)
														
							//Group of interest
							forvalues j=1(1)`ngroup`g'' {
								sum `rid' if `subsetid`g'' == `j'
								local index = r(min)
								
								//Replace ones/zeros if seperated				
								if ((`sumones' == `stratasize') | (`sumzeros' == `stratasize')) & `p'== 1 & strpos("`model'", "bayes") == 0  {
									replace `phat`index'' = `strataevents'/`stratatotal'	
								}
								
								if `j' == 1{
									gen `sumphat`g'' = `phat`index''
								}
								else {
									replace `sumphat`g'' = `sumphat`g'' + `phat`index''
								}
							}
							
							gen `meanphat`g'' = `sumphat`g''/`ngroup`g''
							
							//Generate R 
							gen `meanrrhat`g'' = `meanphat`g'' / `meanphat`baselevel''
							gen `meanrdhat`g'' = - `meanphat`g'' + `meanphat`baselevel''
							gen `meanlrrhat`g'' = ln(`meanrrhat`g'')
							gen `meanorhat`g'' = (`meanphat`g''/(1 - `meanphat`g'')) / (`meanphat`baselevel''/(1 - `meanphat`baselevel''))
							gen `meanlorhat`g'' = ln(`meanorhat`g'')

							//Obtain mean of modelled estimates
							if strpos("`model'", "bayes") == 0 {
								sum `modelp' if `vari' == `g' & `insample' == 1
								local meanmodelp`g' = r(mean)
								
								local meanmodelrr`g' = `meanmodelp`g'' / `meanmodelp`baselevel''
								local meanmodelrd`g' = - `sign'`meanmodelp`g'' + `sign'`meanmodelp`baselevel''
								local meanmodellrr`g' = ln(`meanmodelp`g'' / `meanmodelp`baselevel'')
								local meanmodelor`g' = (`meanmodelp`g''/(1 - `meanmodelp`g'')) / (`meanmodelp`baselevel''/(1 - `meanmodelp`baselevel''))
								local meanmodellor`g' = ln((`meanmodelp`g''/(1 - `meanmodelp`g'')) / (`meanmodelp`baselevel''/(1 - `meanmodelp`baselevel'')))
							}
							
							//Standard error
							sum `meanrrhat`g''
							local postserr = r(sd)
							local postmeanrr = r(mean)
							
							sum `meanrdhat`g''
							local postserd = r(sd)
							local postmeanrd = r(mean)
							
							sum `meanlrrhat`g''
							local postselrr = r(sd)
							local postmeanlrr = r(mean)
							
							sum `meanorhat`g''
							local postseor = r(sd)
							local postmeanor = r(mean)
							
							sum `meanlorhat`g''
							local postselor = r(sd)
							local postmeanlor = r(mean)
							
							//Obtain the quantiles
							centile `meanrrhat`g'', centile(50 `=(100-`level')/2' `=100 - (100-`level')/2')
							local medianrr = r(c_1) //Median
							local lowerprr = r(c_2) //Lower centile
							local upperprr = r(c_3) //Upper centile
							local nrrreps = r(N)
							
							centile `meanrdhat`g'', centile(50 `=(100-`level')/2' `=100 - (100-`level')/2')
							local medianrd = r(c_1) //Median
							local lowerprd = r(c_2) //Lower centile
							local upperprd = r(c_3) //Upper centile
							local nrdreps = r(N)
							
							centile `meanlrrhat`g'', centile(50 `=(100-`level')/2' `=100 - (100-`level')/2')
							local medianlrr = r(c_1) //Median
							local lowerplrr = r(c_2) //Lower centile
							local upperplrr = r(c_3) //Upper centile
							local nlrrreps = r(N)
							
							centile `meanorhat`g'', centile(50 `=(100-`level')/2' `=100 - (100-`level')/2')
							local medianor = r(c_1) //Median
							local lowerpor = r(c_2) //Lower centile
							local upperpor = r(c_3) //Upper centile
							local norreps = r(N)
							
							centile `meanlorhat`g'', centile(50 `=(100-`level')/2' `=100 - (100-`level')/2')
							local medianlor = r(c_1) //Median
							local lowerplor = r(c_2) //Lower centile
							local upperplor = r(c_3) //Upper centile
							local nlorreps = r(N)
							
							if strpos("`model'", "bayes") == 0 { 						
								mat `poprrouti`g'' = (`meanmodelrr`g'', `postserr', `medianrr', `lowerprr', `upperprr', `nrrreps')
								mat `poprdouti`g'' = (`meanmodelrd`g'', `postserd', `medianrd', `lowerprd', `upperprd', `nrdreps')
								mat `poplrrouti`g'' = (`meanmodellrr`g'', `postselrr', `medianlrr', `lowerplrr', `upperplrr', `nlrrreps')
								mat `poporouti`g'' = (`meanmodelor`g'', `postseor', `medianor', `lowerpor', `upperpor', `norreps')
								mat `poplorouti`g'' = (`meanmodellor`g'', `postselor', `medianlor', `lowerplor', `upperplor', `nlorreps')
							}
							else {
								mat `poprrouti`g'' = (`postmeanrr', `postserr', `medianrr', `lowerprr', `upperprr',  `nrrreps')
								mat `poprdouti`g'' = (`postmeanrd', `postserd', `medianrd', `lowerprd', `upperprd',  `nrdreps')
								mat `poplrrouti`g'' = (`postmeanlrr', `postselrr', `medianlrr', `lowerplrr', `upperplrr', `nlrrreps')
								mat `poporouti`g'' = (`postmeanor', `postseor', `medianor', `lowerpor', `upperpor', `norreps')
								mat `poplorouti`g'' = (`postmeanlor', `postselor', `medianlor', `lowerplor', `upperplor', `nlorreps')
							}
							
							local glab = ustrregexra("`glab'", " ", "_")
							mat rownames `poprrouti`g'' = `vari':`glab'
							mat rownames `poprdouti`g'' = `vari':`glab'
							mat rownames `poplrrouti`g'' = `vari':`glab'
							mat rownames `poporouti`g'' = `vari':`glab'
							mat rownames `poplorouti`g'' = `vari':`glab'
						}
						if `g' == 1 {
							mat `poprrouti' = `poprrouti`g''
							mat `poprdouti' = `poprdouti`g''
							mat `poplrrouti' = `poplrrouti`g''
							mat `poporouti' = `poporouti`g''
							mat `poplorouti' = `poplorouti`g''
						}
						else {
							//Stack the matrices
							mat `poprrouti' = `poprrouti'	\  `poprrouti`g''
							mat `poprdouti' = `poprdouti'	\  `poprdouti`g''
							mat `poplrrouti' = `poplrrouti'	\  `poplrrouti`g''
							mat `poporouti' = `poporouti'	\  `poporouti`g''
							mat `poplorouti' = `poplorouti'	\  `poplorouti`g''
						}
					}
					//Stack the matrices
					local ++newnrows
					if `newnrows' == 1 {
						mat `poprrout' = `poprrouti'
						mat `poprdout' = `poprdouti'
						mat `poplrrout' = `poplrrouti'
						mat `poporout' = `poporouti'
						mat `poplorout' = `poplorouti'
					}
					else {
						mat `poprrout' = `poprrout'	\  `poprrouti'
						mat `poprdout' = `poprdout'	\  `poprdouti'
						mat `poplrrout' = `poplrrout'	\  `poplrrouti'
						mat `poporout' = `poporout'	\  `poporouti'
						mat `poplorout' = `poplorout'	\  `poplorouti'
					}
				}
			}
			
			if "`comparative'`mcbnetwork'`aliasdesign'`mpair'" != "" {
				//Comparative R
				local mindex 0
				local newnrows 0
				
				if strpos("`model'", "bayes") == 1 {
					//Add overall
					if `nrows' > 1 {
						local eqnames = "`eqnames' _"
						local rnames = "`rnames' Overall"
					}
				}
								
				foreach vari of local eqnames {		
					local ++mindex
					local group : word `mindex' of `rnames'
					
					cap drop `subset'
					if "`group'" != "Overall" {
						cap drop `hold'
						decode `vari', gen(`hold')
						
						local latentgroup = ustrregexra("`group'", "_", " ")
						gen `subset' = 1 if `hold' == "`latentgroup'"  & `insample' == 1
					}
					else {
						//All
						gen `subset' = 1  if `insample' == 1
					}
					
					cap drop `subsetid'
					egen `subsetid' = seq()	 if `subset' == 1
					
					count if `subset' == 1 & `idpair' == 2
					local nsubset = r(N)
					
					//Compute mean of simulated values
					/*
					cap drop `sumrrhat' `sumorhat' `sumlorhat'
					gen `sumrrhat' = 0
					gen `sumorhat' = 0
					gen `sumlorhat' = 0
					*/
					
					local rrlistvar
					local rdlistvar
					local lrrlistvar
					local orlistvar
					local lorlistvar
		
					forvalues j=1(1)`nobs' { 
						sum `gid' if `subsetid' == `j'
						local index = r(min)
						
						sum `idpair' if `subsetid' == `j'
						local pair = r(min)
						
						if `pair' == 2 {
							local rrlistvar = `"`rrlistvar' `rrhat`index''"'
							local rdlistvar = `"`rdlistvar' `rdhat`index''"'
							local lrrlistvar = `"`lrrlistvar' `lrrhat`index''"'
							local orlistvar = `"`orlistvar' `orhat`index''"'
							local lorlistvar = `"`lorlistvar' `lorhat`index''"'
																
							/*
							replace `sumrrhat' = `sumrrhat' + `rrhat`index''
							replace `sumorhat' = `sumorhat' + `orhat`index''
							replace `sumlorhat' = `sumlorhat' + `lorhat`index''
							*/
						}
					}
										
					//Obtain mean of modelled estimates
					if strpos("`model'", "bayes") == 0 {

						sum `modelrr' if `subset' == 1
						local meanmodelrr = r(mean)
						
						sum `modelrd' if `subset' == 1
						local meanmodelrd = r(mean)
						
						sum `modellrr' if `subset' == 1
						local meanmodellrr = r(mean)
						
						sum `modelor' if `subset' == 1
						local meanmodelor = r(mean)
						
						sum `modellor' if `subset' == 1
						local meanmodellor = r(mean)
					}
					
					cap drop `meanrrhat'  `meanrdhat' `meanlrrhat' `meanorhat' `meanlorhat'
					/*
					gen `meanrrhat' = `sumrrhat'/`nsubset'
					gen `meanorhat' = `sumorhat'/`nsubset'
					gen `meanlorhat' = `sumlorhat'/`nsubset'
					*/
					
					egen `meanrrhat' = rowmean(`rrlistvar')
					egen `meanrdhat' = rowmean(`rdlistvar')
					egen `meanlrrhat' = rowmean(`lrrlistvar')
					egen `meanorhat' = rowmean(`orlistvar')
					egen `meanlorhat' = rowmean(`lorlistvar')

					//Standard error
					sum `meanrrhat'	
					local postserr = r(sd)
					local postmeanrr = r(mean)
					
					sum `meanrdhat'	
					local postserd = r(sd)
					local postmeanrd = r(mean)
					
					sum `meanlrrhat'	
					local postselrr = r(sd)
					local postmeanlrr = r(mean)
					
					sum `meanorhat'	
					local postseor = r(sd)
					local postmeanor = r(mean)
					
					sum `meanlorhat'	
					local postselor = r(sd)
					local postmeanlor = r(mean)
					
					//Obtain the quantiles
					centile `meanrrhat' , centile(50 `=(100-`level')/2' `=100 - (100-`level')/2')
					local medianrr = r(c_1) //Median
					local lowerprr = r(c_2) //Lower centile
					local upperprr = r(c_3) //Upper centile
					local nrrreps = r(N)
					
					centile `meanrdhat' , centile(50 `=(100-`level')/2' `=100 - (100-`level')/2')
					local medianrd = r(c_1) //Median
					local lowerprd = r(c_2) //Lower centile
					local upperprd = r(c_3) //Upper centile
					local nrdreps = r(N)
					
					centile `meanlrrhat' , centile(50 `=(100-`level')/2' `=100 - (100-`level')/2')
					local medianlrr = r(c_1) //Median
					local lowerplrr = r(c_2) //Lower centile
					local upperplrr = r(c_3) //Upper centile
					local nlrrreps = r(N)
					
					centile `meanorhat' , centile(50 `=(100-`level')/2' `=100 - (100-`level')/2')
					local medianor = r(c_1) //Median
					local lowerpor = r(c_2) //Lower centile
					local upperpor = r(c_3) //Upper centile
					local norreps = r(N)
					
					centile `meanlorhat' , centile(50 `=(100-`level')/2' `=100 - (100-`level')/2')
					local medianlor = r(c_1) //Median
					local lowerplor = r(c_2) //Lower centile
					local upperplor = r(c_3) //Upper centile
					local nlorreps = r(N)
					
					if strpos("`model'", "bayes") == 0 {
						mat `poprrouti' = (`meanmodelrr', `postserr', `medianrr', `lowerprr', `upperprr', `nrrreps')
						mat `poprdouti' = (`meanmodelrd', `postserd', `medianrd', `lowerprd', `upperprd', `nrdreps')
						mat `poplrrouti' = (`meanmodellrr', `postselrr', `medianlrr', `lowerplrr', `upperplrr', `nlrrreps')
						mat `poporouti' = (`meanmodelor', `postseor', `medianor', `lowerpor', `upperpor', `norreps')
						mat `poplorouti' = (`meanmodellor', `postselor', `medianlor', `lowerplor', `upperplor', `nlorreps')
					}
					else {
						mat `poprrouti' = (`postmeanrr', `postserr', `medianrr', `lowerprr', `upperprr', `nrrreps')
						mat `poprdouti' = (`postmeanrd', `postserd', `medianrd', `lowerprd', `upperprd', `nrdreps')
						mat `poplrrouti' = (`postmeanlrr', `postselrr', `medianlrr', `lowerplrr', `upperplrr', `nlrrreps')
						mat `poporouti' = (`postmeanor', `postseor', `medianor', `lowerpor', `upperpor', `norreps')
						mat `poplorouti' = (`postmeanlor', `postselor', `medianlor', `lowerplor', `upperplor', `nlorreps')
					}
					
					mat rownames `poprrouti' = `vari':`group'
					mat rownames `poprdouti' = `vari':`group'
					mat rownames `poplrrouti' = `vari':`group'
					mat rownames `poporouti' = `vari':`group'
					mat rownames `poplorouti' = `vari':`group'
					
					//Stack the matrices
					local ++newnrows
					if `newnrows' == 1 {
						mat `poprrout' = `poprrouti'
						mat `poprdout' = `poprdouti'
						mat `poplrrout' = `poplrrouti'
						mat `poporout' = `poporouti'
						mat `poplorout' = `poplorouti'
					}
					else {
						mat `poprrout' = `poprrout'	\  `poprrouti'
						mat `poprdout' = `poprdout'	\  `poprdouti'
						mat `poplrrout' = `poplrrout'	\  `poplrrouti'
						mat `poporout' = `poporout'	\  `poporouti'
						mat `poplorout' = `poplorout'	\  `poplorouti'
					}
				}
			}
		}
		
		if "`todo'" == "smooth" {
			
			foreach metric of local outplot {
				if strpos("`model'", "bayes") == 1 {
					if "`metric'" == "abs" {
						//# of obs
						count if `insample' == 1
						local nobs = r(N)
						
						gen `subsetid' = 1
					}
					else {
						count if `insample' == 1 & `idpair' == 2
						local nobs = r(N)
						
						egen `subsetid' = seq()	if `insample' == 1 & `idpair' == 2
					}
					
					forvalues j=1(1)`nobs' {			
						if "`metric'" == "abs" {
						
							sum `phat`j''
							local muphat = r(mean)
							
							centile `phat`j'', centile(50 `=(100-`level')/2' `=100 - (100-`level')/2')
						}
						else {
							sum `gid' if `subsetid' == `j'
							local index = r(min)
								
							sum ``metric'hat`index''
							local murhat = r(mean)
							
							centile ``metric'hat`index'', centile(50 `=(100-`level')/2' `=100 - (100-`level')/2')
						}

						local median = r(c_1) //Median
						local lowerp = r(c_2) //Lower centile
						local upperp = r(c_3) //Upper centile	
						
						if "`metric'" == "abs" {							
							if "`stat'" == "Median" {
								replace `modelp' = `median' if `insample' == 1 & `rid' == `j'
							}
							else {
								replace `modelp' = `muphat' if `insample' == 1 & `rid' == `j'
							}
							replace `modelplci' = `lowerp' if `insample' == 1 & `rid' == `j'
							replace `modelpuci' = `upperp' if `insample' == 1 & `rid' == `j'
						}
						else {
							
							if "`stat'" == "Median" {
								replace `model`metric'' = `median' if (`gid' == `index') & (`idpair' == 2) &  `insample' == 1
							}
							else {
								replace `model`metric'' = `murhat' if (`gid' == `index') & (`idpair' == 2) &  `insample' == 1
							}
							
							replace `model`metric'lci' = `lowerp' if (`gid' == `index') & (`idpair' == 2) &  `insample' == 1
							replace `model`metric'uci' = `upperp' if (`gid' == `index') & (`idpair' == 2) &  `insample' == 1
						}
					
					}
					
					drop `subsetid'
				}
				else {
					if "`metric'" == "abs" {
						
							//# of obs
							count if `insample' == 1
							local nobs = r(N)
														
							forvalues j=1(1)`nobs' {	
								
								*sum `phat`j''
								*local muphat = r(mean)
									
								centile `phat`j'', centile(50 `=(100-`level')/2' `=100 - (100-`level')/2')
								
								local median = r(c_1) //Median
								local lowerp = r(c_2) //Lower centile
								local upperp = r(c_3) //Upper centile
								
								if "`stat'" == "Median" {
									replace `modelp' = `median' if `insample' == 1 & `rid' == `j'
								}
								//This introduces sampling error
								*replace `modelplci' = `lowerp' if `insample' == 1 & `rid' == `j'  
								*replace `modelpuci' = `upperp' if `insample' == 1 & `rid' == `j'
							}
						
						//obtain the CI's -- quick way
						local critvalue -invnorm((100-`level')/200)
						replace `modelplci' = `invfn'(`eta' - `sign' `critvalue'*`modelse')`closebracket' if  `insample' == 1 //lower
						replace `modelpuci' = `invfn'(`eta' +  `sign' `critvalue'*`modelse')`closebracket' if  `insample' == 1 //upper
					}
					else {
						*sum `gid' if `insample' == 1 
						count if `insample' == 1 & `idpair' == 2
						local nstudies = r(N)
						
						egen `subsetid' = seq()	if `insample' == 1 & `idpair' == 2
					
						forvalues j=1(1)`nstudies' {							
							sum `gid' if `subsetid' == `j'
							local index = r(min)
							
							//Obtain the quantiles
							
							centile ``metric'hat`index'', centile(50 `=(100-`level')/2' `=100 - (100-`level')/2')
							
							local median = r(c_1) //Median
							local lowerp = r(c_2) //Lower centile
							local upperp = r(c_3) //Upper centile
							if "`stat'" == "Median" {
								replace `model`metric'' = `median' if (`gid' == `index') & (`idpair' == 2) &  `insample' == 1							
							}
							replace `model`metric'lci' = `lowerp' if (`gid' == `index') & (`idpair' == 2) &  `insample' == 1
							replace `model`metric'uci' = `upperp' if (`gid' == `index') & (`idpair' == 2) &  `insample' == 1
						}
						
						drop `subsetid'
					}
				}
			}	
		}
		drop if `present' != 1
		
		//drop the extra variables from the simulations
		if strpos("`model'", "bayes") == 1 {
			drop _ysim1_* _mu* _frequency _chain _index
		}	
	}
		
	//Return matrices
	if "`todo'" =="p" {
		mat colnames `popabsout' = Mean SE Median Lower Upper Sample_size
		return matrix outmatrix = `popabsout'
	}
	if "`todo'" == "r" {
		mat colnames `poprdout' = Mean SE Median Lower Upper Sample_size
		mat colnames `poprrout' = Mean SE Median Lower Upper Sample_size
		mat colnames `poplrrout' = Mean SE Median Lower Upper Sample_size
		mat colnames `poporout' = Mean SE Median Lower Upper Sample_size
		mat colnames `poplorout' = Mean SE Median Lower Upper Sample_size
		
		return matrix rdoutmatrix = `poprdout'
		return matrix rroutmatrix = `poprrout'
		return matrix lrroutmatrix = `poplrrout'
		return matrix oroutmatrix = `poporout'
		return matrix loroutmatrix = `poplorout'
	}
end

/*+++++++++++++++++++++++++	SUPPORTING FUNCTIONS: PRINTMAT +++++++++++++++++++++++++
							Print the outplot matrix micely
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
cap program drop printmat
program define printmat
	#delimit ;
	syntax, matrixout(name) type(string) [sumstat(string) dp(integer 2) power(integer 0) stratify 
		mpair mcbnetwork pcbnetwork abnetwork general comparative continuous p(integer 0) model(string) 
		nsims(string) link(string) inference(string) power(integer 0)]
	;
	#delimit cr
		local nrows = rowsof(`matrixout')
		local ncols = colsof(`matrixout')
		local rnames : rownames `matrixout'
		local eqnames : roweq `matrixout'
		local rspec "--`="&"*`=`nrows' - 1''-"
		
		local rownames = ""
		local rownamesmaxlen = 10 /*Default*/
		forvalues r = 1(1)`nrows' {
			local rname : word `r' of `rnames'
			local nlen : strlen local rname
			local rownamesmaxlen = max(`rownamesmaxlen', min(`nlen', 32)) //Check if there is a longer name
		}
		
		if "`eqnames'" != "" {
			local neqnames = wordcount("`eqnames'")
			forvalues r = 1(1)`neqnames ' {
				local eqname : word `r' of `eqnames'
				local nlen : strlen local eqname
				local rownamesmaxlen = max(`rownamesmaxlen', min(`nlen', 32)) //Check if there is a longer name
			}
		}
		
		local nlensstat : strlen local sumstat
		local nlensstat = max(10, `nlensstat')
		if "`type'" == "rde" | "`type'" == "rre" | "`type'" == "ore" {
			di as res _n "****************************************************************************************"
			if "`inference'" == "frequentist" {
				di as txt _n "Wald-type test for nonlinear hypothesis"
				if "`type'" == "rre" { 
					di as txt _n "{phang}H0: All (log)RR equal vs. H1: Some (log)RR different {p_end}"
				}
				else if "`type'" == "rde" {
					di as txt _n "{phang}H0: All RD equal vs. H1: Some RD different {p_end}"
				}
				else if "`type'" == "ore" {
					di as txt _n "{phang}H0: All (log)OR equal vs. H1: Some (log)OR different {p_end}"
				}
			}
			#delimit ;
			noi matlist `matrixout', rowtitle(Parameter) 
						cspec(& %`rownamesmaxlen's |  %8.`=`dp''f &  %8.0f &  %8.`=`dp''f o2&) 
						rspec(`rspec') underscore nodotz
			;
			#delimit cr			
		}
		if "`type'" == "popabs" | "`type'" == "poprr" | "`type'" == "popor" | "`type'" == "poprd" {
			local patho 0
			if "`type'" == "popabs" {
				local parm "Proportion"
			}
			else if "`type'" == "poprr" {
				local parm "Proportion Ratio"
			}
			else if "`type'" == "poprd" {
				local parm "Proportion Difference"
			}
			else if "`type'" == "popor" {
				local parm "Odds Ratio"
			}
		
			di as res _n "****************************************************************************************"
			di as res _n "Population-averaged estimates: `parm' "
			
			tempname mat2print
			mat `mat2print' = `matrixout'
			local nrows = rowsof(`mat2print')
			
			forvalues r = 1(1)`nrows' {
				mat `mat2print'[`r', 1] = `mat2print'[`r', 1]*10^`power'
				mat `mat2print'[`r', 3] = `mat2print'[`r', 3]*10^`power'
				mat `mat2print'[`r', 4] = `mat2print'[`r', 4]*10^`power'
				mat `mat2print'[`r', 5] = `mat2print'[`r', 5]*10^`power'
						
				forvalues c = 1(1)6 {
					local cell = `mat2print'[`r', `c'] 
					if "`cell'" == "." {
						mat `mat2print'[`r', `c'] == .z
					}
				}
				
				//Diagnose the simulation
				local cellreps = `mat2print'[`r', 6] 
				if `cellreps' < `nsims' {
					local patho 1
				}
			}

			#delimit ;
			noi matlist `mat2print', rowtitle(Parameter) 
						cspec(& %`rownamesmaxlen's |  %8.`=`dp''f &  %8.`=`dp''f & %8.`=`dp''f & %8.`=`dp''f & %8.`=`dp''f & %15.0f o2&) 
						rspec(`rspec') underscore nodotz
			;
			#delimit cr				
		}
		if ("`type'" == "exactor")  {
			local typeinf "Exact"
			
			di as res _n "****************************************************************************************"
			if ("`type'" == "exactor") {
				di as res "{pmore2} `typeinf' summary: Odds Ratio {p_end}"
			}
			di as res    "****************************************************************************************" 
			tempname mat2print

			mat `mat2print' = `matrixout'
			local nrows = rowsof(`mat2print')
			forvalues r = 1(1)`nrows' {
				mat `mat2print'[`r', 1] = `mat2print'[`r', 1]*10^`power'
				mat `mat2print'[`r', 3] = `mat2print'[`r', 3]*10^`power'
				mat `mat2print'[`r', 4] = `mat2print'[`r', 4]*10^`power'
						
				forvalues c = 1(1)4 {
					local cell = `mat2print'[`r', `c'] 
					if "`cell'" == "." {
						mat `mat2print'[`r', `c'] == .z
					}
				}
			}

			#delimit ;
			noi matlist `mat2print', rowtitle(Parameter) 
						cspec(& %`rownamesmaxlen's |  %`nlensstat'.`=`dp''f &  %9.`=`dp''f &  %9.`=`dp''f &  %9.`=`dp''f o2&) 
						rspec(`rspec') underscore  nodotz
			;
			#delimit cr
		}
		if ("`type'" == "raw") | ("`type'" == "abs") | ("`type'" == "exactabs") | ("`type'" == "rr")| ("`type'" == "or") |("`type'" == "rd")   {
			if strpos("`model'", "random") != 0 {
				local typeinf "Conditional"
			}
			else if strpos("`model'",  "fixed") != 0 | strpos("`model'", "betabin") != 0  {
				local typeinf "Marginal"
			}
			else if ("`model'" == "hexact")  {
				local typeinf "Exact"
			}
			di as res _n "****************************************************************************************"
			if ("`type'" == "raw") { 
				if "`link'" == "cloglog" {
					local expression "complementary log-log estimates"
				}
				else if "`link'" == "loglog" {
					local expression "log-log estimates"
				}
				else {
					local expression "log-odds estimates"
				}
				di as res "{pmore2} `typeinf' Summary: `expression' {p_end}"
			}
			if `power' > 0 {
				local multiplier "*10^`power'"
			}
			if strpos("`type'", "abs") != 0 { 
				di as res "{pmore2} `typeinf' summary: Proportion`multiplier' {p_end}"
			}
			if ("`type'" == "rr") {
				di as res "{pmore2} `typeinf' summary: Proportion Ratio`multiplier'  {p_end}"
			}
			if ("`type'" == "rd") {
				di as res "{pmore2} `typeinf' summary: Proportion Difference`multiplier'  {p_end}"
			}
			if ("`type'" == "or") {
				di as res "{pmore2} `typeinf' summary: Odds Ratio`multiplier'  {p_end}"
			}
			di as res    "****************************************************************************************" 
			tempname mat2print
			if "`model'" == "hexact" {
				mat `matrixout' = `matrixout'[1..., 1..6]
			}
			
			mat `mat2print' = `matrixout'
			local nrows = rowsof(`mat2print')
			forvalues r = 1(1)`nrows' {
				mat `mat2print'[`r', 1] = `mat2print'[`r', 1]*10^`power'
				mat `mat2print'[`r', 5] = `mat2print'[`r', 5]*10^`power'
				mat `mat2print'[`r', 6] = `mat2print'[`r', 6]*10^`power'
						
				forvalues c = 1(1)6 {
					local cell = `mat2print'[`r', `c'] 
					if "`cell'" == "." {
						mat `mat2print'[`r', `c'] == .z
					}
				}
			}
			
			#delimit ;
			noi matlist `mat2print', rowtitle(Parameter) 
						cspec(& %`rownamesmaxlen's |  %`nlensstat'.`=`dp''f &  %9.`=`dp''f &  %8.`=`dp''f &  %9.`=`dp''f &  %9.`=`dp''f &  %9.`=`dp''f o2&) 
						rspec(`rspec') underscore  nodotz
			;
			#delimit cr
			
			/*if "`inference'" == "bayesian" {
				#delimit ;
				noi matlist `mat2print', rowtitle(Parameter) 
							cspec(& %`rownamesmaxlen's |  %`nlensstat'.`=`dp''f &  %9.`=`dp''f &  %8.`=`dp''f & %9.`=`dp''f & %9.`=`dp''f & %9.`=`dp''f &  %9.`=`dp''f &  %9.`=`dp''f o2&) 
							rspec(`rspec') underscore  nodotz
				;
				#delimit cr
			}
			else {
				#delimit ;
				noi matlist `mat2print', rowtitle(Parameter) 
							cspec(& %`rownamesmaxlen's |  %`nlensstat'.`=`dp''f &  %9.`=`dp''f &  %8.`=`dp''f & %9.`=`dp''f & %9.`=`dp''f & %9.`=`dp''f & %9.`=`dp''f &  %9.`=`dp''f &  %9.`=`dp''f o2&) 
							rspec(`rspec') underscore  nodotz
				;
				#delimit cr
			}*/
		}
		if ("`type'" == "het") {
			di as res _n "****************************************************************************************"
			if strpos("`model'", "betabin")== 1 {
				di as txt _n "Test of heterogeneity - LR Test: beta-binomial vs binomial model"
			}
			else {
				di as txt _n "Test of heterogeneity - LR Test: RE model vs FE model"
			}
			
			tempname mat2print
			mat `mat2print' = `matrixout'
			forvalues r = 1(1)`nrows' {
				forvalues c = 1(1)`ncols' {
					local cell = `mat2print'[`r', `c'] 
					if "`cell'" == "." {
						mat `mat2print'[`r', `c'] == .z
					}
				}
			}
				
			#delimit ;
			noi matlist `mat2print', 
						cspec(& %`rownamesmaxlen's |  %8.0f `="&  %10.`=`dp''f "*`=`ncols'-1'' o2&) 
						rspec(`rspec') underscore nodotz
			;
			#delimit cr	
			
			if strpos("`model'", "betabin")!= 0 {
				di as txt "NOTE: H0: phi = 0 vs. H1: phi > 0"
			}
			else {
				di as txt "NOTE: H0: Between-study variance(s) = 0  vs. H1: Between-study variance(s) > 0"
			}
			
		}
		if ("`type'" == "mc") {
			di as res _n "****************************************************************************************"
			di as txt _n "Model comparison(s): Leave-one-out LR Test(s)"
			local rownamesmaxlen = max(`rownamesmaxlen', 17) //Check if there is a longer name
			
			tempname mat2print
			mat `mat2print' = `matrixout'
			local nrows = rowsof(`mat2print')
			local flag 0
			forvalues r = 1(1)`nrows' {					
				local pcell = `mat2print'[`r', 3] 
				if "`pcell'" == "." {
					local flag 1
				}
				if `flag' {
					continue, break
				}
			}
			
			#delimit ;
			noi matlist `matrixout', rowtitle(Omitted Parameter) 
				cspec(& %`=`rownamesmaxlen' + 2's |  %8.`=`dp''f &  %8.0f &  %8.`=`dp''f &  %15.`=`dp''f o2&) 
				rspec(`rspec') underscore nodotz
			;
			#delimit cr
			if "`interaction'" !="" {
				di as txt "*NOTE: Model with and without interaction parameter(s)"
			}
			else {
				di as txt "*NOTE: Model with and without main parameter(s)"
			}
			if `flag' {
				di as txt "*NOTE: Some p-value are missing because one or more assumptions of the LR test were violated"
			}
			if 	"`inference'" == "frequentist" {
				di as txt "*NOTE: Delta BIC = BIC (specified model) - BIC (reduced model) "
			}
			else {
				di as txt "*NOTE: Delta DIC = DIC (specified model) - DIC (reduced model) "
			}
		}
		
		if ("`continuous'" != "") {
			di as txt "NOTE: For continuous variable margins are computed at their respective mean"
		}
		if 	"`inference'" == "frequentist" {	
			if ("`type'" == "abs") | ("`type'" == "exactabs") {
				di as txt "NOTE: H0: Est = 0.5 vs. H1: Est != 0.5"
			}
			if ("`type'" == "rr") {
				di as txt "NOTE: H0: Est = 1 vs. H1: Est != 1"
			}
			if ("`type'" == "raw") {
				di as txt "NOTE: H0: Est = 0 vs. H1: Est != 0"
			}
		}
		
		if ("`type'" == "popabs") | ("`type'" == "poprr") | ("`type'" == "poprd") | ("`type'" == "popor")  {
			di as txt "NOTE: `level'% centiles obtained from `nsims' simulations of the posterior distribution"
		}		
end	
/*+++++++++++++++++++++++++	SUPPORTING FUNCTIONS: BAYESESTRCORE +++++++++++++++++++++++++
							Obtain the RR after modelling
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
cap program drop bayesestrcore
program define bayesestrcore, rclass

syntax, event(varname) [confounders(varlist) varx(varname) scimethod(string) ///
			level(integer 95)  baselevel(integer 1) link(string) model(string) MUOFF interaction ] 
		
	tempname RRoutmatrix ORoutmatrix RDoutmatrix RRoutmatrixi ORoutmatrixi RDoutmatrixi etimat hpdmat
	
	//Inverse link function
	if "`link'" == "logit" {
		local invfn "invlogit"
	}
	else if "`link'" == "log" {
		local invfn "exp"
	}
	
	local runs 0
	foreach c of local confounders {
		local EstRDexpression //RD
		local EstRRexpression //RR
		local EstORexpression //OR
		
		*qui label list `c'
		*local nlevels = r(max)
		
		qui levelsof `c', local(codelevels)
		local nlevels = r(r)
		
		local test_`c'
		
		if "`muoff'" != "" {
			if "`varx'" != "" {
				foreach l of local codelevels {
				*forvalues l = 1/`nlevels' {
					
					if "`interaction'" != "" & `l' > 1 {
						local xterm = "+ {`event':`l'.`c'#2.`varx'}"
					}
					
					if "`link'" == "logit" {
						local EstORexpression = "`EstORexpression' (`c'_`l':exp({`event':2.`varx'} `xterm'))"  //OR
					}
					else if "`link'" == "log" {
						local EstRRexpression = "`EstRRexpression' (`c'_`l':exp({`event':2.`varx'} `xterm'))"  //RR
					}
				}
			}
			else {
				*forvalues l = 2/`nlevels' {
				foreach l of local codelevels {
					if (`l' != `baselevel') {
						if "`link'" == "logit" {
							local EstORexpression = "`EstORexpression' (`c'_`l':exp({`event':`l'.`c'} `xterm'))"  //OR
						}
						else if "`link'" == "log" {
							local EstRRexpression = "`EstRRexpression' (`c'_`l':exp({`event':`l'.`c'} `xterm'))"  //RR
						}
					}
				}
			}
			
			if "`link'" == "logit" {
				local rr "norr"
			}
			else if "`link'" == "log" {
				local or "noor"
			}
		}
		else {
			if "`varx'" != "" {
				forvalues l = 1/`nlevels' {
				
					if `l' == 1 {
						*local p_`c'_`1'_`varx'_1 "`invfn'({`event':mu})"   //p for each level of c and first level of varx
						*local p_`c'_`1'_`varx'_2 "`invfn'({`event':mu} + {`event':2.`varx'})"   //p for each level of c and second level of varx
						
						
						
						local EstRDexpression = "`EstRDexpression' (`c'_`l':-`invfn'({`event':mu} + {`event':2.`varx'}) + `invfn'({`event':mu}))"
						
						if "`link'" == "logit" {
							local EstORexpression = "`EstORexpression' (`c'_`l':exp({`event':2.`varx'}))"
							local EstRRexpression = "`EstRRexpression' (`c'_`l':`invfn'({`event':mu} + {`event':2.`varx'}) / `invfn'({`event':mu}))"
						}
						else if "`link'" == "log" {
							local EstRRexpression = "`EstRRexpression' (`c'_`l':exp({`event':2.`varx'}))"
							local EstORexpression = "`EstORexpression' (`c'_`l':exp(logit(`invfn'({`event':mu} + {`event':2.`varx'})) - logit(`invfn'({`event':mu}))))"
						}
					}
					else {
					
						*local p_`c'_`1'_`varx'_1 "`invfn'({`event':mu} + {`event':`l'.`c'})"   //p for each level of c and first level of varx
						*local p_`c'_`1'_`varx'_2 "`invfn'({`event':mu} + {`event':`l'.`c'} + {`event':2.`varx'})"   //p for each level of c and second level of varx
						
						if "`interaction'" != "" {
							local xterm = "+ {`event':`l'.`c'#2.`varx'}"
						}
						
						local EstRDexpression = "`EstRDexpression' (`c'_`l':-`invfn'({`event':mu} + {`event':2.`varx'} + {`event':`l'.`c'} `xterm' ) + `invfn'({`event':mu} + {`event':`l'.`c'}))" 
						
						if "`link'" == "logit" {
							local EstORexpression = "`EstORexpression' (`c'_`l':exp({`event':2.`varx'} `xterm'))"
							local EstRRexpression = "`EstRRexpression' (`c'_`l':`invfn'({`event':mu} + {`event':2.`varx'} + {`event':`l'.`c'} `xterm' ) / `invfn'({`event':mu} + {`event':`l'.`c'}))" 
						}
						else if "`link'" == "log" {
							local EstRRexpression = "`EstRRexpression' (`c'_`l':exp({`event':2.`varx'} `xterm'))"
							local EstORexpression = "`EstORexpression' (`c'_`l':exp(logit(`invfn'({`event':mu} + {`event':2.`varx'} + {`event':`l'.`c'} `xterm' )) -logit(`invfn'({`event':mu} + {`event':`l'.`c'}))))" 
						}
					}
					
					*local EstRDexpression = "`EstRDexpression' (`c'_`l':`p_`c'_`1'_`varx'_2' - `p_`c'_`1'_`varx'_1')"
				}
			}
			else {	
				if `baselevel' == 1 {
					local basep = "`invfn'({`event':mu})"
					
					if "`link'" == "logit" {
						local baseodds = "exp({`event':mu})"
					}
					else if "`link'" == "log" {
						local baseodds = "exp(logit(exp({`event':mu})))"
					}
				}
				else {
					local basep = "`invfn'({`event':`baselevel'.`c'} + {`event':mu})"
					
					if "`link'" == "logit" {
						local baseodds = "exp({`event':`baselevel'.`c'} + {`event':mu})"
					}
					else if "`link'" == "log" {
						local baseodds = "exp(logit(exp({`event':`baselevel'.`c'} + {`event':mu})))"
					}
				}
				
				forvalues l = 1/`nlevels' {
					if `l' != `baselevel' {
						
						local EstRDexpression = "`EstRDexpression' (`c'_`l':-`invfn'({`event':`l'.`c'} + {`event':mu}) + `basep')"
						local EstRRexpression = "`EstRRexpression' (`c'_`l':`invfn'({`event':`l'.`c'} + {`event':mu})/`basep')"
						if "`link'" == "logit" {	
							local EstORexpression = "`EstORexpression' (`c'_`l':exp({`event':`l'.`c'} + {`event':mu})/`baseodds')"
						}
						else if "`link'" == "log" {
							local EstORexpression = "`EstORexpression' (`c'_`l':exp(invlogit(exp({`event':`l'.`c'} + {`event':mu})))/`baseodds')"
						}					
					}
					else {
						local EstRRexpression = "`EstRRexpression' (`c'_`l':`basep'/`basep')"
						local EstRDexpression = "`EstRDexpression' (`c'_`l':`basep' -`basep')"
						local EstORexpression = "`EstORexpression' (`c'_`l':`baseodds'/`baseodds')"
					}	
				}
			}
		}
		
		//OR
		//eti
		bayesstats summary `EstORexpression', clevel(`level') 
		mat `etimat' = r(summary)	
		
		//hpd
		bayesstats summary `EstORexpression', clevel(`level') hpd
		mat `hpdmat' = r(summary)	
		mat `hpdmat' = `hpdmat'[1..., 5..6]
		
		mat `ORoutmatrixi' = (`etimat', `hpdmat')
		
		//Nice labels
		local rnames :rownames `ORoutmatrixi'	
		local rownames = ""
		local rownamesmaxlen = 10 /*Default*/
		
		local nrows = rowsof(`ORoutmatrixi')
		forvalues r = 1(1)`nrows' {
			local rname`r':word `r' of `rnames'
			tokenize `rname`r'', parse("_")					
			local left = "`1'"
			local right = "`3'"
			if "`3'" != "" {
				local lab:label `left' `right'
				local lab = ustrregexra("`lab'", " ", "_")
				local nlen : strlen local lab
				local rownamesmaxlen = max(`rownamesmaxlen', min(`nlen', 32)) //Check if there is a longer name
				local rownames = "`rownames' `left':`lab'" 
			}
		}
		
		if "`muoff'" != "" & "`varx'" == "" { 
			local lab:label `c' 1
			local lab = ustrregexra("`lab'", " ", "_")
			local nlen : strlen local lab
			local rownamesmaxlen = max(`rownamesmaxlen', min(`nlen', 32)) //Check if there is a longer name
			local rownames = "`c':`lab' `rownames'" 
			
			mat `ORoutmatrixi' = ( 1, 0, 0, 1, 1, 1, 1, 1  \ `ORoutmatrixi') 
			mat colnames `ORoutmatrixi' = Mean SD MCSE Median eti_Lower eti_Upper hpd_Lower hpd_Upper
		} 
			
		if "`rr'`or'" == "" {
			*-------------------------RR
			//eti
			bayesstats summary `EstRRexpression', clevel(`level')  
			mat `etimat' = r(summary)	
			
			//hpd
			bayesstats summary `EstRRexpression', clevel(`level') hpd 
			mat `hpdmat' = r(summary)	
			mat `hpdmat' = `hpdmat'[1..., 5..6]
			
			mat `RRoutmatrixi' = (`etimat', `hpdmat')
			
			*---------------------RD
			//eti
			bayesstats summary `EstRDexpression', clevel(`level')  
			mat `etimat' = r(summary)	
			
			//hpd
			bayesstats summary `EstRDexpression', clevel(`level') hpd
			mat `hpdmat' = r(summary)
			mat `hpdmat' = `hpdmat'[1..., 5..6]
			
			mat `RDoutmatrixi' = (`etimat', `hpdmat')
		}
		else {
			mat `RDoutmatrixi' = J(`nlevels', 8, .)
			mat colnames `RDoutmatrixi' = Mean SD MCSE Median eti_Lower eti_Upper hpd_Lower hpd_Upper
			
			if "`rr'" != "" {
				mat `RRoutmatrixi' = J(`nlevels', 8, .)
				mat colnames `RRoutmatrixi' = Mean SD MCSE Median eti_Lower eti_Upper hpd_Lower hpd_Upper
			}
			
			if "`or'" != "" {
				mat `ORoutmatrixi' = J(`nlevels', 8, .)
				mat colnames `ORoutmatrixi' = Mean SD MCSE Median eti_Lower eti_Upper hpd_Lower hpd_Upper
			}
		}
		
		mat rownames `ORoutmatrixi' = `rownames'
		mat rownames `RRoutmatrixi' = `rownames'
		mat rownames `RDoutmatrixi' = `rownames'
		
		mat colnames `RRoutmatrixi' = Mean SD MCSE Median eti_Lower eti_Upper hpd_Lower hpd_Upper
		mat colnames `RDoutmatrixi' = Mean SD MCSE Median eti_Lower eti_Upper hpd_Lower hpd_Upper
		mat colnames `ORoutmatrixi' = Mean SD MCSE Median eti_Lower eti_Upper hpd_Lower hpd_Upper
		
		//Stack the matrices
		local ++runs
		if `runs' == 1 {
			mat `ORoutmatrix' = `ORoutmatrixi'
			mat `RRoutmatrix' = `RRoutmatrixi'
			mat `RDoutmatrix' = `RDoutmatrixi'
		}
		else {
			mat `ORoutmatrix' = `ORoutmatrix' \ `ORoutmatrixi'
			mat `RRoutmatrix' = `RRoutmatrix' \ `RRoutmatrixi'
			mat `RDoutmatrix' = `RDoutmatrix' \`RDoutmatrixi'
		}
	}
	
	return matrix rroutmatrix = `RRoutmatrix'
	return matrix rdoutmatrix = `RDoutmatrix'
	return matrix oroutmatrix = `ORoutmatrix'

end	

/*+++++++++++++++++++++++++	SUPPORTING FUNCTIONS: BAYESESTR +++++++++++++++++++++++++
							Estimate RR after modelling
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
	cap program drop bayesestr
	program define bayesestr, rclass
		syntax, event(varname)[ regexpression(string) catreg(varlist) typevarx(string) varx(varname) comparator(varname) scimethod(string) ///
			level(integer 95) mpair mcbnetwork pcbnetwork abnetwork general comparative stratify power(integer 0) by(varname) sid(varname)   ///
			baselevel(integer 1)  interaction link(string) model(string) cov(string)]
				
		if "`comparative'`mcbnetwork'`pcbnetwork'`mpair'" != "" {
			local idpairconcat "#`varx'"
		}
		
		if strpos("`regexpression'", "`sid'") != 0 | "`cov'" == "freeint" {
			local mu "muoff"
		}
		
		if "`mcbnetwork'`pcbnetwork'" != "" {
			tokenize `regexpression'
			if "`interaction'" != "" {
				tokenize `2', parse(".")
			 }
			 else {
				tokenize `3', parse(".")
			 }
			local index "`3'"
			if "`by'" ! = "`index'" {
				local catreg = "`3' `catreg'"
			}
			local stratify //nullify
		}
		
		local confounders "`catreg'"
		
		tempname RRoutmatrix RDoutmatrix ORoutmatrix row ///
				outmatrixr overallRR overallRD overallOR  bymatRR bymatRD bymatOR ///
				compmatRR compmatRD compmatOR  ///
				catregmatRD catregmatRR catregmatOR  ///
				exactlorout exactorout exactlorouti exactorouti exactrrouti ///
				coefor coeflor lorci bymatRD
				 		
		local nrowsout 0
		local nrowsnl 0
		local nby 0
		local ncomp 0
		local ncatreg 0
		
		if ("`by'" != "") & ("`typevarx'" == "i") & ("`stratify'" == "") & ("`mcbnetwork'" == "")  {		
			bayesestrcore, event(`event') varx(`varx')  confounders(`by')  scimethod(`scimethod') link(`link') `mu' baselevel(`baselevel')
			
			matrix `bymatRD' = r(rdoutmatrix)
			matrix `bymatRR' = r(rroutmatrix)
			matrix `bymatOR' = r(oroutmatrix)
			local nby = rowsof(`bymatRR')
			
			mat `RRoutmatrix' = `bymatRR'
			mat `RDoutmatrix' = `bymatRD'
			mat `ORoutmatrix' = `bymatOR'
			
			local nrowsout = rowsof(`RRoutmatrix')
		}
		
		if ("`by'" != "`comparator'") & ("`comparator'" != ""){
			*qui label list `comparator'
			*local nc = r(max)
			
			qui levelsof `comparator'
			local nc = r(r)
			
			if (`nc' > 1) {	
		
				bayesestrcore, event(`event')  varx(`varx') confounders(`comparator')  scimethod(`scimethod') link(`link') `mu' baselevel(`baselevel')
				
				matrix `compmatRR' = r(rroutmatrix)
				matrix `compmatRD' = r(rdoutmatrix)
				matrix `compmatOR' = r(oroutmatrix)
				
				local ncomp = rowsof(`compmatRR')
				
				if `nrowsout' > 0 {
					matrix `RRoutmatrix' = `RRoutmatrix' \ `compmatRR'
					matrix `RDoutmatrix' = `RDoutmatrix' \ `compmatRD'
					matrix `ORoutmatrix' = `ORoutmatrix' \ `compmatOR'
				}
				else {
					matrix `RRoutmatrix' = `compmatRR'	
					matrix `RDoutmatrix' = `compmatRD'
					matrix `ORoutmatrix' = `compmatOR'
				}
				local nrowsout = rowsof(`RRoutmatrix')
			}
		}		
			
		if "`catreg'" != "" {			
			if "`mcbnetwork'`pcbnetwork'`comparative'" != "" { 
				bayesestrcore, event(`event') varx(`varx') confounders(`catreg') baselevel(`baselevel') scimethod(`scimethod') link(`link') `interaction' `mu' 
			}
			else {
				bayesestrcore, event(`event') confounders(`catreg') baselevel(`baselevel')  scimethod(`scimethod') link(`link') `mu'
			}
			
			matrix `catregmatRR' = r(rroutmatrix)
			matrix `catregmatRD' = r(rdoutmatrix)
			matrix `catregmatOR' = r(oroutmatrix)
			
			local ncatreg = rowsof(`catregmatRR')

			if `nrowsout' > 0 {
				matrix `RRoutmatrix' = `RRoutmatrix' \ `catregmatRR'
				matrix `RDoutmatrix' = `RDoutmatrix' \ `catregmatRD'
				matrix `ORoutmatrix' = `ORoutmatrix' \ `catregmatOR'
			}
			else {
				matrix `RRoutmatrix' = `catregmatRR'
				matrix `RDoutmatrix' = `catregmatRD'
				matrix `ORoutmatrix' = `catregmatOR'
			}
			
			local nrowsout = rowsof(`RRoutmatrix')
		}
		
		//Overall when no confounders
		if ("`comparative'`mcbnetwork'`pcbnetwork'`mpair'" != "")  & "`catreg'" == "" {
						
			bayesestrcore, event(`event') confounders(`varx') baselevel(`baselevel') scimethod(`scimethod') link(`link') `mu'
			
			mat `overallRR' = r(rroutmatrix)
			mat `overallRD' = r(rdoutmatrix)
			mat `overallOR' = r(oroutmatrix)
			
			mat `overallRR' = `overallRR'[2, 1...]  //The first row is rendundant
			mat `overallRD' = `overallRD'[2, 1...]
			mat `overallOR' = `overallOR'[2, 1...]
						
			mat rownames `overallRR' = :Overall
			mat rownames `overallRD' = :Overall
			mat rownames `overallOR' = :Overall
			
			if `nrowsout' > 0 {
				matrix `RRoutmatrix' = `RRoutmatrix' \ `overallRR'
				matrix `RDoutmatrix' = `RDoutmatrix' \ `overallRD'
				matrix `ORoutmatrix' = `ORoutmatrix' \ `overallOR'
			}
			else {
				matrix `RRoutmatrix' = `overallRR'
				matrix `RDoutmatrix' = `overallRD'
				matrix `ORoutmatrix' = `overallOR'
			}
			local nrowsout = rowsof(`RRoutmatrix')
		}
		
		local inltest = "no"
		
		return local inltest = "`inltest'"
		return matrix rroutmatrix = `RRoutmatrix'
		return matrix rdoutmatrix = `RDoutmatrix'
		return matrix oroutmatrix = `ORoutmatrix'
	end	

/*+++++++++++++++++++++++++	SUPPORTING FUNCTIONS: FREQESTRCORE +++++++++++++++++++++++++
							Obtain the RR after modelling
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
cap program drop freqestrcore
program define freqestrcore, rclass

syntax, estimates(string) [marginlist(string) scimethod(string) varx(varname) by(varname) confounders(varlist) level(integer 95) ///
	baselevel(integer 1) link(string) model(string) total(varname)]
		
	tempname lRRcoef lRRV RRoutmatrix RDcoef RDV  RDoutmatrix nltestRR nltestRD lORcoef lORV ORoutmatrix nltestOR 
	

	//Expression for logodds prediction
	if "`link'" == "cloglog" {
		local expression "expression(logit(invcloglog(predict(xb))))"  //logit
		local expressionp "expression(invcloglog(predict(xb)))"  //p
	}
	else if "`link'" == "loglog" {
		if "`model'" == "crbetabin" {
			local expression "expression(-logit(exp(-exp(-(predict(xb))))))"  //logit
			local expressionp "expression(exp(-exp(-(predict(xb)))))"  //p
		}
		else {
			local expression "expression(logit(1-invcloglog(predict(xb))))"  //logit
			local expressionp "expression(1-invcloglog(predict(xb)))"  //p
		}
	}
	else if "`link'" == "log" {
		local expression "expression(logit(predict(ir)))"
	}
	else {
		local expression "expression(predict(xb))"     //logit
		local expressionp "expression(invlogit(predict(xb)))"  //p
	}
	if "`model'" == "cbbetabin" {
		local expression "expression(xb() - _b[_cons]) at(mu==1)"
	}
		
	qui {
		//Approximate sampling distribution critical value
		estimates restore `estimates'
		if "`model'" == "random" {
			local df = e(N) -  e(k_f) - e(k_r)
		}
		else {
			local df = e(N) -  e(k)
		}	
		local predcmd = e(predict)
		
		if "`predcmd'" == "mepoisson_p" {
				local expression "expression(logit(exp(predict(xb))/`total'))"
		}
					
		local crittvalue invttail(`df', `=(100-`level')/200')
		local critzvalue -invnorm((100-`level')/200)
		
		local EstRRlnexpression //log RR
		local EstORlnexpression //log OR
		local EstRDexpression // RD
		foreach c of local confounders {	
			*qui label list `c'
			*local nlevels = r(max)
			
			qui levelsof `c', local(codelevels)
			local nlevels = r(r)
					
			local test_`c'
			
			if "`varx'" != "" {
				foreach l of local codelevels {
				*forvalues l = 1/`nlevels' {
					if `l' == 1 {
						local test_`c' = "_b[`c'_`l']"
					}
					else {
						local test_`c' = "_b[`c'_`l'] = `test_`c''"
					}
					local EstRDexpression = "`EstRDexpression' (`c'_`l': - (_b[`l'.`c'#2.`varx']) + (_b[`l'.`c'#1.`varx']))"
					local EstRRlnexpression = "`EstRRlnexpression' (`c'_`l': ln(invlogit(_b[`l'.`c'#2.`varx'])) - ln(invlogit(_b[`l'.`c'#1.`varx'])))"
					local EstORlnexpression = "`EstORlnexpression' (`c'_`l': _b[`l'.`c'#2.`varx'] - _b[`l'.`c'#1.`varx'])"
				}
			}
			else {					
				local test_`c' = "_b[`c'_`baselevel']"
				local init 1
				
				foreach l of local codelevels {
				*forvalues l = 1/`nlevels' {
					if `l' != `baselevel' {
						local test_`c' = "_b[`c'_`l'] = `test_`c''"
						*local EstRDexpression = "`EstRDexpression' (`c'_`l': (invlogit(_b[`l'.`c'])) - (invlogit(_b[`baselevel'.`c'])))"
					}
					local EstRDexpression = "`EstRDexpression' (`c'_`l': - (_b[`l'.`c']) + (_b[`baselevel'.`c']))"
					local EstRRlnexpression = "`EstRRlnexpression' (`c'_`l': ln(invlogit(_b[`l'.`c'])) - ln(invlogit(_b[`baselevel'.`c'])))"
					local EstORlnexpression = "`EstORlnexpression' (`c'_`l': _b[`l'.`c'] - _b[`baselevel'.`c'])"	
				}
			}
		}
		
		//RD
		
		estimates restore `estimates'
		*margins `marginlist', `expressionrd' over(`by') post level(`level')
		margins `marginlist', `expressionp' over(`by') post level(`level') //work at p level
		nlcom `EstRDexpression', post level(`level')
		mat `RDcoef' = e(b)
		mat `RDV' = e(V)
		mat `RDV' = vecdiag(`RDV')	
		local ncols = colsof(`RDcoef') //length of the vector
		local rnames :colnames `RDcoef'
		
			
		local i = 1

		foreach c of local confounders {
			*qui label list `c'
			*local nlevels = r(max)
			
			qui levelsof `c'
			local nlevels = r(r)
			
			if (`nlevels' > 2 & "`varx'" == "") | (`nlevels' > 1 & "`varx'" != "" ){
				qui testnl (`test_`c'')
				local testnl_`c'_chi2 = r(chi2)				
				local testnl_`c'_df = r(df)
				local testnl_`c'_p = r(p)

				if `i'==1 {
					mat `nltestRD' =  [`testnl_`c'_chi2', `testnl_`c'_df', `testnl_`c'_p']
				}
				else {
					mat `nltestRD' = `nltestRD' \ [`testnl_`c'_chi2', `testnl_`c'_df', `testnl_`c'_p']
				}
				 
				local ++i
			}
		}
		
		
		//RR
		estimates restore `estimates'
		margins `marginlist', `expression' over(`by') post level(`level') //work at logit level
		nlcom `EstRRlnexpression', post level(`level')
		mat `lRRcoef' = e(b)
		mat `lRRV' = e(V)
		mat `lRRV' = vecdiag(`lRRV')	
		local ncols = colsof(`lRRcoef') //length of the vector
		local rnames :colnames `lRRcoef'

		local rowtestnl			
		local i = 1

		foreach c of local confounders {
			*qui label list `c'
			*local nlevels = r(max)
			
			qui levelsof `c'
			local nlevels = r(r)
			
			if (`nlevels' > 2 & "`varx'" == "") | (`nlevels' > 1 & "`varx'" != "" ){
				qui testnl (`test_`c'')
				local testnl_`c'_chi2 = r(chi2)				
				local testnl_`c'_df = r(df)
				local testnl_`c'_p = r(p)

				if `i'==1 {
					mat `nltestRR' =  [`testnl_`c'_chi2', `testnl_`c'_df', `testnl_`c'_p']
				}
				else {
					mat `nltestRR' = `nltestRR' \ [`testnl_`c'_chi2', `testnl_`c'_df', `testnl_`c'_p']
				}
				 
				local ++i
				local rowtestnl = "`rowtestnl' `c' "
			}
		}
		
		//OR
		estimates restore `estimates'
		margins `marginlist', `expression' over(`by') post level(`level')  //work at logit level
		nlcom `EstORlnexpression', post level(`level')
		mat `lORcoef' = e(b)
		mat `lORV' = e(V)
		mat `lORV' = vecdiag(`lORV')	
		local ncols = colsof(`lORcoef') //length of the vector
		local rnames :colnames `lORcoef'
					
		local i = 1
		foreach c of local confounders {
			*qui label list `c'
			*local nlevels = r(max)
			
			qui levelsof `c'
			local nlevels = r(r)
			
			if (`nlevels' > 2 & "`varx'" == "") | (`nlevels' > 1 & "`varx'" != "" ){
				testnl (`test_`c'')
				local testnl_`c'_chi2 = r(chi2)				
				local testnl_`c'_df = r(df)
				local testnl_`c'_p = r(p)

				if `i'==1 {
					mat `nltestOR' =  [`testnl_`c'_chi2', `testnl_`c'_df', `testnl_`c'_p']
				}
				else {
					mat `nltestOR' = `nltestOR' \ [`testnl_`c'_chi2', `testnl_`c'_df', `testnl_`c'_p']
				}
				local ++i
			}
		}
	}
	
	mat `RDoutmatrix' = J(`ncols', 9, .)
	mat `RRoutmatrix' = J(`ncols', 9, .)
	mat `ORoutmatrix' = J(`ncols', 9, .)
	
	*mat `tstats' = J(`ncols', 3, .)
	
	forvalues r = 1(1)`ncols' {
		mat `RDoutmatrix'[`r', 1] = (`RDcoef'[1,`r']) /*Estimate*/
		mat `RDoutmatrix'[`r', 2] = sqrt(`RDV'[1, `r']) /*se */
		mat `RDoutmatrix'[`r', 3] = `RDcoef'[1,`r']/sqrt(`RDV'[1, `r']) /*Z */
		
		mat `RRoutmatrix'[`r', 1] = exp(`lRRcoef'[1,`r']) /*Estimate*/
		mat `RRoutmatrix'[`r', 2] = sqrt(`lRRV'[1, `r']) /*se in log scale, power 1*/
		mat `RRoutmatrix'[`r', 3] = `lRRcoef'[1,`r']/sqrt(`lRRV'[1, `r']) /*Z in log scale*/
		
		mat `ORoutmatrix'[`r', 1] = exp(`lORcoef'[1,`r']) /*Estimate*/
		mat `ORoutmatrix'[`r', 2] = sqrt(`lORV'[1, `r']) /*se in log scale, power 1*/
		mat `ORoutmatrix'[`r', 3] = `lORcoef'[1,`r']/sqrt(`lORV'[1, `r']) /*Z in log scale*/
		
		mat `RDoutmatrix'[`r', 4] =  normprob(-abs(`RDoutmatrix'[`r', 3]))*2  /*z p-value*/
		mat `RRoutmatrix'[`r', 4] =  normprob(-abs(`RRoutmatrix'[`r', 3]))*2  /*z p-value*/
		mat `ORoutmatrix'[`r', 4] =  normprob(-abs(`ORoutmatrix'[`r', 3]))*2  /*z p-value*/
		
		//z
		mat `RDoutmatrix'[`r', 5] = (`RDcoef'[1, `r'] - `critzvalue' * sqrt(`RDV'[1, `r'])) /*lower*/
		mat `RDoutmatrix'[`r', 6] = (`RDcoef'[1, `r'] + `critzvalue' * sqrt(`RDV'[1, `r'])) /*upper*/
		
		mat `RRoutmatrix'[`r', 5] = exp(`lRRcoef'[1, `r'] - `critzvalue' * sqrt(`lRRV'[1, `r'])) /*lower*/
		mat `RRoutmatrix'[`r', 6] = exp(`lRRcoef'[1, `r'] + `critzvalue' * sqrt(`lRRV'[1, `r'])) /*upper*/
		
		mat `ORoutmatrix'[`r', 5] = exp(`lORcoef'[1, `r'] - `critzvalue' * sqrt(`lORV'[1, `r'])) /*lower*/
		mat `ORoutmatrix'[`r', 6] = exp(`lORcoef'[1, `r'] + `critzvalue' * sqrt(`lORV'[1, `r'])) /*upper*/
		
		//t
		mat `RDoutmatrix'[`r', 7] = ttail(`df', abs(`RDoutmatrix'[`r', 3]))*2   /*t p-value*/
		mat `RRoutmatrix'[`r', 7] = ttail(`df', abs(`RRoutmatrix'[`r', 3]))*2   /*t p-value*/
		mat `ORoutmatrix'[`r', 7] = ttail(`df', abs(`ORoutmatrix'[`r', 3]))*2   /*t p-value*/
			
		mat `RDoutmatrix'[`r', 8] = (`RDcoef'[1, `r'] - `crittvalue' * sqrt(`RDV'[1, `r'])) /*lower*/
		mat `RDoutmatrix'[`r', 9] = (`RDcoef'[1, `r'] + `crittvalue' * sqrt(`RDV'[1, `r'])) /*upper*/
		
		mat `RRoutmatrix'[`r', 8] = exp(`lRRcoef'[1, `r'] - `crittvalue' * sqrt(`lRRV'[1, `r'])) /*lower*/
		mat `RRoutmatrix'[`r', 9] = exp(`lRRcoef'[1, `r'] + `crittvalue' * sqrt(`lRRV'[1, `r'])) /*upper*/
		
		mat `ORoutmatrix'[`r', 8] = exp(`lORcoef'[1, `r'] - `crittvalue' * sqrt(`lORV'[1, `r'])) /*lower*/
		mat `ORoutmatrix'[`r', 9] = exp(`lORcoef'[1, `r'] + `crittvalue' * sqrt(`lORV'[1, `r'])) /*upper*/
	}
	
	local rownames = ""
	local rownamesmaxlen = 10 /*Default*/
	
	local nrows = rowsof(`RRoutmatrix')
	forvalues r = 1(1)`nrows' {
		local rname`r':word `r' of `rnames'
		tokenize `rname`r'', parse("_")					
		local left = "`1'"
		local right = "`3'"
		if "`3'" != "" {
			local lab:label `left' `right'
			local lab = ustrregexra("`lab'", " ", "_")
			local nlen : strlen local lab
			local rownamesmaxlen = max(`rownamesmaxlen', min(`nlen', 32)) //Check if there is a longer name
			local rownames = "`rownames' `left':`lab'" 
		}
	}
	mat rownames `RDoutmatrix' = `rownames'
	mat rownames `RRoutmatrix' = `rownames'
	mat rownames `ORoutmatrix' = `rownames'
	
	if `i' > 1 {
		mat rownames `nltestRD' = `rowtestnl'
		mat rownames `nltestRR' = `rowtestnl'
		mat rownames `nltestOR' = `rowtestnl'
		
		return matrix nltestRD = `nltestRD'
		return matrix nltestRR = `nltestRR'	
		return matrix nltestOR = `nltestOR'
	}
	return local i = "`i'"
	
	return matrix rdoutmatrix = `RDoutmatrix'
	return matrix rroutmatrix = `RRoutmatrix'
	return matrix oroutmatrix = `ORoutmatrix'

end	
/*+++++++++++++++++++++++++	SUPPORTING FUNCTIONS: FREQESTR +++++++++++++++++++++++++
							Estimate RR after modelling
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
	cap program drop freqestr
	program define freqestr, rclass

		syntax, estimates(string) studyid(varname) [event(varname) total(varname) catreg(varlist) typevarx(string) varx(varname) comparator(varname) scimethod(string) ///
			level(integer 95) DP(integer 2) mpair mcbnetwork pcbnetwork abnetwork general comparative aliasdesign(string) stratify power(integer 0) by(varname) ///
			regexpression(string) baselevel(integer 1)  interaction link(string) model(string) inference(string) total(varname) cov(string)]
		
		//Expression for logodds prediction		
		if "`link'" == "cloglog" {
			local expression "expression(logit(invcloglog(predict(xb))))"  //logit
			local expressionrd "expression(invcloglog(predict(xb)))"  //p
		}
		else if "`link'" == "loglog" {
			if "`model'" == "crbetabin" {
				local expression "expression(logit(exp(-exp(-(predict(xb))))))"  //logit
				local expressionrd "expression(exp(-exp(-(predict(xb)))))"  //p
			}
			else {
				*local expression "expression(-logit(invcloglog(predict(xb))))"  //logit
				local expression "expression(logit(1-invcloglog(predict(xb))))"  //logit
				local expressionrd "expression(1-invcloglog(predict(xb)))"  //p
			}
		}
		else if "`link'" == "log" {
			local expression "expression(logit(predict(ir)))"  //logit
			local expressionrd "expression(predict(ir))"  //p
		}
		else {
			local expression "expression(predict(xb))"  //logit
			local expressionrd "expression(invlogit(predict(xb)))"  //p
		}
		
		if "`model'" == "cbbetabin" {
			local expression "expression(xb() - _b[_cons]) at(mu==1)"  //logit
		}
		
		local invfn "invlogit"
		//Approximate sampling distribution critical value
		
		estimates restore `estimates'
		local df = e(N) -  e(k)
		local crittvalue invttail(`df', `=(100-`level')/200')
		local critzvalue -invnorm((100-`level')/200)
		
		if "`comparative'`mcbnetwork'`pcbnetwork'`mpair'" != "" {
			local idpairconcat "#`varx'"
		}
		
		if "`mcbnetwork'`pcbnetwork'" != "" {
			tokenize `regexpression'
			if "`interaction'" != "" {
				tokenize `2', parse(".")
			 }
			 else {
				tokenize `3', parse(".")
			 }
			local index "`3'"
			if "`by'" ! = "`index'" {
				local catreg = "`3' `catreg'"
			}
			local stratify //nullify
		}
		//CI option
		if "`scimethod'"== "t" {
			local statistic "t"
		}
		else {
			local statistic "z"
		}
		
		if "`aliasdesign'" == "comparative" {
			gettoken varx catreg:catreg
		}
		
		local confounders "`catreg'"
		local marginlist
		while "`catreg'" != "" {
			tokenize `catreg'
			local first "`1'"
			macro shift 
			local catreg `*'
			if "`first'" != "`studyid'" {
				local marginlist = `"`marginlist' `first'`idpairconcat'"'
			}
		}
		
		tempname lRRcoef lRRV RRoutmatrix lORcoef lORV ORoutmatrix row ///
				outmatrixr overallRR overallRD overallOR  nltestRR nltestRD nltestOR rowtestnl testmat2print bymatRR bymatOR ///
				bynltestRR bynltestRD bynltestOR compmatRR compmatRD compmatOR compnltestRR ///
				compnltestOR catregmatRR catregmatRD catregmatOR catregnltestRR catregnltestRD catregnltestOR varxcoef ///
				exactlorout exactorout exactlorouti exactorouti exactrrouti ///
				coefor coeflor lorci RDcoef RDV RDoutmatrix bymatRD bynltestRD ///
				compmatRD compnltestRD catregmatRD catregnltestRD
				 		
		local nrowsout 0
		local nrowsnl 0
		local nby 0
		local ncomp 0
		local ncatreg 0
		
		if "`by'" != "" & "`typevarx'" == "i" & "`stratify'" == "" /*& ("`comparator'" == "")*/{		
			freqestrcore, marginlist(`varx') varx(`varx') by(`by') confounders(`by')  estimates(`estimates') scimethod(`scimethod') link(`link') model(`model') total(`total')
			
			matrix `bymatRD' = r(rdoutmatrix)
			matrix `bymatRR' = r(rroutmatrix)
			matrix `bymatOR' = r(oroutmatrix)
			
			local nby = rowsof(`bymatRR')
			local iby = r(i)
			if `iby' > 1 {
				matrix `bynltestRD' = r(nltestRD)
				matrix `bynltestRR' = r(nltestRR)
				matrix `bynltestOR' = r(nltestOR)
				matrix `nltestRR' = `bynltestRR'
				matrix `nltestRD' = `bynltestRD'
				matrix `nltestOR' = `bynltestOR'
				local nrowsnl = rowsof(`nltestRR')
			}
			
			mat `RDoutmatrix' = `bymatRD'
			mat `RRoutmatrix' = `bymatRR'
			mat `ORoutmatrix' = `bymatOR'
			local nrowsout = rowsof(`RRoutmatrix')
		}
		/*
		//  :when does this apply?
		if ("`by'" != "`comparator'") & ("`comparator'" != ""){
			*label list `comparator'
			*local nc = r(max)
				
			qui levelsof `comparator'
			local nc = r(r)
			
			if (`nc' > 1) {	
		 
				freqestrcore, marginlist(`varx') varx(`varx') by(`comparator') confounders(`comparator') estimates(`estimates') scimethod(`scimethod') link(`link') model(`model') total(`total')
				
				matrix `compmatRD' = r(rdoutmatrix)
				matrix `compmatRR' = r(rroutmatrix)
				matrix `compmatOR' = r(oroutmatrix)
				
				local ncomp = rowsof(`compmatRR')
				local icomp = r(i)
				if `icomp' > 1 {
					matrix `compnltestRR' = r(nltestRR)
					matrix `compnltestRD' = r(nltestRD)
					matrix `compnltestOR' = r(nltestOR)
					
					if `nrowsnl' > 0 {
						matrix `nltestRR' = `nltestRR' \ `compnltestRR'
						matrix `nltestRD' = `nltestRD' \ `compnltestRD'
						matrix `nltestOR' = `nltestOR' \ `compnltestOR'
					}
					else {
						matrix `nltestRR' = `compnltestRR'
						matrix `nltestRD' = `compnltestRD'
						matrix `nltestOR' = `compnltestOR'
					}
					local nrowsnl = rowsof(`nltestRR')
				}
				
				if `nrowsout' > 0 {
					matrix `RDoutmatrix' = `RDoutmatrix' \ `compmatRD'
					matrix `RRoutmatrix' = `RRoutmatrix' \ `compmatRR'
					matrix `ORoutmatrix' = `ORoutmatrix' \ `compmatOR'
				}
				else {
					matrix `RDoutmatrix' = `compmatRD'	
					matrix `RRoutmatrix' = `compmatRR'
					matrix `ORoutmatrix' = `compmatOR'
				}
				local nrowsout = rowsof(`RRoutmatrix')
			}
		}*/	
			
			
		if "`marginlist'" != "" & "`aliasdesign'" == "" {
			
			if "`comparative'`mcbnetwork'`pcbnetwork'`mpair'" != "" { 
				freqestrcore, marginlist(`marginlist') varx(`varx') confounders(`confounders') baselevel(`baselevel') estimates(`estimates') scimethod(`scimethod') link(`link') model(`model') total(`total')
			}
			else {
				freqestrcore, marginlist(`marginlist') confounders(`confounders') baselevel(`baselevel') estimates(`estimates') scimethod(`scimethod') link(`link') model(`model') total(`total')
			}
			
			matrix `catregmatRD' = r(rdoutmatrix)
			matrix `catregmatRR' = r(rroutmatrix)
			matrix `catregmatOR' = r(oroutmatrix)
			
			local ncatreg = rowsof(`catregmatRR')
			local icatreg = r(i)
			if `icatreg' > 1 {
				matrix `catregnltestRD' = r(nltestRD)
				matrix `catregnltestRR' = r(nltestRR)
				matrix `catregnltestOR' = r(nltestOR)
				
				if `nrowsnl' > 0 {
					matrix `nltestRD' = `nltestRD' \ `catregnltestRD'
					matrix `nltestRR' = `nltestRR' \ `catregnltestRR'
					matrix `nltestOR' = `nltestOR' \ `catregnltestOR'
				}
				else {
					matrix `nltestRD' = `catregnltestRD'
					matrix `nltestRR' = `catregnltestRR'
					matrix `nltestOR' = `catregnltestOR'					
				}
				local nrowsnl = rowsof(`nltestRR')
			}
			if `nrowsout' > 0 {
				matrix `RDoutmatrix' = `RDoutmatrix' \ `catregmatRD'
				matrix `RRoutmatrix' = `RRoutmatrix' \ `catregmatRR'
				matrix `ORoutmatrix' = `ORoutmatrix' \ `catregmatOR'
			}
			else {
				matrix `RRoutmatrix' = `catregmatRR'
				matrix `RDoutmatrix' = `catregmatRD'
				matrix `ORoutmatrix' = `catregmatOR'
			}
			
			local nrowsout = rowsof(`RRoutmatrix')
		}

		if ("`comparative'`mcbnetwork'`pcbnetwork'`aliasdesign'`mpair'" != "") {			
			mat `overallRR' = J(1, 9, .)
			mat `overallRD' = J(1, 9, .)
			mat `overallOR' = J(1, 9, .)			
		
			//RR and OR
			estimates restore `estimates'
			local df = e(N) -  e(k)
			local predcmd = e(predict)
			
			if "`predcmd'" == "mepoisson_p" {
				local expression "expression(logit(exp(predict(xb))/`total'))"  //logit
				local expressionrd "expression(exp(predict(xb))/`total')"  //p
			}
			
			margins `varx', `expression' post level(`level')
			
			mat `varxcoef' = e(b)'
			local varxcats:rownames `varxcoef'
			
			forvalues r = 1(1)2 {
				local rname`r':word `r' of `varxcats'
				tokenize `rname`r'', parse(.)
				local coef`r' = "`1'"
				if strpos("`coef`r''", "bn") != 0 {
					local coef`r' = ustrregexra("`coef`r''", "bn", "")
				}
			}
					
			//log rr metric
			nlcom (Overall: ln(`invfn'(_b[`coef2'.`varx'])) - ln(`invfn'(_b[`coef1'.`varx']))) 		  
			mat `lRRcoef' = r(b)
			mat `lRRV' = r(V)
			mat `lRRV' = vecdiag(`lRRV')
			
			//log or metric
			nlcom (Overall: _b[`coef2'.`varx'] - _b[`coef1'.`varx']) 	  
			mat `lORcoef' = r(b)
			mat `lORV' = r(V)
			mat `lORV' = vecdiag(`lORV')
			
			//RD
			estimates restore `estimates'
			margins `varx', `expressionrd' post level(`level')
			nlcom (Overall: -(_b[`coef2'.`varx']) + (_b[`coef1'.`varx'])) 		  
			mat `RDcoef' = r(b)
			mat `RDV' = r(V)
			mat `RDV' = vecdiag(`RDV')
			
			mat `overallRD'[1, 1] = (`RDcoef'[1,1])  //rr
			mat `overallRD'[1, 2] = sqrt(`RDV'[1, 1]) //se
			mat `overallRD'[1, 3] = `RDcoef'[1, 1]/sqrt(`RDV'[1, 1]) //zvalue
			
			mat `overallRR'[1, 1] = exp(`lRRcoef'[1,1])  //rr
			mat `overallRR'[1, 2] = sqrt(`lRRV'[1, 1]) //se
			mat `overallRR'[1, 3] = `lRRcoef'[1, 1]/sqrt(`lRRV'[1, 1]) //zvalue
			
			mat `overallOR'[1, 1] = exp(`lORcoef'[1,1])  //or
			mat `overallOR'[1, 2] = sqrt(`lORV'[1, 1]) //se
			mat `overallOR'[1, 3] = `lORcoef'[1, 1]/sqrt(`lORV'[1, 1]) //zvalue
						
			//z			
			mat `overallRD'[1, 4] = normprob(-abs(`overallRD'[1, 3]))*2 //z pvalue
			mat `overallRR'[1, 4] = normprob(-abs(`overallRR'[1, 3]))*2 //z pvalue
			mat `overallOR'[1, 4] = normprob(-abs(`overallOR'[1, 3]))*2 //z pvalue

			mat `overallRD'[1, 5] = (`RDcoef'[1, 1] - `critzvalue'*sqrt(`RDV'[1, 1])) //ll
			mat `overallRD'[1, 6] = (`RDcoef'[1, 1] + `critzvalue'*sqrt(`RDV'[1, 1])) //ul
			
			mat `overallRR'[1, 5] = exp(`lRRcoef'[1, 1] - `critzvalue'*sqrt(`lRRV'[1, 1])) //ll
			mat `overallRR'[1, 6] = exp(`lRRcoef'[1, 1] + `critzvalue'*sqrt(`lRRV'[1, 1])) //ul
			
			mat `overallOR'[1, 5] = exp(`lORcoef'[1, 1] - `critzvalue'*sqrt(`lORV'[1, 1])) //ll
			mat `overallOR'[1, 6] = exp(`lORcoef'[1, 1] + `critzvalue'*sqrt(`lORV'[1, 1])) //ul
			
			//t
			mat `overallRR'[1, 7] = ttail(`df', abs(`overallRR'[1, 3]))*2 //t pvalue
			mat `overallRD'[1, 7] = ttail(`df', abs(`overallRD'[1, 3]))*2 //t pvalue
			mat `overallOR'[1, 7] = ttail(`df', abs(`overallOR'[1, 3]))*2 //t pvalue
				
			mat `overallRD'[1, 8] = (`RDcoef'[1, 1] - `crittvalue'*sqrt(`RDV'[1, 1])) //ll
			mat `overallRD'[1, 9] = (`RDcoef'[1, 1] + `crittvalue'*sqrt(`RDV'[1, 1])) //ul
			
			mat `overallRR'[1, 8] = exp(`lRRcoef'[1, 1] - `crittvalue'*sqrt(`lRRV'[1, 1])) //ll
			mat `overallRR'[1, 9] = exp(`lRRcoef'[1, 1] + `crittvalue'*sqrt(`lRRV'[1, 1])) //ul
			
			mat `overallOR'[1, 8] = exp(`lORcoef'[1, 1] - `crittvalue'*sqrt(`lORV'[1, 1])) //ll
			mat `overallOR'[1, 9] = exp(`lORcoef'[1, 1] + `crittvalue'*sqrt(`lORV'[1, 1])) //ul
			
			mat rownames `overallRR' = :Overall
			mat rownames `overallRD' = :Overall
			mat rownames `overallOR' = :Overall
			
			if `nrowsout' > 0 {
				matrix `RDoutmatrix' = `RDoutmatrix' \ `overallRD'
				matrix `RRoutmatrix' = `RRoutmatrix' \ `overallRR'
				matrix `ORoutmatrix' = `ORoutmatrix' \ `overallOR'
			}
			else {
				matrix `RRoutmatrix' = `overallRR'
				matrix `RDoutmatrix' = `overallRD'
				matrix `ORoutmatrix' = `overallOR'
			}
			local nrowsout = rowsof(`RRoutmatrix')
		}
		
		mat colnames `RDoutmatrix' = Mean SE `statistic' P>|z| z_Lower z_Upper P>|t| t_Lower t_Upper	
		
		mat colnames `RRoutmatrix' = Mean SE(lrr) `statistic'(lrr) P>|z| z_Lower z_Upper P>|t| t_Lower t_Upper
		mat colnames `ORoutmatrix' = Mean SE(lor) `statistic'(lor) P>|z| z_Lower z_Upper P>|t| t_Lower t_Upper
		
		if "`model'" == "hexact" {
			tempvar subset insample hold holdleft holdright
			
			gen `insample' = e(sample)
			
			//Summarize OR
			local nrows = rowsof(`ORoutmatrix') //length of the vector
			local rnames :rownames `ORoutmatrix'
			local eqnames :roweq  `ORoutmatrix'
			local newnrows 0
			
			if "`comparative'`mcbnetwork'" == "" {
				local catvars : list uniq eqnames	
				foreach vari of local catvars {
					
					cap drop `hold'	
					decode `vari', gen(`hold')
					label list `vari'
					local ngroups = r(max)
					local baselab:label `vari' `baselevel'
					
					//count in basegroup
					tempvar meanphat`baselevel' meanrrhat`baselevel' meanorhat`baselevel' meanlorhat`baselevel' gid`baselevel' sumphat`baselevel' subsetid`baselevel'
					tempname exactrrouti`baselevel' exactorouti`baselevel' exactlorouti`baselevel'
										
					cap gen `varx' = 0 if `vari' == `baselevel' & `insample' == 1
										
	
					mat `exactorouti`baselevel'' = (1, 0, 1, 1)
					mat `exactlorouti`baselevel'' = (1, 0, 1, 1)
					
					local baselab = ustrregexra("`baselab'", " ", "_")

					mat rownames `exactorouti`baselevel'' = `vari':`baselab'
					mat rownames `exactlorouti`baselevel'' = `vari':`baselab'
					
					//Other groups
					forvalues g=1(1)`ngroups' {
					
						if `g' != `baselevel' {
							tempvar meanphat`g' meanrrhat`g' meanorhat`g' meanlorhat`g' gid`g' sumphat`g' subsetid`g'
							tempname exactrrouti`g' exactorouti`g' exactlorouti`g'
							
							local glab:label `vari' `g'
							*count if `vari' == `g' & `insample' == 1
							*local ngroup`g' = r(N)	
							
							replace `varx' = 1 if `vari' == `g' & `insample' == 1
							
							cap noisily exlogistic `event' `varx' if (`vari' == `g' | `vari' == `baselevel') & `insample' == 1, binomial(`total') level(`level') `progress'
							
							if _rc == 0 {
								mat `lorci' = e(ci)	

								estat se, coef
								mat `coeflor' = r(estimates)
								
								estat se
								mat `coefor' = r(estimates)
								
								local lorlci = `lorci'[1, 1]
								if `lorlci' == . {
									local orlci = 0
								}
								else {
									local orlci = exp(`lorlci')
								}																			
								
								mat `exactorouti`g'' = (`coefor'[1, 1], `coefor'[2, 1], `orlci', exp(`lorci'[2, 1]))
								mat `exactlorouti`g'' = (`coeflor'[1, 1], `coeflor'[2, 1],  `lorlci', `lorci'[2, 1])
							}
							else {
								mat `exactorouti`g'' = J(1, 4, .)
								mat `exactlorouti`g'' = J(1, 4, .)
							}
							
							local glab = ustrregexra("`glab'", " ", "_")

							mat rownames `exactorouti`g'' = `vari':`glab'
							mat rownames `exactlorouti`g'' = `vari':`glab'
						}
						if `g' == 1 {
							mat `exactorouti' = `exactorouti`g''
							mat `exactlorouti' = `exactlorouti`g''
						}
						else {
							//Stack the matrices

							mat `exactorouti' = `exactorouti'	\  `exactorouti`g''
							mat `exactlorouti' = `exactlorouti'	\  `exactlorouti`g''
						}
					}
					//Stack the matrices
					local ++newnrows
					if `newnrows' == 1 {

						mat `exactorout' = `exactorouti'
						mat `exactlorout' = `exactlorouti'
					}
					else {

						mat `exactorout' = `exactorout'	\  `exactorouti'
						mat `exactlorout' = `exactlorout'	\  `exactlorouti'
					}
				}
			}
			
			if "`comparative'" != "" | "`mcbnetwork'" != "" | "`mpair'" != "" {
				//Comparative R
				local mindex 0
				local newnrows 0
				foreach vari of local eqnames {		
					local ++mindex
					local group : word `mindex' of `rnames'
					
					cap drop `subset'
					if "`group'" != "Overall" {
						cap drop `hold'
						decode `vari', gen(`hold')
						
						local latentgroup = ustrregexra("`group'", "_", " ")
						gen `subset' = 1 if `hold' == "`latentgroup'"  & `insample' == 1
					}
					else {
						//All
						gen `subset' = 1  if `insample' == 1
					}

					cap exlogistic `event' `varx' if `subset' == 1, binomial(`total') level(`level') `progress'
					
					if _rc == 0 {	
						mat `lorci' = e(ci)

						local lorlci = `lorci'[1, 1]
						if `lorlci' == . {
							local orlci = 0
						}
						else {
							local orlci = exp(`lorlci')
						}					

						estat se, coef
						mat `coeflor' = r(estimates)
						
						estat se
						mat `coefor' = r(estimates)
						
						mat `exactorouti' = (`coefor'[1, 1], `coefor'[2, 1], `orlci', exp(`lorci'[2, 1]))
						mat `exactlorouti' = (`coeflor'[1, 1], `coeflor'[2, 1],  `lorlci', `lorci'[2, 1])
					}
					else {
						mat `exactorouti' = J(1, 4, .)
						mat `exactlorouti' = J(1, 4, .)						
					}

					mat rownames `exactorouti' = `vari':`group'
					mat rownames `exactlorouti' = `vari':`group'
					
					//Stack the matrices
					local ++newnrows
					if `newnrows' == 1 {
						mat `exactorout' = `exactorouti'
						mat `exactlorout' = `exactlorouti'
					}
					else {
						mat `exactorout' = `exactorout'	\  `exactorouti'
						mat `exactlorout' = `exactlorout'	\  `exactlorouti'
					}
				}
			}
			
			mat colnames `exactlorout' = Mean SE Lower Upper
			mat colnames `exactorout' = Mean SE Lower Upper
			
			return matrix exactorout = `exactorout'
			return matrix exactlorout = `exactlorout'
		}
			
		if `nrowsnl' > 0 {
			local inltest = "yes"
			mat colnames `nltestRR' = chi2 df p
			mat colnames `nltestRD' = chi2 df p
			mat colnames `nltestOR' = chi2 df p
			
			return matrix nltestRD = `nltestRD'
			return matrix nltestRR = `nltestRR'
			return matrix nltestOR = `nltestOR'
		}
		else {
			local inltest = "no"
		}
		return local inltest = "`inltest'"
		return matrix rroutmatrix = `RRoutmatrix'
		return matrix rdoutmatrix = `RDoutmatrix'
		return matrix oroutmatrix = `ORoutmatrix'
	end	
	
	/*+++++++++++++++++++++++++	SUPPORTING FUNCTIONS: 	RDCI +++++++++++++++++++++++++
								CI for RD 
	++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
	cap program drop rdci
	program define rdci
	
		syntax varlist, R(name) lowerci(name) upperci(name) icimethod(string) [alpha(real 0.05)]
		
		qui {	
			tokenize `varlist'
			gen `r' = . 
			gen `lowerci' = .
			gen `upperci' = .
			
			local zstar =  -invnorm(`alpha'/2)
			local chisq = invchi2(1, `=1-`alpha'')
			
			if "`transform'" != "" {
				local transform "ln"
			}
			
			count
			forvalues i = 1/`r(N)' {
				local n1 = `1'[`i']  //1 t
				local N1 = `2'[`i'] //1 t
				local n2 = `3'[`i'] //0 c
				local N2 = `4'[`i'] //0 c
				
				forvalues t = 1/2 {
					cii proportions `N`t'' `n`t'', wilson level(`level')
					local p`t' = r(proportion)
					local lo`t' = r(lb)
					local up`t' = r(ub)
				}
				
				replace `r' = `p2' - `p1' in `i'
				
				//Newcombe hybrid CI
				replace `lowerci' = (`p2' - `p1') - sqrt((`p2' - `lo2')^2 + (`up1' - `p1')^2) in `i'
				replace `upperci' = (`p2' - `p1') + sqrt((`p1' - `lo1')^2 + (`up2' - `p2')^2) in `i'
			}
		}
	end
	
	/*+++++++++++++++++++++++++	SUPPORTING FUNCTIONS: 	RRCI +++++++++++++++++++++++++
								CI for RR
	++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
	cap program drop rrci
	program define rrci
	
		syntax varlist, R(name) lowerci(name) upperci(name) icimethod(string) [alpha(real 0.05) transform]
		
		qui {	
			tokenize `varlist'
			gen `r' = . 
			gen `lowerci' = .
			gen `upperci' = .
			
			local zstar =  -invnorm(`alpha'/2)
			local chisq = invchi2(1, `=1-`alpha'')
			
			if "`transform'" != "" {
				local transform "ln"
			}
			
			count
			forvalues i = 1/`r(N)' {
				local n1 = `1'[`i']
				local N1 = `2'[`i']
				local n2 = `3'[`i']
				local N2 = `4'[`i']
				
				//Get the intervals
				if "`icimethod'" == "koopman" {
					koopmancii `n1' `N1' `n2' `N2', alpha(`alpha')
					
					mat ci = r(ci)
					local lorr  = ci[1, 1]
					local uprr = ci[1, 2]
					
					if (`n1' == 0) &(`n2'==0) {
						local rr  = 0 in `i'
					}
					else {
						local rr = (`n1'/`N1')/(`n2'/`N2')	
					}
				}
				if ("`icimethod'" == "adlog") {
					if ((`n1' == `N1') & (`n2' == `N2')) {
						local rr = (`n1'/`N1')/(`n2'/`N2')
						local n1 = `N1' - 0.5
						local n2 = `N2' - 0.5
						local nrr = ((`n1' + 0.5)/(`N1' + 0.5))/((`n2' + 0.5)/(`N2' +  0.5))
						local varhat =(1/(`n1' + 0.5)) - (1/(`N1' + 0.5)) + (1/(`n2' + 0.5)) - (1/(`N2' + 0.5))
						local lorr = `nrr' * exp(-1 * `zstar' * sqrt(`varhat'))
						local uprr = `nrr' * exp((`zstar' * sqrt(`varhat'))
					}
					else if (`n1' == 0 & `n2' == 0) {
						local lorr  = 0
						local uprr = .
						local rr = 0
						local varhat = (1/(`n1' + 0.5)) - (1/(`N1' + 0.5)) + (1/(`n2' + 0.5)) - (1/(`N2' + 0.5))
					}
					else {
						rat = (`n1'/`N1')/(`n2'/`N2')
						local nrr  = ((`n1' + 0.5)/(`N1' + 0.5))/((`n2' + 0.5)/(`N2' + 0.5))
						local varhat = (1/(`n1' + 0.5)) - (1/(`N1' + 0.5)) + (1/(`n2' + 0.5)) - (1/(`N2' + 0.5))
						local lorr = `nrr' * exp(-1 * `zstar' * sqrt(`varhat'))
						local uprr = `nrr' * exp((`zstar' * sqrt(`varhat'))
					}
				}
				if ("`icimethod'" == "bailey") {
				
					local rr = (`n1'/`N1')/(`n2'/`N2')
					
					if (`n1' == `N1') & (`n2' == `N2'){
						local varhat = (1/(`N1' - 0.5)) - (1/(`N1')) + (1/(`N2' - 0.5)) - (1/(`N2'))
					} 
					else {
						local varhat = (1/(`n1')) - (1/(`N1')) + (1/(`n2')) - (1/(`N2'))
					}

					local phat1 = `n1'/`N1'
					local phat2 = `n2'/`N2'
					local qhat1 = 1 - `phat1'
					local qhat2 = 1 - `phat2'
					
					if (`n1' == 0 | `n2' == 0) {
						
						if `n1' == 0 {
							local xn = 0.5
						}
						else {
							local xn =`n1'
						}
										
						if `n2' == 0 { 
							local yn =  0.5
						}					
						else {
							local yn =  `n2'
						}
						
						local nrr = (`xn'/`N1')/(`yn'/`N2')
						local phat1 = `xn'/`N1'
						local phat2 = `yn'/`N2'
						local qhat1 = 1 - `phat1'
						local qhat2 = 1 - `phat2'
						if (`xn' == `N1' | `yn' == `N2') {
							if  `xn' == `N1' { 
								local xn = `N1' - 0.5
							} 
							else {
								local xn = `xn'
							 }
						
							if `yn' == `N2'{ 
								local yn = `N2' - 0.5
							}
							else {
								local yn = `yn'
							}
							local nrr = (`xn'/`N1')/(`yn'/`N2')
							local phat1 = `xn'/`N1'
							local phat2 = `yn'/`N2'
							local qhat1 = 1 - `phat1'
							local qhat2 = 1 - `phat2'
						}
					}
					if (`n1' == 0 | `n2' == 0) {
						if (`n1' == 0 & `n2' == 0) {
						  local nrr = .
						  local lorr = 0
						  local uprr = .
						}
						if (`n1' == 0 & `n2' != 0) {
						  local lorr = 0
						  local uprr = `nrr' * ((1 + `zstar' * sqrt((`qhat1'/`xn') + (`qhat2'/`yn') - (`zstar'^2 * `qhat1' * `qhat2')/(9 * `xn' * `yn'))/3)/((1 - (`zstar'^2 * `qhat2')/(9 * `yn'))))^3
						}
						if (`n2' == 0 & `n1' != 0) {
						  local uprr = .
						  local lorr =`nrr' * ((1 - `zstar' * sqrt((`qhat1'/`xn') + (`qhat2'/`yn') - (`zstar'^2 * `qhat1' * `qhat2')/(9 * `xn' * `yn'))/3)/((1 - (`zstar'^2 * `qhat2')/(9 * `yn'))))^3
						}
					}
					else if (`n1' == `N1' | `n2' == `N2') {
						
						if `n1' == `N1'{
							local xn = `N1' - 0.5
						}
						else {					
							local xn = `n1'
						}
											
						if `n2' == `N2' {
							local yn = `N2' - 0.5
						}
						else {					
							local yn = `n2'
						}
						
						local nrr = (`xn'/`N1')/(`yn'/`N2')
						local phat1 = `xn'/`N1'
						local phat2 = `yn'/`N2'
						local qhat1 = 1 - `phat1'
						local qhat2 = 1 - `phat2'
						local lorr = `nrr' * ((1 - `zstar' * sqrt((`qhat1'/`xn') + (`qhat2'/`yn') - (`zstar'^2 * `qhat1' * `qhat2')/(9 * `xn' * `yn'))/3)/((1 - (`zstar'^2 * `qhat2')/(9 * `yn'))))^3
						local uprr = `nrr' * ((1 + `zstar' * sqrt((`qhat1'/`xn') + (`qhat2'/`yn') - (`zstar'^2 * `qhat1' * `qhat2')/(9 * `xn' * `yn'))/3)/((1 - (`zstar'^2 * `qhat2')/(9 * `yn'))))^3
					}
					else {
						local lorr = `nrr' * ((1 - `zstar' * sqrt((`qhat1'/`n1') + (`qhat2'/`n2') - (`zstar'^2 * `qhat1' * `qhat2')/(9 * `n1' * `n2'))/3)/((1 - (`zstar'^2 * `qhat2')/(9 * `n2'))))^3
						local uprr = `nrr' * ((1 + `zstar' * sqrt((`qhat1'/`n1') + (`qhat2'/`n2') - (`zstar'^2 * `qhat1' * `qhat2')/(9 * `n1' * `n2'))/3)/((1 - (`zstar'^2 * `qhat2')/(9 * `n2'))))^3
					}
				}
				if ("`icimethod'" == "katz") {
					if ((`n1' == 0 & `n2' == 0) | (`n1' == 0 & `n2' != 0) | (`n1' != 0 & `n2' == 0) | (`n1' == `N1' & `n2' == `N2')) {
						if (`n1' == 0 & `n2' == 0) {
						  local lorr = 0
						  local uprr = .
						  local rr = 0
						  local varhat = .
						}
						if (`n1' == 0 & `n2' != 0) {
						  local lorr = 0
						  local rr = (`n1'/`N1')/(`n2'/`N2')
						  local n1 = 0.5
						  local nrr = (`n1'/`N1')/(`n2'/`N2')
						  local varhat = (1/`n1') - (1/`N1') + (1/`n2') - (1/`N2')
						  local uprr = `nrr' * exp(`zstar' * sqrt(`varhat'))
						}
						if (`n1' != 0 & `n2' == 0) {
						  local uprr = .
						  local rr = (`n1'/`N1')/(`n2'/`N2')
						  
						  local n2 = 0.5
						  local nrr = (`n1'/`N1')/(`n2'/`N2')
						  local varhat = (1/`n1') - (1/`N1') + (1/`n2') - (1/`N2')
						  local lorr = `nrr' * exp(-1 * `zstar' * sqrt(`varhat'))
						}
						if (`n1' == `N1' & `n2' == `N2') {
						  local rr = (`n1'/`N1')/(`n2'/`N2')
						  
						  local n1 = `N1' - 0.5
						  local n2 = `N2' - 0.5
						  local nrr = (`n1'/`N1')/(`n2'/`N2')
						  local varhat = (1/`n1') - (1/`N1') + (1/`n2') - (1/`N2')
						  local lorr = `nrr' * exp(-1 * `zstar' * sqrt(`varhat'))
						  
						  local n1 = `N1' - 0.5
						  local n2 = `N2' - 0.5
						  local nrr = (`n1'/`N1')/(`n2'/`N2')
						  local varhat = (1/`n1') - (1/`N1') + (1/`n2') - (1/`N2')
						  local uprr = `nrr' * exp(`zstar' * sqrt(`varhat'))
						}
					}
					else {
						local rr = (`n1'/`N1')/(`n2'/`N2')
						local varhat = (1/`n1') - (1/`N1') + (1/`n2') - (1/`N2')
						local lorr = `rr' * exp(-1 * `zstar' * sqrt(`varhat'))
						local uprr = `rr' * exp(`zstar' * sqrt(`varhat'))
					}
				}
				if ("`icimethod'" == "asinh") {
					if ((`n1' == 0 & `n2' == 0) | (`n1' == 0 & `n2' != 0) | (`n1' != 0 & `n2' == 0) | (`n1' == `N1' & `n2' == `N2')) {
						if (`n1' == 0 & `n2' == 0) {
						  local lorr = 0
						  local uprr = .
						  local rr = 0
						  local varhat = .
						}
						if (`n1' = 0 & `n2' != 0) {
						  local rr = (`n1'/`N1')/(`n2'/`N2')
						  local lorr = 0
						  local n1 = `zstar'
						  local nrr = (`n1'/`N1')/(`n2'/`N2')
						  local varhat = 2 * asinh((`zstar'/2) * sqrt(1/`n1' + 1/`n2' - 1/`N1' - 1/`N2'))
						  local uprr = exp(log(`nrr') + `varhat')
						}
						if (`n1' != 0 & `n2' == 0) {
						  local rr = .
						  local uprr = .
						  local n2 = `zstar'
						  local nrr = (`n1'/`N1')/(`n2'/`N2')
						  local varhat = 2 * asinh((`zstar'/2) * sqrt(1/`n1' + 1/`n2' - 1/`N1' - 1/`N2'))
						  local lorr = exp(log(`nrr') - `varhat')
						}
						if (`n1' = `N1' & `n2' == `N2') {
						  local rr = (`n1'/`N1')/(`n2'/`N2')
						  `n1' = `N1' - 0.5
						  `n2' = `N2' - 0.5
						  local nrr = (`n1'/`N1')/(`n2'/`N2')
						  local varhat = 2 * asinh((`zstar'/2) * sqrt(1/`n1' + 1/`n2' - 1/`N1' - 1/`N2'))
						  local lorr = exp(log(`nrr') - `varhat')
						  local uprr = exp(log(`nrr') + `varhat')
						}
					}
					else {
						local rr = (`n1'/`N1')/(`n2'/`N2')
						local varhat = 2 * asinh((`zstar'/2) * sqrt(1/`n1' + 1/`n2' - 1/`N1' - 1/`N2'))
						local lorr = exp(log(`rr') - `varhat')
						local uprr = exp(log(`rr') + `varhat')
					}
				}
				if ("`icimethod'" == "noether") {
					if ((`n1' == 0 & `n2' == 0) | (`n1' == 0 & `n2' != 0) | (`n1' != 0 & `n2' == 0) | (`n1' == `N1' & `n2' == `N2')) {
						if (`n1' == 0 & `n2' == 0) {
						  local lorr = 0
						  local uprr = .
						  local rr = 0
						  local sehat = .
						  local varhat = .
						}
						if (`n1' == 0 & `n2' != 0) {
						  local rr = (`n1'/`N1')/(`n2'/`N2')
						  local lorr = 0
						  local n1 = 0.5
						  local nrr = (`n1'/`N1')/(`n2'/`N2')
						  local sehat = `nrr' * sqrt((1/`n1') - (1/`N1') + (1/`n2') - (1/`N2'))
						  local uprr = `nrr' + `zstar' * `sehat'
						}
						if (`n1' != 0 & `n2' == 0) {
						  local rr = Inf
						  local uprr = Inf
						  local n2 = 0.5
						  local nrr = (`n1'/`N1')/(`n2'/`N2')
						  local  sehat = `nrr' * sqrt((1/`n1') - (1/`N1') + (1/`n2') - (1/`N2'))
						  local lorr = `nrr' - `zstar' * `sehat'
						}
						if (`n1' == `N1' & `n2' == `N2') {
						  local rr = (`n1'/`N1')/(`n2'/`N2')
						  local n1 = `N1' - 0.5
						  local n2 = `N2' - 0.5
						  local nrr = (`n1'/`N1')/(`n2'/`N2')
						  local sehat = `nrr' * sqrt((1/`n1') - (1/`N1') + (1/`n2') - (1/`N2'))
						  local uprr = `nrr' + `zstar' * `sehat'
						  local lorr = `nrr' - `zstar' * `sehat'
						}
					}
					else {
						local rr = (`n1'/`N1')/(`n2'/`N2')
						local sehat = `rr'* sqrt((1/`n1') - (1/`N1') + (1/`n2') - (1/`N2'))
						local lorr = `rr' - `zstar' * `sehat'
						local uprr = `rr' + `zstar' * `sehat'
					}
					local lorr = max(0, `lorr')
				}
				
				replace `r' = `transform'(`rr') in `i'
				replace `lowerci' = `transform'(`lorr') in `i'
				replace `upperci' = `transform'(`uprr') in `i'
			}
		}
	end

	/*+++++++++++++++++++++++++	SUPPORTING FUNCTIONS: KOOPMANCII +++++++++++++++++++++++++
								CI for RR
	++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
	cap program drop koopmancii
	program define koopmancii, rclass

		syntax anything(name=data id="data"), [alpha(real 0.05)]
		
		local len: word count `data'
		if `len' != 4 {
			di as error "Specify full data: n1 N1 n2 N2"
			exit
		}
		
		foreach num of local data {
			cap confirm integer number `num'
			if _rc != 0 {
				di as error "`num' found where integer expected"
				exit
			}
		}
		
		tokenize `data'
		cap assert ((`1' <= `2') & (`3' <= `4'))
		if _rc != 0{
			di as err "Order should be n1 N1 n2 N2"
			exit _rc
		}
		
		mata: koopman_ci((`1', `2', `3', `4'), `alpha')
		
		return matrix ci = ci
		return scalar alpha = `alpha'	

	end
/*+++++++++++++++++++++++++	SUPPORTING FUNCTIONS: 	cmlCI +++++++++++++++++++++++++
								CI for RR
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
	cap program drop cmlci
	program define cmlci

		syntax varlist, r(name) lowerci(name) upperci(name) [alpha(real 0.05) transform]
		
		if "`transform'" != "" {
			local transform "ln"
		}
		
		qui {	
			tokenize `varlist'
			gen `r' = . 
			gen `lowerci' = .
			gen `upperci' = .
			
			tempname matci
			count
			forvalues i = 1/`r(N)' {
				local a = `1'[`i']
				local b = `2'[`i']
				local c = `3'[`i']
				local d = `4'[`i']

				cmlcii `a' `b' `c' `d', alpha(`alpha')
				mat `matci' = r(ci)
				
				local n = `a' + `b' + `c' + `d'
	
				local p1 = (`a' + `b')/`n'
				local p0 = (`a' + `c')/`n'
				
				local RR = `p1'/`p0'
				
				replace `r' = `transform'(`RR') in `i'
				replace `lowerci' = `transform'(`matci'[1, 1]) in `i'
				replace `upperci' = `transform'(`matci'[1, 2]) in `i'
			}
		}
	end
	
	/*+++++++++++++++++++++++++	SUPPORTING FUNCTIONS: cmlCII +++++++++++++++++++++++++
								CI for RR
	++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
	cap program drop cmlcii
	program define cmlcii, rclass
		syntax anything(name=data id="data"), [alpha(real 0.05)]
		
		local len: word count `data'
		if `len' != 4 {
			di as error "Specify full data: a b c d"
			exit
		}
		
		foreach num of local data {
			cap confirm integer number `num'
			if _rc != 0 {
				di as error "`num' found where integer expected"
				exit
			}
		}
		
		tokenize `data'
		mata: cml_ci((`1', `2', `3', `4'), `alpha')
		
		return matrix ci = ci
		return scalar alpha = `alpha'	

	end
/*+++++++++++++++++++++++++	SUPPORTING FUNCTIONS: 	ORCCCI +++++++++++++++++++++++++
								CI for OR
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
	cap program drop orccci
	program define orccci

		syntax varlist, r(name) lowerci(name) upperci(name) [level(real 95) mpair mcbnetwork icimethod(string) transform]
		
		if "`transform'" != "" {
			local transform "ln"
		}
		
		qui {	
			tokenize `varlist'
			gen `r' = . 
			gen `lowerci' = .
			gen `upperci' = .
			
			count
			forvalues i = 1/`r(N)' {
				//matched data
				if "`mcbnetwork'`mpair'" != ""  {
					local a = `1'[`i']
					local b = `2'[`i']
					local c = `3'[`i']
					local d = `4'[`i']

					mcci `a' `b' `c' `d', l(`level') 
				}
				else {
					//unmatched data
					local n1 = `1'[`i']
					local N1 = `2'[`i']
					local n2 = `3'[`i']
					local N2 = `4'[`i']

					cci `n1' `=`N1'-`n1'' `n2' `=`N2'-`n2'', l(`level') `icimethod'
				}
				local OR =  r(or)
				local lor = r(lb_or)
				local uor = r(ub_or)
				
				replace `r' = `transform'(`OR') in `i'
				replace `lowerci' = `transform'(`lor') in `i'
				replace `upperci' = `transform'(`uor') in `i'
			}
		}
	end

	/*+++++++++++++++++++++++++	SUPPORTING FUNCTIONS: absexactci +++++++++++++++++++++++++
								CI for proportions
	++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
	cap program drop absexactci
	program define absexactci, rclass

		syntax anything(name=data id="data"), [level(real 95) icimethod(string)]
		
		tempname absexact
		mat `absexact' = J(1, 6, .)
		
		local len: word count `data'
		if `len' != 2 {
			di as error "Specify full data: N n"
			exit
		}
		
		foreach num of local data {
			cap confirm integer number `num'
			if _rc != 0 {
				di as error "`num' found where integer expected"
				exit
			}
		}
		tokenize `data'
		cap assert (`2' <= `1') 
		if _rc != 0 {
			di as err "Order should be N n"
			exit _rc
		}
		
		cii proportions `1' `2', `icimethod' level(`level')
		
		mat `absexact'[1, 1] = r(proportion) 
		mat `absexact'[1, 2] = r(se)
		mat `absexact'[1, 5] = r(lb) 
		mat `absexact'[1, 6] = r(ub)
		
		local zvalue = (`absexact'[1, 1] - 0.5)/sqrt(0.25/`1')
		mat `absexact'[1, 3] = `zvalue'
		
		local pvalue = normprob(-abs(`zvalue'))*2
		mat `absexact'[1, 4] = `pvalue'
		
		return matrix absexact = `absexact'
	end	
/*==================================== GETWIDTH  ================================================*/
/*===============================================================================================*/
capture program drop getlen
program define getlen
//From metaprop

qui{

	gen `2' = 0
	count
	local N = r(N)
	forvalues i = 1/`N'{
		local this = `1'[`i']
		local width: _length "`this'"
		replace `2' =  `width' in `i'
	}
} 
end

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
METABAYESOPTCHECK : Options for bayesian inference
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
capture program drop metabayesoptscheck
program define metabayesoptscheck, rclass
	#delimit ;
	syntax [,
		nchains(integer 3)  /*3*/
		thinning(integer 5) /*5*/
		burnin(integer 5000) /*5000*/
		mcmcsize(integer 2500) /*2500*/
		rseed(integer 1)
		refsampling(integer 1)
		varprior(passthru)
		feprior(passthru)
		*
	]
	;
	#delimit cr
	
	local bayesopts = `"nchains(`nchains') thinning(`thinning') burnin(`burnin') mcmcsize(`mcmcsize') rseed(`rseed') `feprior' `varprior'  `options'"'
	
	return local modelopts = `"`bayesopts'"'
	return local mcmcsize = "`=`mcmcsize'*`nchains''"
	return local refsampling = "`refsampling'"
end



/*	SUPPORTING FUNCTIONS: 	METAPLOTCHECK ++++++++++++++++++++++++++++++++++++++++++
			Advance housekeeping for the metapplot
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
	capture program drop metapplotcheck
	program define metapplotcheck, rclass
	#delimit ;
	syntax  [,
		/*Passed from top*/
		SUMMARYonly
		
		/*passed via foptions/coptions*/
		noFPlot
		CATPplot
		noOVerall 
		noOVLine 
		noSTats 
		noBox
		DOUBLE 
		AStext(integer 50) 
		CIOpts(passthru) 
		DIAMopts(passthru) 
		OLineopts(passthru) 
		POINTopts(passthru) 
		BOXopts(passthru) 
		PREDciOpts(passthru)
		PREDIction  //prediction
		SORtby(varlist) //varlist
		LCols(varlist) 
		RCols(varlist) 		
		SUBLine
		TEXts(real 1.0) 
		XLAbel(passthru)
		XLIne(passthru)	/*silent option*/	
		XTick(passthru)  
		graphsave(passthru)
		logscale	
		grid
		*
	  ];
	#delimit cr
	
	if `astext' > 90 | `astext' < 10 {
		di as error "Percentage of graph as text (ASTEXT) must be within 10-90%"
		di as error "Must have some space for text and graph"
		exit
	}
	if `texts' < 0 {
		di as res "Warning: Negative text size (TEXTSize) are ignored"
		local texts 1
	}	
	
	if "`summaryonly'" != ""  {
		local box "nobox"
	}
	
	foreach var of local rcols {
		cap confirm var `var'
		if _rc!=0  {
			di in re "Variable `var' not in the dataset"
			exit _rc
		}
	}
	foreach var of local lcols {
		cap confirm var `var'
		if _rc!=0  {
			di in re "Variable `var' not in the dataset"
			exit _rc
		}
	}
	if "`lcols'" =="" {
		local lcols " "
	}
	if "`rcols'" =="" {
		local rcols " "
	}
	if "`astext'" != "" {
		local astext "astext(`astext')"
	}
	if "`texts'" != "" {
		local texts "texts(`texts')"
	}
	local plotopts `"`overall' `ovline' `stats' `box' `double' `astext' `ciopts' `diamopts' `olineopts' `pointopts' `boxopts' `predciopts' `prediction' `subline' `texts' `xlabel' `xline' `xtick' `graphsave' `logscale' `grid' `options'"'
	return local lcols ="`lcols'"
	return local rcols ="`rcols'"
	return local plotopts = `"`plotopts'"'
end

/*	SUPPORTING FUNCTIONS: 	METAPLOT ++++++++++++++++++++++++++++++++++++++++++++++++
			The forest/catterpillar plot
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
// Some re-used code from metaprop, metadta

	capture program drop metapplot
	program define metapplot

	#delimit ;
	syntax varlist [if] [in] [,
		STudyid(varname)
		POWer(integer 0)
		DP(integer 2) 
		Level(integer 95)
		Groupvar(varname)
		design(string)
		aliasdesign(string)
		smooth
		model(string)
		varxlabs(string)
		type(string)
		noWT
		
		/*Passed through*/
		logscale
		AStext(integer 50)
		ARRowopt(string) 		
		CIOpts(string) 
		DIAMopts(string) 
		DOUble 
 		LCols(varlist)
		RCols(varlist) 		
		noOVerall 
		noOVLine 		
		noSTATS
		noBox
		OLineopts(string) 
		outplot(string) 
		SUMStat(string asis) 
		POINTopts(string) 
		BOXopts(string) 
		predciopts(string)
		SUBLine 
		TEXts(real 1.0) 
		XLAbel(string) 
		XLIne(string) 
		XTick(string)
		GRID
		GRAphsave(string asis)
		prediction

		*
	  ];
	#delimit cr
	
	preserve
	
	local plotopts `"`options'"'
	
	if strpos(`"`plotopts'"', "graphregion") == 0 {
			local plotopts `"graphregion(color(white)) `plotopts'"'
	}
	
	tempvar es modeles lci modellci uci modeluci lpi upi ilci iuci predid use label tlabel id newid rid gid df expand expanded order orig flag ///
	
	tokenize "`varlist'", parse(" ")
	
	qui {
		gen `es'		=`1'*(10^`power')
		gen `lci'   	=`2'*(10^`power')
		gen `uci'   	=`3'*(10^`power')
		gen byte `use'	=`4'
		gen str `label'	=`5'
		gen `df' 		= `6'
		gen `id' 		= `7'
		
		if "`smooth'" !="" {
			gen `modeles' 	= `8'*(10^`power')
			gen `modellci' 	= `9'*(10^`power')
			gen `modeluci' 	= `10'*(10^`power')
		}
		
		//Toggle for comparative
		if ("`design'" == "comparative" | "`aliasdesign'" == "comparative" | "`design'" == "mpair") & "`outplot'" == "abs" & "`groupvar'" != "" {
			local compabs "compabs"
		}		
		
		gen  `newid' = `id'
		
		//Add five spaces on top of the dataset and 1 space below
		*qui summ `id'
		gen `expand' = 1
		replace `expand' = 1 + 5*(_n==1)  + 1*(_n==_N) 
		expand `expand'
		sort `newid' `use'

		replace `newid' = _n in 1/6
		replace `newid' = `newid' + 5 if _n>6
		replace `label' = "" in 1/5
		replace `use' = -2 in 1/4
		replace `use' = 3 in 5
		replace `newid' = _N  if _N==_n
		replace `use' = 3  if _N==_n
		replace `label' = "" if _N==_n
		
		gen `flag' = 1
		replace `flag' = 0 in 1/4
						
		//studylables
		if "`abnetwork'" != "" & "`outplot'" != "abs" {
			local studylb: variable label `groupvar'
			if "`studylb'" == "" {
				label var `label' "`groupvar'"
			}
			else {
				label var `label' "`studylb'"
			}
		}
		else {
			local studylb: variable label `studyid'
			if "`studylb'" == "" {
				label var `label' "`studyid'"
			}
			else {
				label var `label' "`studylb'"
			}
		}
		
		if "`lcols'" == "" {
			local lcols = "`label'"
		}
		else {
			local lcols "`label' `lcols'"
		}
		
		if "`compabs'" != "" { 
			replace `id' = `newid'
			drop `newid'
		}
		else{		
			egen `rid' = group(`newid')
			replace `id' = `rid'
			drop `rid' `newid'
		}
	
		tempvar estText index predText predLabel wtText modelestText
		
		gen str `estText' = string(`es', "%10.`=`dp''f") + " (" + string(`lci', "%10.`=`dp''f") + ", " + string(`uci', "%10.`=`dp''f") + ")"  if (`use' == 1 | `use' == 2 | `use' == 5)
		
		if "`smooth'" !="" {
			gen str `modelestText' = string(`modeles', "%10.`=`dp''f") + " (" + string(`modellci', "%10.`=`dp''f") + ", " + string(`modeluci', "%10.`=`dp''f") + ")"  if (`use' == 1 )
			replace `modelestText' = string(`es', "%10.`=`dp''f") + " (" + string(`lci', "%10.`=`dp''f") + ", " + string(`uci', "%10.`=`dp''f") + ")"  if (`use' == 2 | `use' == 5)
			replace `estText' = " " if (`use' == 2 | `use' == 5)
		}
		if "`wt'" == "" {
			gen str `wtText' = string(_WT, "%10.`=`dp''f") if (`use' == 1 | `use' == 2 | `use' == 5) & _WT !=.
		}
		
		if "`prediction'" != "" {
			tempvar lenestext

			replace `estText' =  " (" + string(`lci', "%10.`=`dp''f") + ", " + string(`uci', "%10.`=`dp''f") + ")"  if (`use' == 4)
			qui gen `lenestext' = length(`estText')
			qui summ `lenestext' if `use' == 1
			local lentext = r(max)
			qui summ `lenestext' if `use' == 4
			local lenic = r(min)
			local lenwhite = `lentext' - `lenic' 
			
			replace `estText' = ".  `=`lenwhite'*" "'" + `estText'  if (`use' == 4)
			
		}
		// GET MIN AND MAX DISPLAY
		// SORT OUT TICKS- CODE PINCHED FROM MIKE AND FIRandomED. TURNS OUT I'VE BEEN USING SIMILAR NAMES...
		// AS SUGGESTED BY JS JUST ACCEPT ANYTHING AS TICKS AND RESPONSIBILITY IS TO USER!
	
		if "`logscale'" != "" {
			tempvar tes tlci tuci 
			gen `tes' = ln(`es') if `es' != . & `es' > 0
			gen `tlci' = ln(`lci') if `lci' != . & `lci' > 0
			gen `tuci' = ln(`uci') if `uci' != . & `uci' > 0
			
			replace `tes' = ln(0.00001) if `es' != . & `es' == 0
			replace `tlci' = ln(0.00001) if `lci' != . & `lci' == 0
			replace `tuci' = ln(0.00001) if `uci' != . & `uci' == 0

			replace `es' = `tes'
            replace `lci' = `tlci'
			replace `uci' = `tuci'
				
			if "`smooth'" !="" { 
				tempvar tmodeles tmodellci tmodeluci
				gen `tmodeles'  = ln(`modeles') if `modeles' != . & `modeles' > 0
				gen `tmodellci' = ln(`modellci') if `modellci' != . & `modellci' > 0
				gen `tmodeluci' = ln(`modeluci') if `modeluci' != . & `modeluci' > 0
				
				replace `tmodeles'  = ln(0.00001) if `modeles' != . & `modeles' == 0
				replace `tmodellci' = ln(0.00001) if `modellci' != . & `modellci' == 0
				replace `tmodeluci' = ln(0.00001) if `modeluci' != . & `modeluci' == 0
				 
				replace `modeles' = `tmodeles'
				replace `modellci' = `tmodellci'
				replace `modeluci' = `tmodeluci'
			}
		}
		qui summ `lci', detail
		local DXmin = r(min)
		
		qui summ `uci', detail
		local DXmax = r(max)
		
		if "`DXmax'" == "." & "`DXmin'" != "" {
			if `DXmin' < 0 {
				local DXmax = - `DXmin'
			}
		}  
		
		/*
		if "`DXmax'" != "" & "`DXmin'" == "" {
			if `DXmax' > 1 {
				local DXmin = 0
			}
			else {
			}
		}
		*/
		
		if "`xlabel'" != "" {
			if "`logscale'" != "" {
				*local DXmin = ln(max(min(`xlabel'), 0.00001))
				local DXmin = ln(max(`xlabel'))
				local DXmax = ln(max(`xlabel'))
			}
			else{
				local DXmin = min(`xlabel')
				local DXmax = max(`xlabel')
			}
		}
		if "`xlabel'"=="" {
			*local xlabel "`DXmin', `DXmax'"
			
			if `DXmin' < 0 {
				local xlabel "`DXmin', 0, `DXmax'"
			}
			else {
				if `DXmax' > 1 {
					local xlabel "0, 1, `DXmax'"
				}
				else {
					local xlabel "0, `DXmax'"
				}
			}
		}

		local lblcmd ""
		tokenize "`xlabel'", parse(",")
		while "`1'" != ""{
			if "`1'" != ","{
				local lbl = string(`1',"%7.3g")
				if "`logscale'" != "" {
					*local val = ln(max(`1', 0.00001))
					
					if "`1'" == "0" {
						local val = ln(`=10^(-`dp')')
					}
					else {
						local val = ln(`1')
					}
				}
				else {
					local val = `1'
				}

				local lblcmd `lblcmd' `val' "`lbl'"
			}
			mac shift
		}
		
		if "`xtick'" == ""{
			local xtick = "`xlabel'"
		}

		local xtick2 = ""
		tokenize "`xtick'", parse(",")
		while "`1'" != ""{
			if "`1'" != ","{
				if "`logscale'" != "" {
					*local val = ln(max(`1', 0.00001))
					
					if "`1'" == "0" {
						local val = ln(`=10^(-`dp')')
					}
					else {
						local val = ln(`1')
					}
				}
				else {
					local val = `1'
				}
				local xtick2 = "`xtick2' " + string(`val')
			}
			if "`1'" == ","{
				local xtick2 = "`xtick2'`1'"
			}
			mac shift
		}
		local xtick = "`xtick2'"
		
		local DXmin = (min(`xtick',`DXmin'))
		local DXmax = (max(`xtick',`DXmax'))
		
		*local DXmin= (min(`xlabel',`xtick',`DXmin'))
		*local DXmax= (max(`xlabel',`xtick',`DXmax'))

		local DXwidth = `DXmax'-`DXmin'
	} // END QUI

	/*===============================================================================================*/
	/*==================================== COLUMNS   ================================================*/
	/*===============================================================================================*/
	qui {	// KEEP QUIET UNTIL AFTER DIAMONDS
			
		// DOUBLE LINE OPTION
		if "`double'" != "" & ("`lcols'" != "" | "`rcols'" != ""){
			replace `expand' = 1
			replace `expand' = 2 if `use' == 1
			expand `expand'
			sort `id' `use'
			bys `id' : gen `index' = _n
			sort  `id' `use' `index'
			egen `newid' = group(`id' `index')
			replace `id' = `newid'
			drop `newid'
			
			replace `use' = 1 if `index' == 2
			replace `es' = . if `index' == 2
			replace `lci' = . if `index' == 2
			replace `uci' = . if `index' == 2
			replace `estText' = "" if `index' == 2			

			foreach var of varlist `lcols' `rcols' {
			   cap confirm string var `var'
			   if _rc == 0 {				
					tempvar length words tosplit splitwhere best
					gen `splitwhere' = 0
					gen `best' = .
					gen `length' = length(`var')
					summ `length', det
					gen `words' = wordcount(`var')
					gen `tosplit' = 1 if `length' > r(max)/2+1 & `words' >= 2
					summ `words', det
					local max = r(max)
					forvalues i = 1/`max'{
						replace `splitwhere' = strpos(`var', word(`var',`i')) ///
						 if abs( strpos(`var',word(`var',`i')) - length(`var')/2 ) < `best' ///
						 & `tosplit' == 1
						replace `best' = abs(strpos(`var',word(`var',`i')) - length(`var')/2) ///
						 if abs(strpos(`var',word(`var',`i')) - length(`var')/2) < `best' 
					}

					replace `var' = substr(`var',1,(`splitwhere'-1)) if (`tosplit' == 1) & (`index' == 1)
					replace `var' = substr(`var',`splitwhere',length(`var')) if (`tosplit' == 1) & (`index' == 2)
					replace `var' = "" if (`tosplit' != 1) & (`index' == 2) & (`use' == 1)
					drop `length' `words' `tosplit' `splitwhere' `best'
			   }
			   if _rc != 0{
				replace `var' = . if (`index' == 2) & (`use' == 1)
			   }
			}
		}
				
		local maxline = 1

		if "`lcols'" != "" {
			tokenize "`lcols'"
			local lcolsN = 0

			while "`1'" != "" {
				cap confirm var `1'
				if _rc!=0  {
					di in re "Variable `1' not defined"
					exit _rc
				}
				local lcolsN = `lcolsN' + 1
				tempvar left`lcolsN' leftLB`lcolsN' leftWD`lcolsN'
				cap confirm string var `1'
				if _rc == 0{
					gen str `leftLB`lcolsN'' = `1'
				}
				if _rc != 0{
					cap decode `1', gen(`leftLB`lcolsN'')
					if _rc != 0{
						local f: format `1'
						gen str `leftLB`lcolsN'' = string(`1', "`f'")
						replace `leftLB`lcolsN'' = "" if `leftLB`lcolsN'' == "."
					}
				}
				replace `leftLB`lcolsN'' = "" if (`use' != 1) & (`lcolsN' != 1)
				local colName: variable label `1'
				if "`colName'"==""{
					local colName = "`1'"
				}

				// WORK OUT IF TITLE IS BIGGER THAN THE VARIABLE
				// SPREAD OVER UP TO FOUR LINES IF NECESSARY
				local titleln = length("`colName'")
				tempvar tmpln
				gen `tmpln' = length(`leftLB`lcolsN'')
				qui summ `tmpln' if `use' == 1
				local otherln = r(max)
				drop `tmpln'
				// NOW HAVE LENGTH OF TITLE AND MAX LENGTH OF VARIABLE
				local spread = int(`titleln'/`otherln') + 1
				if `spread' > 4{
					local spread = 4
				}
				local line = 1
				local end = 0
				gettoken now remain : colName

				while `end' == 0 {
					replace `leftLB`lcolsN'' =  `leftLB`lcolsN'' + " " + "`now'" in `line' 
					
					gettoken now remain : remain
					if ("`now'" == "") | (`line' == 4) {
						local end = 1
					}
					if length("`remain'") > `titleln'/`spread' {
						if `end' == 0 {
							local line = `line' + 1
						}
					}
				}
				if `line' > `maxline' {
					local maxline = `line'
				}
				mac shift
			}
		}
		if "`wt'" == "" {
			local rcols = "`wtText' " + "`rcols'"
			label var `wtText' "% Weight"
		}
		if "`smooth'" !=""  {
			local rcols = "`modelestText' " + "`rcols'"
			label var `modelestText' "Smooth Est (`level'% CI)"
		}
		if "`stats'" == "" {
			local rcols = "`estText' " + "`rcols'"
			label var `estText' "`sumstat' (`level'% CI)"
		}

		tempvar extra
		gen `extra' = " "
		label var `extra' " "
		local rcols = "`rcols' `extra'"

		local rcolsN = 0
		if "`rcols'" != "" {
			tokenize "`rcols'"
			local rcolsN = 0
			
			while "`1'" != ""{
				cap confirm var `1'
				if _rc!=0  {
					di in re "Variable `1' not defined"
					exit _rc
				}
				local rcolsN = `rcolsN' + 1
				tempvar right`rcolsN' rightLB`rcolsN' rightWD`rcolsN'
				cap confirm string var `1'
				if _rc == 0{
					gen str `rightLB`rcolsN'' = `1'
				}
				if _rc != 0{
					local f: format `1'
					gen str `rightLB`rcolsN'' = string(`1', "`f'")
					replace `rightLB`rcolsN'' = "" if `rightLB`rcolsN'' == "."
				}
				
				replace `rightLB`rcolsN'' = "" if (`use' == 3 |  `use' == -2 )
				
				local colName: variable label `1'
				if "`colName'"==""{
					local colName = "`1'"
				}

				// WORK OUT IF TITLE IS BIGGER THAN THE VARIABLE
				// SPREAD OVER UP TO FOUR LINES IF NECESSARY
				local titleln = length("`colName'")
				tempvar tmpln
				gen `tmpln' = length(`rightLB`rcolsN'')
				qui summ `tmpln' if `use' == 1
				local otherln = r(max)
				drop `tmpln'
				
				// NOW HAVE LENGTH OF TITLE AND MAX LENGTH OF VARIABLE
				local spread = int(`titleln'/`otherln')+1
				if `spread' > 4{
					local spread = 4
				}

				local line = 1
				local end = 0

				gettoken now remain : colName
				while `end' == 0 {
					replace `rightLB`rcolsN'' = `rightLB`rcolsN'' + " " + "`now'" in `line'
					gettoken now remain : remain

					if ("`now'" == "") | (`line' == 4) {
						local end = 1
					}
					if  length("`remain'") > `titleln'/`spread' {
						if `end' == 0 {
							local line = `line' + 1
						}
					}
				}
				if `line' > `maxline'{
					local maxline = `line'
				}
				mac shift
			}
		}

		// now get rid of extra title rows if they weren't used
		if `maxline'==3 {
			drop in 4	
		}
		if `maxline'==2 {
			drop in 3/4 
		}
		if `maxline'==1 {
			drop in 2/4 
		}
				
		egen `newid' = group(`id')
		replace `id' = `newid'
		drop `newid'
		
		
		if "`compabs'" != "" {
			summ `id' if `use'==1 & `groupvar'==1

			local startg1 = r(min)
			local stopg1 = r(max)
			drop if (`use' == 0 | `use' == -2) & (_n > `startg1' & _n!=_N)

			egen `newid' = group(`id')
			replace `id' = `newid'
			drop `newid'

			bys `groupvar' : egen `gid' =seq() if `use'==1 | `use'==2
			replace `gid' = `gid' + `startg1' - 1
			replace `id' = `gid' if `gid'!=.
			sort `id' `groupvar'
			
			cap drop `expand'
			cap gen `expand' = 1 + 2*(`id'[_n]==`id'[_n-1])

			expand `expand', gen(`expanded')
			gsort `id' `expanded' `groupvar'

			replace `id' = _n
			
			replace `use' = 3 if `expanded' //blanks
			
			sort `expanded' `gid' `id'
			
			by `expanded' `gid': egen `order' = seq() if !`expanded' & `gid' != .
			replace `id' = `id' + 0.75 if `order' == 2
			
			sort `id'
			replace `id' = `id' + 0.75 if _n == 2
		
		}
		
		if "`compabs'" != "" {
			local borderline = `maxline' + 1.5
		}
		else {
			local borderline = `maxline' + 0.75
		}
		 
		local leftWDtot = 0
		local rightWDtot = 0
		local leftWDtotNoTi = 0

		forvalues i = 1/`lcolsN'{
			getlen `leftLB`i'' `leftWD`i''
			qui summ `leftWD`i'' if `use' != 5 	// DON'T INCLUDE OVERALL STATS AT THIS POINT
			local maxL = r(max)
			local leftWDtotNoTi = `leftWDtotNoTi' + `maxL'
			replace `leftWD`i'' = `maxL'
		}
		tempvar titleLN				// CHECK IF OVERALL LENGTH BIGGER THAN REST OF LCOLS
		getlen `leftLB1' `titleLN'	
		qui summ `titleLN' if `use' == 5
		local leftWDtot = max(`leftWDtotNoTi', r(max))

		forvalues i = 1/`rcolsN'{
			getlen `rightLB`i'' `rightWD`i''
			qui summ `rightWD`i'' if  `use' != 5
			
			replace `rightWD`i'' = r(max)
			local rightWDtot = `rightWDtot' + r(max)
		}
		// CHECK IF NOT WIDE ENOUGH (I.E., OVERALL INFO TOO WIDE)
		// LOOK FOR EDGE OF DIAMOND summ `lci' if `use' == ...

		tempvar maxLeft
		getlen `leftLB1' `maxLeft'
		qui count if `use' == 2 | `use' == 5 
		if r(N) > 0 {
			summ `maxLeft' if `use' == 2 | `use' == 5 	// NOT TITLES THOUGH!
			local max = r(max)
			if `max' > `leftWDtotNoTi'{
				// WORK OUT HOW FAR INTO PLOT CAN EXTEND
				// WIDTH OF LEFT COLUMNS AS FRACTION OF WHOLE GRAPH
				local x = `leftWDtot'*(`astext'/100)/(`leftWDtot'+`rightWDtot')
				tempvar y
				// SPACE TO LEFT OF DIAMOND WITHIN PLOT (FRAC OF GRAPH)
				gen `y' = ((100-`astext')/100)*(`lci'-`DXmin') / (`DXmax'-`DXmin') 
				qui summ `y' if `use' == 2 | `use' == 5
				local extend = 1*(r(min)+`x')/`x'
				local leftWDtot = max(`leftWDtot'/`extend',`leftWDtotNoTi') // TRIM TO KEEP ON SAFE SIDE
													// ALSO MAKE SURE NOT LESS THAN BEFORE!
			}
		}
		local LEFT_WD = `leftWDtot'
		local RIGHT_WD = `rightWDtot'
		
		local textWD = (`DXwidth'*(`astext'/(100-`astext'))) /(`leftWDtot' + `rightWDtot')
		
		local AXmin = `DXmin' - 0.05*(`DXmax' - `DXmin') - `leftWDtot'*`textWD'
		forvalues i = 1/`lcolsN'{
			gen `left`i'' = `DXmin' - 0.05*(`DXmax' - `DXmin') - `leftWDtot'*`textWD'
			local leftWDtot = `leftWDtot'-`leftWD`i''
		}

		gen `right1' = `DXmax' + 0.05*(`DXmax' - `DXmin')
		forvalues i = 2/`rcolsN'{
			local r2 = `i' - 1
			gen `right`i'' = `right`r2'' + `rightWD`r2''*`textWD'
		}

		*local AXmin = `left1'
		local AXmax = `DXmax' + `rightWDtot'*`textWD'

		// DIAMONDS 
		tempvar DIAMleftX DIAMrightX DIAMbottomX DIAMtopX DIAMleftY1 DIAMrightY1 DIAMleftY2 DIAMrightY2 DIAMbottomY DIAMtopY
		
		//Complete diamond
		gen `DIAMleftX'   = `lci' if `use' == 2 | `use' == 5 
		gen `DIAMleftY1'  = `id' if (`use' == 2 | `use' == 5) 
		gen `DIAMleftY2'  = `id' if (`use' == 2 | `use' == 5) 
		
		gen `DIAMrightX'  = `uci' if (`use' == 2 | `use' == 5)
		gen `DIAMrightY1' = `id' if (`use' == 2 | `use' == 5)
		gen `DIAMrightY2' = `id' if (`use' == 2 | `use' == 5)
		
		gen `DIAMbottomY' = `id' - 0.4 if (`use' == 2 | `use' == 5)
		gen `DIAMtopY' 	  = `id' + 0.4 if (`use' == 2 | `use' == 5)
		gen `DIAMtopX'    = `es' if (`use' == 2 | `use' == 5)
		
		//Incomplete diamonds
		replace `DIAMleftX' = `DXmin' if (`lci' < `DXmin' ) & (`use' == 2 | `use' == 5) //cut the left side to the left limit
		replace `DIAMleftX' = . if (`es' < `DXmin' ) & (`use' == 2 | `use' == 5) //miss it if outside limit
		replace `DIAMleftX' = . if (`df' < 2) & (`use' == 2 | `use' == 5)  //miss it if one study
		
		replace `DIAMleftY1' = `id' + 0.4*(abs((`DXmin' -`lci')/(`es'-`lci'))) if (`lci' < `DXmin' ) & (`use' == 2 | `use' == 5) 
		replace `DIAMleftY1' = . if (`es' < `DXmin' ) & (`use' == 2 | `use' == 5) //miss it if outside limit
	
		replace `DIAMleftY2' = `id' - 0.4*( abs((`DXmin' -`lci')/(`es'-`lci')) ) if (`lci' < `DXmin' ) & (`use' == 2 | `use' == 5) 
		replace `DIAMleftY2' = . if (`es' < `DXmin' ) & (`use' == 2 | `use' == 5) 
		
		//Cutting the right side 
		replace `DIAMrightX' = `DXmax' if (`uci' > `DXmax' ) & (`use' == 2 | `use' == 5) 
		replace `DIAMrightX' = . if (`es' > `DXmax' ) & (`use' == 2 | `use' == 5) 
		
		//If one study, no diamond
		replace `DIAMrightX' = . if (`df' == 1) & (`use' == 2 | `use' == 5) 
	
		replace `DIAMrightY1' = `id' + 0.4*( abs((`uci'-`DXmax' )/(`uci'-`es')) ) if (`uci' > `DXmax' ) & (`use' == 2 | `use' == 5) 
		replace `DIAMrightY1' = . if (`es' > `DXmax' ) & (`use' == 2 | `use' == 5) 

		replace `DIAMrightY2' = `id' - 0.4*( abs((`uci'-`DXmax' )/(`uci'-`es')) ) if (`uci' > `DXmax' ) & (`use' == 2 | `use' == 5) 
		replace `DIAMrightY2' = . if (`es' > `DXmax' ) & (`use' == 2 | `use' == 5) 
			
		
		replace `DIAMbottomY' = `id' - 0.4*( abs((`uci'-`DXmin' )/(`uci'-`es')) ) if (`es' < `DXmin' ) & (`use' == 2 | `use' == 5) & (abs((`uci'-`DXmin' )/(`uci'-`es')) < 1)
		replace `DIAMbottomY' = `id' - 0.4*( abs((`DXmax' -`lci')/(`es'-`lci')) ) if (`es' > `DXmax' ) & (`use' == 2 | `use' == 5) & abs((`DXmax' -`lci')/(`es'-`lci')) < 1

		replace `DIAMtopY' = `id' + 0.4*( abs((`uci'-`DXmin' )/(`uci'-`es')) ) if (`es' < `DXmin' ) & (`use' == 2 | `use' == 5) & (abs((`uci'-`DXmin' )/(`uci'-`es')) < 1)
		replace `DIAMtopY' = `id' + 0.4*( abs((`DXmax' -`lci')/(`es'-`lci')) ) if (`es' > `DXmax' ) & (`use' == 2 | `use' == 5) & (abs((`DXmax' -`lci')/(`es'-`lci')) < 1)
		
		
		replace `DIAMtopX' = `DXmin'  if (`es' < `DXmin' ) & (`use' == 2 | `use' == 5) 
		replace `DIAMtopX' = `DXmax'  if (`es' > `DXmax' ) & (`use' == 2 | `use' == 5) 
		replace `DIAMtopX' = . if ((`uci' < `DXmin' ) | (`lci' > `DXmax' )) & (`use' == 2 | `use' == 5) //miss it if outside limit
		
		gen `DIAMbottomX' = `DIAMtopX'
	} // END QUI
	
	if "`compabs'" != "" {
		tempvar strgroupvar 
		decode `groupvar', gen(`strgroupvar')
	}
	forvalues i = 1/`lcolsN'{
		if "`compabs'" != "" {
			qui replace `leftLB`i'' = "" if (`expanded')
			if `i'==1 {
				qui replace `leftLB`i'' = "" if (`order' == 2 & `use'==1) | (`use' == -2 & `id' == `=`startg1'-1')
				*qui replace `leftLB`i'' = "Mean: " + "`groupvar'" + " = " + `strgroupvar' if `use'==2
				qui replace `leftLB`i'' = "Summary: " + "`groupvar'" + " = " + `strgroupvar' if `use'==2
			}
		}
		local lcolCommands`i' "(scatter `id' `left`i'', msymbol(none) mlabel(`leftLB`i'') mlabcolor(black) mlabpos(3) mlabsize(`texts'))"
	}

	forvalues i = 1/`rcolsN' {
		if "`compabs'" != "" {
			qui replace `rightLB`i'' = "" if (`expanded') 
		}
		local rcolCommands`i' "(scatter `id' `right`i'', msymbol(none) mlabel(`rightLB`i'') mlabcolor(black) mlabpos(3) mlabsize(`texts'))"
	}
	
	if `"`diamopts'"' == "" {
		local diamopts "lcolor(red)"
	}
	else {
		if strpos(`"`diamopts'"',"hor") != 0 | strpos(`"`diamopts'"',"vert") != 0 {
			di as error "Options horizontal/vertical not allowed in diamopts()"
			exit
		}
		if strpos(`"`diamopts'"',"con") != 0{
			di as error "Option connect() not allowed in diamopts()"
			exit
		}
		if strpos(`"`diamopts'"',"lp") != 0{
			di as error "Option lpattern() not allowed in diamopts()"
			exit
		}
		local diamopts `"`diamopts'"'
	}
	
	if "`compabs'" != "" {
		/*
		local diamopts0 "lcolor("0 0 0")"
		local diamopts1 "lcolor("255 127 0")"
		*/
		
		local diamopts0 "lcolor("40 40 85")"
		local diamopts1 "lcolor("163 100 249")"
	}
	
	//Box options
	if "`box'" == "" {
		local iw = "[aw = _WT]"
		if `"`boxopts'"' != "" & strpos(`"`boxopts'"',"msy") == 0{
			local boxopts = `"`boxopts' msymbol(square)"' 
		}
		if `"`boxopts'"' != "" & strpos(`"`boxopts'"',"msi") == 0{
			local boxopts = `"`boxopts' msize(0.5)"' 
		}
		if `"`boxopts'"' != "" & strpos(`"`boxopts'"',"mco") == 0{
			local boxopts = `"`boxopts' mcolor("180 180 180")"' 
		}
		if `"`boxopts'"' == "" {
			local boxopts "msymbol(square) msize(.5) mcolor("180 180 180")"
		}
		else{
			local boxopts `"`boxopts'"'
		}
	}
	if ("`box'" != "") {
		local boxopts "msymbol(none)"
		local iw
	}
	
	//Point options
	if "`smooth'" != "" {
		local pointsymbol "msymbol(Oh)"
		local pointcolor "mcolor("128 128 128")"
		local pointsize "msize(small)"
	}
	else {
		local pointsymbol "msymbol(O)"
		local pointcolor "mcolor("0 0 0")"
		local pointsize "msize(vsmall)"
	}
	
	if `"`pointopts'"' != "" & strpos(`"`pointopts'"',"msy") == 0 {
		local pointopts = `"`pointopts' `pointsymbol'"' 
	}
	if `"`pointopts'"' != "" & strpos(`"`pointopts'"',"ms") == 0 {
		local pointopts = `"`pointopts' `pointsize'"' 
	}
	if `"`pointopts'"' != "" & strpos(`"`pointopts'"',"mc") == 0 {
		local pointopts = `"`pointopts' `pointcolor'"' 
	}
	if `"`pointopts'"' == "" {
		local pointopts `"`pointsymbol' `pointsize' `pointcolor'"'
	}
	else{
		local pointopts `"`pointopts'"'
	}
	
	//Smooth Point options

	if `"`smoothpointopts'"' == "" {
		local smoothpointopts "msymbol(D) msize(vsmall) mcolor("0 0 0")"
		if "`compabs'" != "" {
			/*
			local smoothpointopts0 "msymbol(D) msize(vsmall) mcolor("0 0 0")"
			local smoothpointopts1 "msymbol(D) msize(vsmall) mcolor("255 127 0")"
			*/
			local smoothpointopts0 "msymbol(D) msize(vsmall) mcolor("40 40 85")"
			local smoothpointopts1 "msymbol(D) msize(vsmall) mcolor("163 100 249")"
		}
	}
	else {
		local smoothpointopts `"`smoothpointopts'"'
	}
	
	// CI options
	if `"`ciopts'"' == "" {
		
		if "`compabs'" != "" {
				/*
				local ciopts0 = `"lcolor("0 0 0")"' 
				local ciopts1 = `" lcolor("255 127 0")"'
				*/
				local ciopts0 = `"lcolor("40 40 85")"' 
				local ciopts1 = `" lcolor("163 100 249")"'
		}
		
		if "`smooth'" != "" {
			local ciopts = `"`ciopts' lwidth(1.25) lcolor(gs13)"'
		}
		else {
			local ciopts = `"`ciopts' lcolor("0 0 0")"' 	
		}
	}
	else {
		if strpos(`"`ciopts'"',"hor") != 0 | strpos(`"`ciopts'"',"vert") != 0{
			di as error "Options horizontal/vertical not allowed in ciopts()"
			exit
		}
		if strpos(`"`ciopts'"',"con") != 0{
			di as error "Option connect() not allowed in ciopts()"
			exit
		}
		if strpos(`"`ciopts'"',"lp") != 0 {
			di as error "Option lpattern() not allowed in ciopts()"
			exit
		}
		/*
		if "`smooth'" != "" {
			if strpos(`"`ciopts'"',"lc") == 0 {
				if "`compabs'" != "" {
					local ciopts0 = `"`ciopts' lcolor("0 0 0")"' 
					local ciopts1 = `"`ciopts' lcolor("255 127 0")"'
				}
				else {
					local ciopts = `"`ciopts' lcolor(red)"' 
				}
			}
			if strpos(`"`ciopts'"',"lw") == 0 {
				local ciopts = `"`ciopts' lwidth(2)"' 
			}
		}
		*/		
		local ciopts `"`ciopts'"'
	}
	//Smooth ci
	if `"`smoothciopts'"' == "" {
		local smoothciopts "lcolor("0 0 0")"
		
		if "`compabs'" != "" {
			*local smoothciopts0 "lcolor("0 0 0")"
			*local smoothciopts1 "lcolor("255 127 0")"
			*local smoothciopts0 "lcolor("9 44 72")"
			local smoothciopts0 "lcolor("40 40 85")"
			local smoothciopts1 "lcolor("163 100 249")"
			*local smoothciopts1 "lcolor("40 40 85")"
			*local smoothciopts1 "lcolor("255 127 0")"
		}
	}
	else {
		if strpos(`"`smoothciopts'"',"hor") != 0 | strpos(`"`smoothciopts'"',"vert") != 0{
			di as error "Options horizontal/vertical not allowed in ciopts()"
			exit
		}
		if strpos(`"`smoothciopts'"',"con") != 0{
			di as error "Option connect() not allowed in ciopts()"
			exit
		}
		if strpos(`"`smoothciopts'"',"lp") != 0 {
			di as error "Option lpattern() not allowed in ciopts()"
			exit
		}
		if strpos(`"`smoothciopts'"',"lw") == 0 {
				local smoothciopts = `"`smoothciopts' lwidth(.5)"' 
		}

		local smoothciopts `"`smoothciopts'"'
	}
	
	/*
	// PREDCI options
	if `"`predciopts'"' == "" {
		local predciopts "lcolor(red) lpattern(solid)"
	}
	else {
		if strpos(`"`predciopts'"',"hor") != 0 | strpos(`"`predciopts'"',"vert") != 0{
			di as error "Options horizontal/vertical not allowed in predciopts()"
			exit
		}
		if strpos(`"`predciopts'"',"con") != 0{
			di as error "Option connect() not allowed in predciopts()"
			exit
		}
		if `"`predciopts'"' != "" & strpos(`"`predciopts'"',"lp") == 0 {
			local predciopts = `"`predciopts' lpattern(solid)"' 
		}
		if `"`predciopts'"' != "" & strpos(`"`predciopts'"',"lc") == 0{
			local predciopts = `"`predciopts' lcolor(red)"' 
		}
		local predciopts `"`predciopts'"'
	}
	*/
	// Arrow options
	if `"`arrowopts'"' == "" {
		if "`smooth'" != "" {
			local arrowopts "mcolor(red) lstyle(none)"
			
			if "`compabs'" != "" {
				local arrowopts0 "mcolor("0 0 0") lstyle(none)"
				local arrowopts1 "mcolor("255 127 0") lstyle(none)"
			}
		}
		else {
			local arrowopts "mcolor("0 0 0") lstyle(none)"
		}
	}
	else {
		local forbidden "connect horizontal vertical lpattern lwidth lcolor lsytle"
		foreach option of local forbidden {
			if strpos(`"`arrowopts'"',"`option'")  != 0 {
				di as error "Option `option'() not allowed in arrowopts()"
				exit
			}
		}
		if `"`arrowopts'"' != "" & strpos(`"`arrowopts'"',"mc") == 0{
			local arrowopts = `"`arrowopts' mcolor("0 0 0")"' 
		}
		local arrowopts `"`arrowopts' lstyle(none)"'
	}

	// END GRAPH OPTS

	tempvar tempOv overrallLine ovMin ovMax h0Line
	
	if `"`olineopts'"' == "" {
		local olineopts "lwidth(thin) lcolor(red) lpattern(shortdash)"
	}
	qui summ `id'
	local DYmin = r(min)
	local DYmax = r(max) + 2
	
	qui summ `es' if `use' == 5 
	local overall = r(max)
	if `overall' > `DXmax' | `overall' < `DXmin' | "`ovline'" != "" {	// ditch if not on graph
		local overallCommand ""
	}
	else {
		local overallCommand `" (pci `=`DYmax'-2' `overall' `borderline' `overall', `olineopts') "'
	
	}
	if "`ovline'" != "" {
		local overallCommand ""
	}
	if "`subline'" != "" & "`groupvar'" != "" {
		local sublineCommand ""		
		*qui label list `groupvar'
		*local nlevels = r(max)
		
		qui levelsof `groupvar', local(codelevels)
		local nlevels = r(r)
		
		foreach l of local codelevels {			
		*forvalues l = 1/`nlevels' {
			qui summ `es' if `use' == 2  & `groupvar' == `l' 
			local tempSub`l' = r(mean)
			qui summ `id' if `use' == 1 & `groupvar' == `l'
			local subMax`l' = r(max) + 1
			local subMin`l' = r(min) - 2
			qui count if `use' == 1 & `groupvar' == `l' 
			if r(N) > 1 {
				local sublineCommand `" `sublineCommand' (pci `subMin`l'' `tempSub`l'' `subMax`l'' `tempSub`l'', `olineopts')"'
			}
		}
	}
	else {
		local sublineCommand ""
	}

	if `"`xline'"' != "" {
		tokenize "`xline'", parse(",")
		if "`logscale'" != "" {
			if "`1'" == "0" {
				local xlineval = ln(`=10^(-`dp')')
			}
			else {
				local xlineval = ln(`1')
			}
		}
		else {
			local xlineval = `1'
		}
		if "`3'" == "" {
			local xlineopts = "`3'"
		}
		else {
			local xlineopts = "lcolor(black)"
		}
		local xlineCommand `" (pci `=`DYmax'-2' `xlineval' `borderline' `xlineval', `xlineopts') "'
	}

	qui {
		//Generate indicator on direction of the off-scale arro
		tempvar rightarrow leftarrow biarrow noarrow rightlimit leftlimit offRhiY offRhiX offRloY offRloX offLloY offLloX offLhiY offLhiX
		gen `rightarrow' = 0
		gen `leftarrow' = 0
		gen `biarrow' = 0
		gen `noarrow' = 0
		
		replace `rightarrow' = 1 if ///
			(round(`uci', 0.001) > round(`DXmax' , 0.001)) & ///
			(round(`lci', 0.001) >= round(`DXmin' , 0.001))  & ///
			(`use' == 1 | `use' == 4) & (`uci' != .) & (`lci' != .)
			
		replace `leftarrow' = 1 if ///
			(round(`lci', 0.001) < round(`DXmin' , 0.001)) & ///
			(round(`uci', 0.001) <= round(`DXmax' , 0.001)) & ///
			(`use' == 1 | `use' == 4) & (`uci' != .) & (`lci' != .)
		
		replace `biarrow' = 1 if ///
			(round(`lci', 0.001) < round(`DXmin' , 0.001)) & ///
			(round(`uci', 0.001) > round(`DXmax' , 0.001)) & ///
			(`use' == 1 | `use' == 4) & (`uci' != .) & (`lci' != .)
			
		replace `noarrow' = 1 if ///
			(`leftarrow' != 1) & (`rightarrow' != 1) & (`biarrow' != 1) & ///
			(`use' == 1 | `use' == 4) & (`uci' != .) & (`lci' != .)	

		replace `lci' = `DXmin'  if (round(`lci', 0.001) < round(`DXmin' , 0.001)) & (`use' == 1 | `use' == 4) 
		replace `uci' = `DXmax'  if (round(`uci', 0.001) > round(`DXmax' , 0.001)) & (`uci' !=.) & (`use' == 1 | `use' == 4) 
		
		//Upper bound not estimable
		replace `rightarrow' = 1 if  (`lci' !=.) & (`uci' ==.) & (`use' == 1 | `use' == 4) 
		replace `uci' = `DXmax'  if  (`lci' !=.) & (`uci' ==.) & (`use' == 1 | `use' == 4) 
		
		replace `lci' = . if (round(`uci', 0.001) < round(`DXmin' , 0.001)) & (`uci' !=. ) & (`use' == 1 | `use' == 4) 
		replace `uci' = . if (round(`lci', 0.001) > round(`DXmax' , 0.001)) & (`lci' !=. ) & (`use' == 1 | `use' == 4)
		replace `es' = . if (round(`es', 0.001) < round(`DXmin' , 0.001)) & (`use' == 1 | `use' == 4) 
		replace `es' = . if (round(`es', 0.001) > round(`DXmax' , 0.001)) & (`use' == 1 | `use' == 4) 

		if "`smooth'" != "" {		
			replace `modellci' = `DXmin'  if (round(`modellci', 0.001) < round(`DXmin' , 0.001)) & (`use' == 1) 
			replace `modeluci' = `DXmax'  if (round(`modeluci', 0.001) > round(`DXmax' , 0.001)) & (`modeluci' !=.) & (`use' == 1 ) 
			
			replace `modellci' = . if (round(`modeluci', 0.001) < round(`DXmin' , 0.001)) & (`modeluci' !=. ) & (`use' == 1 ) 
			replace `modeluci' = . if (round(`modellci', 0.001) > round(`DXmax' , 0.001)) & (`modellci' !=. ) & (`use' == 1 )
			replace `modeles' = . if (round(`modeles', 0.001) < round(`DXmin' , 0.001)) & (`use' == 1 ) 
			replace `modeles' = . if (round(`modeles', 0.001) > round(`DXmax' , 0.001)) & (`use' == 1 )
		}
		
		summ `id'
		local xaxislineposition = r(max)

		local xaxis "(pci `xaxislineposition' `DXmin' `xaxislineposition' `DXmax', lwidth(thin) lcolor(black))"
		/*Xaxis 1 title */
		local xaxistitlex `=(`DXmax' + `DXmin')*0.5'
		
		/*if "`design'" == "comparative" & "`outplot'" != "abs"  {
			local varx :word 1 of `varxlabs'
			local indexlab :word 2 of `varxlabs'
			local baselab :word 3 of `varxlabs'
			
			local xaxistitle  (scatteri `=`xaxislineposition' + 2.25' `xaxistitlex' "`sumstat' (p_i (`varx' = `indexlab') / p_i(`varx' = `baselab'))", msymbol(i) mlabcolor(black) mlabpos(0) mlabsize(`texts'))
		}
		else {*/
			local xaxistitle  (scatteri `=`xaxislineposition' + 2.25' `xaxistitlex' "`sumstat'", msymbol(i) mlabcolor(black) mlabpos(0) mlabsize(`texts'))
		*}
		/*xticks*/
		local ticksx
		tokenize "`xtick'", parse(",")	
		while "`1'" != "" {
			if "`1'" != "," {
				local ticksx "`ticksx' (pci `xaxislineposition'  `1' 	`=`xaxislineposition'+.25' 	`1' , lwidth(thin) lcolor(black)) "
			}
			macro shift 
		}
		/*labels*/
		local xaxislabels
		tokenize `lblcmd'
		while "`1'" != ""{			
			local xaxislabels "`xaxislabels' (scatteri `=`xaxislineposition'+1' `1' "`2'", msymbol(i) mlabcolor(black) mlabpos(0) mlabsize(`texts'))"
			macro shift 2
		}
		if "`grid'" != "" {
			tempvar gridy gridxmax gridxmin
			
			gen `gridy' = `id' + 0.5
			gen `gridxmax' = `AXmax'
			gen `gridxmin' = `left1'
			local betweengrids "(pcspike `gridy' `gridxmin' `gridy' `gridxmax'  if `use' == 1 , lwidth(vvthin) lcolor(gs12))"	
		}
		/*
		//prediction
		if "`prediction'" != "" {
			gen `ilci' = `lci'[_n-1]
			gen `iuci' = `uci'[_n-1]
			replace `ilci' = . if `use' != 4 
			replace `iuci' = . if `use' != 4
			gen `predid' = `id'[_n-1]
			*replace `id' = `predid' if `use' == 4
			
			local cipred0 "(pcspike `predid' `lci' `predid' `ilci' if `use' == 4 , `predciopts') (pcspike `predid' `uci' `predid' `iuci' if `use' == 4 , `predciopts')"
			local cipred1 "(pcarrow `predid' `ilci' `predid' `lci' if `leftarrow' == 1  & `use' == 4, `arrowopts')	(pcarrow `predid' `iuci' `predid' `uci' if `rightarrow' == 1 & `use' == 4, `arrowopts')"
		}	
		*/		
	}	// end qui	
	/*===============================================================================================*/
	/*====================================  GRAPH    ================================================*/
	/*===============================================================================================*/
	//Smooth stats
	if "`smooth'" != "" {
		local xboxcenter "`modeles'"
		local smoothcommands1 "(pcspike `id' `modellci' `id' `modeluci' if `use' == 1 , `smoothciopts')"
		local smoothcommands2 "(scatter `id' `modeles' if `use' == 1 , `smoothpointopts')"
		
		if "`compabs'" != "" {
			local smoothcommands10 "(pcspike `id' `modellci' `id' `modeluci' if `use' == 1 & `groupvar'==1  , `smoothciopts0')"
			local smoothcommands11 "(pcspike `id' `modellci' `id' `modeluci' if `use' == 1 & `groupvar'==2 , `smoothciopts1')"
			
			local smoothcommands20 "(scatter `id' `modeles' if `use' == 1 & `groupvar'==1, `smoothpointopts0')"
			local smoothcommands21 "(scatter `id' `modeles' if `use' == 1 & `groupvar'==2, `smoothpointopts1')"
		}
	}
	else {
		local xboxcenter "`es'"
	}
	
	//Observed CI
	local cicommand "(pcspike `id' `lci' `id' `uci' if `use' == 1 , `ciopts')"
	if "`compabs'" != "" {
			local cicommand1 "(pcspike `id' `lci' `id' `uci' if `use' == 1 & `groupvar'==1, `ciopts0')"
			local cicommand2 "(pcspike `id' `lci' `id' `uci' if `use' == 1 & `groupvar'==2, `ciopts1')"
	} 
	
	//Diamonds
	if "`compabs'" != "" {
		local diamondcommand11 "(pcspike `DIAMleftY1' `DIAMleftX' `DIAMtopY' `DIAMtopX' if (`use' == 2 & `groupvar'==1) , `diamopts0')"
		local diamondcommand12 " (pcspike `DIAMtopY' `DIAMtopX' `DIAMrightY1' `DIAMrightX' if (`use' == 2 & `groupvar'==1) , `diamopts0')"
		local diamondcommand13 " (pcspike `DIAMrightY2' `DIAMrightX' `DIAMbottomY' `DIAMbottomX' if (`use' == 2 & `groupvar'==1) , `diamopts0')"
		local diamondcommand14 " (pcspike `DIAMbottomY' `DIAMbottomX' `DIAMleftY2' `DIAMleftX' if (`use' == 2 & `groupvar'==1) , `diamopts0')" 
		
		local diamondcommand21 " (pcspike `DIAMleftY1' `DIAMleftX' `DIAMtopY' `DIAMtopX' if (`use' == 2 & `groupvar'==2) , `diamopts1')"
		local diamondcommand22 " (pcspike `DIAMtopY' `DIAMtopX' `DIAMrightY1' `DIAMrightX' if (`use' == 2 & `groupvar'==2) , `diamopts1')"
		local diamondcommand23 " (pcspike `DIAMrightY2' `DIAMrightX' `DIAMbottomY' `DIAMbottomX' if (`use' == 2 & `groupvar'==2) , `diamopts1')"
		local diamondcommand24 " (pcspike `DIAMbottomY' `DIAMbottomX' `DIAMleftY2' `DIAMleftX' if (`use' == 2 & `groupvar'==2) , `diamopts1')"
		
		local diamondcommand1 "(pcspike `DIAMleftY1' `DIAMleftX' `DIAMtopY' `DIAMtopX' if (`use' == 5) , `diamopts')"
		local diamondcommand2 "(pcspike `DIAMtopY' `DIAMtopX' `DIAMrightY1' `DIAMrightX' if (`use' == 5) , `diamopts')"
		local diamondcommand3 "(pcspike `DIAMrightY2' `DIAMrightX' `DIAMbottomY' `DIAMbottomX' if (`use' == 5) , `diamopts')"
		local diamondcommand4 "(pcspike `DIAMbottomY' `DIAMbottomX' `DIAMleftY2' `DIAMleftX' if (`use' == 5) , `diamopts')"

	}
	else {
		local diamondcommand1 "(pcspike `DIAMleftY1' `DIAMleftX' `DIAMtopY' `DIAMtopX' if ( `use' == 2 |`use' == 5) , `diamopts')"
		local diamondcommand2 "(pcspike `DIAMtopY' `DIAMtopX' `DIAMrightY1' `DIAMrightX' if (`use' == 2 |`use' == 5) , `diamopts')"
		local diamondcommand3 "(pcspike `DIAMrightY2' `DIAMrightX' `DIAMbottomY' `DIAMbottomX' if (`use' == 2 |`use' == 5) , `diamopts')"
		local diamondcommand4 "(pcspike `DIAMbottomY' `DIAMbottomX' `DIAMleftY2' `DIAMleftX' if (`use' == 2 |`use' == 5) , `diamopts')"
	}
	
	//legend
	if "`compabs'" != "" {
		local lab0:label `groupvar' 1
		local lab1:label `groupvar' 2
		local legendon "legend(order(1 2) lab(1 "`groupvar' = `lab0'") lab(2 "`groupvar' = `lab1'") bmargin(zero) size(`texts') cols(2) ring(2) position(5))"
	}
	else {
		local legendoff "legend(off)"
	}
	
	//Give name if none
	if "`type'" == "fplot" {
		local fullname "forest"
	}
	else {
		local fullname "caterpillar"
	}
	
	if strpos(`"`plotopts'"',"name") == 0 {
		local plotname = "name(`outplot'`type', replace)"
	}
	if "$by_index_" != "" {
		local plotname = "name(`outplot'`type'" + "$by_index_" + ", replace)"
		noi di as res _n  "NOTE: `fullname' plot name -> `outplot'`type'$by_index_"
	}

	#delimit ;
	twoway
	 // Draw diamond to make it to construct the legend
		`diamondcommand11' `diamondcommand21'
		
		`xlineCommand' `xaxis' `xaxistitle' 
		`ticksx' `xaxislabels' 
		
	 /*COLUMN VARIABLES */	 
		`lcolCommands1' `lcolCommands2' `lcolCommands3' `lcolCommands4'  `lcolCommands5'  `lcolCommands6'
		`lcolCommands7' `lcolCommands8' `lcolCommands9' `lcolCommands10' `lcolCommands11' `lcolCommands12'
		`rcolCommands1' `rcolCommands2' `rcolCommands3' `rcolCommands4'  `rcolCommands5'  `rcolCommands6' 
		`rcolCommands7' `rcolCommands8' `rcolCommands9' `rcolCommands10' `rcolCommands11' `rcolCommands12' 
		
	 /*PLOT BOXES AND PUT ALL THE GRAPH OPTIONS IN THERE */ 
		(scatter `id' `xboxcenter' `iw' if `use' == 1, 
			`boxopts'		
			yscale(range(`DYmin' `DYmax') noline reverse)
			ylabel(none) ytitle("")
			xscale(range(`AXmin' `AXmax') noline)
			xlabel(none)
			yline(`borderline', lwidth(thin) lcolor(gs12))
			xtitle("") `legendoff' xtick(""))
			
	 /*HERE ARE GRIDS */
		`betweengrids'			
	 /*HERE ARE THE CONFIDENCE INTERVALS */
	
		`cicommand' /*`cicommand1' `cicommand2'*/
		`smoothcommands1'	`smoothcommands10' `smoothcommands11' `smoothcommands2' `smoothcommands21' `smoothcommands20'	
	 /*ADD ARROWS */
		(pcarrow `id' `uci' `id' `lci' if `leftarrow' == 1 &  `use' == 1 , `arrowopts')	
		(pcarrow `id' `lci' `id' `uci' if `rightarrow' == 1 &  `use' == 1, `arrowopts')	
		(pcbarrow `id' `lci' `id' `uci' if `biarrow' == 1 &  `use' == 1, `arrowopts')
		
	 /*DIAMONDS FOR SUMMARY ESTIMATES  */
		`diamondcommand1' `diamondcommand11' `diamondcommand21'
		`diamondcommand2' `diamondcommand12' `diamondcommand22'
		`diamondcommand3' `diamondcommand13' `diamondcommand23'
		`diamondcommand4' `diamondcommand14' `diamondcommand24'
		
	 /*HERE ARE THE PREDICTION INTERVALS */
		`cipred0'		
	 /*ADD ARROWS */
		`cipred1'
	 /*POINTS */
	 
		(scatter `id' `es' if `use' == 1 , `pointopts')
		
	//overall & sublines	
		`overallCommand' `sublineCommand'
		
	//Others	
		/*`smoothcommands2' `overallCommand'*/	
		,`plotopts' `legendon' `plotname' 
		;
		#delimit cr	
		
		/*	
		if "$by_index_" != "" {
			qui graph dir
			local gnames = r(list)
			local gname: word 1 of `gnames'
			tokenize `gname', parse(".")
			local gname `1'
			if "`3'" != "" {
				local ext =".`3'"
			}
			
			qui graph rename `gname'`ext' `gname'_$by_index_`ext', replace
			if "`graphsave'" != "" {
				graph save `graphsave'_$by_index, replace
			}
		}
		else {
			if "`graphsave'" != "" {
				di _n
				graph save `graphsave', replace
			}			
		}*/
		if `"`graphsave'"' != `""' {
			di _n
			noi graph save `graphsave', replace
		}
		restore
end