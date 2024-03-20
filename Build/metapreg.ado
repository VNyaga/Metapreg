/*
CREATED:	8 Sep 2017
AUTHOR:		Victoria N Nyaga
PURPOSE: 	Generalized linear fixed, mixed & random effects modelling of binomial data.
VERSION: 	3.0.2
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
31.10.2023					clolog ink (might be better if p > .9)	lolog ink (might be better if p < .1)	
10.01.2024					Fit FE/Hexact if RE fails.
							If isq = ., suppress the text in the graph
01.02.2024					Introduce beta-binomial regression :cbbetabin - common beta beta-binomial, crbetabin - common rho beta-binomial	
05.02.2024					Introduce catterpillar plot	
							seperate forest and catterpillar plot options
26.02.2024					Introduce Bayesian in version >16.1
							inference(frequentist|bayesian)
							More options for rr ci's
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
		ALphasort
		CImethod(string) //i=[wald, exact, score], o=[wald, exact, score, t]
		GOF //Goodness of fit
		DOWNload(string asis) 
		DP(integer 2)
		POwer(integer 0)		
		Level(integer 95) 
		INTeraction
		LABEL(string) 
		Model(string) //model(random|mixed|fixed|hexact|cbbetabin|crbetabin, options)
		INFerence(string) //inference(FREQuentist|BAYesian)
		noGRaph //Synonym with nofplot
		noFPlot
		CATPplot
		noSUBgroup 
		noOVerall 
		noITAble
		noWT
		outplot(string) //abs|rr|or
		SUMTable(string) //none|logit|abs|rr|all
		DESign(string asis) //design(general|matched-mcbnetwork|paired-pcbnetwork|comparative|network-abnetwork, baselevel(string) | cov(inde|unstr))
		MC /*Model comparison - Saves time*/
		PROGress /*See the model fitting*/
		by(varname)
		STRatify  /*Stratified analysis, requires byvar()*/
		link(string) /*logit|cloglog|loglog*/
		SMooth nsims(integer 800) //max dim in stata ic
		FOptions(string asis) /*Options specific to the forest plot*/
		COptions(string asis) /*Options specific to the catterpilar plot*/

		/*passthrough*/
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
		PREDciOpt(passthru)
		RCols(varlist) 
		PREDIction  //prediction
		SORtby(varlist) //varlist
		LCols(varlist) 		
		SUBLine
		SUMMARYonly
		SUMStat(string)
		TEXts(real 1.0) 
		XLAbel(passthru)
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
			modeles modellci modeluci holder uniqstudyid clone clones strata numsid
			
	tempname nltestRR nltestOR mctest samtrix rawest rawesti logodds rrout orout absout logoddsi orouti rrouti absouti  exactabsouti exactabsout absexact ///
			coefmat coefvar BVar WVar  omat isq2 bghet bshet lrtestp V dftestnl ptestnl lrtest matgof ///
			outr absoutp absoutpi hetout hetouti popabsout popabsouti poprrout poprrouti poporout ///
			poporouti poplorout poplorouti exactorouti  exactlorouti exactorout  exactlorout
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
		cap gen `modeles' = .
		cap gen `modellci' = .
		cap gen `modeluci' = .
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
			cap assert ("`design'" == "comparative")
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
			else {
				di as error "`cov' not allowed in specifying the design()"
				exit
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
	local fopts `"`options'"'
	
	/*Check if variables exist*/
	foreach var of local varlist {
		cap confirm var `var'
		if _rc!=0  {
			di in re "Variable `var' not in the dataset"
			exit _rc
		}
	}
	
	//General housekeeping	
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
		
	//Check if crbetabin is installed
	if "`model'" == "crbetabin" {
		capture which betabin
		if _rc != 0 {
			di as res "The user-package betabin is required"
			di  `"{stata "search betabin": Click to search and install the package}"'
			exit
		}
	}
	//Logscale only in or/rr 
	if "`logscale'" != "" {
		cap assert "`outplot'" != "abs"
		if _rc != 0 {
			di as res "Option `logscale' not allowed"
			local logscale
		}
	}
	//smooth is redundant in betabin, we cannot recover the individual estimates
	/*if ("`model'" == "betabin" ) & "`smooth'" != "" {
		local smooth
		di as res "The option -smooth- is ignored. The model-based study estimates are irrecoverable."
	}*/
	
	//Weights only for the logit link
	/*if "`link'" != "logit" {
		local wt nowt
	}*/

	//Avoid Incosistencies & Redundancies
	if "`stratify'" != "" & "`summaryonly'" != "" {
		local wt nowt
	}
	
	/*
	if "`outplot'" == "abs" & "`design'" == "comparative" & "`model'" == "fixed" {
		local design "general" 	
	}
	*/
	
	qui count
	if `=r(N)' < 2 {
		di as err "Insufficient data to perform meta-analysis"
		exit 
	}
	if `=r(N)' < 3 & "`model'" != "hexact" {
		local model hexact //If less than 3 studies, use exact
		di as res _n  "Note: Homo-exact model imposed whenever number of studies is less than 3."
		if "`modelopts'" != "" {
			local modelopts
			di as res _n  "Warning: Model options ignored."
			di as res _n  "Warning: Consider re-specifying options for the fixed-effects model should the model not converge."
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
	if _rc!=0 {
		di as err "The study identifier variable is not specified"
		di as err "Specify it with STUDYID(varname) "
		exit _rc
	}
	
	tokenize `varlist'
	if "`design'" == "general" | "`design'" == "comparative" | "`design'" == "abnetwork"  {
		local event = "`1'"
		local total = "`2'"
		
		*gen `total' = `2'
		*gen `event' = `1'
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
			cap confirm string variable `first'
			if _rc != 0 {
				di as error "The first covariate in `design' analysis should be a string"
				exit _rc
			}
		}
		if ("`by'" != "") & ("`outplot'" != "rr") {
			cap assert ("`first'" != "`by'") 
			if _rc != 0 { 
					di as error "Remove the option by(`by') or specify a different by-variable"
					exit _rc
			}
		}
	}
	
	if "`outplot'" == "" {
		local outplot = "abs"
	}
	else {
		if "`outplot'" == "rr" | "`outplot'" == "or"  {
			cap assert "`design'" != "general"
			if _rc != 0 {
				di as error "Option outplot(`outplot') only avaialable for comparative/mcbnetwork/pcbnetwork/abnetwork designs with first covariate as string"
				di as error "Specify the first string covariate and the appropriate design(comparative/mcbnetwork/pcbnetwork/abnetwork)"
				exit _rc
			}
		}
	}
	
	if ("`design'" == "mcbnetwork" | "`design'" == "pcbnetwork") & ("`outplot'" != "or" ) {
		local outplot = "rr"  //default
	}
	
	//
	if ("`model'" == "hexact") & ("`outplot'" == "rr") {
		local outplot = "or"  //only option
	}
	
	cap assert ("`outplot'" == "rr") | ("`outplot'" == "or") | ("`outplot'" == "abs") 
	if _rc != 0  {
		di as error "Invalid option in outplot(`outplot')"
		exit _rc
	}	

	if "`sumstat'" == "" {
		if "`outplot'" == "abs" {
			local sumstat = "Proportion"
		}
		else if "`outplot'" == "rr"  {
			local sumstat = "Proportion Ratio"
		}
		else if "`outplot'" == "or"  {
			local sumstat = "Odds Ratio"
		}
	}
	
	if "`outplot'" != "abs"  {
		if "`cimethod'" != "" & (strpos("`cimethod'", ",") != 1) { 
			tokenize "`cimethod'", parse(",")
			local icimethod "`1'"
			
			if "`3'" != "" {
				local ocimethod = strltrim("`3'")
			}
		}
		if "`cimethod'" != "" & (strpos("`cimethod'", ",") == 1) { 
			tokenize "`cimethod'", parse(",")		
			if "`2'" != "" {
				local ocimethod = strltrim("`2'")
			}
		}
		
		//iOR
		if ("`icimethod'" != "") & "`outplot'" == "or" {
			if (strpos("`icimethod'", "ex") != 1) & (strpos("`icimethod'", "wo") != 1) & (strpos("`icimethod'", "co") != 1)  {
				di as error "Option `icimethod' not allowed in cimethod(`cimethod')"
				exit	
			}
		}
		if ("`icimethod'" == "") & "`outplot'" == "or"  {
			local icimethod "woolf"
		}
		
		//iRR
		if "`design'" == "mcbnetwork" & "`outplot'" == "rr"  {
			local icimethod "CML"
		}
		if ("`design'" == "pcbnetwork" | "`design'" == "comparative") & "`outplot'" == "rr"   {
			if ("`icimethod'" != "") {
				if (strpos("`icimethod'", "koo") != 1) & (strpos("`icimethod'", "ka") != 1) & (strpos("`icimethod'", "bail") != 1) & (strpos("`icimethod'", "adlo") != 1) {
					di as error "Option `icimethod' not allowed in cimethod(`cimethod')"
					exit	
				}
			}
			if ("`icimethod'" == "")  {
				local icimethod "koopman"
			}
		}
		
		//Summary
		if "`ocimethod'" != "" {
			if (strpos("`ocimethod'", "t") != 1) &  (strpos("`ocimethod'", "w") != 1){
				di as error "Option `ocimethod' not allowed in cimethod(`cimethod')"
				exit	
			}
		}
		else {
			local ocimethod "wald"
		}
	}
	if "`outplot'" == "abs"  {
		if "`cimethod'" != "" { 
			tokenize "`cimethod'", parse(",")
			local icimethod "`1'"
			if "`3'" != "" {
				local ocimethod = strltrim("`3'")
			}
		}
		if "`icimethod'" != "" {
			if (strpos("`icimethod'", "ex") != 1) & (strpos("`icimethod'", "wi") != 1) &  (strpos("`icimethod'", "wa") != 1) & (strpos("`icimethod'", "e") != 1) & (strpos("`icimethod'", "ag") != 1) & (strpos("`icimethod'", "je") != 1)   {
				di as error "Option `icimethod' not allowed in cimethod(`cimethod')"
				exit	
			}
		}
		else {
			local icimethod "wilson"
		}
		if "`inference'" == "frequenstist" {
			if "`ocimethod'" != "" {
				if "`model'" == "random" | "`model'" == "fixed" {
					if (strpos("`ocimethod'", "t") != 1) &  (strpos("`ocimethod'", "w") != 1){
						di as error "Option `ocimethod' not allowed in cimethod(`cimethod')"
						exit	
					}
				}
				else {
					if (strpos("`ocimethod'", "ex") != 1) & (strpos("`ocimethod'", "wi") != 1) &  (strpos("`ocimethod'", "wa") != 1) & (strpos("`icimethod'", "e") != 1) & (strpos("`ocimethod'", "ag") != 1) & (strpos("`ocimethod'", "je") != 1)   {
						di as error "Option `ocimethod' not allowed in cimethod(`cimethod')"
						exit	
					}
				}
			}
			else {
				if "`model'" == "random" | "`model'" == "fixed" { 
					local ocimethod "wald" 
				}
				else {
					local ocimethod "wilson"
				}
			}
		}
		else {
			//Bayesian credible intervals; default is equal-tailed interval
			if "`ocimethod'" != "" {
				if "`ocimethod'" != "hpd" | "`ocimethod'" != "eti" {
					di as error "Option `ocimethod' not allowed in cimethod(`cimethod')"
					exit
				}
			}
			if "`ocimethod'" == "eti" {
				local ocimethod
			}
		}
	}
		
	if "`prediction'" != "" & "`outplot'" != "abs" {
		local prediction 
		di as res "NOTE: Predictions only computed for absolute measures. The option _prediction_ will be ignored"
	}
		
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
	if "`design'" == "mcbnetwork" | "`design'" =="pcbnetwork" {
		longsetup `varlist', rid(`rid') assignment(`assignment') event(`event') total(`total') idpair(`idpair') `design'

		qui gen `ipair' = "Yes"
		qui replace `ipair' = "No" if `idpair'
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
		local baselab:label `varx' `basecode'
		
		if `basecode' == 1 {
			local indexcode "2"
		}
		else {
			local indexcode "1"
		}
		local indexlab:label `varx' `indexcode'
		
		if "`outplot'" != "abs" {
			local varxlabs "`varx' `indexlab' `baselab'"
		}
	}
		
	if "`design'" == "pcbnetwork" | "`design'" == "mcbnetwork" { 
		local varx = "`ipair'"
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
			local catregs = "`catreg' `Index'"
		}

		if "`design'" == "comparative" {
			local catregs = "`catreg' `varx'" 
		}
		if "`design'" == "abnetwork" {
			tokenize `catreg'
			macro shift
			local catregs "`*'"
		}
		if "`design'" == "general" {
			local catregs "`catreg'"
		}
	}
	
	if "`subgroup'" == "" & ("`catreg'" != "" | "`typevarx'" =="i" ) {
		if "`outplot'" == "abs" {
			if "`typevarx'" =="i" {
				local groupvar = "`varx'"
			}
			else {
				local groupvar : word 1 of `catreg'
			}
		}
		if "`outplot'" != "abs" & "`varx'" != "" {
			local groupvar : word 1 of `catreg'
		}
	}
	
	if "`by'" != "" {
		local groupvar "`by'"
		local byvar "`by'"
		*How many times to loop
		qui label list `by'
		local nlevels = r(max)
	}
	if "`design'" == "abnetwork" {
		local groupvar "`first'"
		local overall "nooverall"
		if "`outplot'" != "abs" {
			local itable "noitable"
		}
	} 
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
	//nullify groupvar if its the studyid
	if "`groupvar'" == "`studyid'" {
		local groupvar
	}	
	if "`groupvar'" == "" {
		local subgroup nosubgroup
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
			local hetdim 4
		}
	}
	
	if ("`model'" == "cbbetabin" ) | ("`design'" == "general") {
		egen `numsid' = group(`rid')
	}
	else {
		qui egen `numsid' = group(`studyid')
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
		if (`i' > `nlevels') & ("`stratify'" != "") & ("`design'" == "comparative") {
			local overall "nooverall"
			continue, break
		}
		
		*Stratify except the last loop for the overall
		if (`i' < `=`nlevels' + 1') & ("`stratify'" != "") {
			local strataif `"if `by' == `i'"'
			local ilab:label `by' `i'
			local stratalab `":`by' = `ilab'"'
			local ilab = ustrregexra("`ilab'", " ", "_")
			local byrownames = "`byrownames' `by':`ilab'"
			local byrowname = "`by'|`ilab'"
			if "`design'" == "comparative" & "`stratify'" != "" {
				local bybirownames = "`bybirownames' `ilab':`baselab' `ilab':`ilab'"
			}

			*Check if there is enough data in each strata
			//Number of obs in the analysis
			qui egen `obsid' = group(`rid') if `by' == `i'
			qui summ `obsid'
			local Nobs= r(max)
			drop `obsid'	

			//Number of studies in the analysis
			qui egen `uniq' = group(`studyid') if `by' == `i'
			qui summ `uniq'
			local Nuniq = r(max)
			drop `uniq'	
		}
		else {
			//Skip if overall not needed
			if "`overall'" != "" & (`i' > `=`nlevels'+1')  & ("`stratify'" != "" | "`design'" == "comparative")   {
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
		*if (`Nobs' < 3 & "`modeli'" != "hexact" & "`design'" != "comparative") | ((`Nobs' < 5 ) & ("`modeli'" == "random") & ("`design'" == "comparative")) {
		if `Nuniq' < 3 & "`modeli'" == "random"  {
			local modeli fixed //If less than 3 studies, use exact model
			if "`modeloptsi'" != "" {
				local modeloptsi
				noi di as res _n  "Warning: `model'-effects model options ignored."
				noi di as res _n  "Warning: Homo-exact model fitted instead."
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
			preg `event' `nonevent' `total' `strataif', rid(`rid') sid(`numsid') studyid(`studyid') use(`use') regexpression(`regexpression') regexpression2(`regexpression2') nu(`nu')  ///
				regressors(`regressors')  catreg(`catreg') contreg(`contreg') level(`level') varx(`varx') typevarx(`typevarx')  /// 
				`progress' model(`modeli') modelopts(`modeloptsi') `mc' `interaction' `design' by(`by') `stratify' baselevel(`basecode') ///
				comparator(`Comparator') cimethod(`ocimethod') `gof' nsims(`nsims') link(`link') bayesrepsfilename(`metapregbayesreps') ///
				modeles(`modeles')  modellci(`modellci') modeluci(`modeluci') outplot(`outplot') `smoothi' cov(`cov') `computewti' ///
				inference(`inference') refsampling(`refsampling')
	
			mat `rawesti' = r(rawest)
			mat `popabsouti' = r(popabsout)
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
					mat `poprrouti' = r(poprrout)
					mat `orouti' = r(orout)
					mat `poporouti' = r(poporout)
					mat `poplorouti' = r(poplorout)
					local inltest = r(inltest)
					if "`inltest'" == "yes" & "`stratify'" == "" {
						mat `nltestRR' = r(nltestRR) 
						mat `nltestOR' = r(nltestOR) 
					}
				}
			}
			else {
				local rr "norr"
			}
			/*
			if (`p' > 0) & ("`mc'" =="") { 
				mat `mctest' = r(mctest) 
			}*/
			
			mat `absouti' = r(absout)
			mat `absoutpi' = r(absoutp)
			if strpos("`modeli'", "random") !=0 | strpos("`model'", "betabin") != 0 { 
				mat `hetouti' = r(hetout)
			}
			else {
				mat `hetouti' = J(1, `hetdim', .)
			}			
		}
		*if 1 study or exact inference
		else {
			mat `rawesti' = J(1, 6, .)
			mat `popabsouti' = J(1, 6, .)
			mat `exactabsouti' = J(1, 11, .)
			mat `absouti' = J(1, 6, .)
			mat `absoutpi' = J(1, 2, .)
			mat `hetouti' = J(1, `hetdim', .)		
			mat `rrouti' = J(1, 6, 1)
			mat `poprrouti' = J(1, 6, 1)
			mat `orouti' = J(1, 6, 1)
			mat `poporouti' = J(1, 6, 1)
			mat `poplorouti' = J(1, 6, 1)

			mat `exactorouti' = J(1, 4, .)
			mat `exactlorouti' = J(1, 4, .)
			
			mat rownames `rawesti' = Overall
			mat rownames `popabsouti' = Overall
			mat rownames `exactabsouti' = Overall
			mat rownames `absouti' = Overall
			mat rownames `absoutpi' = Overall
			mat rownames `hetouti' = Overall		
			mat rownames `rrouti' = Overall
			mat rownames `poprrouti' = Overall
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
			mat roweq `absoutpi' = `byrowname'
			mat roweq `rawesti' = `byrowname'
			if "`rr'" == "" {
				if "`model'" == "hexact" {
					mat roweq `exactorouti' = `byrowname'
					mat roweq `exactlorouti' = `byrowname'
				}
				else {
					mat roweq `rrouti' = `byrowname'
					mat roweq `poprrouti' = `byrowname'
					mat roweq `orouti' = `byrowname'
					mat roweq `poporouti' = `byrowname'
					mat roweq `poplorouti' = `byrowname'
				}
			}
		}
		
		*Stack the matrices
		if `i' == 1 {
			mat `absout' =	`absouti'
			if "`rr'" == "" {
				if "`model'" != "hexact" {
					mat `rrout' =	`rrouti'
					mat `poprrout' = `poprrouti'
					
					mat `orout' =	`orouti'
					mat `poporout' = `poporouti'
					mat `poplorout' = `poplorouti'
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
		}
		else {
			mat `absout' = `absout' \ `absouti'
			mat `popabsout' = `popabsout' \ `popabsouti'
			mat `exactabsout' = `exactabsout' \ `exactabsouti'			
			if "`rr'" == "" {
				if "`model'" != "hexact" {
					mat `rrout' = `rrout' \ `rrouti'
					mat `poprrout' = `poprrout' \ `poprrouti'
					
					mat `orout' = `orout' \ `orouti'
					mat `poporout' = `poporout' \ `poporouti'
					mat `poplorout' = `poplorout' \ `poplorouti'
				}
				else {
					mat `exactorout' = `exactorout' \ `exactorouti'
					mat `exactlorout' = `exactlorout' \ `exactlorouti'
				}
			}
			mat `rawest' = `rawest' \ `rawesti'
			mat `absoutp' = `absoutp' \ `absoutpi'
			mat `hetout' = `hetout' \ `hetouti'
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
		else if "`design'" == "mcbnetwork"  {
			di "{phang} `1' + `2'  ~ binomial(p, `1' + `2' + `3' + `4'){p_end}"
			di "{phang} `1' + `3' ~ binomial(p, `1' + `2' + `3' + `4'){p_end}"
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
				di "{phang} `link'(p) = `nu' + `first'.`studyid' + `studyid'{p_end}"	
			}
			if "`cov'" == "" {
				di "{phang}`studyid' ~ N(0, tau2){p_end}"
			}
			if "`cov'" == "independent" {
				di "{phang}`studyid' ~ N(0, tau2){p_end}"
				di "{phang}`first'.`studyid' ~ N(0, sigma2){p_end}"
			}
			if "`cov'" == "unstructured" {
				di "{phang}`studyid',`first'.`studyid'  ~ biv.normal(0, Sigma){p_end}"
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
			di "{phang} `link'(p) = `nu'{p_end}"		
		}
		if "`design'" == "pcbnetwork" | "`design'" == "mcbnetwork" {
			di "{phang} Ipair = 0 if 1st pair{p_end}"
			di "{phang} Ipair = 1 if 2nd pair{p_end}"
		}
		if "`design'" == "abnetwork" {
			di "{phang}`first' ~ N(0, sigma2){p_end}"
			qui label list `first'
			local nfirst = r(max)
		}		
		if ("`catreg'" != " " | "`typevarx'" =="i" | ("`design'" == "comparative" | "`design'" == "mcbnetwork" | "`design'" == "pcbnetwork"))  {
			di _n "{phang}Base levels{p_end}"
			di _n as txt "{pmore} Variable  -- Base Level{p_end}"
		}
		foreach fv of local catregs  {			
			local lab:label `fv' 1
			if "`fv'" != "`studyid'" {
				di "{pmore} `fv'  -- `lab'{p_end}"
			}			
		}
		if "`design'" == "abnetwork" {
			local lab:label `first' `basecode'
			di "{pmore} `first'  -- `lab'{p_end}"
		}
			
		di _n
		di "{phang}" as txt "Number of observations = " as res "`Nobs'{p_end}"
		di "{phang}" as txt "Number of studies = " as res "`Nuniq'{p_end}"
		if "`design'" == "abnetwork" {
			di "{phang}" as txt "Number of `first's = " as res "`nfirst'{p_end}"
		}
				
		if `Nobs' > 1 & "`model'" != "hexact" {
		
			//Display GOF
			if ("`gof'" != "") {
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
					mat rownames `matgof' = Value
					local widthc = 15
				}
				
				mat rownames `matgof' = Value
				
				di _n
				di as text "Goodness of Fit Criterion"
				#delimit ;
				noi matlist `matgof',  
							cspec(& %7s |   %8.`=`dp''f &  %`=`widthc''.`=`dp''f o2&) 
							rspec(&-&) underscore  nodotz
				;
				#delimit cr 
			}
			
			if "`stratify'" != "" {
				qui estimates restore metapreg_modest
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
					qui estimates restore metapreg_modest
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
				capture noisily mcpreg `event' `nonevent' `total' `strataif', sid(`numsid')  regexpression(`regexpression') nu(`nu')  ///
					regressors(`regressors') level(`level') varx(`varx') typevarx(`typevarx')  `progress' /// 
					model(`modeli') modelopts(`modeloptsi')  `interaction' `design' by(`by') `stratify' baselevel(`basecode') ///
					comparator(`Comparator')  link(`link')   ///
					 cov(`cov') inference(`inference') refsampling(`refsampling')
					 
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
				mat colnames `hetout' = DF Chisq p tau2 sigma2 rho I2tau I2sigma 
			}
			else {
				mat colnames `hetout' = DF Chisq p tau2 sigma2 I2tau I2sigma 
			}
		}
		else {
			if (`p' == 0) & ("`model'" == "random") & "`pcbnetwork'`mcbnetwork'" == "" {
				mat colnames `hetout' = DF Chisq p tau2 I2tau 
			}
			else {
				mat colnames `hetout' = DF Chisq p tau2 
			}
		}
	}
	//If stratify, no overall for comparative
	if "`stratify'" != "" & "`design'" == "comparative" {
		local overall "nooverall"
	}
	
	//CI
	if "`outplot'" != "abs" {
		local se
	}
	if "`outplot'" != "abs" & "`design'" == "abnetwork" {
		local summaryonly "summaryonly"
		local smooth
	}
	
	if "`design'" == "mcbnetwork" | "`design'" == "pcbnetwork" {
		*widesetup `event' `total', sid(`rid') idpair(`assignment')  jvar(`comparator')

		sort `rid'
		cap drop `assignment' `ipair' `nonevent'
		qui reshape wide `event' `total' _WT `modeles' `modellci' `modeluci', i(`rid') j(`idpair')

		//Add the weights
		qui gen _WT = _WT0 + _WT1
		qui drop _WT0  _WT1
		
		qui gen `modeles' = `modeles'1
		qui drop `modeles'0 `modeles'1
		
		qui gen `modellci' = `modellci'1
		qui drop `modellci'0 `modellci'1
		
		qui gen `modeluci' = `modeluci'1
		qui drop `modeluci'0 `modeluci'1
	}
		
	metapregci `depvars', studyid(`studyid') first(`first') es(`es') se(`se') uci(`uci') lci(`lci') `design' ///
		id(`id') rid(`rid') regressors(`regressors') outplot(`outplot') level(`level') ///
		cimethod(`icimethod') lcols(`lcols') rcols(`rcols')  sortby(`sortby') by(`by')  ///
		modeles(`modeles') modellci(`modellci') modeluci(`modeluci') `smooth'
	
	local depvars = r(depvars)
	local rcols = r(rcols)
	local lcols = r(lcols)
	local sortby = r(sortby)
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
		
	prep4show `id' `cid' `use' `neolabel' `es' `lci' `uci' `modeles' `modellci' `modeluci', `design' ///
		sortby(`sortby') groupvar(`groupvar') grptotal(`grptotal') se(`se') 	///
		outplot(`outplot') rrout(`rrout') poprrout(`poprrout') orout(`orout') poporout(`poporout') exactorout(`exactorout') ///
		absout(`absout') popabsout(`popabsout') exactabsout(`exactabsout') absoutp(`absoutp') hetout(`hetout')	///
	    `subgroup' `summaryonly' dp(`dp') pcont(`pcont') model(`model') `prediction'	///
		`overall' download(`download') indvars(`indvars') depvars(`depvars') `stratify' level(`level') `enhance'
		
	
	if "`model'" == "hexact"  {
			//exactabs
		if  (("`sumtable'" == "all") |(strpos("`sumtable'", "abs") != 0)) {
			printmat, matrixout(`exactabsout') type(exactabs) p(`p') dp(`dp') power(`power') `continuous'  model(`model')
		}
		
		//or
		if (("`sumtable'" == "all") | (strpos("`sumtable'", "or") != 0)) & (("`catreg'" != " ") | ("`typevarx'" == "i"))   {
			
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
		if strpos("`model'", "random") != 0 | strpos("`model'", "betabin")!= 0  {			
			printmat, matrixout(`hetout') type(het) dp(`dp') `design' model(`model')
		}
		//Raw estimates
		if  (("`sumtable'" == "all") |(strpos("`sumtable'", "logit") != 0)) {
			printmat, matrixout(`rawest') type(raw) p(`p') dp(`dp') power(`power') `continuous' model(`model') `stratify' link(`link')
		}
		//abs
		if  (("`sumtable'" == "all") |(strpos("`sumtable'", "abs") != 0)) {
			printmat, matrixout(`absout') type(abs) p(`p') dp(`dp') power(`power') `continuous'  model(`model')
		}
		//Pop p 
		if  (("`sumtable'" == "all") | (strpos("`sumtable'", "abs") != 0))  {
			printmat, matrixout(`popabsout') type(popabs) dp(`dp') power(`power') nsims(`nsims') model(`model')
		}
		//rr
		if (("`sumtable'" == "all") | (strpos("`sumtable'", "rr") != 0)) & (("`catreg'" != " ") | ("`typevarx'" == "i"))   {
			//rr
			cap confirm matrix `rrout'
			if _rc == 0 {
				printmat, matrixout(`rrout') type(rr) p(`p') dp(`dp') power(`power')  model(`model')
			}
			
			//rr equal
			if "`inltest'" == "yes" {
				cap confirm matrix `nltestRR'
				if _rc == 0 {
					printmat, matrixout(`nltestRR') type(rre) dp(`dp') inference(`inference')
				}
			}
			cap confirm matrix `poprrout'
			if _rc == 0 {
				printmat, matrixout(`poprrout') type(poprr) p(`p') dp(`dp') power(`power')  model(`model') nsims(`nsims')
			}
		}
		
		//or
		if (("`sumtable'" == "all") | (strpos("`sumtable'", "or") != 0)) & (("`catreg'" != " ") | ("`typevarx'" == "i"))   {
			
			//rr
			cap confirm matrix `orout'
			if _rc == 0 {
				printmat, matrixout(`orout') type(or) p(`p') dp(`dp') power(`power')  model(`model')
			}
			
			//or equal
			if "`inltest'" == "yes" {
				cap confirm matrix `nltestOR'
				if _rc == 0 {
					printmat, matrixout(`nltestOR') type(ore) dp(`dp') inference(`inference')
				}
			}
			cap confirm matrix `poprrout'
			if _rc == 0 {
				printmat, matrixout(`poporout') type(popor) p(`p') dp(`dp') power(`power')  model(`model') nsims(`nsims')
			}
		}
	}
	
	if "`itable'" == "" {
		sort `id'
		
		disptab `id'  `use' `neolabel' `es' `lci' `uci' `grptotal' `modeles' `modellci' `modeluci', ///
			`itable' dp(`dp') power(`power') design(`design') `summaryonly' ///
			`subgroup' sumstat(`sumstat') level(`level') `wt' `smooth' ///
			ocimethod(`ocimethod') icimethod(`icimethod') model(`model') groupvar(`groupvar') outplot(`outplot')
    }
	
	if "`graph'" == "" {
		//Gather graph options
		if "`astext'" != "" {
			local astext "astext(`astext')"
		}
		if "`texts'" != "" {
			local texts "texts(`texts')"
		}
		if "`lcols'" != "" {
			local lcols "lcols(`lcols')"
		}
		if "`rcols'" != "" {
			local rcols "rcols(`rcols')"
		}
		local goptions `"`lcols' `rcols' `overall' `ovline' `stats' `box' `double' `astext' `ciopts' `diamopts' `olineopts' `pointopts' `boxopts' `predciopts' `prediction' `subline' `texts' `xlabel' `xline' `xtick'  `logscale' `options'"'
		
		metapplotcheck,`summaryonly' `goptions'  //plots advance housekeeping
		local goptions = r(plotopts)
		local glcols = r(lcols)
		local grcols = r(rcols)
		if "`glcols'" == " " { //if empty
			local glcols
		}
		if "`grcols'" == " " { //if empty
			local grcols
		}
	}
	
	//Draw the forestplot
	if "`graph'`fplot'" == "" {
		if `"`foptions'"' != "" {
			metapplotcheck,`summaryonly' `goptions' `foptions'  //Forest plot advance housekeeping
			local foptions = r(plotopts)
			local flcols = r(lcols)
			local frcols = r(rcols)
			if "`flcols'" == " " { //if empty
				local flcols
			}
			if "`frcols'" == " " { //if empty
				local frcols
			}
		}
		else /*`"`goptions'"' != ""*/ {
			local foptions `"`goptions'"'
			if "`glcols'" != "" {
				local flcols = "`flcols' `glcols'"
			}
			if "`grcols'" != "" {
				local frcols = "`frcols' `grcols'"
			}
		}
		
		sort `id' 
				
		metapplot `es' `lci' `uci' `use' `neolabel' `grptotal' `id' `modeles' `modellci' `modeluci', model(`model') ///	
			studyid(`studyid') power(`power') dp(`dp') level(`level') ///
			groupvar(`groupvar')  type(fplot) sumstat(`sumstat') ///
			outplot(`outplot') lcols(`flcols') rcols(`frcols') ///
			`foptions'  design(`design') `wt' `smooth' varxlabs(`varxlabs')
	}
	
	//Draw the catterpillar plot
	if "`graph'" == "" & "`catpplot'" != "" {
		if "`coptions'" != "" {
			metapplotcheck,`summaryonly' `goptions' `coptions'  //Catterpillar plot advance housekeeping
			local coptions = r(plotopts)
			local clcols = r(lcols)
			local crcols = r(rcols)
			if "`clcols'" == " " { //if empty
				local clcols
			}
			if "`crcols'" == " " { //if empty
				local crcols
			}
		}
		else {
			local coptions = `" `goptions'"'
			if "`glcols'" != "" {
				local clcols = "`clcols' `glcols'"
			}
			if "`grcols'" != "" {
				local crcols = "`crcols' `grcols'"
			}
		}
		
		sort `cid'
		
		metapplot `es' `lci' `uci' `use' `neolabel' `grptotal' `cid' `modeles' `modellci' `modeluci', model(`model') ///	
			studyid(`studyid') power(`power') dp(`dp') level(`level') ///
			groupvar(`groupvar')  type(catpplot) sumstat(`sumstat') ///
			outplot(`outplot') lcols(`clcols') rcols(`crcols') ///
			`coptions' design(`design') `wt' `smooth' varxlabs(`varxlabs')
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

	cap confirm matrix `mctest'
	if _rc == 0 {
		ereturn matrix mctest = `mctest'
	}
	cap confirm matrix `hetout'
	if _rc == 0 {
		ereturn matrix hetout = `hetout'
	}
	cap confirm matrix `nltestRR'
	if _rc == 0 {
		ereturn matrix rrtest = `nltestRR'
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
	cap confirm matrix `absoutp'
	if _rc == 0 {
		ereturn matrix absoutp = `absoutp'
	}
	cap confirm matrix `rrout'
	if _rc == 0 {
		ereturn matrix rrout = `rrout'
		ereturn matrix poprrout = `poprrout'
		ereturn matrix orout = `orout'
		ereturn matrix poporout = `poporout'
		ereturn matrix poplorout = `poplorout'
	}
		
	restore	
end

/*+++++++++++++++++++++++++	SUPPORTING FUNCTIONS: ESTCOVAR +++++++++++++++++++++++++
							Compose the var-cov matrix
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
cap program drop estcovar
program define estcovar, rclass

	syntax, matrix(name) cov(string) [predcmd(string)]
	*matrix is colvector
	tempname matcoef rosevar rawvar
	mat `matcoef' = `matrix''
	local nrows = rowsof(`matcoef')
	*Initialize - Default
	mat	`rosevar' = (0, 0\ ///
				0, 0)
	mat	`rawvar' = (0, 0\ ///
				0, 0)			
	if "`predcmd'" == "meqrlogit_p" {
		if strpos("`cov'", "uns") != 0 {
			mat	`rosevar' = (exp(`matcoef'[`nrows' - 1 , 1])^2, exp(`matcoef'[ `nrows' - 1, 1])*exp(`matcoef'[`nrows' - 2, 1])*tanh(`matcoef'[ `nrows', 1])\ ///
						exp(`matcoef'[ `nrows' - 1, 1])*exp(`matcoef'[`nrows' - 2, 1])*tanh(`matcoef'[ `nrows', 1]), exp(`matcoef'[ `nrows' - 2, 1])^2)
						
			mat	`rawvar' = (`matcoef'[`nrows' - 1 , 1], `matcoef'[ `nrows', 1]\ ///
						`matcoef'[ `nrows', 1], `matcoef'[ `nrows' - 2, 1])			
			local k = 3
		}		
		else if strpos("`cov'", "ind") != 0 {
			mat	`rawvar' = (`matcoef'[ `nrows', 1], 0\ ///
						0, `matcoef'[ `nrows' - 1, 1])
			mat	`rosevar' = (exp(`matcoef'[ `nrows', 1])^2, 0\ ///
						0, exp(`matcoef'[ `nrows'-1, 1])^2)			
			local k = 2
		}
	}
	else {
		if strpos("`cov'", "uns") != 0 {
			mat	`rosevar' = (`matcoef'[`nrows' - 1 , 1], `matcoef'[ `nrows', 1] \ ///
						`matcoef'[ `nrows', 1], `matcoef'[ `nrows' - 2, 1])
						
			mat	`rawvar' = (`matcoef'[`nrows' - 1 , 1], `matcoef'[ `nrows', 1] \ ///
						`matcoef'[ `nrows', 1], `matcoef'[ `nrows' - 2, 1])			
			local k = 3
		}		
		else if strpos("`cov'", "ind") != 0 {
			mat	`rawvar' = (`matcoef'[ `nrows', 1], 0 \ ///
						0, `matcoef'[ `nrows' - 1, 1])
						
			mat	`rosevar' = (`matcoef'[ `nrows', 1], 0 \ ///
						0, `matcoef'[ `nrows'-1, 1])			
			local k = 2
		}

	}
	
	return matrix rosevar = `rosevar' 
	return matrix rawvar = `rawvar' 
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
		cimethod(string) lcols(varlist) rcols(varlist) mcbnetwork pcbnetwork sortby(varlist) 
		comparative abnetwork general modeles(varname) modellci(varname) modeluci(varname) smooth
		];
	#delimit cr
	tempvar uniq event event1 event2 total total1 total2 a b c d idpair
	*gettoken idpair confounders : regressors
	
	tokenize `varlist'
	if "`mcbnetwork'`pcbnetwork'" == "" {
		generate `event' = `1'
		generate `total' = `2'
		local depvars "`1' `2'"
	}
	else if "`mcbnetwork'" != ""  {
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
	
	if "`outplot'" != "abs" {
		if "`abnetwork'" != "" {
			gen `id' = _n
			gen `es' = .
			gen `lci' = .
			gen `uci' = .
		}
		if "`mcbnetwork'" != "" { 
			if "`outplot'" == "rr" {
				//constrained maximum likelihood estimation
				cmlci `a' `b' `c' `d', r(`es') upperci(`uci') lowerci(`lci') alpha(`=1 - `level'*0.01')
			}
			if "`outplot'" == "or" {
				orccci `a' `b' `c' `d', r(`es') upperci(`uci') lowerci(`lci') level(`level') `mcbnetwork' cimethod(`cimethod')
			}
		}
		if "`pcbnetwork'" !="" {
			if "`outplot'" == "rr" {
				rrci `event1' `total1' `event2' `total2', r(`es') upperci(`uci') lowerci(`lci') alpha(`=1 - `level'*0.01') cimethod(`cimethod')
			}
			if "`outplot'" == "or" {
				orccci `event1' `total1' `event2' `total2', r(`es') upperci(`uci') lowerci(`lci') level(`level')  cimethod(`cimethod')
			}
		}
		if "`comparative'" != "" {
			egen `id' = group(`studyid' `by')
			sort `id' `rid'
			by `id': egen `idpair' = seq()

			qui count
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
			
			if "`outplot'" == "rr" {
				rrci `event'1 `total'1 `event'0 `total'0, r(`es') upperci(`uci') lowerci(`lci') alpha(`=1 - `level'*0.01') cimethod(`cimethod')
			}
			if "`outplot'" == "or" {
				orccci `event'1 `total'1 `event'0 `total'0, r(`es') upperci(`uci') lowerci(`lci') level(`level')  cimethod(`cimethod')
			}
						
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
			qui {
				gen _WT = 0 
				replace _WT = _WT  +  _WT_1 if _WT_1 != .
				replace _WT = _WT  +  _WT_2 if _WT_2 != .
			}
			
			*qui gen _WT = _WT_1 + _WT_2

			qui drop _WT_1  _WT_2
			
			//Remove unnecessary columns
			if "`smooth'" !="" {
				qui drop `modeles'_1
				rename `modeles'_2 `modeles'
				
				qui drop `modellci'_1
				rename `modellci'_2 `modellci'
				
				qui drop `modeluci'_1
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
		metapreg_propci `total' `event', p(`es') se(`se') lowerci(`lci') upperci(`uci') cimethod(`cimethod') level(`level')
		gen `id' = _n
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
		cimethod(string) cov(string)
		level(integer 95)
		DP(integer 2)
		progress
		model(string) modelopts(string asis) outplot(string)
		noMC noCONstant
		interaction	
		comparative mcbnetwork pcbnetwork abnetwork
		by(varname) stratify
		GOF nsims(string) link(string) bayesrepsfilename(string asis)
		modeles(varname) modelse(varname) modellci(varname) modeluci(varname) 
		smooth computewt inference(string) refsampling(string)
			*];

	#delimit cr
	marksample touse, strok 
	
	tempfile bayesest bayesreps	
	local bayesreps = subinstr("`bayesreps'", ".tmp", ".dta", 1)
	local bayesest = subinstr("`bayesest'", ".tmp", ".dta", 1)
	
	tempvar event nonevent total invtotal predevent ill iw
	tempname coefmat coefvar testlr V logodds absout absoutp rrout orout nltestRR nltestOR  ///
			 hetout mctest absexact newobs matgof popabsout poprrout poporout poplorout ///
			 rosevar rawvar rawest rawestp coeflor coefor lorci exactabsout exactorout ///
			 exactlorout lnrho bayestest  bayesstats
	
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
	if "`link'" == "loglog" {
		local outcome "`nonevent'"
	}
	else {
		local outcome "`event'"
	}
	
	`echo' fitmodel `outcome' `total' if `touse', `modelopts' model(`model') regexpression(`regexpression') ///
		sid(`sid') level(`level') nested(`first') `abnetwork' cov(`cov') link(`link') p(`p')  ///
		bayesreps(`bayesreps') bayesest(`bayesest') inference(`inference') refsampling(`refsampling')
		
	//Returned model	
	local getmodel = r(model)
	estimates store metapreg_modest
	local lnvar = r(lnvar)
	
	qui {
		replace _ESAMPLE = e(sample) 
		replace `use' = 1 if (_ESAMPLE == 1)
	}
	
	local mdf = .

	if "`inference'" == "bayesian" {
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
				bayespredict `predevent', mean  rseed(1) 			
			}
			else {
				//fixed
				predict `predevent', mu
			}
	
			if "`model'" == "random" {
				//Revert to original data if filler data was generated
				if (`nobs' < 7)  {
					keep if !`newobs'
				}
			}
					
			//compute the weight
			gen `iw' = `total'*(`predevent'/`total')*(1 - `predevent'/`total') if  `predevent' <`total'
			replace `iw' = `total' if  `predevent'==`total'
		
			//compute the relative weight
			sum `iw' if (_ESAMPLE == 1 & mu == 1)
			local W = r(sum)
			
			//compute the weights
			replace _WT = (`iw'/`W')*100 if (_ESAMPLE == 1 & mu == 1) & (_WT == .) & (`iw' != .)
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
		if "`cov'" != "" {
			estcovar, matrix(`coefmat') cov(`cov') predcmd(`predcmd')
			local DF_BHET = r(k)
			mat `rosevar' = r(rosevar)  //var-cov matrix
			mat `rawvar' = r(rawvar)  //raw var-cov matrix
			mat colnames `rosevar' = intercept slope
			mat rownames `rosevar' = intercept slope
			mat colnames `rawvar' = intercept slope
			mat rownames `rawvar' = intercept slope
		}
	}
	else if "`getmodel'" == "crbetabin" {
		local BHET = e(chi2_c)
		local P_BHET = `=chi2tail(1, e(chi2_c))*0.5'
		local DF_BHET = 1
		local SIGMA = e(sigma)
	}
	else if "`getmodel'" == "cbbetabin" {
		//obtain ln rho  = ln (1/(1 + alpha + beta))
		qui margins , expression(-xb() - _b[_cons]) at(mu==1) 
		
		qui margins , expression(ln(1/(1 + exp(xb()) + exp(_b[_cons])))) at(mu==1) 

		mat `lnrho' = r(table)

		//fit the fixed effects model, to test if overdispersion
		qui fitmodel `outcome' `total' if `touse' & mu == 1, `modelopts' model(fixed) regexpression(`regexpression2') ///
			sid(`sid') level(`level') nested(`first') `abnetwork' cov(`cov') link(`link') p(`p') inference(`inference') refsampling(`refsampling')
			
		qui estimates store metapreg_Null
		
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
	else if "`getmodel'" == "bayesrandom" {	
		tempfile bayesfixedest
		local bayesfixedest = subinstr("`bayesfixedest'", ".tmp", ".dta", 1)
		
		qui fitmodel `outcome' `total' if `touse', `modelopts' model(bayesfixed) regexpression(`regexpression') ///
		sid(`sid') level(`level') nested(`first') `abnetwork' cov(`cov') link(`link') p(`p')  ///
		bayesreps(`bayesreps') bayesest(`bayesfixedest') inference(`inference') refsampling(`refsampling')
		
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
		
		if _rc == 0 {
			mat `bayestest' = r(test)
			local BHET = `bayestest'[2,2] - `bayestest'[1,2] //Difference in marginal likelihood
			local P_BHET = `bayestest'[1,4]  // posterior prob of the model
			local DF_BHET = 1
		}

		estimates drop metapreg_Null
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
			local TAU21 = `scalefn'(`coefmat'[1, `npar'])^`scalepow' //Between study variance	1
			local TAU22 = 0
		}
		else if "`abnetwork'" != "" {
			local TAU21 = `scalefn'(`coefmat'[1, `=`npar'-1'])^`scalepow' //Between study variance	1
			local TAU22 = `scalefn'(`coefmat'[1, `npar'])^`scalepow' //Between study variance	2
		}
		else if "`cov'" != "" {
			local TAU21 = `rosevar'[1, 1]
			local TAU22 = `rosevar'[2, 2]
			if "`cov'" == "unstructured" {
				local rho = tanh(`rawvar'[1, 2])
			}  
		}
	}
	
	if "`getmodel'" == "bayesrandom" {
		qui estimate restore metapreg_modest
		
		if `lnvar' {
			local parmsigma2 = `"(exp({lnsigmasq}))"'
			local parmtau2 = `"(exp({lntausq}))"'
		}
		else {
			local parmsigma2 = `"{sigmasq}"'
			local parmtau2 = `"{tausq}"'
		}
		
		if "`abnetwork'`cov'" == "" {
			qui bayesstats summary `parmtau2'
			mat `rosevar' = r(summary)
		
			local TAU21 = `rosevar'[1, 1] //Between study variance	1
			local TAU22 = 0
		}
		else if "`abnetwork'" != "" |  "`cov'" == "independent" {
				#delimit ;
				qui bayesstats summary `parmtau2' `parmsigma2'
					(isq1:`parmtau2'/(`parmtau2' + `parmsigma2'))
					(isq2:`parmtau2'/(`parmtau2' + `parmsigma2'))
				;
				#delimit cr	
				mat `rosevar' = r(summary)
				local TAU21 = `rosevar'[1, 1]
				local TAU22 = `rosevar'[2, 1]
				local ISQ1 = `rosevar'[3, 1]*100
				local ISQ2 = `rosevar'[4, 1]*100 
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

			mat `rosevar' = r(summary)				
			local TAU21 = `rosevar'[1, 1]
			local TAU22 = `rosevar'[2, 1]			
			local rho = `rosevar'[3, 1]
			local ISQ1 = `rosevar'[4, 1]*100
			local ISQ2 = `rosevar'[5, 1]*100  
		}
	}

	if (`p' == 0) & ("`getmodel'" == "random") & ("`pcbnetwork'`mcbnetwork'" == "") & "`link'" == "logit" {
		/*Compute I2*/				
		qui gen `invtotal' = 1/`total'
		qui summ `invtotal' if `touse'
		local invtotal= r(sum)
		local K = r(N)
		local Esigma = (exp(`TAU21'*0.5 + `coefmat'[1, 1]) + exp(`TAU21'*0.5 - `coefmat'[1, 1]) + 2)*(1/(`K'))*`invtotal'
		local ISQ1 = `TAU21'/(`Esigma' + `TAU21')*100	
	}
	else if ("`abnetwork'`cov'" != "") & ("`getmodel'" == "random") {
		local ISQ1 = `TAU21'/(`TAU21' + `TAU22')*100
		local ISQ2 = `TAU22'/(`TAU21' + `TAU22')*100		
	}

	if "`inference'" == "frequentist" {
		//Raw estimates in logit/cloglog scale
		cap freqsummary, event(`event') total(`total') studyid(`studyid') estimates(metapreg_modest) ///
			`interaction' catreg(`catreg') contreg(`contreg') level(`level') model(`getmodel') cimethod(`cimethod')  ///
			varx(`varx') typevarx(`typevarx') by(`by') regexpression(`regexpression') ///
			`mcbnetwork' `comparative' `pcbnetwork' `abnetwork' `stratify'  ///
			comparator(`comparator') link(`link') 
			
		mat `rawestp' = r(outmatrixp)
		mat `rawest' = r(outmatrix)
		mat `exactabsout' = r(exactabsout)
	
		estp, rawestmat(`rawestp') link(`link') cimethod(`cimethod') 
		mat `absoutp' = r(outmatrix)
		
		//Conditional ABS
		estp, rawestmat(`rawest') link(`link') cimethod(`cimethod') 
		mat `absout' = r(outmatrix)
	}
	else {
		//Conditional odds/abs
		cap bayessummary, event(`event') total(`total') studyid(`studyid') estimates(metapreg_modest) ///
			`interaction' catreg(`catreg') contreg(`contreg') level(`level') model(`getmodel') cimethod(`cimethod')  ///
			varx(`varx') typevarx(`typevarx') by(`by') regexpression(`regexpression') ///
			`mcbnetwork' `comparative' `pcbnetwork' `abnetwork' `stratify'  ///
			comparator(`comparator') link(`link')
		
		mat `rawest' = r(loddsout)
		mat `exactabsout' = r(exactabsout)
		mat `absout' = r(absout)
	}
		
	//Population abs
	if "`getmodel'" != "hexact" {
		//simulations ABS
		postsim, event(`event') total(`total') orderid(`rid') studyid(`studyid') todo(p) estimates(metapreg_modest) rawest(`rawest') ///
				level(`level')  model(`getmodel')  by(`by') `comparative' `interaction' `abnetwork' ///
				`mcbnetwork' varx(`varx') cov(`cov') nsims(`nsims') link(`link') p(`p')  bayesreps(`bayesreps')
				
		mat `popabsout' = r(outmatrix)
	}

	//R
	local rrsuccess 0
	if "`catreg'" != "" | "`typevarx'" == "i" {
		if "`inference'" == "frequentist" {
			capture  freqestr, event(`event') total(`total') studyid(`studyid') estimates(metapreg_modest)  catreg(`catreg') ///
				level(`level') comparator(`comparator') `interaction' cimethod(`cimethod') ///
				varx(`varx') typevarx(`typevarx') by(`by') `mcbnetwork' `pcbnetwork' ///
				`comparative' `abnetwork' `stratify' model(`getmodel') ///
				regexpression(`regexpression') `baselevel' link(`link') inference(`inference')
		}
		else {
			capture  bayesestr, event(`event') catreg(`catreg') ///
				level(`level') comparator(`comparator') `interaction' cimethod(`cimethod') ///
				varx(`varx') typevarx(`typevarx') by(`by') `mcbnetwork' `pcbnetwork' ///
				`comparative' `abnetwork' `stratify' model(`getmodel') ///
				regexpression(`regexpression') `baselevel' link(`link') 
		}	
			
		if _rc == 0 {
			if "`getmodel'" != "hexact" {
				mat `rrout' = r(rroutmatrix)
				mat `orout' = r(oroutmatrix)
				
				if "`inference'" == "frequentist" {
					local inltest = r(inltest)
					if "`inltest'" == "yes" {
						mat `nltestRR' = r(nltestRR) //if RR by groups are equal
						mat `nltestOR' = r(nltestOR) //if RR by groups are equal
					}
				}
				//simulations
				postsim,  event(`event') total(`total') orderid(`rid') studyid(`studyid') todo(r) estimates(metapreg_modest) rrout(`rrout') orout(`orout') ///
						 level(`level')  model(`getmodel')  by(`by') `comparative' `interaction' `abnetwork' ///
						 `baselevel' `mcbnetwork' varx(`varx')  cov(`cov') nsims(`nsims') link(`link') p(`p') bayesreps(`bayesreps')
				
				mat `poprrout' = r(rroutmatrix)
				mat `poporout' = r(oroutmatrix)
				mat `poplorout' = r(loroutmatrix)
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
			level(`level')  model(`getmodel')  by(`by') `comparative'  ///
			modeles(`modeles') modellci(`modellci') modeluci(`modeluci') outplot(`outplot')	///
			`interaction' `abnetwork'  `mcbnetwork' varx(`varx')  cov(`cov') nsims(`nsims') link(`link')  bayesreps(`bayesreps')
	}
	//===================================================================================
	//Return the matrices
	if "`abnetwork'`cov'" != "" {
		if "`cov'" == "unstructured" {
			mat `hetout' = (`DF_BHET', `BHET' ,`P_BHET', `TAU21', `TAU22', `rho', `ISQ1', `ISQ2')
			mat colnames `hetout' = DF Chisq p tau2 sigma2 rho I2tau I2sigma 
		}
		else {
			mat `hetout' = (`DF_BHET', `BHET' ,`P_BHET', `TAU21', `TAU22', `ISQ1', `ISQ2')
			mat colnames `hetout' = DF Chisq p tau2 sigma2 I2tau I2sigma 
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
		return matrix  exactabsout = `exactabsout'
	}
	cap confirm matrix `rrout'
	if _rc == 0 {
		return matrix rrout = `rrout'
		return matrix poprrout = `poprrout'
	}
	cap confirm matrix `orout'
	if _rc == 0 {
		return matrix orout = `orout'
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
		return matrix nltestOR = `nltestOR'
	}
	
	return scalar mdf = `mdf'
	return local model = "`getmodel'"
	return local rrsuccess ="`rrsuccess'"
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
		cimethod(string) cov(string)
		level(integer 95)
		DP(integer 2)
		progress
		model(string) modelopts(string) 
		noMC noCONstant
		interaction	
		comparative mcbnetwork pcbnetwork abnetwork
		link(string) inference(string) refsampling(string)
			*];

	#delimit cr
	marksample touse, strok 
	
	tempvar event nonevent total 
	tempfile bayesnullest
	local bayesnullest = subinstr("`bayesnullest'", ".tmp", ".dta", 1)
		
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
	if "`link'" == "loglog" {
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
		
		`echo' fitmodel `outcome' `total' if `touse',  `modelopts' model(`model') ///
			regexpression(`nureduced') sid(`sid') level(`level')  nested(`first') `abnetwork' ///
			cov(`newcov') link(`link') p(`p') inference(`inference') ///
			refsampling(`refsampling') bayesest(`bayesnullest')
		
		if "`inference'" == "frequentist" {
			qui estat ic
			mat `matgof' = r(S)
			local BICmc = `matgof'[1, 6]
		}
		estimates store metapreg_Null
		
		if "`inference'" == "frequentist" {
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
			sid(`sid') level(`level')  nested(`first') `abnetwork' link(`link') p(`p')  ///
			inference(`inference') refsampling(`refsampling') bayesest(`bayesnullest')
		
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
		[model(string) regexpression(string) sid(varname) p(string) 
		level(integer 95) mcbnetwork pcbnetwork abnetwork general comparative nested(string) cov(string) link(string) 
		bayesest(string asis) bayesreps(string asis) inference(string) 
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
	
	tokenize `varlist'
	local events = "`1'"
	local total = "`2'"
	
	if "`cov'" != "" {
		local varx "`nested'"
	}
	else{
		local varx
	}

	if "`inference'" == "frequentist" {	
		if "`cov'" != "" {
			if "`cov'" == "unstructured" {
				local cov "cov(`cov')"
			}
			else {
				local cov
			}
		}
		else {
			local cov
		}
	
		if "`abnetwork'" != "" & "`model'" == "random" {
			local nested = `"|| (`nested': )"'
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
		else {
			local fitcommand "mecloglog"
		}
	}
	
	if "`inference'" == "bayesian" {
		if  "`bayesest'" != "" {
			local saving = `"saving(`bayesest', replace)"'
		}
		
		if "`varprior'" == "" {
			if "`cov'" == "unstructured" {
				local varprior = "iwishart(2,3,I(2))"
				*local varprior = "jeffreys(2)"
			}
			else {
				*local varprior = "jeffreys"
				local varprior = "normal(0, 20)"
			}
		}
		
		if "`feprior'" == "" {
			local feprior = "normal(0, 50)"
		}
		
		//Get initial values from binreg
		noi di as txt _n "Obtaining  initial values"
		capture binreg `events' `regexpression' if `touse', noconstant n(`total') ml  
		if _rc == 0 {
			tempname coefmat
			mat `coefmat' = e(b)
			local parmnames : colnames `coefmat'

			local nfeparm : word count `parmnames'

			local scales "0.5 1 1.5 2 0.25 1.75"
			
			forvalues c=1(1)`nchains' {
				local init`c'
				local scale : word `c' of `scales'
				forvalues t=1(1)`nfeparm' {
					local term : word `t' of `parmnames'
					local mle = `coefmat'[1, `t']
					if strpos("`model'", "random") != 0 {
						if strpos("`term'", "1b.`varx'") != 0 & "`cov'" =="independent" {
							continue
						}
						local init`c' = "`init`c'' {`events':`term'} `mle' "
					}
					else {
						local init`c' = "`init`c'' {`events':`term'} `=`scale'*`mle'' "
					}
				}
			}
		}
	}
	
	//fit bayesian fixed
	if "`model'" == "bayesfixed" {
		/*local nterms: word count `regexpression'
		if `nterms' > 1 {
			gettoken first variableterms : regexpression
		}

		if "`variableterms'" != "" {
			local priorfe = `" prior({`events':`regexpression'}, `feprior')"'
			local blockfe = `"block({`events':mu `variableterms'})"'
		}
		else {
			local priormu = `"prior({`events':mu}, `feprior')"'
			local blockmu = `"block({`events':mu})"'
		}*/
		//Initial values 
		local initialvalues 
		forvalues c=1(1)`nchains' {
			local initialvalues = "init`c'(`init`c'') `initialvalues'"
		}
				
		local priorfe = `" prior({`events':`regexpression'}, `feprior')"'
		local blockfe = `"block({`events':`regexpression'})"'
			
		#delim ;
		capture noisily bayesmh `events' `regexpression', noconstant likelihood(binomial(`total')) 
			`priormu' `priorfe' 
			`blockfe' `blockmu' 
			nchains(`nchains') thinning(`thinning') burnin(`burnin') mcmcsize(`mcmcsize') rseed(`rseed')
			/*`initialvalues'*/
			`saving'
			;
		#delimit cr 
		local success = _rc			
	}
	
	//fit bayesian re
	local lnvar = 0
	if "`model'" == "bayesrandom" {
		fvset base none `sid'
			
		gettoken mu variableterms: regexpression
		
		if strpos("`varprior'", "normal") != 0 {
			local expsigma2 = "sigmasq:exp({lnsigmasq})"
			local parmsigma2 = "{lnsigmasq}"
			
			local exptau2 = "tausq:exp({lntausq})"
			local parmtau2 = "{lntausq}"
			local lnvar = 1
			
			local tauscales "-4 1 -3 2 -6 4"
			local sigmascales "-6 -5 -4 -3 1 2 3"
		}
		else {
			local expsigma2 = `"{sigmasq}"'
			local parmsigma2 = `"{sigmasq}"'
			
			local exptau2 = `"{tausq}"'
			local parmtau2 = `"{tausq}"'
			
			local tauscales "0 2 0.5 1 2.5 0.25 3.5"
			local sigmascales "0 0.05 0.1 0.15 0.2 2 3"
		}
		
		if "`cov'" != "" {
			tokenize `variableterms'
			macro shift
			local variableterms `*'
						
			local neoregexpression = `"i.`sid' i.`sid'#2.`varx' `variableterms'"'
			
			//Assign re priors		
			if "`cov'" == "independent" {
				local priorslopes = `"prior({`events':i.`sid'#2.`varx'}, normal({`events':2.`varx'}, `expsigma2'))"'
				local priorsid = `"prior({`events':i.`sid'}, normal({`events':mu}, `exptau2'))"'	
				
				local priorsigmasq = `"prior(`parmsigma2', `varprior')"'
				local priortausq = `"prior(`parmtau2', `varprior')"'
				
				local initialvalues
				forvalues c=1(1)`nchains' {
					local tauinit : word `c' of `tauscales'
					local sigmainit : word `c' of `tauscales'
					local init`c' = `"`init`c'' `parmtau2' `tauinit' `parmsigma2' `sigmainit' "'
					local initialvalues = "init`c'(`init`c'') `initialvalues'"
				}
			}
			else if "`cov'" == "unstructured" {
				local priorre = `"prior({`events':i.`sid' i.`sid'#2.`varx'}, mvnormal(2, {`events':mu}, {`events':2.`varx'}, {Sigma, matrix}))"'
				local priorvarcov = `"prior({Sigma, matrix}, `varprior')"'
				
				local initialvalues
				forvalues c=1(1)`nchains' {
					local tauinit : word `c' of `tauscales'
					local sigmainit : word `c' of `tauscales'
					local init`c' = `"`init`c'' {Sigma, matrix} 0.1 0 0.05 "'
					local initialvalues = "init`c'(`init`c'') `initialvalues'"
				}
			}
			
			local priorvarx = `"prior({`events':2.`varx'}, `feprior')"'
			
			
			//block			
			local blockvarx = `"block({`events':2.`varx'})"'
			
			if "`cov'" == "independent" {
				if	`refsampling' == 1 {
					local blockslopes = `"block({`events':i.`sid'#2.`varx'}, split)"'
				}
				else if `refsampling' == 2  {
					local blockslopes = `"block({`events':i.`sid'#2.`varx'}, reffects)"'
				}
				if strpos("`varprior'", "gamma") != 0 { 
					local blocksigmasq = `"block(`parmsigma2', gibbs)"'
					local blocktausq = `"block(`parmtau2', gibbs)"'
				}
				else {
					local blocksigmasq = `"block(`parmsigma2')"'
					local blocktausq = `"block(`parmtau2')"'
				}
			}
			else if "`cov'" == "unstructured" {
				local blockvarcov = `"block({Sigma, matrix}, gibbs)"'
			}
		}
		else {
			local neoregexpression = `"i.`sid' `variableterms'"'
			
			//Assign re priors
			local priorsid = `"prior({`events':i.`sid'}, normal({`events':mu}, `exptau2'))"'			
			local priortausq = `"prior(`parmtau2', `varprior')"'
			
			//Block re
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
				local blocktausq = `"block(`parmtau2', gibbs)"'
			}
			else {
				local blocktausq = `"block(`parmtau2')"'
			}
			
			local initialvalues
			forvalues c=1(1)`nchains' {
				local tauinit : word `c' of `tauscales'
				local init`c' = `"`init`c'' `parmtau2' `tauinit' "'
				local initialvalues = "init`c'(`init`c'') `initialvalues'"
			}
		}
		
		//assign fe priors			
		if "`variableterms'" != "" {
			local priorfe = `"prior({`events':`variableterms'}, `feprior')"'
		}
		
		local priormu = `"prior({`events':mu}, `feprior')"'
		
		//block
		if "`cov'" == "unstructured" {
			local blockmu = `"block({`events':mu})"'
		}
		else {
			local blockmu = `"block({`events':mu}, gibbs)"'
		}
		
		if "`variableterms'" != "" {
			local blockfe = `"block({`events':`variableterms'})"'
		}
		
		//Initial values 
		local initialvalues 
		forvalues c=1(1)`nchains' {
			local initialvalues = "init`c'(`init`c'') `initialvalues'"
		}
		
		#delim ;
		capture noisily bayesmh `events' `neoregexpression', noconstant likelihood(binomial(`total'))  
			`priorre' `priorvarcov' `priorslopes' `priorsid'  `priormu' `priorvarx' `priortausq' `priorsigmasq' `priorfe'
			`blockfe' `blockvarx' `blockmu' `blockslopes' `blocksid' `blocksigmasq' `blocktausq'
			nchains(`nchains') thinning(`thinning') burnin(`burnin') mcmcsize(`mcmcsize') rseed(`rseed')
			/*`initialvalues'*/
			`saving'
			;
		#delimit cr 
		local success = _rc	
	}

	/*===========================Frequenstist models===========================================*/
	//Fit cbbetabin - common beta beta-binomial
	if "`model'" == "cbbetabin" {
		//Default iterations
		if strpos(`"`modelopts'"', "iterate") == 0  {
			local iterate = `"iterate(100)"'
		}
			
		qui xtset `sid'
		
		capture noisily xtnbreg `events' `regexpression' if `touse', fe `modelopts' `iterate' l(`level')	
		local success = _rc	
	}
	
	//Fit crbetabin - common rho beta-binomial
	if "`model'" == "crbetabin" {
		//Default iterations
		if strpos(`"`modelopts'"', "iterate") == 0  {
			local iterate = `"iterate(100)"'
		}
		
		capture noisily betabin `events' `regexpression' if `touse', noconstant n(`total') link(`link') `modelopts' `iterate' l(`level')	
		local success = _rc
	}
			
	//Fit the FE model
	if ("`model'" == "fixed") |("`model'" == "hexact"){
		if "`link'" == "logit" { 
			capture noisily binreg `events' `regexpression' if `touse', noconstant n(`total') ml `modelopts' l(`level')
		}
		else {
			capture noisily glm `events' `regexpression' if `touse', noconstant family(binomial `total') link(cloglog) ml `modelopts' l(`level')	
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
		#delim ;
		capture noisily `fitcommand' (`events' `regexpression' if `touse', noconstant)|| 
		  (`sid': `varx' , `cov') `nested' ,
		  binomial(`total') `ipoints' `modelopts' l(`level') `iterate';
		#delimit cr 
		
		local success = _rc
		local converged = e(converged)
		
		//Got to meqrlogit if melogit fails
		if  ("`fitcommand'" == "melogit") & (`success' != 0)  {
			local fitcommand = "meqrlogit"
			local iterate = `"iterate(30)"'
			
			local ++try
			#delim ;
			capture noisily  `fitcommand' (`events' `regexpression' if `touse', noconstant)|| 
			  (`sid': `varx' , `cov') `nested' ,
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
			capture noisily  `fitcommand' (`events' `regexpression' if `touse', noconstant)|| 
				(`sid': `varx' , `cov') `nested' ,
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
				capture noisily  `fitcommand' (`events' `regexpression' if `touse', noconstant)|| 
				  (`sid': `varx', `cov') `nested' ,
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
			capture noisily  `fitcommand' (`events' `regexpression' if `touse', noconstant )|| 
				(`sid': `varx' , `cov') `nested' ,
				binomial(`total') `ipoints'  l(`level') `refineopts' matlog `iterate';
			#delimit cr
			
			local success = _rc 
			
			local converged = e(converged)
		}
		
		*Try laplace if not for other commands
		if (`success' != 0) & ("`fitcommand'" != "meqrlogit") & (strpos(`"`modelopts'"', "laplac") == 0) {
			#delim ;
			capture noisily  `fitcommand' (`events' `regexpression' if `touse', noconstant )|| 
				(`sid': `varx' , `cov') `nested' ,
				binomial(`total') `modelopts' l(`level') intmethod(laplace) `iterate';
			#delimit cr
			
			local success = _rc 
			local converged = e(converged)		
		}
	}
	//Revert to FE if ME fails
	if (`success' != 0) & ("`model'" == "random") {
		
		if "`link'" == "logit" { 
			capture noisily binreg `events' `regexpression' if `touse', noconstant n(`total') ml `modelopts' l(`level')
		}
		else {
			capture noisily glm `events' `regexpression' if `touse', noconstant family(binomial `total') link(cloglog) ml `modelopts' l(`level')	
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
	/*if "`model'" == "hexact" {
		return matrix absexact = `absexact'
	}*/
end

	/*+++++++++++++++++++++++++	SUPPORTING FUNCTIONS: metadta_PROPCI +++++++++++++++++++++++++
								CI for proportions
	++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
	cap program drop metapreg_propci
	program define metapreg_propci

		syntax varlist [if] [in], p(name) se(name)lowerci(name) upperci(name) [cimethod(string) level(real 95)]
		
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

				cii proportions `N' `n', `cimethod' level(`level')
				
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

syntax varlist, rid(name) [assignment(name) event(name) total(name) idpair(name)  mcbnetwork pcbnetwork abnetwork general comparative panelize]

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

	syntax varlist, sid(varlist) idpair(varname) [sortby(varlist) jvar(varname) mcbnetwork pcbnetwork abnetwork general comparative]

		qui{
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
	syntax varlist, [exactorout(name) poprrout(name) rrout(name) poporout(name) orout(name) popabsout(name) exactabsout(name) absout(name) absoutp(name) sortby(varlist) 
		by(varname) hetout(name) model(string) prediction 
		groupvar(varname) se(varname) summaryonly nooverall nosubgroup outplot(string) grptotal(name) download(string asis) 
		indvars(varlist) depvars(varlist) dp(integer 2) stratify pcont(integer 0) level(integer 95)
		comparative abnetwork general pcbnetwork mcbnetwork enhance
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
	
	if "`outplot'" != "abs" & "`design'" == "abnetwork" {
		local summaryonly "summaryonly"
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
			replace `label' = "Group Mean" if `use' == 2 
			replace _WT = . if `use' == 2 
						
			qui label list `groupvar'
			local nlevels = r(max)
			forvalues l = 1/`nlevels' {
				if "`outplot'" == "abs" {
					if "`model'" == "hexact" {
						local S_1 = `exactabsout'[`l', 1]
						local S_3 = `exactabsout'[`l', 5]
						local S_4 = `exactabsout'[`l', 6]
					}
					else { 					
						local S_1 = `popabsout'[`l', 1]
						local S_3 = `popabsout'[`l', 4]
						local S_4 = `popabsout'[`l', 5]
						
						//Get Conditional CI
						local C_1 = `absout'[`=`pcont' +`l'', 1]
						local C_3 = `absout'[`=`pcont' +`l'', 5]
						local C_4 = `absout'[`=`pcont' +`l'', 6]
						
						//Get exact estimates
						local E_1 = `exactabsout'[`l', 1]
						local E_3 = `exactabsout'[`l', 5]
						local E_4 = `exactabsout'[`l', 6]
						if "`enhance'" != "" {
							//if simulated RE more than 5 times larger, replace with conditional stats 
							if (`=(`S_4' - `S_3')/(`C_4' - `C_3')' > 5) & (`C_4' == .) & (`C_3' == .)  {
								local S_3 = `C_3'
								local S_4 = `C_4'
							}
							
							//if simulated FE variance 5 times larger replace with exact
							if (`=(`S_4' - `S_3')/(`E_4' - `E_3')' > 5) & (`E_4' != .) & (`E_3' != .)  {
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
						replace `label' = "Group Mean" + " (Isq = " + string(`isq', "%10.`=`dp''f") + "%, p = " + string(`phet', "%10.`=`dp''f") + ")" if (`use' == 2) & (`groupvar' == `l') & (`grptotal' > 2)  & (`isq' != .)	
						
					}	 
				}
				else {
					if "`model'" == "hexact" {
						local S_1 = `exactorout'[`l', 1]
						local S_3 = `exactorout'[`l', 3]
						local S_4 = `exactorout'[`l', 4]
					}
					else {
						if "`outplot'" == "rr" {
							local S_1 = `poprrout'[`l', 1]
							local S_3 = `poprrout'[`l', 4]
							local S_4 = `poprrout'[`l', 5]
							
							//Conditional stats
							local C_1 = `rrout'[`l', 1]
							local C_3 = `rrout'[`l', 5]
							local C_4 = `rrout'[`l', 6]
						}
						if "`outplot'" == "or" {
							local S_1 = `poporout'[`l', 1]
							local S_3 = `poporout'[`l', 4]
							local S_4 = `poporout'[`l', 5]
							
							//Conditional stats
							local C_1 = `orout'[`l', 1]
							local C_3 = `orout'[`l', 5]
							local C_4 = `orout'[`l', 6]
						}
					}
					
					if "`enhance'" != "" {
					//if simulated more than 5 times larger, replace with conditional stats 
						if `=(`S_4' - `S_3')/(`C_4' - `C_3')' > 5 {
							local S_3 = `C_3'
							local S_4 = `C_4'
						}
					}
				}
				local lab:label `groupvar' `l'
				replace `label' = "`groupvar' = `lab'" if `use' == -2 & `groupvar' == `l' & (("`abnetwork'" == "") |("`outplot'" == "abs" & "`abnetwork'" != ""))	
				replace `label' = "`lab'" if `use' == 2 & `groupvar' == `l' & "`outplot'" != "abs" & "`abnetwork'" != ""		
				replace `es'  = `S_1' if `use' == 2 & `groupvar' == `l'	
				replace `lci' = `S_3' if `use' == 2 & `groupvar' == `l'	
				replace `uci' = `S_4' if `use' == 2 & `groupvar' == `l'	
				//Predictions
				if "`outplot'" == "abs" & "`prediction'" != "" {
					replace `lci' = `S_5' if `use' == 4 & `groupvar' == `l'	
					replace `uci' = `S_6' if `use' == 4 & `groupvar' == `l'	
				}
				//Weights
				sum _WT if `use' == 1 & `groupvar' == `l'
				local groupwt = r(sum)
				replace _WT = `groupwt' if `use' == 2 & `groupvar' == `l'	
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
				else {
					local nrows = rowsof(`popabsout')
					local S_1 = `popabsout'[`nrows', 1]
					local S_3 = `popabsout'[`nrows', 4]
					local S_4 = `popabsout'[`nrows', 5]
				}
				//predictions
				if "`prediction'" != "" {
					local nrows = rowsof(`absoutp')
					local S_5 = `absoutp'[`nrows', 1]
					local S_6 = `absoutp'[`nrows', 2]
				}			
			}
			else {
				if "`model'" == "hexact" {
					local nrows = rowsof(`exactorout')
					local S_1 = `exactorout'[`nrows', 1]
					local S_3 = `exactorout'[`nrows', 3]
					local S_4 = `exactorout'[`nrows', 4]
				}
				else {
					if "`outplot'" == "rr" {
						local nrows = rowsof(`poprrout')
						local S_1 = `poprrout'[`nrows', 1]
						local S_3 = `poprrout'[`nrows', 4]
						local S_4 = `poprrout'[`nrows', 5]
					}
					if "`outplot'" == "or" {
						local nrows = rowsof(`poporout')
						local S_1 = `poporout'[`nrows', 1]
						local S_3 = `poporout'[`nrows', 4]
						local S_4 = `poporout'[`nrows', 5]
					}
				}
			}
			replace `label' = "Population Mean" if `use' == 5
			
			if "`model'" == "random" & "`indvars'" == ""  & "`outplot'" == "abs" {
				local nrows = rowsof(`hetout')
				local isq = `hetout'[`nrows', 5]
				local phet = `hetout'[`nrows', 3]
				*replace `label' = "Overall (Isq = " + string(`isq', "%10.`=`dp''f") + "%, p = " + string(`phet', "%10.`=`dp''f") + ")" if `use' == 3
				replace `label' = "Population Mean (Isq = " + string(`isq', "%10.`=`dp''f") + "%, p = " + string(`phet', "%10.`=`dp''f") + ")" if `use' == 5 & (`isq' != .)	
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
				
		replace `label' = "Predictive t Interval" if `use' == 4 & "`model'" == "random"
		replace `label' = "t Interval" if `use' == 4 & "`model'" != "random"
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
		
		if "`abnetwork'" != "" & "`outplot'" != "abs" {
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
	dp(integer 2) power(integer 0) nowt smooth icimethod(string) ocimethod(string) groupvar(varname) 
	design(string) outplot(string) summaryonly]
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
	
	if "`outplot'" != "abs" & "`design'" == "abnetwork" {
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
					local citext "- `icimethod' (Wald) CI - "
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
		
		if "`outplot'" =="abs" & "`design'" == "comparative" {
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
						local mes "`modeles'[`i']*(10^`power')"
						local mlci "`modellci'[`i']*(10^`power')"
						local muci "`modeluci'[`i']*(10^`power')"
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
						local sumtext = "Group Mean"
					}
					else {
						local sumtext = "Population Mean"			
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
		
		syntax varlist, [interaction alphasort mcbnetwork pcbnetwork abnetwork general comparative ipair(varname) baselevel(string) studyid(varname) model(string)]
		
		tempvar holder
		tokenize `varlist'

		if "`mcbnetwork'`pcbnetwork'" == "" {
			macro shift 2
			local regressors "`*'"
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
		else if "`mcbnetwork'`pcbnetwork'" != "" {
			if "`interaction'" != "" {				
				local regexpression = "ibn.`ipair'#ibn.`comparator' i.`index'"
				//nulllify
				local interaction
			}
			else {
				local regexpression = "mu i.`ipair' i.`index'"	
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
			capture confirm numeric var ``i''
			if _rc != 0 {
				if "`alphasort'" != "" {
					sort ``i''
				}
				my_ncod `holder', oldvar(``i'')
				drop ``i''
				rename `holder' ``i''
				local prefix_`i' "i"
			}
			else {
				local prefix_`i' "c"
			}
			if "`baselevel'" != "" & `i'==1 {
				//Find the base level
				qui label list ``i''
				local nlevels = r(max)
				local found = 0
				local level 1
				while !`found' & `level' < `=`nlevels'+1' {
					local lab:label ``i'' `level'
					if "`lab'" == "`baselevel'" {
						local found = 1
						local basecode `level'
					}
					local ++level
				}
				
				if "`abnetwork'" != "" {
					local prefix_`i' "ibn"
				}
				if "`general'`comparative'" != "" {
					local prefix_`i' "ib`basecode'"
				}
			}
			/*Add the proper expression for regression*/
			local regexpression2 = "`regexpression2' `prefix_`i''.``i''#c.mu"
			local regexpression = "`regexpression' `prefix_`i''.``i''"
				
			if `i' > 1 & "`interaction'" != "" {
				local regexpression = "`regexpression' `prefix_`i''.``i''#`prefix_1'.`1'"
				local regexpression2 = "`regexpression2' `prefix_`i''.``i''#`prefix_1'.`1'#c.mu"				
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
	
	syntax, rawestmat(name)  [link(string) cimethod(string)]
	
	if "`link'" == "cloglog" {
		local invfn "invcloglog"
	}
	else if "`link'" == "loglog" {
		local invfn "1 - invcloglog"
	}
	else {
		local invfn "invlogit"
	}
	tempname matrixout
	mat `matrixout' = `rawestmat'
	
	local nrows = rowsof(`matrixout')
	local ncols = colsof(`matrixout')
			
	if `ncols' > 2 {	
		forvalues r = 1(1)`nrows' {
			mat `matrixout'[`r', 1] = `invfn'(`matrixout'[`r', 1]) //p
			mat `matrixout'[`r', 5] = `invfn'(`matrixout'[`r', 5]) //lower
			mat `matrixout'[`r', 6] = `invfn'(`matrixout'[`r', 6]) //upper
		}
		
		if "`cimethod'"== "t" {
			mat colnames `matrixout' = Mean SE(`link') t(`link') P>|t| Lower Upper
		}
		else{
			mat colnames `matrixout' = Mean SE(`link') z(`link') P>|z| Lower Upper
		}	
	}
	else {
		forvalues r = 1(1)`nrows' {
			mat `matrixout'[`r', 1] = `invfn'(`matrixout'[`r', 1])  //lower
			mat `matrixout'[`r', 2] = `invfn'(`matrixout'[`r', 2])  //upper
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
			comparator(varname) cimethod(string) mcbnetwork pcbnetwork abnetwork general comparative stratify interaction ///
			catreg(varlist) contreg(varlist) power(integer 0) level(integer 95) by(varname) link(string)]
		
		tempname absout loddsout exactabsout exactabsouti absexact 
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
		}
		else {
			if "`contreg'" == "" {
				local contreg = "`varx'"
			}
		}
		
		local invfn "invlogit"
						
		local ncatreg 0
		local parmlodds
		local parmp
		
		local catvars = "`catreg'"
		if "`catvars'" != "" {
			foreach c of local catvars {
				qui label list `c'
				local nlevels = r(max)
				forvalues l = 1/`nlevels' {
					if `l' == 1 {
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
		else {
			//mu
			local parmlodds = "(Overall:{`event':mu})"
			local parmp = "(Overall:`invfn'({`event':mu}))"
		}
		
		if "`parmlodds'" != "" | (  "`parmlodds'" == "" & ("`abnetwork'`mcbnetwork'`pcbnetwork'" == "" )) {
			bayesstats summary `parmlodds', clevel(`level') /*`cimethod'*/
			
			mat `loddsout' = r(summary)			
			local rnames :rownames `loddsout'	
			local ncatreg = rowsof(`loddsout')
			
			bayesstats summary `parmp', clevel(`level') /*`cimethod'*/
			mat `absout' = r(summary)	
		}
														
		local catrownames = ""
		if "`mcbnetwork'`pcbnetwork'" != "" {
			local rownamesmaxlen : strlen local Index
			local rownamesmaxlen = max(`rownamesmaxlen', 10)
		}
		else {
			local rownamesmaxlen = 10 /*Default*/
		}
		
		//Nice labels
		forvalues r = 1(1)`ncatreg' {
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
			else {
				local rownames = "`rownames' `rname`r''" 
			}
		}

		mat rownames `loddsout' = `rownames'
		mat rownames `absout' = `rownames'
										
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
				if (strpos("`vari'", "_") == 1) & ("`group'" != "Overall"){
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
			comparator(varname) cimethod(string) mcbnetwork pcbnetwork abnetwork general comparative stratify interaction ///
			catreg(varlist) contreg(varlist) power(integer 0) level(integer 95) by(varname) link(string)]
		
		tempname coefmat outmatrix outmatrixp matrixout bycatregmatrixout catregmatrixout contregmatrixout row ///
		outmatrixr overall Vmatrix byVmatrix exactabsouti exactabsout absexact
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
		else {
			local expression "exp(predict(xb))"
		}
		
		if "`link'" == "cloglog" {
			local invfn "invcloglog"
		}
		else if "`link'" == "loglog" {
			local invfn "1 - invcloglog"
		}
		else {
			local invfn "invlogit"
		}
		
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
		
		qui estimates restore `estimates'
		local df = e(N) -  e(k) 
		mat `coefmat' = e(b)
		local predcmd = e(predict)
		
		if "`model'" == "random" {
			local npar = colsof(`coefmat')
			if "`predcmd'" == "meqrlogit_p" {
				local scalefn "exp"
				local scalepow "2"
			}
			else {
				local scalepow "1"
			}
			
			if "`abnetwork'`cov'" == "" {
				local TAU21 = `scalefn'(`coefmat'[1, `npar'])^`scalepow' //Between study variance	1
				local TAU22 = 0
			}
			else if "`abnetwork'" != "" {
				local TAU21 = `scalefn'(`coefmat'[1, `=`npar'-1'])^`scalepow' //Between study variance	1
				local TAU22 = `scalefn'(`coefmat'[1, `npar'])^`scalepow' //Between study variance	2
			}
			else if "`cov'" != "" {
				local TAU21 = `rosevar'[1, 1]
				local TAU22 = `rosevar'[2, 2]
				if "`cov'" == "unstructured" {
					local rho = tanh(`rawvar'[1, 2])
				}  
			}
		}
		else {
			local TAU21 = 0
			local TAU22 = 0
		}
		
		local byncatreg 0
		if ("`by'" != "" & "`stratify'"  == "")  {
			qui margin , `expression' `at' over(`by') level(`level')
			
			
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
		if "`marginlist'" != "" | (  "`marginlist'" == "" & ("`abnetwork'`mcbnetwork'`pcbnetwork'" == "" )) {
			qui margin `marginlist', `expression' `at' `grand' level(`level')
						
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
				
		mat `outmatrixp' = J(`=`byncatreg' + `ncatreg'', 2, .)
		forvalues r = 1(1)`byncatreg' {
			mat `outmatrixp'[`r', 1] = (`bycatregmatrixout'[`r',1] - invttail((`df'), 0.5-`level'/200) * sqrt(`bycatregmatrixout'[`r',2]^2 + `TAU21'^2  + `TAU22'^2))
			mat `outmatrixp'[`r', 2] = (`bycatregmatrixout'[`r',1] + invttail((`df'), 0.5-`level'/200)* sqrt(`bycatregmatrixout'[`r',2]^2 + `TAU21'^2 + `TAU22'^2))
		}
		forvalues r = `=`byncatreg' + 1'(1)`=`byncatreg' + `ncatreg''{
			mat `outmatrixp'[`r', 1] = (`catregmatrixout'[`=`r' - `byncatreg'', 1] - invttail((`df'), 0.5-`level'/200) * sqrt(`catregmatrixout'[`=`r' - `byncatreg'', 2]^2 + `TAU21'^2  + `TAU22'^2))
			mat `outmatrixp'[`r', 2] = (`catregmatrixout'[`=`r' - `byncatreg'', 1] + invttail((`df'), 0.5-`level'/200)* sqrt(`catregmatrixout'[`=`r' - `byncatreg'', 2]^2 + `TAU21'^2  + `TAU22'^2))
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
		
		//t distribution
		if "`cimethod'" == "t" {
			forvalues r = 1(1)`=`byncatreg' + `ncatreg' + `ncontreg''  {
					local tstat = `matrixout'[`r', 3]
					mat `matrixout'[`r', 4] = ttail(`df', abs(`tstat'))*2
					mat `matrixout'[`r', 5] = `matrixout'[`r', 1] - invttail((`df'), 0.5-`level'/200) * `matrixout'[`r', 2]
					mat `matrixout'[`r', 6] = `matrixout'[`r', 1] + invttail((`df'), 0.5-`level'/200) * `matrixout'[`r', 2]
			}
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
		mat colnames `outmatrixp' = Lower Upper
		mat rownames `outmatrixp' = `catrownames' `Overall'
					
		if "`cimethod'" == "t" { 
			mat colnames `matrixout' = Mean SE t P>|t| Lower Upper
		}
		else {
			mat colnames `matrixout' = Mean SE z P>|z| Lower Upper
		}

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
		return matrix outmatrixp = `outmatrixp'	
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
	modeles(varname) modellci(varname) modeluci(varname) outplot(string) baselevel(integer 1) cov(string)
	model(string) comparative  by(varname) level(real 95) interaction abnetwork mcbnetwork  varx(varname) nsims(string) p(integer 0)]
	;
	#delimit cr 
	
	marksample touse, strok
	
	tempname betacoef rawcoef varrawcoef ///
				fullrawcoef fullvarrawcoef X beta sims ///
				popabsout popabsouti poprrout poprrouti ///
				poporout poporouti poplorout poplorouti simvar absexact
	
	tempvar feff sfeff reff sreff reff1 sreff1 reff2 sreff2 eta insample ///
			newobs idpair gid rid hold holdleft holdright ///
			simmu sumphat meanphat subset subsetid subsetid1 sumphat1 ///
			meanphat1 gid1 modelp modelrr modelor modellor modelse sumrrhat ///
			meanrrhat sumorhat meanorhat sumlorhat meanlorhat varint varslope fisherrho ///
			simvarint simvarslope simfisherrho simrho covar simcovar lnsigma simlnsigma
			
	if "`link'" == "cloglog" {
		local invfn "invcloglog"
	}
	else if "`link'" == "loglog" {
		local invfn "1-invcloglog"
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
		gen `insample' = e(sample) * mu
		
		if strpos("`model'", "bayes") == 0 {
			//Restore estimates
			estimates restore `estimates'
			local predcmd = e(predict)
			
			//Coefficients estimates and varcov
			mat `fullrawcoef' = e(b)
			mat `fullvarrawcoef' = e(V)		
			
			local ncoef = colsof(`fullrawcoef')
			local rho = 0
			if "`model'" == "random" {
				if "`abnetwork'`cov'" == "" {
					local nfeff = `=`ncoef' - 1'
					local varnames "`varint'"
					local simvarnames "`simvarint'"
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
			
			//log-odds
			if "`model'" == "cbbetabin" {
				predictnl `feff' = xb() - _b[_cons], se(`sfeff')
			}
			else {
				//depends on the link used
				predict `feff', xb  //FE
				predict `sfeff', stdp //se of FE
			}
			
			if "`model'" == "random" {
				if "`abnetwork'`cov'" == "" {
					predict `reff', reffects reses(`sreff')
					gen `eta' = `feff' + `reff' //linear predictor
					gen `modelse' = sqrt(`sreff'^2 + `sfeff'^2) if `insample'==1
				}
				else if "`abnetwork'" !="" | "`cov'" !="" {
					predict `reff1' `reff2', reffects reses(`sreff1' `sreff2')  //slope=1  int=2 
					
					if "`abnetwork'"  == "" {
						gen `reff' = `reff1'*`varx' + `reff2' 
						gen `sreff' = sqrt(`varx'^2 * `sreff1'^2 + `sreff2'^2)	
					}
					else {
						gen `reff' = `reff1' + `reff2' 
						gen `sreff' = sqrt(`sreff1'^2 + `sreff2'^2)	
					}
								
					gen `eta' = `feff' + `reff' //linear predictor				
					gen `modelse' = sqrt(`sreff1'^2 + `sreff2'^2 + `sfeff'^2) if `insample'==1
				}
			}
			else {
				gen `eta' = `feff' 
				gen `modelse' = `sfeff' if `insample'==1
			}
			
			//Revert to original data if filler data was generated
			if (("`model'" == "random") & (`nobs' < 7))  {
				keep if !`newobs'
			}
			
			//Smooth p estimates
			gen `modelp' = `invfn'(`eta') if `insample'==1
			
			//identifiers
			sort `insample' `orderid'
			*egen `rid' = seq() if `insample'==1  //rowid
			
			if "`comparative'`mcbnetwork'`abnetwork'" != ""  {
				egen `gid' = group(`studyid' `by') if `insample'==1  
				sort `gid' `orderid' `varx'
				by `gid': egen `idpair' = seq()
				egen `rid' = seq() if `insample'==1  //rowid
				
				if "`abnetwork'" == "" {
					gen `modelrr' = `modelp'[_n] / `modelp'[_n-1] if (`gid'[_n]==`gid'[_n-1]) & (`idpair' == 2)
					gen `modelor' = (`modelp'[_n]/(1 - `modelp'[_n])) / (`modelp'[_n-1]/(1 - `modelp'[_n-1])) if (`gid'[_n]==`gid'[_n-1]) & (`idpair' == 2)
					gen `modellor' = ln(`modelor')
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
						replace `sreff1' = exp(`simvarslope') //Marginal se
						replace `sreff2' = sqrt((1 - (`simrho')^2)*(exp(`simvarint')^2)) //Conditional se
					}
					else {
						replace `sreff1' = sqrt(`simvarslope') //Marginal se
						replace `sreff2' = sqrt((1 - (`simrho')^2)*(`simvarint')) //Conditional se
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
					
				if "`comparative'`mcbnetwork'`abnetwork'" == "" {
					tempvar phat`j' 
				
					if "`model'" == "random" {
						tempvar restudy`j'
						//EB re 
						sum `reff' if `rid' == `j' 
						local reff_`j' = r(mean)
						
						if "`predcmd'" == "meqrlogit_p" {
							qui gen `restudy`j'' = rnormal(0, exp(`simvarint'))
						}
						else {
							qui gen `restudy`j'' = rnormal(0, sqrt(`simvarint'))
						}
						
						gen `phat`j'' = `invfn'(`reff_`j'' + `restudy`j'' + `festudy`j'')
					}
					else {
						gen `phat`j'' = `invfn'(`festudy`j'')
					}
				}
				if "`comparative'" != "" | "`mcbnetwork'" != "" | "`abnetwork'" != "" {
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
								replace `reff2' = rnormal(`simrho'*`sreff2'*(`reff1'/`sreff1'), `sreff2')
								
								gen `restudy`index'' = `reff1' + `reff2'
							}
						}					
						gen `phat`j'' = `invfn'(`reff_`j'' + `restudy`index'' +  `festudy`j'')
					}
					else {
						gen `phat`j'' = `invfn'(`festudy`j'')
					}
					
					
					if "`abnetwork'" == "" {
						//Create the pairs
						gen `phat_`pair'`index'' = `phat`j''
						
						if `pair' == 2 {
							tempvar rrhat`index' orhat`index' lorhat`index'
							gen `rrhat`index'' = `phat_2`index'' / `phat_1`index''
							gen `orhat`index'' = (`phat_2`index''/(1 - `phat_2`index'')) / (`phat_1`index'' /(1 - `phat_1`index''))
							gen `lorhat`index'' = ln(`orhat`index'')
						}
					}
				}
			}
		}	
		else {
			//identifiers
			sort `insample' `orderid'
			egen `rid' = seq() if `insample'==1  //rowid
			
			if "`comparative'`mcbnetwork'`abnetwork'" != ""  {
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
			merge 1:1 _n  using `bayesreps', nogenerate
		
			//# of obs
			count if `insample' == 1
			local nobs = r(N)
			
			//Generate the p's and r's
			forvalues j=1(1)`nobs' { 
				tempvar phat`j'  
			
				//toal 
				sum `total' if `rid' == `j' 
				local total_`j' = r(mean)
				
				gen `phat`j'' = _ysim1_`j'/`total_`j''
					
				if "`comparative'" != "" | "`mcbnetwork'" != "" | "`abnetwork'" != "" {
					sum `gid' if `rid' == `j'
					local index = r(min)
					
					sum `idpair' if `rid' == `j'
					local pair = r(min)
					
					//Generate the variables
					tempvar phat_`pair'`index'
					
					if "`abnetwork'" == "" {
						//Create the pairs
						gen `phat_`pair'`index'' = `phat`j''
						
						if `pair' == 2 {
							tempvar rrhat`index' orhat`index' lorhat`index'
							gen `rrhat`index'' = `phat_2`index'' / `phat_1`index''
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
			
			//Add overall
			if strpos("`model'", "bayes") == 1 {
				local eqnames = "`eqnames' _"
				local rnames = "`rnames' Overall"
			}
	
			foreach vari of local eqnames {		
				local ++mindex
				local group : word `mindex' of `rnames'
				
				//Skip if continous variable
				if (strpos("`vari'", "_") == 1) & ("`group'" != "Overall"){
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
				
				forvalues j=1(1)`nsubset' { 
				
					sum `rid' if `subsetid' == `j'
					local index = r(min)
					
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
				
				gen `meanphat' = `sumphat'/`nsubset'
								
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
					tempname poprrouti`baselevel' poporouti`baselevel' poplorouti`baselevel'
					
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
					mat `poporouti`baselevel'' = (1, 0, 1, 1, 1, .)
					mat `poplorouti`baselevel'' = (1, 0, 1, 1, 1, .)
					
					local baselab = ustrregexra("`baselab'", " ", "_")
					mat rownames `poprrouti`baselevel'' = `vari':`baselab'
					mat rownames `poporouti`baselevel'' = `vari':`baselab'
					mat rownames `poplorouti`baselevel'' = `vari':`baselab'
					
					//Other groups
					forvalues g=1(1)`ngroups' {
						if `g' != `baselevel' {
							tempvar meanphat`g' meanrrhat`g' meanorhat`g' meanlorhat`g' gid`g' sumphat`g' subsetid`g'
							tempname poprrouti`g' poporouti`g' poplorouti`g'
							
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
							gen `meanorhat`g'' = (`meanphat`g''/(1 - `meanphat`g'')) / (`meanphat`baselevel''/(1 - `meanphat`baselevel''))
							gen `meanlorhat`g'' = ln(`meanorhat`g'')

							//Obtain mean of modelled estimates
							if strpos("`model'", "bayes") == 0 {
								sum `modelp' if `vari' == `g' & `insample' == 1
								local meanmodelp`g' = r(mean)
								
								local meanmodelrr`g' = `meanmodelp`g'' / `meanmodelp`baselevel''
								local meanmodelor`g' = (`meanmodelp`g''/(1 - `meanmodelp`g'')) / (`meanmodelp`baselevel''/(1 - `meanmodelp`baselevel''))
								local meanmodellor`g' = ln((`meanmodelp`g''/(1 - `meanmodelp`g'')) / (`meanmodelp`baselevel''/(1 - `meanmodelp`baselevel'')))
							}
							
							//Standard error
							sum `meanrrhat`g''
							local postserr = r(sd)
							local postmeanrr = r(mean)
							
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
								mat `poporouti`g'' = (`meanmodelor`g'', `postseor', `medianor', `lowerpor', `upperpor', `norreps')
								mat `poplorouti`g'' = (`meanmodellor`g'', `postselor', `medianlor', `lowerplor', `upperplor', `nlorreps')
							}
							else {
								mat `poprrouti`g'' = (`postmeanrr', `postserr', `medianrr', `lowerprr', `upperprr',  `nrrreps')
								mat `poporouti`g'' = (`postmeanor', `postseor', `medianor', `lowerpor', `upperpor', `norreps')
								mat `poplorouti`g'' = (`postmeanlor', `postselor', `medianlor', `lowerplor', `upperplor', `nlorreps')

							}
							
							local glab = ustrregexra("`glab'", " ", "_")
							mat rownames `poprrouti`g'' = `vari':`glab'
							mat rownames `poporouti`g'' = `vari':`glab'
							mat rownames `poplorouti`g'' = `vari':`glab'
						}
						if `g' == 1 {
							mat `poprrouti' = `poprrouti`g''
							mat `poporouti' = `poporouti`g''
							mat `poplorouti' = `poplorouti`g''
						}
						else {
							//Stack the matrices
							mat `poprrouti' = `poprrouti'	\  `poprrouti`g''
							mat `poporouti' = `poporouti'	\  `poporouti`g''
							mat `poplorouti' = `poplorouti'	\  `poplorouti`g''
						}
					}
					//Stack the matrices
					local ++newnrows
					if `newnrows' == 1 {
						mat `poprrout' = `poprrouti'
						mat `poporout' = `poporouti'
						mat `poplorout' = `poplorouti'
					}
					else {
						mat `poprrout' = `poprrout'	\  `poprrouti'
						mat `poporout' = `poporout'	\  `poporouti'
						mat `poplorout' = `poplorout'	\  `poplorouti'
					}
				}
			}
			
			if "`comparative'" != "" | "`mcbnetwork'" != "" {
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
					local orlistvar
					local lorlistvar
		
					forvalues j=1(1)`nobs' { 
						sum `gid' if `subsetid' == `j'
						local index = r(min)
						
						sum `idpair' if `subsetid' == `j'
						local pair = r(min)
						
						if `pair' == 2 {
							local rrlistvar = `"`rrlistvar' `rrhat`index''"'
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
						
						sum `modelor' if `subset' == 1
						local meanmodelor = r(mean)
						
						sum `modellor' if `subset' == 1
						local meanmodellor = r(mean)
					}
					
					cap drop `meanrrhat' `meanorhat' `meanlorhat'
					/*
					gen `meanrrhat' = `sumrrhat'/`nsubset'
					gen `meanorhat' = `sumorhat'/`nsubset'
					gen `meanlorhat' = `sumlorhat'/`nsubset'
					*/
					
					egen `meanrrhat' = rowmean(`rrlistvar')
					egen `meanorhat' = rowmean(`orlistvar')
					egen `meanlorhat' = rowmean(`lorlistvar')

					
					//Standard error
					sum `meanrrhat'	
					local postserr = r(sd)
					local postmeanrr = r(mean)
					
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
						mat `poporouti' = (`meanmodelor', `postseor', `medianor', `lowerpor', `upperpor', `norreps')
						mat `poplorouti' = (`meanmodellor', `postselor', `medianlor', `lowerplor', `upperplor', `nlorreps')
					}
					else {
						mat `poprrouti' = (`postmeanrr', `postserr', `medianrr', `lowerprr', `upperprr', `nrrreps')
						mat `poporouti' = (`postmeanor', `postseor', `medianor', `lowerpor', `upperpor', `norreps')
						mat `poplorouti' = (`postmeanlor', `postselor', `medianlor', `lowerplor', `upperplor', `nlorreps')
					}
					
					mat rownames `poprrouti' = `vari':`group'
					mat rownames `poporouti' = `vari':`group'
					mat rownames `poplorouti' = `vari':`group'
					
					//Stack the matrices
					local ++newnrows
					if `newnrows' == 1 {
						mat `poprrout' = `poprrouti'
						mat `poporout' = `poporouti'
						mat `poplorout' = `poplorouti'
					}
					else {
						mat `poprrout' = `poprrout'	\  `poprrouti'
						mat `poporout' = `poporout'	\  `poporouti'
						mat `poplorout' = `poplorout'	\  `poplorouti'
					}
				}
			}
		}
		
		if "`todo'" == "smooth" {
			if strpos("`model'", "bayes") == 1 {
				//Smooth estimates		
				if "`outplot'" == "abs" {
					//# of obs
					count if `insample' == 1
					local nobs = r(N)
				}
				else {
					count if `insample' == 1 & `idpair' == 2
					local nobs = r(N)
					
					egen `subsetid' = seq()	if `insample' == 1 & `idpair' == 2
				}
								
				forvalues j=1(1)`nobs' {			
					if "`outplot'" == "abs" {
						centile `phat`j'', centile(50 `=(100-`level')/2' `=100 - (100-`level')/2')
					}
					else {
						sum `gid' if `subsetid' == `j'
						local index = r(min)
							
						if "`outplot'" == "rr" {
							centile `rrhat`index'', centile(50 `=(100-`level')/2' `=100 - (100-`level')/2')
						}
						else {
							centile `orhat`index'', centile(50 `=(100-`level')/2' `=100 - (100-`level')/2')
						}
					}
					
					local median = r(c_1) //Median
					local lowerp = r(c_2) //Lower centile
					local upperp = r(c_3) //Upper centile
					
					if "`outplot'" == "abs" {		
						replace `modeles' = `median' if `insample' == 1 & `rid' == `j'
						replace `modellci' = `lowerp' if `insample' == 1 & `rid' == `j'
						replace `modeluci' = `upperp' if `insample' == 1 & `rid' == `j'
					}
					else {
						replace `modeles' = `median' if (`gid' == `index') & (`idpair' == 2) &  `insample' == 1
						replace `modellci' = `lowerp' if (`gid' == `index') & (`idpair' == 2) &  `insample' == 1
						replace `modeluci' = `upperp' if (`gid' == `index') & (`idpair' == 2) &  `insample' == 1
					}
				}
			}
			else {			
				replace `modeles' = `modelp' if  `insample' == 1
					
				//Smooth p's
				if "`outplot'" == "abs" {
					//postci
					local critvalue -invnorm((100-`level')/200)
					if "`link'" == "loglog" {
						local sign -
					}
					replace `modellci' = `invfn'(`eta' - `sign' `critvalue'*`modelse') if  `insample' == 1 //lower
					replace `modeluci' = `invfn'(`eta' +  `sign' `critvalue'*`modelse') if  `insample' == 1 //upper
				
				}
							
				//Smooth r's
				if "`outplot'" == "rr" | "`outplot'" == "or" {
					if "`outplot'" == "rr" { 
						replace `modeles' = `modelrr' if  `insample' == 1
					}
					else {
						replace `modeles' = `modelor' if  `insample' == 1
					}
					
					*sum `gid' if `insample' == 1 
					count if `insample' == 1 & `idpair' == 2
					local nstudies = r(N)
					
					egen `subsetid' = seq()	if `insample' == 1 & `idpair' == 2
				
					forvalues j=1(1)`nstudies' {
						
						sum `gid' if `subsetid' == `j'
						local index = r(min)
						
						//Obtain the quantiles
						if "`outplot'" == "rr" {
							centile `rrhat`index'', centile(50 `=(100-`level')/2' `=100 - (100-`level')/2')
						}
						else {
							centile `orhat`index'', centile(50 `=(100-`level')/2' `=100 - (100-`level')/2')
						}
						local median = r(c_1) //Median
						local lowerp = r(c_2) //Lower centile
						local upperp = r(c_3) //Upper centile
						
						replace `modellci' = `lowerp' if (`gid' == `index') & (`idpair' == 2) &  `insample' == 1
						replace `modeluci' = `upperp' if (`gid' == `index') & (`idpair' == 2) &  `insample' == 1
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
		*mat colnames `popabsout' = Mean SE Median Lower Upper Total Events Studies Ones Zeros
		mat colnames `popabsout' = Mean SE Median Lower Upper Sample_size
		return matrix outmatrix = `popabsout'
	}
	if "`todo'" == "r" {
		mat colnames `poprrout' = Mean SE Median Lower Upper Sample_size
		mat colnames `poporout' = Mean SE Median Lower Upper Sample_size
		mat colnames `poplorout' = Mean SE Median Lower Upper Sample_size
		
		return matrix rroutmatrix = `poprrout'
		return matrix oroutmatrix = `poporout'
		return matrix loroutmatrix = `poplorout'
	}
end

/*+++++++++++++++++++++++++	SUPPORTING FUNCTIONS: PRINTMAT +++++++++++++++++++++++++
							Print the outplot matrix beautifully
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
cap program drop printmat
program define printmat
	#delimit ;
	syntax, matrixout(name) type(string) [sumstat(string) dp(integer 2) power(integer 0) stratify 
		mcbnetwork pcbnetwork abnetwork general comparative continuous p(integer 0) model(string) 
		nsims(string) link(string) inference(string)]
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
		if "`type'" == "rre" | "`type'" == "ore" {
			di as res _n "****************************************************************************************"
			if "`inference'" == "frequentist" {
				di as txt _n "Wald-type test for nonlinear hypothesis"
				if "`type'" == "rre" { 
					di as txt _n "{phang}H0: All (log)RR equal vs. H1: Some (log)RR different {p_end}"
				}
				else {
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
		if "`type'" == "popabs" | "`type'" == "poprr" | "`type'" == "popor" {
			local patho 0
			if "`type'" == "popabs" {
				local parm "Proportion"
			}
			else if "`type'" == "poprr" {
				local parm "Proportion Ratio"
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
			
			*if `patho' {
			#delimit ;
			noi matlist `mat2print', rowtitle(Parameter) 
						cspec(& %`rownamesmaxlen's |  %8.`=`dp''f &  %8.`=`dp''f & %8.`=`dp''f & %8.`=`dp''f & %8.`=`dp''f & %15.0f o2&) 
						rspec(`rspec') underscore nodotz
			;
			#delimit cr	
			/*}
			else {
				mat `mat2print' = `mat2print'[1..., 1..5]
				#delimit ;
				noi matlist `mat2print', rowtitle(Parameter) 
							cspec(& %`rownamesmaxlen's |  %8.`=`dp''f &  %8.`=`dp''f & %8.`=`dp''f & %8.`=`dp''f & %8.`=`dp''f o2&) 
							rspec(`rspec') underscore nodotz
				;
				#delimit cr	
			}*/			
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
		if ("`type'" == "raw") | ("`type'" == "abs") | ("`type'" == "exactabs") | ("`type'" == "rr")| ("`type'" == "or")   {
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
			if strpos("`type'", "abs") != 0 { 
				di as res "{pmore2} `typeinf' summary: Proportion {p_end}"
			}
			if ("`type'" == "rr") {
				di as res "{pmore2} `typeinf' summary: Proportion Ratio {p_end}"
			}
			if ("`type'" == "or") {
				di as res "{pmore2} `typeinf' summary: Odds Ratio {p_end}"
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
				di as txt "NOTE: H0: tau = 0 vs. H1: tau > 0"
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
			if ("`type'" == "abs") {
				di as txt "NOTE: H0: Est = 0.5 vs. H1: Est != 0.5"
			}
			if ("`type'" == "rr") {
				di as txt "NOTE: H0: Est = 1 vs. H1: Est != 1"
			}
			if ("`type'" == "logit") {
				di as txt "NOTE: H0: Est = 0 vs. H1: Est != 0"
			}
		}
		
		if ("`type'" == "popabs") | ("`type'" == "poprr") | ("`type'" == "popor")  {
			/*if `patho' {
				di as txt "NOTE: `level'% centiles obtained from Nreps successful replications out of `nsims' simulations of the posterior distribution"
				di as text "NOTE: Estimates not reliable."
			}
			else{*/
				di as txt "NOTE: `level'% centiles obtained from simulations of the posterior distribution"
			*}
			
		}		
end	
/*+++++++++++++++++++++++++	SUPPORTING FUNCTIONS: BAYESESTRCORE +++++++++++++++++++++++++
							Obtain the RR after modelling
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
cap program drop bayesestrcore
program define bayesestrcore, rclass

syntax, event(varname) [confounders(varlist) varx(varname) cimethod(string) ///
			level(integer 95)  baselevel(integer 1) link(string) model(string)] 
		
	tempname RRoutmatrix ORoutmatrix
	
	//Expression for logodds prediction
	if "`link'" == "cloglog" {
		local expression "exp(logit(invcloglog(predict(xb))))"
	}
	else if "`link'" == "loglog" {
		local expression "exp(-logit(invcloglog(predict(xb))))"
	}
	else {
		local expression "exp(predict(xb))"
	}
	if "`model'" == "cbbetabin" {
		local expression "expression(xb() - _b[_cons]) at(mu==1)"
	}
	
	local invfn "invlogit"
	
	local EstRRexpression //RR
	local EstORexpression //OR
	foreach c of local confounders {	
		qui label list `c'
		local nlevels = r(max)
		local test_`c'
		
		if "`varx'" != "" {
			forvalues l = 1/`nlevels' {
			
				if `l' == 1 {
					local EstRRexpression = "`EstRRexpression' (`c'_`l':`invfn'({`event':mu} + {`event':2.`varx'}) / `invfn'({`event':mu}))"
					local EstORexpression = "`EstRRexpression' (`c'_`l':exp({`event':2.`varx'}))"
				}
				else {
					if "`interaction'" != "" {
						local xterm = "+ {`event':2.`varx'#`l'.`c'}"
					}
					local EstRRexpression = "`EstRRexpression' (`c'_`l':`invfn'({`event':mu} + {`event':2.`varx'} + {`event':`l'.`c'} `xterm' ) / `invfn'({`events':mu} + {`event':`l'.`c'}))" 
					local EstORexpression = "`EstRRexpression' (`c'_`l':exp({`event':2.`varx'} `xterm'))"
				}
			}
		}
		else {			
			if `baselevel' == 1 {
				local basep = "`invfn'({`event':mu})"
				local baseodds = "exp({`event':mu})"
			}
			else {
				local basep = "`invfn'({`event':`baselevel'.`c'} + {`event':mu})"
				local baseodds = "exp({`event':`baselevel'.`c'} + {`event':mu})"
			}
			
			forvalues l = 1/`nlevels' {
				if `l' != `baselevel' {
					local EstRRexpression = "`EstRRexpression' (`c'_`l':`invfn'({`event':`l'.`c'} + {`event':mu})/`basep')"
					local EstORexpression = "`EstORexpression' (`c'_`l':exp({`event':`l'.`c'} + {`event':mu})/`baseodds')"
				}
				else {
					local EstRRexpression = "`EstRRexpression' (`c'_`l':`basep'/`basep')"
					local EstORexpression = "`EstORexpression' (`c'_`l':`baseodds'/`baseodds')"
				}	
			}
		}
	}
	
	//RR
	bayesstats summary `EstRRexpression', clevel(`level')  /*`cimethod'*/
	mat `RRoutmatrix' = r(summary)	
	
	//OR
	bayesstats summary `EstORexpression', clevel(`level')  /*`cimethod'*/
	mat `ORoutmatrix' = r(summary)	
	
	//Nice labels
	local rnames :rownames `RRoutmatrix'	
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
	mat rownames `RRoutmatrix' = `rownames'
	mat rownames `ORoutmatrix' = `rownames'
	
	return matrix rroutmatrix = `RRoutmatrix'
	return matrix oroutmatrix = `ORoutmatrix'

end	

/*+++++++++++++++++++++++++	SUPPORTING FUNCTIONS: BAYESESTR +++++++++++++++++++++++++
							Estimate RR after modelling
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
	cap program drop bayesestr
	program define bayesestr, rclass
		syntax, event(varname)[ regexpression(string) catreg(varlist) typevarx(string) varx(varname) comparator(varname) cimethod(string) ///
			level(integer 95)  mcbnetwork pcbnetwork abnetwork general comparative stratify power(integer 0) by(varname) ///
			baselevel(integer 1)  interaction link(string) model(string) ]
		
		//Expression for logodds prediction		
		if "`link'" == "cloglog" {
			local expression "exp(logit(invcloglog(predict(xb))))"
		}
		else if "`link'" == "loglog" {
			local expression "exp(-logit(invcloglog(predict(xb))))"
		}
		else {
			local expression "exp(predict(xb))"
		}
		
		if "`model'" == "cbbetabin" {
			local expression "expression(xb() - _b[_cons]) at(mu==1)"
		}
		
		local invfn "invlogit"
		
		if "`comparative'`mcbnetwork'`pcbnetwork'" != "" {
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
		
		local confounders "`catreg'"
		
		tempname lRRcoef lRRV RRoutmatrix lORcoef lORV ORoutmatrix row ///
				outmatrixr overallRR overallOR  nltestRR nltestOR rowtestnl testmat2print bymatRR bymatOR ///
				bynltestRR bynltestOR compmatRR compmatOR compnltestRR ///
				compnltestOR catregmatRR catregmatOR catregnltestRR catregnltestOR varxcoef ///
				exactlorout exactorout exactlorouti exactorouti exactrrouti ///
				coefor coeflor lorci
				 		
		local nrowsout 0
		local nrowsnl 0
		local nby 0
		local ncomp 0
		local ncatreg 0
		
		if "`by'" != "" & "`typevarx'" == "i" & "`stratify'" == "" {		
			bayesestrcore, event(`event') varx(`varx')  confounders(`by')  cimethod(`cimethod') link(`link') 
			
			matrix `bymatRR' = r(rroutmatrix)
			matrix `bymatOR' = r(oroutmatrix)
			local nby = rowsof(`bymatRR')
			
			mat `RRoutmatrix' = `bymatRR'
			mat `ORoutmatrix' = `bymatOR'
			local nrowsout = rowsof(`RRoutmatrix')
		}
		
		if ("`by'" != "`comparator'") & ("`comparator'" != ""){
			qui label list `comparator'
			local nc = r(max)
			if (`nc' > 1) {	
		
				bayesestrcore, event(`event')  varx(`varx') confounders(`comparator')  cimethod(`cimethod') link(`link')
				
				matrix `compmatRR' = r(rroutmatrix)
				matrix `compmatOR' = r(oroutmatrix)
				local ncomp = rowsof(`compmatRR')
				
				if `nrowsout' > 0 {
					matrix `RRoutmatrix' = `RRoutmatrix' \ `compmatRR'
					matrix `ORoutmatrix' = `ORoutmatrix' \ `compmatOR'
				}
				else {
					matrix `RRoutmatrix' = `compmatRR'	
					matrix `ORoutmatrix' = `compmatOR'
				}
				local nrowsout = rowsof(`RRoutmatrix')
			}
		}		
			
		if "`catreg'" != "" {			
			if "`mcbnetwork'`pcbnetwork'" != "" | ( "`comparative'" != "" & "`interaction'" != "" ) { 
				bayesestrcore, event(`event') varx(`varx') confounders(`catreg') baselevel(`baselevel') cimethod(`cimethod') link(`link')
			}
			else {
				bayesestrcore, event(`event') confounders(`catreg') baselevel(`baselevel')  cimethod(`cimethod') link(`link') 
			}
			
			matrix `catregmatRR' = r(rroutmatrix)
			matrix `catregmatOR' = r(oroutmatrix)
			
			local ncatreg = rowsof(`catregmatRR')

			if `nrowsout' > 0 {
				matrix `RRoutmatrix' = `RRoutmatrix' \ `catregmatRR'
				matrix `ORoutmatrix' = `ORoutmatrix' \ `catregmatOR'
			}
			else {
				matrix `RRoutmatrix' = `catregmatRR'
				matrix `ORoutmatrix' = `catregmatOR'
			}
			
			local nrowsout = rowsof(`RRoutmatrix')
		}
		
		//Overall
		if ("`comparative'`mcbnetwork'`pcbnetwork'" != "") {
			
			bayesestrcore, event(`event') confounders(`varx') baselevel(`baselevel') cimethod(`cimethod') link(`link')
			
			mat `overallRR' = r(rroutmatrix)
			mat `overallOR' = r(oroutmatrix)
			
			mat `overallRR' = `overallRR'[2, 1...]  //The first row is rendundant
			mat `overallOR' = `overallOR'[2, 1...]
						
			mat rownames `overallRR' = :Overall
			mat rownames `overallOR' = :Overall
			
			if `nrowsout' > 0 {
				matrix `RRoutmatrix' = `RRoutmatrix' \ `overallRR'
				matrix `ORoutmatrix' = `ORoutmatrix' \ `overallOR'
			}
			else {
				matrix `RRoutmatrix' = `overallRR'
				matrix `ORoutmatrix' = `overallOR'
			}
			local nrowsout = rowsof(`RRoutmatrix')
		}
		
		local inltest = "no"
		
		return local inltest = "`inltest'"
		return matrix rroutmatrix = `RRoutmatrix'
		return matrix oroutmatrix = `ORoutmatrix'
	end	

/*+++++++++++++++++++++++++	SUPPORTING FUNCTIONS: FREQESTRCORE +++++++++++++++++++++++++
							Obtain the RR after modelling
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
cap program drop freqestrcore
program define freqestrcore, rclass

syntax, estimates(string) [marginlist(string) cimethod(string) varx(varname) by(varname) confounders(varlist) level(integer 95) ///
	baselevel(integer 1) link(string) model(string)]
		
	tempname lRRcoef lRRV RRoutmatrix nltestRR lORcoef lORV ORoutmatrix nltestOR
	

	//Expression for logodds prediction
	if "`link'" == "cloglog" {
		local expression "exp(logit(invcloglog(predict(xb))))"
	}
	else if "`link'" == "loglog" {
		local expression "exp(-logit(invcloglog(predict(xb))))"
	}
	else {
		local expression "exp(predict(xb))"
	}
	if "`model'" == "cbbetabin" {
		local expression "expression(xb() - _b[_cons]) at(mu==1)"
	}
	
	//Approximate sampling distribution critical value
	qui estimates restore `estimates'
	local df = e(N) -  e(k)
			
	if "`cimethod'" != "wald" {
		local critvalue invttail(`df', `=(100-`level')/200')
	}
	else {
		local critvalue -invnorm((100-`level')/200)
	}
	
	local EstRRlnexpression //log RR
	local EstORlnexpression //log OR
	foreach c of local confounders {	
		qui label list `c'
		local nlevels = r(max)
		local test_`c'
		
		if "`varx'" != "" {
			forvalues l = 1/`nlevels' {
				if `l' == 1 {
					local test_`c' = "_b[`c'_`l']"
				}
				else {
					local test_`c' = "_b[`c'_`l'] = `test_`c''"
				}
				local EstRRlnexpression = "`EstRRlnexpression' (`c'_`l': ln(invlogit(_b[`l'.`c'#2.`varx'])) - ln(invlogit(_b[`l'.`c'#1.`varx'])))"
				local EstORlnexpression = "`EstORlnexpression' (`c'_`l': _b[`l'.`c'#2.`varx'] - _b[`l'.`c'#1.`varx'])"
			}
		}
		else {					
			local test_`c' = "_b[`c'_`baselevel']"
			local init 1
			
			forvalues l = 1/`nlevels' {
				if `l' != `baselevel' {
					local test_`c' = "_b[`c'_`l'] = `test_`c''"
				}
				local EstRRlnexpression = "`EstRRlnexpression' (`c'_`l': ln(invlogit(_b[`l'.`c'])) - ln(invlogit(_b[`baselevel'.`c'])))"
				local EstORlnexpression = "`EstORlnexpression' (`c'_`l': _b[`l'.`c'] - _b[`baselevel'.`c'])"	
			}
		}
	}
	
	//RR
	qui estimates restore `estimates'
	qui margins `marginlist', `expression' over(`by') post level(`level')
	qui  nlcom `EstRRlnexpression', post level(`level')
	mat `lRRcoef' = e(b)
	mat `lRRV' = e(V)
	mat `lRRV' = vecdiag(`lRRV')	
	local ncols = colsof(`lRRcoef') //length of the vector
	local rnames :colnames `lRRcoef'

	local rowtestnl			
	local i = 1

	foreach c of local confounders {
		qui label list `c'
		local nlevels = r(max)
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
	qui estimates restore `estimates'
	qui margins `marginlist', `expression' over(`by') post level(`level')
	qui nlcom `EstORlnexpression', post level(`level')
	mat `lORcoef' = e(b)
	mat `lORV' = e(V)
	mat `lORV' = vecdiag(`lORV')	
	local ncols = colsof(`lORcoef') //length of the vector
	local rnames :colnames `lORcoef'
				
	local i = 1
	foreach c of local confounders {
		qui label list `c'
		local nlevels = r(max)
		if (`nlevels' > 2 & "`varx'" == "") | (`nlevels' > 1 & "`varx'" != "" ){
			qui testnl (`test_`c'')
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
	
	mat `RRoutmatrix' = J(`ncols', 6, .)
	mat `ORoutmatrix' = J(`ncols', 6, .)
	
	forvalues r = 1(1)`ncols' {
		mat `RRoutmatrix'[`r', 1] = exp(`lRRcoef'[1,`r']) /*Estimate*/
		mat `RRoutmatrix'[`r', 2] = sqrt(`lRRV'[1, `r']) /*se in log scale, power 1*/
		mat `RRoutmatrix'[`r', 3] = `lRRcoef'[1,`r']/sqrt(`lRRV'[1, `r']) /*Z in log scale*/
		
		mat `ORoutmatrix'[`r', 1] = exp(`lORcoef'[1,`r']) /*Estimate*/
		mat `ORoutmatrix'[`r', 2] = sqrt(`lORV'[1, `r']) /*se in log scale, power 1*/
		mat `ORoutmatrix'[`r', 3] = `lORcoef'[1,`r']/sqrt(`lORV'[1, `r']) /*Z in log scale*/
		
		if "`cimethod'" != "wald" {
			mat `RRoutmatrix'[`r', 4] = ttail(`df', abs(`RRoutmatrix'[`r', 3]))*2   /*p-value*/
			mat `ORoutmatrix'[`r', 4] = ttail(`df', abs(`ORoutmatrix'[`r', 3]))*2   /*p-value*/
		}
		else {
			mat `RRoutmatrix'[`r', 4] =  normprob(-abs(`RRoutmatrix'[`r', 3]))*2  /*p-value*/
			mat `ORoutmatrix'[`r', 4] =  normprob(-abs(`ORoutmatrix'[`r', 3]))*2  /*p-value*/
		}
		mat `RRoutmatrix'[`r', 5] = exp(`lRRcoef'[1, `r'] - `critvalue' * sqrt(`lRRV'[1, `r'])) /*lower*/
		mat `RRoutmatrix'[`r', 6] = exp(`lRRcoef'[1, `r'] + `critvalue' * sqrt(`lRRV'[1, `r'])) /*upper*/
		
		mat `ORoutmatrix'[`r', 5] = exp(`lORcoef'[1, `r'] - `critvalue' * sqrt(`lORV'[1, `r'])) /*lower*/
		mat `ORoutmatrix'[`r', 6] = exp(`lORcoef'[1, `r'] + `critvalue' * sqrt(`lORV'[1, `r'])) /*upper*/
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
	mat rownames `RRoutmatrix' = `rownames'
	mat rownames `ORoutmatrix' = `rownames'
	
	if `i' > 1 {
		mat rownames `nltestRR' = `rowtestnl'
		mat rownames `nltestOR' = `rowtestnl'
		return matrix nltestRR = `nltestRR'	
		return matrix nltestOR = `nltestOR'
	}
	return local i = "`i'"
	return matrix rroutmatrix = `RRoutmatrix'
	return matrix oroutmatrix = `ORoutmatrix'

end	
/*+++++++++++++++++++++++++	SUPPORTING FUNCTIONS: FREQESTR +++++++++++++++++++++++++
							Estimate RR after modelling
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
	cap program drop freqestr
	program define freqestr, rclass

		syntax, estimates(string) studyid(varname) [event(varname) total(varname) catreg(varlist) typevarx(string) varx(varname) comparator(varname) cimethod(string) ///
			level(integer 95) DP(integer 2) mcbnetwork pcbnetwork abnetwork general comparative stratify power(integer 0) by(varname) ///
			regexpression(string) baselevel(integer 1)  interaction link(string) model(string) inference(string)]
		
		//Expression for logodds prediction		
		if "`link'" == "cloglog" {
			local expression "exp(logit(invcloglog(predict(xb))))"
		}
		else if "`link'" == "loglog" {
			local expression "exp(-logit(invcloglog(predict(xb))))"
		}
		else {
			local expression "exp(predict(xb))"
		}
		
		if "`model'" == "cbbetabin" {
			local expression "expression(xb() - _b[_cons]) at(mu==1)"
		}
		
		local invfn "invlogit"
		//Approximate sampling distribution critical value
		if "`cimethod'" != "wald" {
			qui estimates restore `estimates'
			local df = e(N) -  e(k)
			local critvalue invttail(`df', `=(100-`level')/200')
		}
		else {
			local critvalue -invnorm((100-`level')/200)
		}
		
		if "`comparative'`mcbnetwork'`pcbnetwork'" != "" {
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
				outmatrixr overallRR overallOR  nltestRR nltestOR rowtestnl testmat2print bymatRR bymatOR ///
				bynltestRR bynltestOR compmatRR compmatOR compnltestRR ///
				compnltestOR catregmatRR catregmatOR catregnltestRR catregnltestOR varxcoef ///
				exactlorout exactorout exactlorouti exactorouti exactrrouti ///
				coefor coeflor lorci
				 		
		local nrowsout 0
		local nrowsnl 0
		local nby 0
		local ncomp 0
		local ncatreg 0
		
		if "`by'" != "" & "`typevarx'" == "i" & "`stratify'" == "" {		
			freqestrcore, marginlist(`varx') varx(`varx') by(`by') confounders(`by')  estimates(`estimates') cimethod(`cimethod') link(`link') model(`model')
			
			matrix `bymatRR' = r(rroutmatrix)
			matrix `bymatOR' = r(oroutmatrix)
			local nby = rowsof(`bymatRR')
			local iby = r(i)
			if `iby' > 1 {
				matrix `bynltestRR' = r(nltestRR)
				matrix `bynltestOR' = r(nltestOR)
				matrix `nltestRR' = `bynltestRR'
				matrix `nltestOR' = `bynltestOR'
				local nrowsnl = rowsof(`nltestRR')
			}
			mat `RRoutmatrix' = `bymatRR'
			mat `ORoutmatrix' = `bymatOR'
			local nrowsout = rowsof(`RRoutmatrix')
		}
		
		if ("`by'" != "`comparator'") & ("`comparator'" != ""){
			qui label list `comparator'
			local nc = r(max)
			if (`nc' > 1) {	
		
				freqestrcore, marginlist(`varx') varx(`varx') by(`comparator') confounders(`comparator') estimates(`estimates') cimethod(`cimethod') link(`link') model(`model')
				
				matrix `compmatRR' = r(rroutmatrix)
				matrix `compmatOR' = r(oroutmatrix)
				local ncomp = rowsof(`compmatRR')
				local icomp = r(i)
				if `icomp' > 1 {
					matrix `compnltestRR' = r(nltestRR)
					matrix `compnltestOR' = r(nltestOR)
					if `nrowsnl' > 0 {
						matrix `nltestRR' = `nltestRR' \ `compnltestRR'
						matrix `nltestOR' = `nltestOR' \ `compnltestOR'
					}
					else {
						matrix `nltestRR' = `compnltestRR'
						matrix `nltestOR' = `compnltestOR'
					}
					local nrowsnl = rowsof(`nltestRR')
				}
				
				if `nrowsout' > 0 {
					matrix `RRoutmatrix' = `RRoutmatrix' \ `compmatRR'
					matrix `ORoutmatrix' = `ORoutmatrix' \ `compmatOR'
				}
				else {
					matrix `RRoutmatrix' = `compmatRR'	
					matrix `ORoutmatrix' = `compmatOR'
				}
				local nrowsout = rowsof(`RRoutmatrix')
			}
		}		
			
		if "`marginlist'" != "" {
			
			if "`comparative'`mcbnetwork'`pcbnetwork'" != "" { 
				freqestrcore, marginlist(`marginlist') varx(`varx') confounders(`confounders') baselevel(`baselevel') estimates(`estimates') cimethod(`cimethod') link(`link') model(`model')
			}
			else {
				freqestrcore, marginlist(`marginlist') confounders(`confounders') baselevel(`baselevel') estimates(`estimates') cimethod(`cimethod') link(`link') model(`model')
			}
			
			matrix `catregmatRR' = r(rroutmatrix)
			matrix `catregmatOR' = r(oroutmatrix)
			
			local ncatreg = rowsof(`catregmatRR')
			local icatreg = r(i)
			if `icatreg' > 1 {
				matrix `catregnltestRR' = r(nltestRR)
				matrix `catregnltestOR' = r(nltestOR)
				
				if `nrowsnl' > 0 {
					matrix `nltestRR' = `nltestRR' \ `catregnltestRR'
					matrix `nltestOR' = `nltestOR' \ `catregnltestOR'
				}
				else {
					matrix `nltestRR' = `catregnltestRR'
					matrix `nltestOR' = `catregnltestOR'					
				}
				local nrowsnl = rowsof(`nltestRR')
			}
			if `nrowsout' > 0 {
				matrix `RRoutmatrix' = `RRoutmatrix' \ `catregmatRR'
				matrix `ORoutmatrix' = `ORoutmatrix' \ `catregmatOR'
			}
			else {
				matrix `RRoutmatrix' = `catregmatRR'
				matrix `ORoutmatrix' = `catregmatOR'
			}
			
			local nrowsout = rowsof(`RRoutmatrix')
		}

		if ("`comparative'`mcbnetwork'`pcbnetwork'" != "") {			
			mat `overallRR' = J(1, 6, .)
			mat `overallOR' = J(1, 6, .)			
						
			qui estimates restore `estimates'
			local df = e(N) -  e(k)
			
			qui margins `varx', `expression' post level(`level')
			
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
			qui nlcom (Overall: ln(`invfn'(_b[`coef2'.`varx'])) - ln(`invfn'(_b[`coef1'.`varx']))) 		  
			mat `lRRcoef' = r(b)
			mat `lRRV' = r(V)
			mat `lRRV' = vecdiag(`lRRV')
			
			//log or metric
			qui nlcom (Overall: _b[`coef2'.`varx'] - _b[`coef1'.`varx']) 	  
			mat `lORcoef' = r(b)
			mat `lORV' = r(V)
			mat `lORV' = vecdiag(`lORV')
			
			mat `overallRR'[1, 1] = exp(`lRRcoef'[1,1])  //rr
			mat `overallRR'[1, 2] = sqrt(`lRRV'[1, 1]) //se
			mat `overallRR'[1, 3] = `lRRcoef'[1, 1]/sqrt(`lRRV'[1, 1]) //zvalue
			
			mat `overallOR'[1, 1] = exp(`lORcoef'[1,1])  //or
			mat `overallOR'[1, 2] = sqrt(`lORV'[1, 1]) //se
			mat `overallOR'[1, 3] = `lORcoef'[1, 1]/sqrt(`lORV'[1, 1]) //zvalue
			
			
			if "`cimethod'" != "wald" {
				mat `overallRR'[1, 4] = ttail(`df', -abs(`overallRR'[1, 3]))*2 //pvalue
				mat `overallOR'[1, 4] = ttail(`df', -abs(`overallOR'[1, 3]))*2 //pvalue
			}
			else {
				mat `overallRR'[1, 4] = normprob(-abs(`overallRR'[1, 3]))*2 //pvalue
				mat `overallOR'[1, 4] = normprob(-abs(`overallOR'[1, 3]))*2 //pvalue
			}
			
			mat `overallRR'[1, 5] = exp(`lRRcoef'[1, 1] - `critvalue'*sqrt(`lRRV'[1, 1])) //ll
			mat `overallRR'[1, 6] = exp(`lRRcoef'[1, 1] + `critvalue'*sqrt(`lRRV'[1, 1])) //ul
			
			mat `overallOR'[1, 5] = exp(`lORcoef'[1, 1] - `critvalue'*sqrt(`lORV'[1, 1])) //ll
			mat `overallOR'[1, 6] = exp(`lORcoef'[1, 1] + `critvalue'*sqrt(`lORV'[1, 1])) //ul
			
			mat rownames `overallRR' = :Overall
			mat rownames `overallOR' = :Overall
			
			if `nrowsout' > 0 {
				matrix `RRoutmatrix' = `RRoutmatrix' \ `overallRR'
				matrix `ORoutmatrix' = `ORoutmatrix' \ `overallOR'
			}
			else {
				matrix `RRoutmatrix' = `overallRR'
				matrix `ORoutmatrix' = `overallOR'
			}
			local nrowsout = rowsof(`RRoutmatrix')
		}
				
		if "`cimethod'" == "wald" {
			mat colnames `RRoutmatrix' = Mean SE(lrr) z(lrr) P>|z| Lower Upper
			mat colnames `ORoutmatrix' = Mean SE(lor) z(lor) P>|z| Lower Upper
		}
		else {
			mat colnames `RRoutmatrix' = Mean SE(lrr) t(lrr) P>|t| Lower Upper
			mat colnames `ORoutmatrix' = Mean SE(lor) t(lor) P>|t| Lower Upper
		}
		
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
										
	
					mat `exactorouti`baselevel'' = (1, 0, 1, 1, 1)
					mat `exactlorouti`baselevel'' = (1, 0, 1, 1, 1)
					
					local baselab = ustrregexra("`baselab'", " ", "_")

					mat rownames `exactorouti`baselevel'' = `vari':`baselab'
					mat rownames `exactlorouti`baselevel'' = `vari':`baselab'
					
					//Other groups
					forvalues g=1(1)`ngroups' {
						if `g' != `baselevel' {
							tempvar meanphat`g' meanrrhat`g' meanorhat`g' meanlorhat`g' gid`g' sumphat`g' subsetid`g'
							tempname exactrrouti`g' exactorouti`g' exactlorouti`g'
							
							local glab:label `vari' `g'
							count if `vari' == `g' & `insample' == 1
							local ngroup`g' = r(N)	
							
							replace `varx' = 1 if `vari' == `g' & `insample' == 1
							
							cap exlogistic `event' `varx' if (`vari' == `g' | `vari' == `baselevel') & `insample' == 1, binomial(`total') level(`level') `progress'
							
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
			
			if "`comparative'" != "" | "`mcbnetwork'" != "" {
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
			mat colnames `nltestOR' = chi2 df p
			return matrix nltestRR = `nltestRR'
			return matrix nltestOR = `nltestOR'
		}
		else {
			local inltest = "no"
		}
		return local inltest = "`inltest'"
		return matrix rroutmatrix = `RRoutmatrix'
		return matrix oroutmatrix = `ORoutmatrix'
	end	
	/*+++++++++++++++++++++++++	SUPPORTING FUNCTIONS: 	RRCI +++++++++++++++++++++++++
								CI for RR
	++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
	cap program drop rrci
	program define rrci
	
		syntax varlist, R(name) lowerci(name) upperci(name) cimethod(string) [alpha(real 0.05)]
		
		qui {	
			tokenize `varlist'
			gen `r' = . 
			gen `lowerci' = .
			gen `upperci' = .
			
			local zstar =  -invnorm(`alpha'/2)
			local chisq = invchi2(1, `=1-`alpha'')
			
			count
			forvalues i = 1/`r(N)' {
				local n1 = `1'[`i']
				local N1 = `2'[`i']
				local n2 = `3'[`i']
				local N2 = `4'[`i']
				
				//Get the intervals
				if "`cimethod'" == "koopman" {
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
				if ("`cimethod'" == "adlog") {
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
				if ("`cimethod'" == "bailey") {
				
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
				if ("`cimethod'" == "katz") {
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
								if ("`method'" == "asinh") {
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
				if ("`method'" == "noether") {
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
				
				replace `r' = `rr' in `i'
				replace `lowerci' = `lorr' in `i'
				replace `upperci' = `uprr' in `i'
			}
		}
	end

	/*+++++++++++++++++++++++++	SUPPORTING FUNCTIONS: 	KOOPMANCI +++++++++++++++++++++++++
								CI for RR
	++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
	cap program drop koopmanci
	program define koopmanci

		syntax varlist, r(name) lowerci(name) upperci(name) [alpha(real 0.05)]
		
		qui {	
			tokenize `varlist'
			gen `r' = . 
			gen `lowerci' = .
			gen `upperci' = .
			
			tempname matci
			count
			forvalues i = 1/`r(N)' {
				local n1 = `1'[`i']
				local N1 = `2'[`i']
				local n2 = `3'[`i']
				local N2 = `4'[`i']

				koopmancii `n1' `N1' `n2' `N2', alpha(`alpha')
				mat `matci' = r(ci)
				
				if (`n1' == 0) &(`n2'==0) {
					replace `r' = 0 in `i'
				}
				else {
					replace `r' = (`n1'/`N1')/(`n2'/`N2')  in `i'	
				}
				replace `lowerci' = `matci'[1, 1] in `i'
				replace `upperci' = `matci'[1, 2] in `i'
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

		syntax varlist, r(name) lowerci(name) upperci(name) [alpha(real 0.05)]
		
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
				
				replace `r' = `RR' in `i'
				replace `lowerci' = `matci'[1, 1] in `i'
				replace `upperci' = `matci'[1, 2] in `i'
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

		syntax varlist, r(name) lowerci(name) upperci(name) [level(real 95) mcbnetwork cimethod(string)]
		
		qui {	
			tokenize `varlist'
			gen `r' = . 
			gen `lowerci' = .
			gen `upperci' = .
			
			count
			forvalues i = 1/`r(N)' {
				if "`mcbnetwork'" != ""  {
					local a = `1'[`i']
					local b = `2'[`i']
					local c = `3'[`i']
					local d = `4'[`i']

					cci `=`a'+`b'' `=`c'+`d'' `=`a'+`c'' `=`b' +`d'', l(`level') `cimethod'
				}
				else {
					local n1 = `1'[`i']
					local N1 = `2'[`i']
					local n2 = `3'[`i']
					local N2 = `4'[`i']

					cci `n1' `=`N1'-`n1'' `n2' `=`N2'-`n2'', l(`level') `cimethod'
				}
				local OR =  r(or)
				local lor = r(lb_or)
				local uor = r(ub_or)
				
				replace `r' = `OR' in `i'
				replace `lowerci' = `lor' in `i'
				replace `upperci' = `uor' in `i'
			}
		}
	end

	/*+++++++++++++++++++++++++	SUPPORTING FUNCTIONS: absexactci +++++++++++++++++++++++++
								CI for proportions
	++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
	cap program drop absexactci
	program define absexactci, rclass

		syntax anything(name=data id="data"), [level(real 95)]
		
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
		
		cii proportions `1' `2', `cimethod' level(`level')
		
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
		nchains(integer 4)
		thinning(integer 5) /*5*/
		burnin(integer 5000) /*5000*/
		mcmcsize(integer 2500) /*3000*/
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
		if "`design'" == "comparative" & "`outplot'" == "abs" {
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
			replace `es' = ln(max(`es', 0.00001)) if `es' != .
			replace `lci' = ln(max(`lci', 0.00001)) if `lci' != .
			replace `uci' = ln(`uci')
			
			if "`smooth'" !="" { 
				replace `modeles'  = ln(max(`modeles', 0.00001)) if `modeles' != .
				replace `modellci' = ln(max(`modellci', 0.00001)) if `modellci' != .
				replace `modeluci' = ln(`modeluci')
			}
		}
		qui summ `lci', detail
		local DXmin = r(min)
		qui summ `uci', detail
		local DXmax = r(max)
				
		if "`xlabel'" != "" {
			if "`logscale'" != "" {
				local DXmin = ln(max(min(`xlabel'), 0.00001))
				local DXmax = ln(max(`xlabel'))
			}
			else{
				local DXmin = min(`xlabel')
				local DXmax = max(`xlabel')
			}
		}
		if "`xlabel'"=="" {
			local xlabel "0, `DXmax'"
		}

		local lblcmd ""
		tokenize "`xlabel'", parse(",")
		while "`1'" != ""{
			if "`1'" != ","{
				local lbl = string(`1',"%7.3g")
				if "`logscale'" != "" {
					local val = ln(max(`1', 0.00001))
					/*
					if "`1'" == "0" {
						local val = ln(`=10^(-`dp')')
					}
					else {
						local val = ln(`1')
					}
					*/
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
					local val = ln(max(`1', 0.00001))
					/*
					if "`1'" == "0" {
						local val = ln(`=10^(-`dp')')
					}
					else {
						local val = ln(`1')
					}
					*/
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
		
		forvalues i = 1/`lcolsN'{
			gen `left`i'' = `DXmin' - 0.05*(`DXmax' - `DXmin') - `leftWDtot'*`textWD'
			local leftWDtot = `leftWDtot'-`leftWD`i''
		}

		gen `right1' = `DXmax' + 0.05*(`DXmax' - `DXmin')
		forvalues i = 2/`rcolsN'{
			local r2 = `i' - 1
			gen `right`i'' = `right`r2'' + `rightWD`r2''*`textWD'
		}

		local AXmin = `left1'
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
				qui replace `leftLB`i'' = "Mean: " + "`groupvar'" + " = " + `strgroupvar' if `use'==2
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
		local diamopts0 "lcolor("0 0 0")"
		local diamopts1 "lcolor("255 127 0")"
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
	if `"`smoothpointopts'"' != "" & strpos(`"`smoothpointopts'"',"msy") == 0 {
		local smoothpointopts = `"`smoothpointopts' msymbol(D)"' 
	}
	if `"`smoothpointopts'"' != "" & strpos(`"`smoothpointopts'"',"ms") == 0 {
		local smoothpointopts = `"`smoothpointopts' msize(vsmall)"' 
	}
	if `"`smoothpointopts'"' != "" & strpos(`"`smoothpointopts'"',"mc") == 0 {
		local smoothpointopts = `"`smoothpointopts' mcolor("0 0 0")"' 
	}
	if `"`smoothpointopts'"' == ""{
		local smoothpointopts "msymbol(D) msize(vsmall) mcolor("0 0 0")"
		if "`compabs'" != "" {
			local smoothpointopts0 "msymbol(D) msize(vsmall) mcolor("0 0 0")"
			local smoothpointopts1 "msymbol(D) msize(vsmall) mcolor("255 127 0")"
		}
	}
	else{
		local smoothpointopts `"`smoothpointopts'"'
	}
	
	// CI options
	if `"`ciopts'"' == "" {
		if "`compabs'" != "" {
				local ciopts0 = `"lcolor("0 0 0")"' 
				local ciopts1 = `" lcolor("255 127 0")"'
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
		local ciopts `"`ciopts'"'
	}
	//Smooth ci
	if `"`smoothciopts'"' == "" {
		local smoothciopts "lcolor("0 0 0")"
		
		if "`compabs'" != "" {
			local smoothciopts0 "lcolor("0 0 0")"
			local smoothciopts1 "lcolor("255 127 0")"
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
		qui label list `groupvar'
		local nlevels = r(max)
		forvalues l = 1/`nlevels' {
			summ `es' if `use' == 2  & `groupvar' == `l' 
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
	}	// end qui	
	/*===============================================================================================*/
	/*====================================  GRAPH    ================================================*/
	/*===============================================================================================*/
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
		local plotname = "name(`type', replace)"
	}
	if "$by_index_" != "" {
		local plotname = "name(`type'" + "$by_index_" + ", replace)"
		noi di as res _n  "NOTE: `fullname' plot name -> `type'$by_index_"
	}

	#delimit ;
	twoway
	 // Draw diamond to make it to construct the legend
		`diamondcommand11' `diamondcommand21'
		
		`overallCommand' `sublineCommand' `xlineCommand' `xaxis' `xaxistitle' 
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
	
		`cicommand' `cicommand1' `cicommand2'
		`smoothcommands1'	`smoothcommands10' `smoothcommands11' `smoothcommands21' `smoothcommands20'	
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
		`smoothcommands2' `overallCommand'	
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