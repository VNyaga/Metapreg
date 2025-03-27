//Process the variables
//Process the variables
replace stat=ustrltrim(stat)
replace model=ustrltrim(model)

//To numeric if string
ds
local vlist = r(varlist)

foreach v of local vlist {
	if "`v'" != "model" & "`v'" != "stat" & "`v'" != "ic" {
		destring `v', replace force
	} 
}

//New covariances
//zeroint >>commonint
replace model = regexr(model, "zeroint", "commonint")

//zeroslope >>commonslope
replace model = regexr(model, "zeroslope", "commonslope")

//crbetabin1 >>commonslope
replace model = regexr(model, "none", "commonslope") if strpos(model, "crbetabin") != 0 & strpos(model, "1") != 0

//crbetabin2 >>randomslope
replace model = regexr(model, "none", "randomslope") if strpos(model, "crbetabin") != 0 & strpos(model, "2") != 0

//fixed-none >>commonint
replace model = regexr(model, "none", "commonint") if strpos(model, "fixed") & strpos(model, "metapreg") != 0

replace model = ustrregexra(model, "-", " ")
forvalues i=1/10 {
	gen v`i' = word(model, `i')
}

/*
forvalues r=1(1)`=_N' {
	local modeltxt = model[`r']
	
	local i = 1
	
	while "`modeltxt'" != "" {
		gettoken parm modeltxt : modeltxt , parse("-")	
		if "`parm'" != "-" {
			cap confirm var v`i'
			if _rc != 0 {
				gen v`i' = ""
			}
			
			replace v`i' = "`parm'" in `r'
			local ++ i
		}
	}
}
*/
rename v1 package
rename v2 inference

qui pwf 

local wf = "`=r(currentframe)'"
if  "`wf'" == "observedresults" {
	replace ic=ustrltrim(ic)

	gen AIC = ""
	gen BIC = ""
	gen DIC = ""
	gen lrp = ""
	gen logBF = ""
	gen run = ""

	forvalues r=1(1)`=_N' {
		local ic = ic[`r']
			
		if "`ic'" != "." {
			tokenize `ic' , parse("-")
			if "`9'" != "" {
				replace AIC = "`1'" in `r'
				replace BIC = "`3'" in `r'
				replace lrp = "`5'`6'`7'" in `r'
				replace run = "`9'" in `r'
			}
			else if "`7'" != "" {
				replace AIC = "`1'" in `r'
				replace BIC = "`3'" in `r'
				replace lrp = "`5'" in `r'
				replace run = "`7'" in `r'
			}
			else if "`4'" == "" {
				replace AIC = "`1'" in `r'
				replace BIC = "`3'" in `r'
			}
			else if "`6'" == ""  {
				replace DIC = "`1'" in `r'
				replace logBF = "`3'" in `r'
				replace run = "`5'" in `r'
			}
			else if "`6'" != ""  {
				replace DIC = "`1'" in `r'
				replace logBF = "`3'`4'" in `r'
				replace run = "`6'" in `r'
			}
		}
	}

	destring AIC, replace force
	destring BIC, replace force
	destring DIC, replace force
	destring logBF, replace force
	destring run, replace force
	destring lrp, replace force
}
else {
	cap confirm var ic
	if _rc == 0 {

		destring ic, replace force
	}
}

//metapreg
gen Design = v3  if package == "metapreg"

gen link = v7 if  package == "metapreg"
gen CImethod = v8 if  package == "metapreg"
replace CImethod = "eti" if  package == "metapreg" & strpos(stat, "pop") != 0
replace CImethod = "z" if  package == "metapreg" & v8 == "" & inference =="frequentist"
replace CImethod = "hpd" if  package == "metapreg" & v8 == "" & inference =="bayesian"


gen dist = v4 + v5 if  package == "metapreg" & v4 =="crbetabin"
replace dist = "bin-normal" if  package == "metapreg" & v4 =="mixed"
replace dist = "binomial" if  package == "metapreg" & (v4 =="fixed" | v4 =="hexact" )
replace dist = "poisson" if  package == "metapreg" & v4 =="fixed" & link == "log"
replace dist = "poi-normal" if  package == "metapreg" & v4 =="mixed" & link == "log"

gen covariance = v6 if  package == "metapreg" 
gen effect = v4 if  package == "metapreg" 
replace effect = "mixed" if package == "metapreg"  & v4 =="crbetabin"

gen Prior = "IG(0.01, 0.01)" if  package == "metapreg" & inference =="bayesian" & v5 == "1" &  v4 =="mixed"
replace Prior = "IG(0.001, 0.001)" if  package == "metapreg" & inference =="bayesian" & v5 == "2" &  v4 =="mixed"
replace Prior = "Jeffreys" if  package == "metapreg" & inference =="bayesian" & v5 == "3" &  v4 =="mixed"
replace Prior = "IW(2,3,I(2))" if  package == "metapreg" & inference =="bayesian" & v5 == "1" &  v4 =="mixed" & covariance=="unstructured"

//The ci for absout in the loglog are switched up
cap gen switched = esthatlo > esthatup & esthatlo != .
cap confirm var esthatlo
if _rc == 0 {
	gen hold = esthatlo if package=="metapreg" & link=="loglog" & stat=="mean-absout" & package == "metapreg"
	replace esthatlo = esthatup if package=="metapreg" & link=="loglog" & stat=="mean-absout" & package == "metapreg" & switched == 1
}
cap confirm var esthatup
if _rc == 0 {
	replace esthatup = hold if package=="metapreg" & link=="loglog" & stat=="mean-absout" & package == "metapreg" & switched == 1
}
cap drop hold switched

gen Env = "R" 
replace Env = "Stata" if package == "metapreg" | package == "smeta" | package == "metan"

//R meta
gen slope = v3 if package == "meta"
gen Weighting = v4 if package == "meta" 

gen Sigmethod = v5 if package == "meta" & slope == "random"
replace Weighting = "IV" if v4=="MH" & package == "meta" & slope == "random"
replace CImethod = v6 if package == "meta"
replace CImethod = "t" if package == "meta" & slope != "random"

replace dist = "normal" if inlist(package, "bayesmeta", "metabma", "meta", "smeta", "metan", "metasem", "randmeta")  
*replace Weighting = "IV" if inlist(package, "bayesmeta", "metabma", "meta", "smeta", "metan", "metaplus") 
replace Weighting = "IV" if inlist(package,  "metabma", "metaplus", "randmeta")  
 
//bayesmeta
replace slope = "random" if package == "bayesmeta" 
replace CImethod = v4 if package == "bayesmeta" 
replace Prior = "0.5Cauchy(1)" if package == "bayesmeta" & v6 == "1"
replace Prior = "0.5normal(1)" if package == "bayesmeta" & v6 == "2"

//metabma
replace slope = "random" if v3=="DL" 
replace slope = "common" if v3=="iv"

//metabma
gen loglest = v4 if package == "metabma"
gen sumethod = v5 if package == "metabma"
replace Prior = "0.5cauchy(1)"  if package == "metabma" & slope == "random" & v6 == "1"
replace Prior = "0.5normal(1)" if package == "metabma" & slope == "random" & v6 == "2"
replace Prior = "IG(0.01, 0.01)" if package == "metabma" & slope == "random" & v6 == "3"
replace Prior = "IG(0.001, 0.001)"  if package == "metabma" & slope == "random" & v6 == "4"

//metafor
replace slope = "common" if v4 == "EE" 
replace slope = "common" if  v3 == "UM.FS" 
replace slope = "random" if v3 == "UM.RS" & v4 == "ML"
replace dist = "NCHG" if v3 == "CM.EL"
replace dist = "ANCHG" if v3 == "CM.AL"
replace dist = "bin-norm" if v3 == "UM.FS" | v3 == "UM.RS"
replace covariance = "independent" if v3 == "UM.RS" & v4 == "ML" & v6 == "1"
replace covariance = "unstructured" if v3 == "UM.RS" & v4 == "ML" & v6 == "2"
replace covariance = "freeint" if  v3 == "UM.FS" 
replace covariance = "commonslope" if v3 == "UM.RS" & v4 == "EE" 
replace CImethod = v5 if package =="metafor"

//metan
replace slope = "common" if inlist(v3, "mh", "iv", "peto") & package =="metan"
replace slope = "random" if !inlist(v3, "mh", "iv", "peto") & package =="metan"
replace Sigmethod = v3 if !inlist(v3, "mh", "iv", "peto", "ivhet") & package =="metan"
replace Sigmethod = "dl" if v3 == "ivhet" & package =="metan"
replace dist = "quasi-norm" if v3 == "ivhet" & package =="metan"
replace Weighting = v3 if inlist(v3, "mh", "peto") & package =="metan"
replace Weighting = "iv" if !inlist(v3, "mh", "peto") & package =="metan"
replace CImethod = "z" if package =="metan"
replace CImethod = "chi2" if package =="metan" & inlist(v3, "mh", "peto", "pl") 
replace CImethod = "t" if package =="metan" & inlist(v3, "hk", "kr") 

//metaplus
replace CImethod = v3 if package == "metaplus"
replace dist = "normal" if package == "metaplus" & v4 == "normal"
replace dist = "t-dist" if package == "metaplus" & v4 == "t"
replace dist = "normmix" if package == "metaplus" & v4 == "mixture"
replace Sigmethod = "ml" if package == "metaplus" & v4 != "mixture"
replace Sigmethod = "em" if package == "metaplus" & v4 == "mixture"

//metasem
replace slope = "common" if v3 == "IV" & package == "metasem"
replace slope = "random" if v3 == "DL" & package == "metasem"
replace CImethod = v4 if package == "metasem"
replace dist = "normal" if package == "metasem"
replace Sigmethod = "ml" if package == "metasem" & slope == "random" 

//metastan
replace CImethod = v4 if package == "metastan"
replace slope = "common" if package == "metastan" & v3 == "FE"
replace slope = "random" if package == "metastan" & v3 == "RE"
replace Prior = "0.5normal(0.5)" if package == "metastan" & v3 == "RE" & v5 =="normal"
replace Prior = "0.5cauchy(0.5)" if package == "metastan" & v3 == "RE" & v5 =="cauchy"
replace Prior = "uniform(0.5)" if package == "metastan" & v3 == "RE" & v5 =="uniform"
gen param = v6 if package == "metastan"
replace dist = "binomial" if slope == "common" & package == "metastan"
replace dist = "bin-normal" if slope == "random" & package == "metastan"
replace covariance = "freeint" if package == "metastan"

//randmeta
replace CImethod = v3 if package == "randmeta"
replace slope = "random" if package == "randmeta"

//smeta
replace slope = "common" if inlist(v3, "iv", "mh") & package == "smeta"
replace slope = "random" if !inlist(v3, "iv", "mh") & package == "smeta"
replace Sigmethod = v3 if slope == "random" & package == "smeta"
replace CImethod = v5 if  package == "smeta"
replace Weighting = "IV" if  package == "smeta" & !inlist(v3, "mh")
replace Weighting = v3 if  package == "smeta" & inlist(v3, "mh")
replace package = "meta" if package == "smeta" 

//Abbreviations
gen Dist = dist
replace Dist = "QN" if dist=="quasi-norm"
replace Dist = "N" if dist=="normal" & slope=="common" 
replace Dist = "NN" if dist=="normal" & slope=="random" | inlist(package, "metasem", "metaplus")
replace Dist = "BB1" if dist =="crbetabin1"
replace Dist = "BB2" if dist =="crbetabin2"
replace Dist = "BN" if dist =="bin-normal" | dist =="bin-norm" 
replace Dist = "PN" if dist =="poi-normal"
replace Dist ="B1" if dist == "binomial"
replace Dist ="B2" if dist == "binomial" & strpos(model, "hexact") != 0 
replace Dist ="P" if dist == "poisson"
replace Dist = "Nt" if dist =="t-dist"
replace Dist = "Nmix" if dist == "normmix"

gen Link = "LGT" if link == "logit"
replace Link = "LLG" if link == "loglog"
replace Link = "CLLG" if link == "cloglog"
replace Link = "L" if link == "log"

gen Effect = "M" if effect == "mixed"
replace Effect = "F" if effect == "fixed" 
replace Effect = "H" if effect == "hexact" 

gen Slope = "C" if slope == "common"
replace Slope = "R" if slope == "random"

gen Covariance = "FI" if covariance == "freeint"
replace Covariance = "CI" if covariance == "commonint"
replace Covariance = "IND" if covariance == "independent"
replace Covariance = "RS" if covariance == "randomslope"
replace Covariance = "CS" if covariance == "commonslope"
replace Covariance = "UN" if covariance == "unstructured"


gen Inf = "F" if inference == "frequentist"
replace Inf = "B" if inference != "frequentist"

gen Type = "C" if strpos(stat, "pop") == 0 & strpos(model, "betabin") == 0
replace Type = "PA" if strpos(stat, "pop") != 0 | strpos(model, "betabin") != 0 | strpos(model, "hexact") != 0 

gen Stat = "Mean" if strpos(stat, "mean") != 0
replace Stat = "Median" if strpos(stat, "median") != 0

gen Measure = "Base Prop"  if strpos(stat, "absout") != 0
replace Measure = "RR"  if strpos(stat, "rrout") != 0
replace Measure = "OR"  if strpos(stat, "orout") != 0
replace Measure = "RD"  if strpos(stat, "rdout") != 0

replace esthat = . if Covariance == "FI" & package == "metapreg" & strpos(stat, "pop") == 0 & Measure == "RR" 
replace esthatlo = . if Covariance == "FI" & package == "metapreg" & strpos(stat, "pop") == 0 & Measure == "RR"
replace esthatup = . if Covariance == "FI" & package == "metapreg" & strpos(stat, "pop") == 0 & Measure == "RR"

replace CImethod = "z" if CImethod == "wald"
replace CImethod = "exact" if Dist =="B2"
replace CImethod = "eti" if CImethod == "central"
replace CImethod = "hpd" if CImethod == "shortest"

destring v9, gen(Covariates) force
replace Covariates = 1 if Covariates == .

gen Interaction = v10

replace Sigmethod = strlower(Sigmethod)
replace Prior = strlower(Prior)
replace CImethod = strlower(CImethod)

//meta
drop if package == "meta" & CImethod == "kr" & slope =="random" & Sigmethod != "reml"  //kr only for reml
replace Sigmethod = "mp" if Sigmethod =="pm" & package == "meta"
replace CImethod = "tkh" if CImethod == "hk"
replace CImethod = "tkr" if CImethod == "kr"
drop if stat=="median-orout" &  package == "meta"

replace Weighting = "IV" if Weighting == "iv"

replace Design = "general"  if Covariance == "FI" & Dist == "B1"
drop if strpos(stat, "pop") != 0 & strpos(model, "betabin") != 0

*tkht is same as tkh
drop if CImethod == "tkht"

label var tau2hat "tau2"	
label var sigma2hat "sigma2"

