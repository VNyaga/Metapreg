/*
Simulations
*/
*True proportion
local op "0.2"
*local op "0.2 0.5 0.9"
local nop: word count `op'

*max var 1/12 depending on mu
//For 0.2 [1/37.5, 1/50, 1/75, 1/150]
local var "0.026 0.02 0.013 0.0067"

//For 0.5 [1/12, 1/16, 1/24, 1/48]
*local var "0.08 0.06 0.04 0.02"

//For 0.9 [1/123, 1/164, 1/246, 1/492]
*local var "0.008 0.006 0.004 0.002"
 
local nvar: word count `var'

*Number of studies
*local studies  "5 10"
local studies  "5 10 20 30"
local nstudies: word count `studies'
 
*N ~ integer between 10 & 500
*alpha & beta are functions of op & var

*Study proportion
* pi = rbeta(alpha, beta)

*Study n
* n ~ binomial(pi, N) 

//File
*cd "C:\DATA\WIV\Projects\GitHub\Metapreg\Simulation\"
cd "\\sciensano.be\fs\1120_Cancer_Employee\BMH\Victoria\Simulation\"

local s = "sim; model; phat; phatlo; phatup; tauhat; ptrue; sigmatrue; phitrue; nstudies; pvalue; pbias; mse; inci"
file open results using simdatap2-1.txt, text write append //append or replace
file write results "`s'" _n

//logs
*cap log close
*log using "C:\DATA\WIV\Projects\GitHub\Metapreg\Simulation\Logs\simulation.txt", text replace

local nsim = 2000 //I need 2000
set more off
set trace off

forvalues opi = 1(1)`nop' {
	local op_i : word `opi' of `op'
	di "`op_i'"
	forvalues vari = 1(1)`nvar' {
		local var_i : word `vari' of `var'
		di "`var_i'"
		local alpha = (`op_i'^2*(1 - `op_i'))/`var_i' - `op_i'
		local beta = `alpha'*(1/`op_i' - 1)
		local sigma = 1/(`alpha'*`beta')
		local phi = 1/(1 + `alpha'*`beta')
		di "`alpha'"
		di "`beta'"
		di "`sigma'"
		di "`phi'"
		forvalues studiesi = 1(1)`nstudies' {
			local studies_i : word `studiesi' of `studies'
			di "`studies_i'"
			forvalues r = 1(1)`nsim' {
				//Generate the dataset
				clear
				set seed `r'
				
				qui set obs `studies_i'
				
				qui gen studyid = _n
				qui gen p = .
				qui gen total = .
				qui gen events = .
				local cobs = 1
				
				while `cobs' < `=`studies_i' + 1' {
					qui replace p = rbeta(`alpha', `beta') in `cobs'				
					qui replace total = runiformint(10, 500) in `cobs'
					qui replace events = rbinomial(total, p) in `cobs'
					local cobs = `cobs' + 1
				}
				//Replace if p=1
				qui replace events = total if events==. & p==1
				
				//Fit the models
				qui {
					//1----------- betabinomial - wald - logit
					capture betabin events, n(total) link(logit) iterate(100) nolog
					mat coef = r(table)
					local chi = e(chi2_c)
					local pvalue = chi2tail(1, `chi')/2
					local success = e(converged) 
					
					local phat = invlogit(coef[1,1])
					local phatlo = invlogit(coef[5,1])
					local  phatup = invlogit(coef[6,1])
					local tauhat = coef[1,3] //sigma
					
					local pbias = `phat' - `op_i'
					local mse = (`phat' - `op_i')^2
					local inci = 0
					if ((`phatlo' == `op_i') | (`phatlo' < `op_i')) & ((`phatup' > `op_i') | (`phatup' == `op_i'))  {
						local inci = 1
					}
					
					if `success' {					
						local s = "`r'; betabin; `phat'; `phatlo'; `phatup'; `tauhat'; `op_i'; `sigma'; `phi'; `studies_i'; `pvalue'; pbias; mse; inci "
						file write results "`s'" _n
					}
															
					//2-------logistic + random + c	
					capture metapreg events total, model(random, sformat(%8.4f) ) ///
						studyid(study) cimethod(wilson, wald) dp(4) ///
						 noitable nograph

					mat absout = e(absout)
					mat hetout = e(hetout)

					local phat = absout[1,1]
					local phatlo = absout[1,5]
					local  phatup = absout[1,6]
					local tauhat = hetout[1,4]
					local pvalue = hetout[1,3]
					
					local pbias = `phat' - `op_i'
					local mse = (`phat' - `op_i')^2
					local inci = 0
					if ((`phatlo' == `op_i') | (`phatlo' < `op_i')) & ((`phatup' > `op_i') | (`phatup' == `op_i'))  {
						local inci = 1
					}
				
					if _rc == 0 {
						local s = "`r'; logistic-re-c; `phat'; `phatlo'; `phatup'; `tauhat'; `op_i'; `sigma'; `phi';  `studies_i'; `pvalue'; pbias; mse; inci"
						file write results "`s'" _n
					}
					
					//3-------logistic + random + m	
					capture metapreg events total, model(random, sformat(%8.4f) ) ///
						studyid(study) cimethod(wilson, wald) dp(4) ///
						 sumstat(ES) sumtable(abs) nograph 
						 
					mat meanout = e(meanout)
					mat hetout = e(hetout)

					local phat = meanout[1,1]
					local phatlo = meanout[1,3]
					local  phatup = meanout[1,4]
					local tauhat = hetout[1,4]
					local pvalue = hetout[1,3]
					
					local pbias = `phat' - `op_i'
					local mse = (`phat' - `op_i')^2
					local inci = 0
					if ((`phatlo' == `op_i') | (`phatlo' < `op_i')) & ((`phatup' > `op_i') | (`phatup' == `op_i'))  {
						local inci = 1
					}
					
					if _rc == 0 {
						local s = "`r'; logistic-re-m; `phat'; `phatlo'; `phatup'; `tauhat'; `op_i'; `sigma'; `phi';  `studies_i'; `pvalue'; pbias; mse; inci"
						file write results "`s'" _n
					}
					
					//4-------logistic + fixed + t	
					capture metapreg events total, model(fixed, sformat(%8.4f) ) ///
						studyid(study) cimethod(wilson, t) dp(4) ///
						 noitable nograph

					mat absout = e(absout)
					mat hetout = e(hetout)

					local phat = absout[1,1]
					local phatlo = absout[1,5]
					local  phatup = absout[1,6]
					local tauhat = 0
					local pvalue = .
				
					local pbias = `phat' - `op_i'
					local mse = (`phat' - `op_i')^2
					local inci = 0
					if ((`phatlo' == `op_i') | (`phatlo' < `op_i')) & ((`phatup' > `op_i') | (`phatup' == `op_i'))  {
						local inci = 1
					}
					
					if _rc == 0 {
						local s = "`r'; logistic-ce-t; `phat'; `phatlo'; `phatup'; `tauhat'; `op_i'; `sigma'; `phi';  `studies_i'; `pvalue'; pbias; mse; inci"
						file write results "`s'" _n
					}
					//2-------exact ce	 
					capture metapreg events total, model(exact, sformat(%8.4f) ) ///
						studyid(study) cimethod(wilson, exact) dp(4) ///
						 sumstat(ES) sumtable(all) nograph

					mat absout = e(absout)
					mat hetout = e(hetout)	 
					local phat = absout[1,1]
					local phatlo = absout[1,5]
					local phatup = absout[1,6]
					local tauhat = 0
					local pvalue = .
				
					local pbias = `phat' - `op_i'
					local mse = (`phat' - `op_i')^2
					local inci = 0
					if ((`phatlo' == `op_i') | (`phatlo' < `op_i')) & ((`phatup' > `op_i') | (`phatup' == `op_i'))  {
						local inci = 1
					}
					
					if _rc == 0 {					
						local s = "`r'; exact-ce; `phat'; `phatlo'; `phatup'; `tauhat'; `op_i'; `sigma'; `phi';  `studies_i'; `pvalue'; pbias; mse; inci"
						file write results "`s'" _n 
					}

					//3-------IV + random
					 capture metaprop events total, ///
					 random nograph  dp(4)

					local phat = r(ES)
					local phatlo = r(ci_low)
					local phatup = r(ci_upp)
					local tauhat = r(tau2) 
					local pvalue = .
					
					local pbias = `phat' - `op_i'
					local mse = (`phat' - `op_i')^2
					local inci = 0
					if ((`phatlo' == `op_i') | (`phatlo' < `op_i')) & ((`phatup' > `op_i') | (`phatup' == `op_i'))  {
						local inci = 1
					}
					
					local s = "`r'; iv-re; `phat'; `phatlo'; `phatup'; `tauhat'; `op_i'; `sigma'; `phi'; `studies_i'; `pvalue'; pbias; mse; inci"
					file write results "`s'" _n 
					 
					//4-------IV + fixed
					capture metaprop events total, ///
					 lcols(study) fixed nograph  dp(4)
					 
					local phat = r(ES)
					local phatlo = r(ci_low)
					local phatup = r(ci_upp)
					local tauhat = r(tau2) 
					local pvalue = .
					
					local pbias = `phat' - `op_i'
					local mse = (`phat' - `op_i')^2
					local inci = 0
					if ((`phatlo' == `op_i') | (`phatlo' < `op_i')) & ((`phatup' > `op_i') | (`phatup' == `op_i'))  {
						local inci = 1
					}
					
					
					local s = "`r'; iv-ce; `phat'; `phatlo'; `phatup'; `tauhat'; `op_i'; `sigma'; `phi'; `studies_i'; `pvalue'; pbias; mse; inci"
					file write results "`s'" _n 

					//5-------IV + ftt (harmonic) + random
					capture metaprop events total, ftt ///
					 lcols(study) random nograph  dp(4)
					 
					local phat = r(ES)
					local phatlo = r(ci_low)
					local phatup = r(ci_upp)
					local tauhat = r(tau2) 
					local pvalue = .
					
					local pbias = `phat' - `op_i'
					local mse = (`phat' - `op_i')^2
					local inci = 0
					if ((`phatlo' == `op_i') | (`phatlo' < `op_i')) & ((`phatup' > `op_i') | (`phatup' == `op_i'))  {
						local inci = 1
					}
					
					
					local s = "`r'; iv-ftt-harm-re; `phat'; `phatlo'; `phatup'; `tauhat'; `op_i'; `sigma'; `phi'; `studies_i'; `pvalue'; pbias; mse; inci"
					file write results "`s'" _n 
					 
					//6-------IV + ftt(harmonic) + fixed
					capture metaprop events total, ftt ///
					 lcols(study) fixed nograph  dp(4)
					 
					local phat = r(ES)
					local phatlo = r(ci_low)
					local phatup = r(ci_upp)
					local tauhat = r(tau2) 
					local pvalue = .
					
					local pbias = `phat' - `op_i'
					local mse = (`phat' - `op_i')^2
					local inci = 0
					if ((`phatlo' == `op_i') | (`phatlo' < `op_i')) & ((`phatup' > `op_i') | (`phatup' == `op_i'))  {
						local inci = 1
					}
					
					local s = "`r'; iv-ftt-harm-ce; `phat'; `phatlo'; `phatup'; `tauhat'; `op_i'; `sigma'; `phi';  `studies_i'; `pvalue'; pbias; mse; inci"
					file write results "`s'" _n 

					//7------- IVHET + ftt (harmonic)
					capture metan events total,  ivhet pr ftt nograph
					local phat = r(prop_eff)
					local phatlo = r(prop_lci)
					local phatup = r(prop_uci)
					local tauhat = r(tausq)
					local pvalue = .
					
					local pbias = `phat' - `op_i'
					local mse = (`phat' - `op_i')^2
					local inci = 0
					if ((`phatlo' == `op_i') | (`phatlo' < `op_i')) & ((`phatup' > `op_i') | (`phatup' == `op_i'))  {
						local inci = 1
					}
					
					
					local s = "`r'; ivhet-ftt-harm; `phat'; `phatlo'; `phatup'; `tauhat'; `op_i'; `sigma'; `phi';  `studies_i'; `pvalue'; pbias; mse; inci"
					file write results "`s'" _n 
					 

					//8------- IVHET + ftt (iv)
					capture metan events total, ivhet pr  tr(ft, iv) nograph
					local phat = r(prop_eff)
					local phatlo = r(prop_lci)
					local phatup = r(prop_uci)
					local tauhat = r(tausq) 
					local pvalue = .
					
					local pbias = `phat' - `op_i'
					local mse = (`phat' - `op_i')^2
					local inci = 0
					if ((`phatlo' == `op_i') | (`phatlo' < `op_i')) & ((`phatup' > `op_i') | (`phatup' == `op_i'))  {
						local inci = 1
					}
					
					local s = "`r'; ivhet-ftt-iv; `phat'; `phatlo'; `phatup'; `tauhat'; `op_i'; `sigma'; `phi';  `studies_i'; `pvalue'; pbias; mse; inci"
					file write results "`s'" _n 

					//9------- ce + ftt (iv)
					capture metan events total, fixed pr  tr(ft, iv)  nograph
					local phat = r(prop_eff)
					local phatlo = r(prop_lci)
					local phatup = r(prop_uci)
					local tauhat = 0 
					local pvalue = .
					
					local pbias = `phat' - `op_i'
					local mse = (`phat' - `op_i')^2
					local inci = 0
					if ((`phatlo' == `op_i') | (`phatlo' < `op_i')) & ((`phatup' > `op_i') | (`phatup' == `op_i'))  {
						local inci = 1
					}
					
					local s = "`r'; ce-ftt-iv; `phat'; `phatlo'; `phatup'; `tauhat'; `op_i'; `sigma'; `phi';  `studies_i'; `pvalue'; pbias; mse; inci"
					file write results "`s'" _n 

					//10------- re + ftt (iv)
					capture metan events total, random pr  tr(ft, iv)  nograph
					local phat = r(prop_eff)
					local phatlo = r(prop_lci)
					local phatup = r(prop_uci)
					local tauhat = r(tausq) 
					local pvalue = .
					
					local pbias = `phat' - `op_i'
					local mse = (`phat' - `op_i')^2
					local inci = 0
					if ((`phatlo' == `op_i') | (`phatlo' < `op_i')) & ((`phatup' > `op_i') | (`phatup' == `op_i'))  {
						local inci = 1
					}
					
					
					local s = "`r'; re-ftt-iv; `phat'; `phatlo'; `phatup'; `tauhat'; `op_i'; `sigma'; `phi';  `studies_i'; `pvalue'; pbias; mse; inci"
					file write results "`s'" _n 

					//11------- ce + ftt (arithmetic)
					capture metan events total, fixed pr  tr(ft, a)  nograph
					local phat = r(prop_eff)
					local phatlo = r(prop_lci)
					local phatup = r(prop_uci)
					local tauhat = 0 	
					local pvalue = .
					
					local pbias = `phat' - `op_i'
					local mse = (`phat' - `op_i')^2
					local inci = 0
					if ((`phatlo' == `op_i') | (`phatlo' < `op_i')) & ((`phatup' > `op_i') | (`phatup' == `op_i'))  {
						local inci = 1
					}
					
					local s = "`r'; ce-ftt-arith; `phat'; `phatlo'; `phatup'; `tauhat'; `op_i'; `sigma'; `phi';  `studies_i'; `pvalue'; pbias; mse; inci"
					file write results "`s'" _n 

					//12------- IVHET + ftt (arithmetic)
					capture metan events total, ivhet pr  tr(ft, a)  nograph
					local phat = r(prop_eff)
					local phatlo = r(prop_lci)
					local phatup = r(prop_uci)
					local tauhat = r(tausq) 
					local pvalue = .
					
					local pbias = `phat' - `op_i'
					local mse = (`phat' - `op_i')^2
					local inci = 0
					if ((`phatlo' == `op_i') | (`phatlo' < `op_i')) & ((`phatup' > `op_i') | (`phatup' == `op_i'))  {
						local inci = 1
					}
					
					local s = "`r'; ivhet-ftt-arith; `phat'; `phatlo'; `phatup'; `tauhat'; `op_i'; `sigma'; `phi'; `studies_i'; `pvalue'; pbias; mse; inci"
					file write results "`s'" _n 

					//13------- re + ftt (arithmetic)
					capture metan events total, random pr  tr(ft, a)  nograph
					local phat = r(prop_eff)
					local phatlo = r(prop_lci)
					local phatup = r(prop_uci)
					local tauhat = r(tausq) 
					local pvalue = .
					
					local pbias = `phat' - `op_i'
					local mse = (`phat' - `op_i')^2
					local inci = 0
					if ((`phatlo' == `op_i') | (`phatlo' < `op_i')) & ((`phatup' > `op_i') | (`phatup' == `op_i'))  {
						local inci = 1
					}
					
					local s = "`r'; re-ftt-arith; `phat'; `phatlo'; `phatup'; `tauhat'; `op_i'; `sigma'; `phi'; `studies_i'; `pvalue'; pbias; mse; inci"
					file write results "`s'" _n 

					//14------- ce + ftt (geometric)
					capture metan events total, fixed pr  tr(ft, g)  nograph
					local phat = r(prop_eff)
					local phatlo = r(prop_lci)
					local phatup = r(prop_uci)
					local tauhat = 0 	
					local pvalue = .
					
					local pbias = `phat' - `op_i'
					local mse = (`phat' - `op_i')^2
					local inci = 0
					if ((`phatlo' == `op_i') | (`phatlo' < `op_i')) & ((`phatup' > `op_i') | (`phatup' == `op_i'))  {
						local inci = 1
					}
					
					local s = "`r'; ce-ftt-geom; `phat'; `phatlo'; `phatup'; `tauhat'; `op_i'; `sigma'; `phi'; `studies_i'; `pvalue'; pbias; mse; inci"
					file write results "`s'" _n 

					//15------- IVHET + ftt (geometric)
					capture metan events total, ivhet pr  tr(ft, g)  nograph
					local phat = r(prop_eff)
					local phatlo = r(prop_lci)
					local phatup = r(prop_uci)
					local tauhat = r(tausq) 
					local pvalue = .
					
					local pbias = `phat' - `op_i'
					local mse = (`phat' - `op_i')^2
					local inci = 0
					if ((`phatlo' == `op_i') | (`phatlo' < `op_i')) & ((`phatup' > `op_i') | (`phatup' == `op_i'))  {
						local inci = 1
					}
					
					local s = "`r'; ivhet-ftt-geom; `phat'; `phatlo'; `phatup'; `tauhat'; `op_i'; `sigma'; `phi';  `studies_i'; `pvalue'; pbias; mse; inci"
					file write results "`s'" _n 

					//16------- re + ftt (geometric)
					capture metan events total, random pr  tr(ft, g)  nograph
					local phat = r(prop_eff)
					local phatlo = r(prop_lci)
					local phatup = r(prop_uci)
					local tauhat = r(tausq) 
					local pvalue = .
					
					local pbias = `phat' - `op_i'
					local mse = (`phat' - `op_i')^2
					local inci = 0
					if ((`phatlo' == `op_i') | (`phatlo' < `op_i')) & ((`phatup' > `op_i') | (`phatup' == `op_i'))  {
						local inci = 1
					}
					
					local s = "`r'; re-ftt-geom; `phat'; `phatlo'; `phatup'; `tauhat'; `op_i'; `sigma';  `phi'; `studies_i'; `pvalue'; pbias; mse; inci"
					file write results "`s'" _n 

					//17------- IVHET + arcsine
					capture metan events total, ivhet pr  tr(ar)  nograph
					local phat = r(prop_eff)
					local phatlo = r(prop_lci)
					local phatup = r(prop_uci)
					local tauhat = r(tausq) 
					local pvalue = .
					
					local pbias = `phat' - `op_i'
					local mse = (`phat' - `op_i')^2
					local inci = 0
					if ((`phatlo' == `op_i') | (`phatlo' < `op_i')) & ((`phatup' > `op_i') | (`phatup' == `op_i'))  {
						local inci = 1
					}
					
					local s = "`r'; ivhet-arcsine; `phat'; `phatlo'; `phatup'; `tauhat'; `op_i'; `sigma'; `phi'; `studies_i'; `pvalue'; pbias; mse; inci"
					file write results "`s'" _n 

					//18------- re+ arcsine
					capture metan events total, random pr  tr(ar)  nograph
					local phat = r(prop_eff)
					local phatlo = r(prop_lci)
					local phatup = r(prop_uci)
					local tauhat = r(tausq) 
					local pvalue = .
					
					local pbias = `phat' - `op_i'
					local mse = (`phat' - `op_i')^2
					local inci = 0
					if ((`phatlo' == `op_i') | (`phatlo' < `op_i')) & ((`phatup' > `op_i') | (`phatup' == `op_i'))  {
						local inci = 1
					}
					
					local s = "`r'; re-arcsine; `phat'; `phatlo'; `phatup'; `tauhat'; `op_i'; `sigma'; `phi';  `studies_i'; `pvalue'; pbias; mse; inci"
					file write results "`s'" _n 

					//19------- ce + arcsine
					capture metan events total, fixed pr  tr(ar)  nograph
					local phat = r(prop_eff)
					local phatlo = r(prop_lci)
					local phatup = r(prop_uci)
					local tauhat = 0 
					local pvalue = .
					
					local pbias = `phat' - `op_i'
					local mse = (`phat' - `op_i')^2
					local inci = 0
					if ((`phatlo' == `op_i') | (`phatlo' < `op_i')) & ((`phatup' > `op_i') | (`phatup' == `op_i'))  {
						local inci = 1
					}
					
					local s = "`r'; ce-arcsine; `phat'; `phatlo'; `phatup'; `tauhat'; `op_i'; `sigma'; `phi';  `studies_i'; `pvalue'; pbias; mse; inci"
					file write results "`s'" _n 
				
				}				
			}
			di "---"
		}
	}
}

	 
file close results