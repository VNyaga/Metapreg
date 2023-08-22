/*
Simulations
*/
*True overall proportion
local op "0.2"
*local op "0.2 0.5 0.9"
local nop: word count `op'

*Number of studies
local studies  "3"
*local studies  "5 10 20 30"
local nstudies: word count `studies'
 
*N ~ integer between 10 & 500

*Study proportion
* pi = op

*Study n
* n ~ binomial(pi, N) 

//File
cd "C:\DATA\WIV\Projects\GitHub\Metapreg\Simulation\"
cd "\\sciensano.be\fs\1120_Cancer_Employee\BMH\Victoria\Simulation\"

local s = "sim; model; phat; phatlo; phatup; tauhat; ptrue; tautrue; nstudies"
file open results using sim.txt, text write append //append or replace
file write results "`s'" _n

//logs
*cap log close
*log using "C:\DATA\WIV\Projects\GitHub\Metapreg\Simulation\Logs\simulation.txt", text replace

local nsim = 1  //20 takes forever
set more off


forvalues opi = 1(1)`nop' {
	local op_i : word `opi' of `op'
	di "`op_i'"
		forvalues studiesi = 1(1)`nstudies' {
			local studies_i : word `studiesi' of `studies'
			di "`studies_i'"
			forvalues r = 1(1)`nsim' {
				//Generate the dataset
				clear
				set seed 1
				
				qui set obs `studies_i'
				
				qui gen studyid = _n
				qui gen p = .
				qui gen total = .
				qui gen events = .
				local cobs = 1
				
				while `cobs' < `=`studies_i' + 1' {
					qui replace p = `op_i' in `cobs'
					qui replace total = runiformint(10, 500) in `cobs'
					qui replace events = rbinomial(total, p) in `cobs'
					local cobs = `cobs' + 1
				}
				
				//Fit the models
				qui {
				
					//beta-binomial - logit
					capture betabin events, n(total) link(logit) iterate(100)
					mat coef = r(table)
					local chi = e(chi2_c)
					local pvalue = chi2tail(1, `chi')/2
					
					local phat = invlogit(coef[1,1])
					local phatlo = invlogit(coef[5,1])
					local  phatup = invlogit(coef[6,1])
					local tauhat = coef[1,3] //sigma
					
					if _rc == 0 {
						local s = "`r'; betabin; `phat'; `phatlo'; `phatup'; `tauhat'; `op_i'; `tau_i'; `studies_i'"
						file write results "`s'" _n
					}
					//-4-------cloglog + fixed + wald	
					capture metapreg events total, model(fixed, sformat(%8.4f) ) ///
						studyid(study) cimethod(wilson, wald) dp(4) ///
						 noitable nograph sumtable(all) link(cloglog)

					mat absout = e(absout)
					mat hetout = e(hetout)

					local phat = absout[1,1]
					local phatlo = absout[1,5]
					local  phatup = absout[1,6]
					local tauhat = 0
					
					if _rc == 0 {
						local s = "`r'; cloglog-ce; `phat'; `phatlo'; `phatup'; `tauhat'; `op_i'; `tau_i'; `studies_i'"
						file write results "`s'" _n
					}
					
					//-3-------loglog + fixed + w	
					capture metapreg events total, model(fixed, sformat(%8.4f) ) ///
						studyid(study) cimethod(wilson, t) dp(4) ///
						 noitable nograph link(loglog)

					mat absout = e(absout)
					mat hetout = e(hetout)

					local phat = absout[1,1]
					local phatlo = absout[1,5]
					local  phatup = absout[1,6]
					local tauhat = 0
					
					if _rc == 0 {
						local s = "`r'; loglog-ce; `phat'; `phatlo'; `phatup'; `tauhat'; `op_i'; `tau_i'; `studies_i'"
						file write results "`s'" _n
					}
					
					//-2-------cloglog + random + wald	
					capture metapreg events total, model(random, sformat(%8.4f) ) ///
						studyid(study) cimethod(wilson, wald) dp(4) ///
						 noitable nograph link(cloglog)

					mat absout = e(absout)
					mat hetout = e(hetout)

					local phat = absout[1,1]
					local phatlo = absout[1,5]
					local  phatup = absout[1,6]
					local tauhat = hetout[1,4]
					if _rc == 0 {
						local s = "`r'; cloglog-re; `phat'; `phatlo'; `phatup'; `tauhat'; `op_i'; `tau_i'; `studies_i'"
						file write results "`s'" _n
					}
					
					//-1-------loglog + random + wald	
					capture metapreg events total, model(random, sformat(%8.4f) ) ///
						studyid(study) cimethod(wilson, wald) dp(4) ///
						 noitable nograph link(loglog)

					mat absout = e(absout)
					mat hetout = e(hetout)

					local phat = absout[1,1]
					local phatlo = absout[1,5]
					local  phatup = absout[1,6]
					local tauhat = hetout[1,4]
					if _rc == 0 {
						local s = "`r'; loglog-re; `phat'; `phatlo'; `phatup'; `tauhat'; `op_i'; `tau_i'; `studies_i'"
						file write results "`s'" _n
					}

					
					//1-------logistic + random + wald	
					capture metapreg events total, model(random, sformat(%8.4f) ) ///
						studyid(study) cimethod(wilson, wald) dp(4) 
						 
						 

					mat absout = e(absout)
					mat hetout = e(hetout)

					local phat = absout[1,1]
					local phatlo = absout[1,5]
					local  phatup = absout[1,6]
					local tauhat = hetout[1,4]
					if _rc == 0 {
						local s = "`r'; logistic-re; `phat'; `phatlo'; `phatup'; `tauhat'; `op_i'; `tau_i'; `studies_i'"
						file write results "`s'" _n
					}
					
					//1-------logistic + fixed + wald	
					capture metapreg events total, model(fixed, sformat(%8.4f) ) ///
						studyid(study) cimethod(wilson, wald) dp(4) ///
						 sumstat(ES) sumtable(abs) nograph
						 

					mat absout = e(absout)
					mat hetout = e(hetout)

					local phat = absout[1,1]
					local phatlo = absout[1,5]
					local  phatup = absout[1,6]
					local tauhat = hetout[1,4]
					if _rc == 0 {
						local s = "`r'; logistic-ce; `phat'; `phatlo'; `phatup'; `tauhat'; `op_i'; `tau_i'; `studies_i'"
						file write results "`s'" _n
					}

					//2-------exact ce	 
					capture metapreg events total, model(exact, sformat(%8.4f) ) ///
						studyid(study) cimethod(wilson, exact) dp(4) ///
						 sumstat(ES) sumtable(abs) nograph

					mat absout = e(absout)
					mat hetout = e(hetout)	 
					local phat = absout[1,1]
					local phatlo = absout[1,5]
					local phatup = absout[1,6]
					local tauhat = 0
					if _rc == 0 {					
						local s = "`r'; exact-ce; `phat'; `phatlo'; `phatup'; `tauhat'; `op_i'; `tau_i'; `studies_i'"
						file write results "`s'" _n 
					}

					//3-------IV + random
					 capture metaprop events total, ///
					 random nograph logit dp(4)

					local phat = r(ES)
					local phatlo = r(ci_low)
					local phatup = r(ci_upp)
					local tauhat = r(tau2) 	
					local s = "`r'; iv-re; `phat'; `phatlo'; `phatup'; `tauhat'; `op_i'; `tau_i'; `studies_i'"
					file write results "`s'" _n 
					 
					//4-------IV + fixed
					capture metaprop events total, ///
					 lcols(study) fixed nograph logit dp(4)
					 
					local phat = r(ES)
					local phatlo = r(ci_low)
					local phatup = r(ci_upp)
					local tauhat = r(tau2) 	
					local s = "`r'; iv-ce; `phat'; `phatlo'; `phatup'; `tauhat'; `op_i'; `tau_i'; `studies_i'"
					file write results "`s'" _n 

					//5-------IV + ftt (harmonic) + random
					capture metaprop events total, ftt ///
					 lcols(study) random nograph logit dp(4)
					 
					local phat = r(ES)
					local phatlo = r(ci_low)
					local phatup = r(ci_upp)
					local tauhat = r(tau2) 	
					local s = "`r'; iv-ftt-harm-re; `phat'; `phatlo'; `phatup'; `tauhat'; `op_i'; `tau_i'; `studies_i'"
					file write results "`s'" _n 
					 
					//6-------IV + ftt(harmonic) + fixed
					capture metaprop events total, ftt ///
					 lcols(study) fixed nograph logit dp(4)
					 
					local phat = r(ES)
					local phatlo = r(ci_low)
					local phatup = r(ci_upp)
					local tauhat = r(tau2) 	
					local s = "`r'; iv-ftt-harm-ce; `phat'; `phatlo'; `phatup'; `tauhat'; `op_i'; `tau_i'; `studies_i'"
					file write results "`s'" _n 

					//7------- IVHET + ftt (harmonic)
					capture metan events total,  ivhet pr ftt nograph
					local phat = r(prop_eff)
					local phatlo = r(prop_lci)
					local phatup = r(prop_uci)
					local tauhat = r(tausq) 	
					local s = "`r'; ivhet-ftt-harm; `phat'; `phatlo'; `phatup'; `tauhat'; `op_i'; `tau_i'; `studies_i'"
					file write results "`s'" _n 
					 

					//8------- IVHET + ftt (iv)
					capture metan events total, ivhet pr  tr(ft, iv) nograph
					local phat = r(prop_eff)
					local phatlo = r(prop_lci)
					local phatup = r(prop_uci)
					local tauhat = r(tausq) 	
					local s = "`r'; ivhet-ftt-iv; `phat'; `phatlo'; `phatup'; `tauhat'; `op_i'; `tau_i'; `studies_i'"
					file write results "`s'" _n 

					//9------- ce + ftt (iv)
					capture metan events total, fixed pr  tr(ft, iv)  nograph
					local phat = r(prop_eff)
					local phatlo = r(prop_lci)
					local phatup = r(prop_uci)
					local tauhat = 0 	
					local s = "`r'; ce-ftt-iv; `phat'; `phatlo'; `phatup'; `tauhat'; `op_i'; `tau_i'; `studies_i'"
					file write results "`s'" _n 

					//10------- re + ftt (iv)
					capture metan events total, random pr  tr(ft, iv)  nograph
					local phat = r(prop_eff)
					local phatlo = r(prop_lci)
					local phatup = r(prop_uci)
					local tauhat = r(tausq) 	
					local s = "`r'; re-ftt-iv; `phat'; `phatlo'; `phatup'; `tauhat'; `op_i'; `tau_i'; `studies_i'"
					file write results "`s'" _n 

					//11------- ce + ftt (arithmetic)
					capture metan events total, fixed pr  tr(ft, a)  nograph
					local phat = r(prop_eff)
					local phatlo = r(prop_lci)
					local phatup = r(prop_uci)
					local tauhat = 0 	
					local s = "`r'; ce-ftt-arith; `phat'; `phatlo'; `phatup'; `tauhat'; `op_i'; `tau_i'; `studies_i'"
					file write results "`s'" _n 

					//12------- IVHET + ftt (arithmetic)
					capture metan events total, ivhet pr  tr(ft, a)  nograph
					local phat = r(prop_eff)
					local phatlo = r(prop_lci)
					local phatup = r(prop_uci)
					local tauhat = r(tausq) 	
					local s = "`r'; ivhet-ftt-arith; `phat'; `phatlo'; `phatup'; `tauhat'; `op_i'; `tau_i'; `studies_i'"
					file write results "`s'" _n 

					//13------- re + ftt (arithmetic)
					capture metan events total, random pr  tr(ft, a)  nograph
					local phat = r(prop_eff)
					local phatlo = r(prop_lci)
					local phatup = r(prop_uci)
					local tauhat = r(tausq) 	
					local s = "`r'; re-ftt-arith; `phat'; `phatlo'; `phatup'; `tauhat'; `op_i'; `tau_i'; `studies_i'"
					file write results "`s'" _n 

					//14------- ce + ftt (geometric)
					capture metan events total, fixed pr  tr(ft, g)  nograph
					local phat = r(prop_eff)
					local phatlo = r(prop_lci)
					local phatup = r(prop_uci)
					local tauhat = 0 	
					local s = "`r'; ce-ftt-geom; `phat'; `phatlo'; `phatup'; `tauhat'; `op_i'; `tau_i'; `studies_i'"
					file write results "`s'" _n 

					//15------- IVHET + ftt (geometric)
					capture metan events total, ivhet pr  tr(ft, g)  nograph
					local phat = r(prop_eff)
					local phatlo = r(prop_lci)
					local phatup = r(prop_uci)
					local tauhat = r(tausq) 	
					local s = "`r'; ivhet-ftt-geom; `phat'; `phatlo'; `phatup'; `tauhat'; `op_i'; `tau_i'; `studies_i'"
					file write results "`s'" _n 

					//16------- re + ftt (geometric)
					capture metan events total, random pr  tr(ft, g)  nograph
					local phat = r(prop_eff)
					local phatlo = r(prop_lci)
					local phatup = r(prop_uci)
					local tauhat = r(tausq) 	
					local s = "`r'; re-ftt-geom; `phat'; `phatlo'; `phatup'; `tauhat'; `op_i'; `tau_i'; `studies_i'"
					file write results "`s'" _n 

					//17------- IVHET + arcsine
					capture metan events total, ivhet pr  tr(ar)  nograph
					local phat = r(prop_eff)
					local phatlo = r(prop_lci)
					local phatup = r(prop_uci)
					local tauhat = r(tausq) 	
					local s = "`r'; ivhet-arcsine; `phat'; `phatlo'; `phatup'; `tauhat'; `op_i'; `tau_i'; `studies_i'"
					file write results "`s'" _n 

					//18------- re+ arcsine
					capture metan events total, random pr  tr(ar)  nograph
					local phat = r(prop_eff)
					local phatlo = r(prop_lci)
					local phatup = r(prop_uci)
					local tauhat = r(tausq) 	
					local s = "`r'; re-arcsine; `phat'; `phatlo'; `phatup'; `tauhat'; `op_i'; `tau_i'; `studies_i'"
					file write results "`s'" _n 

					//19------- ce + arcsine
					capture metan events total, fixed pr  tr(ar)  nograph
					local phat = r(prop_eff)
					local phatlo = r(prop_lci)
					local phatup = r(prop_uci)
					local tauhat = 0 	
					local s = "`r'; ce-arcsine; `phat'; `phatlo'; `phatup'; `tauhat'; `op_i'; `tau_i'; `studies_i'"
					file write results "`s'" _n 
				
				}				
			}
			di "---"
		}
	
}

	 
file close results