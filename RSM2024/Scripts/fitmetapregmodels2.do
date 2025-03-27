//Fit the metapreg models
local designs "abnetwork comparative"				
local links "log"
local models "fixed mixed"
local varcors "unstructured independent commonslope commonint "
*local inferences "frequentist"
local inferences "frequentist bayesian"
*local inferences "bayesian"

foreach design of local designs {
	foreach inference of local inferences {
		//skip
		if "`design'" == "abnetwork" & "`inference'" == "bayesian" {
			continue
		}
		foreach model of local models {
						
			//skip 
			if "`inference'" == "bayesian" & "`model'" == "crbetabin" {
				continue
			}
				
			foreach link of local links {
			
				//Directory to save MCMC results
				if "`inference'" == "bayesian" {
					local bwd "bwd($wdir)"
				}
								
				local covar
				foreach cov of local varcors {
									
					//Loop once and skip
					if "`design'" == "abnetwork" & "`cov'" != "independent" {
						continue
					}
					
					//Loop once and skip the rest 
					if "`model'" == "crbetabin" & "`cov'"  != "commonslope" {
						continue
					}
					
					//Loop twice and skip
					if "`model'" == "fixed"  & ("`cov'"  == "independent" | "`cov'"  == "unstructured" | "`cov'"  == "commonslope") {
						continue
					}
										
					if "`model'" == "mixed" & "`design'" == "comparative" {
						local covar ", cov(`cov')"
					}
					
					/*
					if "`inference'" == "bayesian" & "`cov'" != "unstructured" & "`model'" == "mixed" {
						local nt = 3
					}
					else
					*/
					if "`model'" == "crbetabin" {
						local nt = 2
					}
					else {
						local nt = 1
					}
					
					local varprior
					forvalues t=1(1)`nt' {
						/*if "`inference'" == "bayesian" {
							if `t' == 2 {
								local varprior ", varprior(igamma(0.001, 0.001))"
							}
							if `t' == 3 {
								local varprior ", varprior(jeffreys)"
							}
						} */
						
						local tsvar 
						if "`model'" == "crbetabin" {
							local cov "commonslope"
							
							if `t' == 2 {
								local tsvar "studyid"
								local cov "randomslope"
							}
						}
										
						cap metapreg event total group `tsvar' ,  ///
							model(`model' `varprior') link(`link') inference(`inference') `bwd' ///
							studyid(studyid) design(`design' `covar') nograph sumtable(none) noitable
						
						/*
						set more off
						
						cap log close
						log using "C:\DATA\WIV\Projects\Stata\Metapreg\Logs\doczone.txt", text replace 
						set trace off

						metapreg event total group, progress ///
						model(mixed) link(logit)  ///
						studyid(studyid) design(comparative, cov(unstructured)) 					
						*/	
															
						if _rc == 0 {
							mat poporout = e(poporout)
							mat orout = e(orout)
							mat poprrout = e(poprrout)
							mat rrout = e(rrout)
							mat poprdout = e(poprdout)
							mat rdout = e(rdout)
							mat absout = e(absout)
							mat popabsout = e(popabsout)
							mat hetout = e(hetout)
							mat covmat = e(covmat)
							mat gof = e(gof)
							
							if "$test"  == "yes" & "`inference'" == "bayesian" {
								qui {
									estimates restore metapreg_modest
									estimates store run$run
									
									if $run != 1 {
										bayesstats ic run1 run$run
										mat ic = r(ic)
										local logBF = ic[2, 4]
									}
									else {
										local logBF = .
									}
								}
							}
							
							//Information criterion
							if "`inference'" == "bayesian" {
								local DIC = gof[1, 1]
								local IC = "`DIC'" + "-" + "`logBF'" + "-"  + "$run"
								global run = $run + 1
							}
							else {
								local BIC = gof[1, 2]
								local AIC = gof[1, 1]
								local IC = "`AIC'" + "-" + "`BIC'"
							}
													
							//Default between-study variances
							local tau2hat = .
							local sigma2hat = .
							
							local tau2hatlo = .
							local sigma2hatlo = .
							
							local tau2hatup = .
							local sigma2hatup = .
																			
							//Get between-study variances								
							if "`model'" == "mixed" {							
								if "`inference'" == "bayesian" {
									local statcol 4 //median var
									local statlo 5
									local statup 6
								}
								else {
									local statcol 1 //mean var
									local statlo 2
									local statup 3
								}
								if "`cov'" != "commonint" {
									local tau2hat = covmat[1, `statcol']
									local tau2hatlo = covmat[1, `statlo']
									local tau2hatup = covmat[1, `statup']
									
								}
								if "`cov'" != "commonslope" { 
									local sigma2hat = covmat[2, `statcol']
									local sigma2hatlo = covmat[2, `statlo']
									local sigma2hatup = covmat[2, `statup']
								}
								
								if "`cov'" == "commonint" { 
									local tau2hat = 0
									local tau2hatlo = 0
									local tau2hatup = 0								
								}
								
								if "`cov'" == "commonslope" { 
									local sigma2hat = 0
									local sigma2hatlo = 0
									local sigma2hatup = 0
								}
							}
							else if "`model'" == "crbetabin" {
								local tau2hat = covmat[1,1]
								local sigma2hat = covmat[1,1]
								
								local tau2hatlo = covmat[1,2]
								local sigma2hatlo = covmat[1,2]
								
								local tau2hatup = covmat[1,3]
								local sigma2hatup = covmat[1,3]	
							}
							else {
								local tau2hat = 0
								local sigma2hat = 0
								
								local tau2hatlo = 0
								local sigma2hatlo = 0
								
								local tau2hatup = 0
								local sigma2hatup = 0
							}
																
							local popmatrices "popabsout poprrout poporout poprdout"
							
							foreach popmat of local popmatrices {
								
								if "`design'" == "abnetwork"  & "`popmat'" != "popabsout" {
									local rowid = 2
								}
								else {
									local rowid = 1
								}
								
								local esthatlo = `popmat'[`rowid', 4]
								local esthatup = `popmat'[`rowid', 5]
															
								forvalues l = 1/2 {
									if `l' == 1 {
										local prefix "mean"
										local esthat = `popmat'[`rowid', 1]
									}
									else {
										local prefix "median"
										local esthat = `popmat'[`rowid', 3]
									}
									
									local s = "$seed; metapreg-`inference'-`design'-`model'-`t'-`cov'-`link'; `prefix'-`popmat'; `esthat'; `esthatlo'; `esthatup'; `tau2hat'; `tau2hatlo'; `tau2hatup';`sigma2hat';`sigma2hatlo';`sigma2hatup'; `IC'; $mu0_i; $or_i; $tausq_i; $sigmasq_i; $studies; $studies; $studysize_i"
									file open results using $outfile , text write append //append 
									file write results "`s'" _n
									file close results
								}
							}
							
							local condmatrices "absout rrout orout rdout"
							
							foreach condmat of local condmatrices {
								if "`design'" == "abnetwork"   & "`condmat'" != "absout" {
									local rowid = 2
								}
								else {
									local rowid = 1
								}

								
								if "`inference'" == "frequentist" {
									local nloops = 1
								}
								else {
									local nloops = 2
								}
								
								forvalues l=1/`nloops' {
									if `l' == 1 {
										local prefix "mean"
										local esthat = `condmat'[`rowid', 1]
									}
									else {
										local prefix "median"
										local esthat = `condmat'[`rowid', 4]
									}
																								
									forvalues k=1/2 {
										if `k' == 1 {
											//eti - z
											local esthatlo = `condmat'[`rowid', 5]
											local esthatup = `condmat'[`rowid', 6]
											if "`inference'" == "bayesian" {
												local citype "eti"
											}
											else {
												local citype "z"
											}
										}
										else {
											if "`inference'" == "bayesian" {
												//hpd
												local esthatlo = `condmat'[`rowid', 7]
												local esthatup = `condmat'[`rowid', 8]
												local citype "hpd"
											}
											else {
												//t
												local esthatlo = `condmat'[`rowid', 8]
												local esthatup = `condmat'[`rowid', 9]
												local citype "t"
											}
										}
										
										local s = "$seed; metapreg-`inference'-`design'-`model'-`t'-`cov'-`link'-`citype'; `prefix'-`condmat'; `esthat'; `esthatlo'; `esthatup'; `tau2hat';`tau2hatlo'; `tau2hatup';`sigma2hat'; `sigma2hatlo';`sigma2hatup'; `IC'; $mu0_i; $or_i; $tausq_i; $sigmasq_i; $studies; $studies; $studysize_i"
										file open results using $outfile , text write append //append 
										file write results "`s'" _n
										file close results	
									}
								}										
							}
						}
					}
				}												
			}
		}
	}
}