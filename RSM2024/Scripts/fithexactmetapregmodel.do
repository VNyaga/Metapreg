//Fit the metapreg models
local designs "comparative"				
local links "logit"
local models "hexact"
local inferences "frequentist"

foreach design of local designs {
	foreach inference of local inferences {
		foreach model of local models {			
			foreach link of local links {			
												
				cap metapreg event total group ,  ///
					model(`model') link(`link') inference(`inference')  ///
					studyid(studyid) design(`design' ) nograph sumtable(none) noitable
				
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
					mat orout = e(exactorout)
					mat absout = e(exactabsout)
					mat gof = .
					
					local BIC = .
					local AIC = .
					local IC = "`AIC'" + "-" + "`BIC'"
											
					//Default between-study variances
					local tau2hat = 0
					local sigma2hat = 0
					
					local tau2hatlo = 0
					local sigma2hatlo = 0
					
					local tau2hatup = 0
					local sigma2hatup = 0
																						
					local condmatrices "absout orout"
					
					foreach condmat of local condmatrices {
						
						local rowid = 1
						if "`condmat'" == "absout" {
							local addcol = 2
						}
						else {
							local addcol = 0
						}
						
						local esthat = `condmat'[`rowid', 1]
						local esthatlo = `condmat'[`rowid', `=3+`addcol'']
						local esthatup = `condmat'[`rowid', `=4+`addcol'']
							
						local prefix "mean"
																																
						local s = "$seed; metapreg-`inference'-`design'-`model'-1-zero-`link'; `prefix'-`condmat'; `esthat'; `esthatlo'; `esthatup'; `tau2hat';`tau2hatlo'; `tau2hatup';`sigma2hat'; `sigma2hatlo';`sigma2hatup'; `IC'; $mu0_i; $or_i; $tausq_i; $sigmasq_i; $studies; $studies; $studysize_i"
						file open results using $outfile , text write append //append 
						file write results "`s'" _n
						file close results
																
					}
				}													
			}
		}
	}
}