local models "peto iv mh dl bdl he mp ml reml hm b0 bp sj2s dk2s hk pl kr bt mu ivhet"
local correction "zeroonly"

foreach model of local models {
	set more off
	capture metan event1 noevent1 event0 noevent0, ///
		or model(`model') study(studyid)  nograph `correction'
	/*	
	metan event1 noevent1 event0 noevent0, ///
		rd model(ivhet) study(studyid)  keepall	 forestplot(ysize(7.5) xsize(10) graphregion(color(white)))
	*/
	
	if _rc == 0 {
		mat ovstats = r(ovstats)
		local IC = .
		local k = r(k)
		
		//Default between-study variances
		local tau2hat = .
		local sigma2hat = .
		
		local tau2hatlo = .
		local sigma2hatlo = .
		
		local tau2hatup = .
		local sigma2hatup = .	
		
		local esthat = exp(ovstats[1, 1])
		local esthatlo = exp(ovstats[3, 1])
		local esthatup = exp(ovstats[4, 1])
		
				
		if "`model'" == "iv" | "`model'" == "mh" |  "`model'" == "peto"  {
			local sigma2hat = 0
		}
		else {
			local sigma2hat = ovstats[9, 1]
		}
		
		local s = "$seed; metan-frequentist-`model'-1-`tcc'-logit; mean-orout; `esthat'; `esthatlo'; `esthatup'; `tau2hat'; `tau2hatlo'; `tau2hatup';`sigma2hat';`sigma2hatlo';`sigma2hatup'; `IC'; $mu0_i; $or_i; $tausq_i; $sigmasq_i; `k' ; $studies; $studysize_i"
								
		file open results using $outfile , text write append //append 
		file write results "`s'" _n
		file close results
	}		
}