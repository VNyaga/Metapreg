	capture meta esize event1 noevent1 event0 noevent0, es(lnor) studylabel(studyid)

	local remethods "reml mle eb dl sj he hs"
	foreach method of local remethods {
	
		//Default between-study variances
		local tau2hat = .
		local sigma2hat = .
		
		local tau2hatlo = .
		local sigma2hatlo = .
		
		local tau2hatup = .
		local sigma2hatup = .
		local IC = .
						
		//====================
		cap meta summarize, random(`method') 
		
		/* 
			capture meta esize event1 noevent1 event0 noevent0, es(lnor) studylabel(studyid)
			meta summarize, random(sj) eform(or)
			meta forestplot, random(sj) eform(or)  xline(1) esrefline se(kh)
			
			*rr
			capture meta esize event1 noevent1 event0 noevent0, es(lnrr) studylabel(studyid)
			meta summarize, random(sj) eform(rr)
			meta forestplot, random(sj) eform(rr)  xline(1) esrefline se(kh)
			
			*rd
			capture meta esize event1 noevent1 event0 noevent0, es(rd) studylabel(studyid)
			meta summarize, random(sj) 
			meta forestplot, random(sj)   xline(0) esrefline se(kh)			
		*/
		if _rc == 0 {
		local esthat = exp(r(theta))
		local esthatlo = exp(r(ci_lb))
		local esthatup = exp(r(ci_ub))
		local sigma2hat = r(tau2)
		local k =r(N)
				
		local s = "$seed; smeta-frequentist-`method'-1-wald-logit; mean-orout; `esthat'; `esthatlo'; `esthatup'; `tau2hat'; `tau2hatlo'; `tau2hatup';`sigma2hat';`sigma2hatlo';`sigma2hatup'; `IC'; $mu0_i; $or_i; $tausq_i; $sigmasq_i; `k'; $studies; $studysize_i"
								
		file open results using $outfile , text write append //append 
		file write results "`s'" _n
		file close results		
		}
		//====================
		cap meta summarize, random(`method') tdist
		if _rc == 0 {
		local esthat = exp(r(theta))
		local esthatlo = exp(r(ci_lb))
		local esthatup = exp(r(ci_ub))
		local sigma2hat = r(tau2)
		local k =r(N)
				
		local s = "$seed; smeta-frequentist-`method'-1-t-logit; mean-orout; `esthat'; `esthatlo'; `esthatup'; `tau2hat'; `tau2hatlo'; `tau2hatup';`sigma2hat';`sigma2hatlo';`sigma2hatup'; `IC'; $mu0_i; $or_i; $tausq_i; $sigmasq_i; `k'; $studies; $studysize_i"
		file open results using $outfile , text write append //append 
		file write results "`s'" _n
		file close results
		}
		//====================
		cap meta summarize, random(`method') se(kh)
		if _rc == 0 {
		local esthat = exp(r(theta))
		local esthatlo = exp(r(ci_lb))
		local esthatup = exp(r(ci_ub))
		local sigma2hat = r(tau2)
		local k =r(N)
				
		local s = "$seed; smeta-frequentist-`method'-1-tkh-logit; mean-orout; `esthat'; `esthatlo'; `esthatup'; `tau2hat'; `tau2hatlo'; `tau2hatup';`sigma2hat';`sigma2hatlo';`sigma2hatup'; `IC';  $mu0_i; $or_i; $tausq_i; $sigmasq_i; `k'; $studies; $studysize_i"
		file open results using $outfile , text write append //append 
		file write results "`s'" _n
		file close results
		}
		cap meta summarize, random(`method') se(kh,trunc)
		if _rc == 0 {
		local esthat = exp(r(theta))
		local esthatlo = exp(r(ci_lb))
		local esthatup = exp(r(ci_ub))
		local sigma2hat = r(tau2)
		local k =r(N)
				
		local s = "$seed; smeta-frequentist-`method'-1-tkht-logit; mean-orout; `esthat'; `esthatlo'; `esthatup'; `tau2hat'; `tau2hatlo'; `tau2hatup';`sigma2hat';`sigma2hatlo';`sigma2hatup'; `IC';  $mu0_i; $or_i; $tausq_i; $sigmasq_i; `k'; $studies; $studysize_i"
		
		file open results using $outfile , text write append //append 
		file write results "`s'" _n
		file close results
		}
	}

	local femethods "mh iv"
	foreach method of local femethods {
		cap meta summarize, fixed(`method')
		if _rc == 0 {
		local esthat = exp(r(theta))
		local esthatlo = exp(r(ci_lb))
		local esthatup = exp(r(ci_ub))
		local sigma2hat = 0
		local k =r(N)
				
		local s = "$seed; smeta-frequentist-`method'-1-wald-logit; mean-orout; `esthat'; `esthatlo'; `esthatup'; `tau2hat'; `tau2hatlo'; `tau2hatup';`sigma2hat';`sigma2hatlo';`sigma2hatup'; `IC';  $mu0_i; $or_i; $tausq_i; $sigmasq_i; `k'; $studies; $studysize_i"
		file open results using $outfile , text write append //append 
		file write results "`s'" _n
		file close results
		}
	}	

	cap meta summarize, fixed(iv) tdist
	if _rc == 0 {
	local esthat = exp(r(theta))
	local esthatlo = exp(r(ci_lb))
	local esthatup = exp(r(ci_ub))
	local sigma2hat = 0
	local k =r(N)
				
	local s = "$seed; smeta-frequentist-iv-1-t-logit; mean-orout; `esthat'; `esthatlo'; `esthatup'; `tau2hat'; `tau2hatlo'; `tau2hatup';`sigma2hat';`sigma2hatlo';`sigma2hatup'; `IC';  $mu0_i; $or_i; $tausq_i; $sigmasq_i; `k'; $studies; $studysize_i"
	file open results using $outfile , text write append //append 
	file write results "`s'" _n
	file close results
	}