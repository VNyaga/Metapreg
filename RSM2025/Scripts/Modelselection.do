//Data management
	if "$review" == "hemkens2016b" {
		global fmt "e"
		global txts "1.75"
	}
	else {
		global fmt "f"
		global txts "2.35"
	}
	tostring tau2hat, gen(tau2) format(%4.2f) force 

	gen model1 = package + " " +  inference + " "  +  Design +  " " +  v4 +  " " +  v5 + " "  +  v6  +  " " + v7 
	gen statid = 1 if stat == "median-poprrout" 

	gsort model1 statid 
	bys model1: gen holder = esthat[1] 	
	bys model1: gen holderlo = esthatlo[1] 	
	bys model1: gen holderup = esthatup[1] 

	gen popstat = holder if ((stat == "median-rrout" & Inf == "B")|(stat == "mean-rrout" & Inf == "F")) 
	gen popstatlo = holderlo if ((stat == "median-rrout" & Inf == "B")|(stat == "mean-rrout" & Inf == "F")) 
	gen popstatup= holderup if ((stat == "median-rrout" & Inf == "B")|(stat == "mean-rrout" & Inf == "F")) 
	
	drop statid holder holderlo holderup
	gen statid = 1 if stat == "median-poporout" 

	gsort model1 statid 
	bys model1: gen holder = esthat[1] 	
	bys model1: gen holderlo = esthatlo[1] 	
	bys model1:gen  holderup = esthatup[1] 	 
	
	replace popstat = holder if (stat == "median-orout" & Inf == "B")|(stat == "mean-orout" & Inf == "F") 
	replace popstatlo = holderlo if (stat == "median-orout" & Inf == "B")|(stat == "mean-orout" & Inf == "F") 
	replace popstatup= holderup if (stat == "median-orout" & Inf == "B")|(stat == "mean-orout" & Inf == "F") 
	

	gen str est =  string(esthat, "%10.2$fmt") +  " (" + string(esthatlo, "%10.2$fmt") + ", " + string(esthatup, "%10.2$fmt") + ")"  
	gen str popest =  string(popstat, "%10.2$fmt") +  " (" + string(popstatlo, "%10.2$fmt") + ", " + string(popstatup, "%10.2$fmt") + ")"  
	
	label var est "SPE (95%CI)"
	label var popest "PA (95%CI)"
	
*=========================================	LOG BF 	=========================================	
	
	gsort -logBF
	 
	blobbogram Dist esthat esthatlo esthatup if ///
		stat=="median-orout" & inference !="frequentist" & link == "logit" & CImethod == "eti"  & ///
		strpos(Dist, "B") != 0  &  package == "metapreg" ,  ///	
		sumstat(Median OR) lcols(Prior  Covariance logBF ) rcols(est popest)    logscale superimpose1(popstat popstatlo popstatup) ///
		xlab($allrrange)	texts($txts)  astext(75) xline(1)  fysize(120) fxsize(175)  grid nostats ///
		name(CORestimatesbw, replace) 		 


	blobbogram tau2 esthat esthatlo esthatup if ///
		stat=="median-rrout" & inference !="frequentist" & link == "logit" & CImethod == "eti"  & ///
		strpos(Dist, "B") != 0 &  package == "metapreg" ,  ///	
		sumstat(Median RR) rcols(est popest ) lcols( sigma2hat)    logscale superimpose1(popstat popstatlo popstatup) ///
		xlab($allrrange)	texts($txts)  astext(75) xline(1)  fysize(120) fxsize(150)  grid nostats   ///
		name(CRRestimatesbw, replace) 		 
					 
		
	gr combine  CORestimatesbw CRRestimatesbw, ///
			graphregion(color(white)) cols(2) ycommon imargin(0 2 0 0) ///
			 name(CORR, replace)	 xsize(17.5) ysize(10)	
		
	
	graph export "$graphs\allfits-B.png", as(png) width(2000) height(1000) replace	
	
	
*=========================================	BIC 	=========================================	 

	gsort BIC

	blobbogram Dist esthat esthatlo esthatup if ///
		stat=="mean-orout" & inference =="frequentist" & link == "logit" & ///
		CImethod == "t" &  inlist(Dist, "B1",  "BN") & inlist(Design,  "comparative", "general") & ///
		strpos(Dist, "B") != 0 /*&  Prior != "IG(0.001, 0.001)" & Prior != "Jeffreys"*/  &  package == "metapreg" ,  ///	
		sumstat(OR) lcols(Covariance BIC ) rcols(est popest)    logscale superimpose1(popstat popstatlo popstatup) ///
		xlab($allrrange)	texts($txts)  astext(75) xline(1)  fysize(120) fxsize(175)  grid nostats ///
		name(CORestimatesbw, replace) 		 


	blobbogram tau2 esthat esthatlo esthatup if ///
		stat=="mean-rrout" & inference =="frequentist" & link == "logit" & ///
		CImethod == "t"  &  inlist(Dist, "B1",  "BN") & inlist(Design,  "comparative", "general") & ///
		strpos(Dist, "B") != 0 /*&  Prior != "IG(0.001, 0.001)" & Prior != "Jeffreys"*/  &  package == "metapreg" ,  ///	
		sumstat(RR) rcols(est popest ) lcols( sigma2hat)    logscale superimpose1(popstat popstatlo popstatup) ///
		xlab($allrrange)	texts($txts)  astext(75) xline(1)  fysize(120) fxsize(150)  grid nostats   ///
		name(CRRestimatesbw, replace) 		 
					 
					 
	gr combine  CORestimatesbw CRRestimatesbw, ///
			graphregion(color(white)) cols(2) ycommon imargin(0 2 0 0) ///
			 name(CORR1, replace) xsize(15) ysize(10)	


	graph export "$graphs\allfits-F.png", as(png) width(2000) height(1000) replace	
	
				