
*gen width = "  " + string(ciwidth, "%4.2f") 

*blobbogram Dist esthat esthatlo esthatup /*testhat testhatlo testhatup*/  if  ///
*	(((stat=="median-orout" & Inf=="B") | ///
*	(stat=="mean-orout" & Inf=="F")) & ///
*	(conditional|pick)),  ///	
*	sumstat(OR) lcols(Inf package Env param Weighting Sigmethod CImethod Prior ) rcols( /*Est*/ tau2hat sigma2hat ciwidth /*WidthCI*/ )  ///
*	xlab($simorrange)	texts(1.75)  astext(80) xline(1)  fysize(200) fxsize(180)  grid /*nostats*/ logscale ///
*	name(bestmodelsbw, replace) 

	
blobbogram Dist esthat1 esthatlo1 esthatup1 /*testhat testhatlo testhatup*/  if  ///
	(((stat=="median-orout" & Inf=="B") | ///
	(stat=="mean-orout" & Inf=="F")) & ///
	(conditional|pick)),  ///	
	sumstat(OR) lcols(Inf package Env param Weighting Covariance Sigmethod Prior  tau2hat sigma2hat)   ///
	xlab($orrange)	texts(1.75)  astext(80) xline(1)  ysize(10) xsize(20)  grid  logscale ///
	name(bestmodels, replace)	superimpose1(esthat2 esthatlo2 esthatup2) rcols(OR bratio fratio CI1 CI2) nostats

	
if "$review" == "hemkens2" {	
	blobbogram Dist esthat1 esthatlo1 esthatup1 /*testhat testhatlo testhatup*/  if  ///
		(((stat=="median-orout" & Inf=="B") | ///
		(stat=="mean-orout" & Inf=="F")) & ///
		(conditional|pick)),  ///	
		sumstat(OR) lcols(Inf package Env param Weighting Covariance Sigmethod Prior  tau2hat sigma2hat)   ///
		xlab($orrange)	texts(1.75)  astext(90) xline(1)  fysize(150) fxsize(200)  grid  logscale ///
		name(bestmodelsbw, replace)	superimpose1(esthat2 esthatlo2 esthatup2) rcols(OR CI1 CI2) nostats	
	}
	else {
	blobbogram Dist esthat1 esthatlo1 esthatup1 /*testhat testhatlo testhatup*/  if  ///
		(((stat=="median-orout" & Inf=="B") | ///
		(stat=="mean-orout" & Inf=="F")) & ///
		(conditional|pick)),  ///	
		sumstat(OR) lcols(Inf package Env param Weighting Covariance Sigmethod Prior  tau2hat sigma2hat)   ///
		xlab($orrange)	texts(1.75)  astext(80) xline(1)  fysize(200) fxsize(180)  grid  logscale ///
		name(bestmodelsbw, replace)	superimpose1(esthat2 esthatlo2 esthatup2) rcols(OR CI1 CI2) nostats
	}
		
	graph export "$graphs\bestmodelsbw.png", as(png) width(2000) height(1250) replace	





