//All fits	
	if "$review" == "hemkens2" {
		global fmt "e"
		global txts "2.15"
	}
	else {
		global fmt "f"
		global txts "2.5"
	}
	gen subset = ((stat=="median-orout" & Inf=="B") | (stat=="mean-orout" & Inf=="F")) 
	replace subset =0 if inlist(Link, "L", "LLG", "CLLG")
	
	gen model1 = package + " " +  inference + " "  +  Design + " " + param + Covariance +  Weighting + " " + Sigmethod +  Prior +  " " +  Env + " "  +  Dist  
	gsort subset model1 CImethod
	bys subset model1: egen nest = seq()
	
	drop model sigma2hatlo sigma2hatup v4 v5 v6 v8 WidthCI
	drop if !subset
	reshape wide esthat esthatlo esthatup CImethod ciwidth , i(model1) j(nest) 
	
	gsort -esthat1 -ciwidth1	
	
	gen str CI1 =  CImethod1 + ":" +  " (" + string(esthatlo1, "%10.2$fmt") + ", " + string(esthatup1, "%10.2$fmt") + ")  " 
	gen str CI2=  CImethod2 + ":" +  " (" + string(esthatlo2, "%10.2$fmt") + ", " + string(esthatup2, "%10.2$fmt") + ")  " if !mi(CImethod2)
	gen str CI3 =  CImethod3 + ":" +  " (" + string(esthatlo3, "%10.2$fmt") + ", " + string(esthatup3, "%10.2$fmt") + ")  " if !mi(CImethod3)
	
	tostring esthat1, gen(OR) format(%4.2f) force 
	
	label var OR "OR"
	label var CI1 "95% CI(black)"
	label var CI2 "95% CI (orange)"
	label var CI3 "95% CI (cyan)"
	
	label var CImethod1 "CImethod"
	
	
//meta
{
blobbogram Dist esthat1 esthatlo1 esthatup1 if  ///
	package == "meta" & Env=="R" & stat=="mean-orout" ,  ///	
	sumstat(Mean OR) lcols( Weighting Slope Sigmethod sigma2hat ) rcols(OR CI1 CI2 bratio )  graphregion(fcolor(white))   ///
	xlab($orrange)	texts(2.5)  astext(80) xline($orline)  ysize(12.5) xsize(20)  grid  ///
	optimalmat(optimalor) addtextleft($text2add) addtextright($text2addright) ///
	 name(metaobs, replace) logscale nostats superimpose1(esthat2 esthatlo2 esthatup2)
	 
	 graph export "$graphs\meta.png", as(png) width(1750) height(1000) replace
 
}	 
//bayesmeta - bayesian
{
blobbogram Dist esthat1 esthatlo1 esthatup1 if stat=="median-orout" &  package == "bayesmeta",  ///	
	sumstat(Median OR) lcols(Prior  sigma2hat )  rcols(OR CI1 CI2 bratio )  grid  ///
	xlab($orrange)	texts(2.5)  astext(80) xline($orline)  dp(2) ysize(12.5) xsize(20)   ///
	optimalmat(optimalor) addtextleft($text2add) addtextright($text2addright) ///
	 name(bayesmetaobs, replace)  logscale nostats superimpose1(esthat2 esthatlo2 esthatup2)
	 
	 graph export "$graphs\bayesmeta.png", as(png) width(1750) height(1000) replace

}	 

//metafor 
{

blobbogram Dist esthat1 esthatlo1 esthatup1 if stat=="mean-orout" &  package == "metafor",  ///	
	sumstat(Mean OR) lcols(Covariance tau2hat sigma2hat)  rcols(OR  CI1 CI2 bratio )   grid  ///
	xlab($orrange)	texts(2.5)  astext(80) xline($orline)  dp(2) ysize(12.5) xsize(20)   ///
	optimalmat(optimalor) addtextleft($text2add) addtextright($text2addright) textopts(place(e) color(red) size(2.5)) ///
	 name(metaforobs, replace)  logscale nostats superimpose1(esthat2 esthatlo2 esthatup2)
	 
	graph export "$graphs\metafor.png", as(png) width(1750) height(1000) replace	 

}	 
//metan 
{
 
blobbogram Dist esthat1 esthatlo1 esthatup1 if  package == "metan" &  stat=="mean-orout",  ///	
	sumstat(Mean OR) lcols(Weighting Slope Sigmethod sigma2hat ) rcol( CImethod1 bratio )  graphregion(fcolor(white))  ///
	xlab($orrange)	texts(2.5)  astext(80) xline($orline)  fxsize(200) fysize(120)  grid  ///
	optimalmat(optimalor) addtextleft($text2add) addtextright($text2addright) ///
	 name(metanobs, replace) logscale
	 
graph export "$graphs\metan.png", as(png) width(1750) height(1000) replace	

	 
}

//metaplus
{

blobbogram Dist esthat1 esthatlo1 esthatup1 if  package == "metaplus" &  stat=="mean-orout",  ///	
	sumstat(Mean OR) lcols(sigma2hat CImethod1 ) rcols(bratio )  graphregion(fcolor(white))  ///
	xlab($orrange)	texts(2.5)  astext(80) xline($orline)  ysize(12.5) xsize(20)  grid ///
	optimalmat(optimalor) addtextleft($text2add) addtextright($text2addright) textopts(place(e) color(red) size(2.5)) ///
	 name(metaplusobs, replace) logscale
	 
graph export "$graphs\metaplus.png", as(png) width(1750) height(1000) replace
}

//metastan - bayesian
{
blobbogram Dist esthat1 esthatlo1 esthatup1 if stat=="median-orout" &  package == "metastan",  ///	
	sumstat(Median OR) lcols(param Prior sigma2hat )  rcols(OR  CI1 CI2 bratio) grid   ///
	xlab($orrange)	texts(2.5)  astext(80) xline($orline) ysize(12.5) xsize(20)   ///
	optimalmat(optimalor) addtextleft($text2add) addtextright($text2addright) textopts(place(e) color(red) size(2.5)) ///
	 name(metastanobs, replace) logscale superimpose1(esthat2 esthatlo2 esthatup2) nostats

	 
	graph export "$graphs\metastan.png", as(png) width(1750) height(1000) replace	

}

	 
//smeta	
{
*tkht and tkh are the same

blobbogram Dist esthat1 esthatlo1 esthatup1 if stat=="mean-orout" &  package == "meta" & Env=="Stata" /*& inlist(CImethod, "z", "tkh", "t")*/  ,  ///	
	sumstat(Mean OR) lcols(Slope Weighting Sigmethod /*CImethod*/ sigma2hat ) rcols(OR  CI1 CI2 CI3 bratio) grid  ///
	xlab($orrange)	texts($txts)  astext(80) xline($orline)  ysize(12.5) xsize(20)   ///
	optimalmat(optimalor) addtextleft($text2add) addtextright($text2addright) textopts(place(e) color(red) size($txts)) nostats ///
	 name(smetaobs, replace) logscale  superimpose1(esthat2 esthatlo2 esthatup2) superimpose2(esthat3 esthatlo3 esthatup3)

	 
	graph export "$graphs\smeta.png", as(png) width(1750) height(1000) replace	

}