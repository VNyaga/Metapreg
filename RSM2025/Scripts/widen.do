*	-------------------------REAL DATA
{
	//All fits	
	if "$review" == "hemkens2" {
		global fmt "e"
		global txts "1.75"
	}
	else {
		global fmt "f"
		global txts "2.35"
	}
	gen subset = ((stat=="median-orout" & Inf=="B") | (stat=="mean-orout" & Inf=="F")) 
	replace subset =0 if inlist(Link, "L", "LLG", "CLLG")
	
	gen model1 = package + " " +  inference + " "  +  Design + " " + param + Covariance +  Weighting + " " + Sigmethod +  Prior +  " " +  Env + " "  +  Dist  
	gsort subset model1 CImethod
	bys subset model1: egen nest = seq()
	
	drop model sigma2hatlo sigma2hatup v4 v5 v6 v8 WidthCI
	drop if !subset
	reshape wide esthat esthatlo esthatup CImethod ciwidth, i(model1) j(nest) 
	
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
}

gsort -esthat1 -ciwidth1	
	

	*gen str Est =  string(esthat, "%10.2e") + " [" + string(esthatlo, "%10.2e") + ", " + string(esthatup, "%10.2e") + "]"  
	

replace tau2hat =. if inlist(Dist, "BB1", "BB2")
replace sigma2hat =. if inlist(Dist, "BB1", "BB2")

replace Covariance = " " if inlist(Dist, "BB1", "BB2")
