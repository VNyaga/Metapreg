/*
CREATED:	18 Apr 2024
AUTHOR:		Victoria N Nyaga
PURPOSE: 	Draw a forest plot given es lci es weights and matrix with overall stats
VERSION: 	1.0.0
NOTES:
Re-uses code from metapplot
*/
capture program drop blobbogram
	program define blobbogram

	#delimit ;
	syntax varlist [if] [in] [,
		POWer(integer 0)
		DP(integer 2) 
		Level(integer 95)
		weights(string)
		smooth
		noWT
		ovmat(name)
		optimalmat(name)
		cicolor(varname)
		
		/*Passed through*/
		logscale
		AStext(integer 50)
		ARRowopt(string) 		
		CIOpts(string) 
		DIAMopts(string) 
		DOUble 
 		LCols(varlist)
		RCols(varlist) 		
		noOVerall 
		noOVLine 		
		noSTATS
		noBox
		OLineopts(string) 
		SUMStat(string asis) 
		POINTopts(string) 
		BOXopts(string) 
		predciopts(string)
		SUBLine 
		TEXts(real 1.0) 
		XLAbel(string) 
		XLIne(string) 
		XTick(string)
		GRID
		GRAphsave(string asis)
		prediction
		band1(string asis)
		band2(string asis)
		band3(string asis)
		bands(integer 0)
		base(real 0)
		addtextleft(string asis)
		addtextright(string asis)
		textopts(string asis)
		groupvar(varname)
		compabs
		rowid(varname)
		ses(varname)
		superimpose1(varlist)
		superimpose2(varlist)
		*
	  ];
	#delimit cr
	
	preserve

	marksample touse, strok  novarlist
	markout `touse', sysmissok
	
	qui keep if `touse'
	
	local plotopts `"`options'"'
	
	if strpos(`"`plotopts'"', "graphregion") == 0 {
			local plotopts `"graphregion(color(white)) `plotopts'"'
	}
	
	tempvar es modeles lci modellci uci modeluci lpi upi ilci iuci predid use label tlabel id newid rid gid df expand expanded order orig flag ///
	
	tokenize "`varlist'", parse(" ")
	
	qui {
		gen str `label'	=`1'
		
		gen `es'		=`2'*(10^`power')
		
		if "`ses'" != "" {
			tempvar es2
			gen `es2' =`ses'*(10^`power')
		}
		if "`3'" != "" {
			gen `lci'   	=`3'*(10^`power')
		}
		if "`4'" != "" {
			gen `uci'   	=`4'*(10^`power')
		}
		
		if "`smooth'" !="" {
			gen `modeles' 	= `5'*(10^`power')
			gen `modellci' 	= `6'*(10^`power')
			gen `modeluci' 	= `7'*(10^`power')
		}
		

		
		local studylb: variable label `label'
		if "`studylb'" == "" {
			label var `label' "`1'"
		}
		if "`rowid'" != "" {
				gen `id' = `rowid'
		}
		else {			
			gen `id' = _n
		}
		
		if "`superimpose1'" != "" {
			tempvar es2 lci2 uci2
			
			tokenize "`superimpose1'", parse(" ")
			gen `es2'		=`1'*(10^`power')
			gen `lci2'   	=`2'*(10^`power')
			gen `uci2'   	=`3'*(10^`power')
		}
		
		if "`superimpose2'" != "" {
			tempvar es3 lci3 uci3
			
			tokenize "`superimpose2'", parse(" ")
			gen `es3'		=`1'*(10^`power')
			gen `lci3'   	=`2'*(10^`power')
			gen `uci3'   	=`3'*(10^`power')
		}
		
		
		gen byte `use'	= 1
		count if `use' == 1
		gen `df' = r(N)
		
		if "`weights'" == "" {
			local wt "nowt"
			local box "nobox"
		}
		
		if "`ovmat'" != "" {
			set obs `=_N + 1'
		
			replace `es' = `ovmat'[1,1] if `id' == .
			replace `lci' = `ovmat'[1,2] if `id' == .
			replace `uci' = `ovmat'[1,3] if `id' == .
			
			/*if "`superimpose1'" != "" {
				replace `lci2' = `ovmat'[1,4] if `id' == .
				replace `uci2' = `ovmat'[1,5] if `id' == .
			}*/

			replace `use' = 5 if `id' == .
			replace `label' = "Summary" if `id' == .
			replace `id' = _N if `id' == .
		}
		
		if "`optimalmat'" != "" {
			local seconddiamond = colsof(`optimalmat') > 3
			set obs `=_N + 1'
		
			replace `es' = `optimalmat'[1,1] if `id' == .
			replace `lci' = `optimalmat'[1,2] if `id' == .
			replace `uci' = `optimalmat'[1,3] if `id' == .
			
			if `seconddiamond' {
				if "`superimpose1'" == "" { 
					tempvar  lci2 uci2
					gen `lci2' = .
					gen `uci2' = .
				}
				replace `lci2' = `optimalmat'[1,4] if `id' == .
				replace `uci2' = `optimalmat'[1,5] if `id' == .
			}

			replace `use' = 5 if `id' == .
			replace `label' = "Optimal" if `id' == .
			replace `id' = _N if `id' == .
		}
		else {
			local seconddiamond 0
		}
				
		gen  `newid' = `id'
		

		//Add five spaces on top of the dataset and 1 space below
		*qui summ `id'
		gen `expand' = 1
		replace `expand' = 1 + 5*(_n==1)  + 1*(_n==_N) 
		expand `expand'
		sort `newid' `groupvar' `use'

		replace `newid' = _n in 1/6
		replace `newid' = `newid' + 5 if _n>6
		replace `label' = "" in 1/5
		replace `use' = -2 in 1/4
		replace `use' = 3 in 5
		replace `newid' = _N  if _N==_n
		replace `use' = 3  if _N==_n
		replace `label' = "" if _N==_n
		
		gen `flag' = 1
		replace `flag' = 0 in 1/4
						
		
		if "`lcols'" == "" {
			local lcols = "`label'"
		}
		else {
			local lcols "`label' `lcols'"
		}
		if "`compabs'" != "" { 
			replace `id' = `newid'
			drop `newid'
		}
		else{		
			egen `rid' = group(`newid')
			replace `id' = `rid'
			drop `rid' `newid'
		}
				
		tempvar estText index predText predLabel wtText modelestText
		cap confirm var `lci'
		if _rc==0 {
			gen str `estText' = string(`es', "%10.`=`dp''f") + " (" + string(`lci', "%10.`=`dp''f") + ", " + string(`uci', "%10.`=`dp''f") + ")"  if (`use' == 1 | `use' == 2 | `use' == 5) 
		}
		else {
			gen str `estText' = string(`es', "%10.`=`dp''f")  if (`use' == 1 | `use' == 2 | `use' == 5)
			gen `lci' = .
			gen `uci' = .
		}
		if "`smooth'" !="" {
			gen str `modelestText' = string(`modeles', "%10.`=`dp''f") + " (" + string(`modellci', "%10.`=`dp''f") + ", " + string(`modeluci', "%10.`=`dp''f") + ")"  if (`use' == 1 )
			replace `modelestText' = string(`es', "%10.`=`dp''f") + " (" + string(`lci', "%10.`=`dp''f") + ", " + string(`uci', "%10.`=`dp''f") + ")"  if (`use' == 2 | `use' == 5)
			replace `estText' = " " if (`use' == 2 | `use' == 5)
		}
		if "`wt'" == "" {
			gen str `wtText' = string(`weights', "%10.`=`dp''f") if (`use' == 1 | `use' == 2 | `use' == 5) & `weights' !=.
		}
		
		// GET MIN AND MAX DISPLAY
		// SORT OUT TICKS- CODE PINCHED FROM MIKE AND FIRandomED. TURNS OUT I'VE BEEN USING SIMILAR NAMES...
		// AS SUGGESTED BY JS JUST ACCEPT ANYTHING AS TICKS AND RESPONSIBILITY IS TO USER!
	
		if "`logscale'" != "" {
			replace `es' = ln(max(`es', 0.00001)) if `es' != .
			replace `lci' = ln(max(`lci', 0.00001)) if `lci' != .
			replace `uci' = ln(`uci')
			
			if "`smooth'" !="" { 
				replace `modeles'  = ln(max(`modeles', 0.00001)) if `modeles' != .
				replace `modellci' = ln(max(`modellci', 0.00001)) if `modellci' != .
				replace `modeluci' = ln(`modeluci')
			}
			
			if "`superimpose1'" != "" {
				replace `es2' = ln(max(`es2', 0.00001)) if `es2' != .
				replace `lci2' = ln(max(`lci2', 0.00001)) if `lci2' != .
				replace `uci2' = ln(`uci2')
			}
			else {
				if `seconddiamond' {
					replace `lci2' = ln(max(`lci2', 0.00001)) if `lci2' != .
					replace `uci2' = ln(`uci2')
				}
			}
			
			if "`superimpose2'" != "" {
				replace `es3' = ln(max(`es3', 0.00001)) if `es3' != .
				replace `lci3' = ln(max(`lci3', 0.00001)) if `lci3' != .
				replace `uci3' = ln(`uci3')
			}
		}
		
		qui summ `lci', detail
		local DXmin = r(min)
		qui summ `uci', detail
		local DXmax = r(max)
				
		if "`xlabel'" != "" {
			if "`logscale'" != "" {
				local DXmin = ln(max(min(`xlabel'), 0.00001))
				local DXmax = ln(max(`xlabel'))
			}
			else{
				local DXmin = min(`xlabel')
				local DXmax = max(`xlabel')
			}
		}
		
		if "`xlabel'"=="" {
			local xlabel "0, `DXmax'"
		}

		local lblcmd ""
		tokenize "`xlabel'", parse(",")
		while "`1'" != ""{
			if "`1'" != ","{
				local lbl = string(`1',"%7.3g")
				if "`logscale'" != "" {
					local val = ln(max(`1', 0.00001))
					/*
					if "`1'" == "0" {
						local val = ln(`=10^(-`dp')')
					}
					else {
						local val = ln(`1')
					}
					*/
				}
				else {
					local val = `1'
				}

				local lblcmd `lblcmd' `val' "`lbl'"
			}
			mac shift
		}
		
		if "`xtick'" == "" {
			local xtick = "`xlabel'"
		}

		local xtick2 = ""
		tokenize "`xtick'", parse(",")
		while "`1'" != ""{
			if "`1'" != ","{
				if "`logscale'" != "" {
					local val = ln(max(`1', 0.00001))
					/*
					if "`1'" == "0" {
						local val = ln(`=10^(-`dp')')
					}
					else {
						local val = ln(`1')
					}
					*/
				}
				else {
					local val = `1'
				}
				local xtick2 = "`xtick2' " + string(`val')
			}
			if "`1'" == ","{
				local xtick2 = "`xtick2'`1'"
			}
			mac shift
		}
		local xtick = "`xtick2'"
		
		local DXmin = (min(`xtick',`DXmin'))
		local DXmax = (max(`xtick',`DXmax'))
		
		*local DXmin= (min(`xlabel',`xtick',`DXmin'))
		*local DXmax= (max(`xlabel',`xtick',`DXmax'))

		local DXwidth = `DXmax'-`DXmin'
	} // END QUI

	/*===============================================================================================*/
	/*==================================== COLUMNS   ================================================*/
	/*===============================================================================================*/
	qui {	// KEEP QUIET UNTIL AFTER DIAMONDS
			
		// DOUBLE LINE OPTION
		if "`double'" != "" & ("`lcols'" != "" | "`rcols'" != ""){
			replace `expand' = 1
			replace `expand' = 2 if `use' == 1
			expand `expand'
			sort `id' `use'
			bys `id' : gen `index' = _n
			sort  `id' `use' `index'
			egen `newid' = group(`id' `index')
			replace `id' = `newid'
			drop `newid'
			
			replace `use' = 1 if `index' == 2
			replace `es' = . if `index' == 2
			replace `lci' = . if `index' == 2
			replace `uci' = . if `index' == 2
			replace `estText' = "" if `index' == 2			

			foreach var of varlist `lcols' `rcols' {
			   cap confirm string var `var'
			   if _rc == 0 {				
					tempvar length words tosplit splitwhere best
					gen `splitwhere' = 0
					gen `best' = .
					gen `length' = length(`var')
					summ `length', det
					gen `words' = wordcount(`var')
					gen `tosplit' = 1 if `length' > r(max)/2+1 & `words' >= 2
					summ `words', det
					local max = r(max)
					forvalues i = 1/`max'{
						replace `splitwhere' = strpos(`var', word(`var',`i')) ///
						 if abs( strpos(`var',word(`var',`i')) - length(`var')/2 ) < `best' ///
						 & `tosplit' == 1
						replace `best' = abs(strpos(`var',word(`var',`i')) - length(`var')/2) ///
						 if abs(strpos(`var',word(`var',`i')) - length(`var')/2) < `best' 
					}

					replace `var' = substr(`var',1,(`splitwhere'-1)) if (`tosplit' == 1) & (`index' == 1)
					replace `var' = substr(`var',`splitwhere',length(`var')) if (`tosplit' == 1) & (`index' == 2)
					replace `var' = "" if (`tosplit' != 1) & (`index' == 2) & (`use' == 1)
					drop `length' `words' `tosplit' `splitwhere' `best'
			   }
			   if _rc != 0{
				replace `var' = . if (`index' == 2) & (`use' == 1)
			   }
			}
		}
				
		local maxline = 1

		if "`lcols'" != "" {
			tokenize "`lcols'"
			local lcolsN = 0

			while "`1'" != "" {
				cap confirm var `1'
				if _rc!=0  {
					di in re "Variable `1' not defined"
					exit _rc
				}
				local lcolsN = `lcolsN' + 1
				tempvar left`lcolsN' leftLB`lcolsN' leftWD`lcolsN'
				cap confirm string var `1'
				if _rc == 0{
					gen str `leftLB`lcolsN'' = `1'
				}
				if _rc != 0{
					cap decode `1', gen(`leftLB`lcolsN'')
					if _rc != 0{
						local f: format `1'
						gen str `leftLB`lcolsN'' = string(`1', "`f'")
						replace `leftLB`lcolsN'' = "" if `leftLB`lcolsN'' == "."
					}
				}
				replace `leftLB`lcolsN'' = "" if (`use' != 1) & (`lcolsN' != 1)
				local colName: variable label `1'
				if "`colName'"==""{
					local colName = "`1'"
				}

				// WORK OUT IF TITLE IS BIGGER THAN THE VARIABLE
				// SPREAD OVER UP TO FOUR LINES IF NECESSARY
				local titleln = length("`colName'")
				tempvar tmpln
				gen `tmpln' = length(`leftLB`lcolsN'')
				qui summ `tmpln' if `use' == 1
				local otherln = r(max)
				drop `tmpln'
				// NOW HAVE LENGTH OF TITLE AND MAX LENGTH OF VARIABLE
				local spread = int(`titleln'/`otherln') + 1
				if `spread' > 4{
					local spread = 4
				}
				local line = 1
				local end = 0
				gettoken now remain : colName

				while `end' == 0 {
					replace `leftLB`lcolsN'' =  `leftLB`lcolsN'' + " " + "`now'" in `line' 
					
					gettoken now remain : remain
					if ("`now'" == "") | (`line' == 4) {
						local end = 1
					}
					if length("`remain'") > `titleln'/`spread' {
						if `end' == 0 {
							local line = `line' + 1
						}
					}
				}
				if `line' > `maxline' {
					local maxline = `line'
				}
				mac shift
			}
		}
		if "`wt'" == "" {
			local rcols = "`wtText' " + "`rcols'"
			label var `wtText' "% Weight"
		}
		if "`smooth'" !=""  {
			local rcols = "`modelestText' " + "`rcols'"
			label var `modelestText' "Smooth Est (`level'% CI)"
		}
		if "`stats'" == "" {
			local rcols = "`estText' " + "`rcols'"
			label var `estText' "`sumstat' (`level'% CI)"
		}

		tempvar extra
		gen `extra' = " "
		label var `extra' " "
		local rcols = "`rcols' `extra'"

		local rcolsN = 0
		if "`rcols'" != "" {
			tokenize "`rcols'"
			local rcolsN = 0
			
			while "`1'" != ""{
				cap confirm var `1'
				if _rc!=0  {
					di in re "Variable `1' not defined"
					exit _rc
				}
				local rcolsN = `rcolsN' + 1
				tempvar right`rcolsN' rightLB`rcolsN' rightWD`rcolsN'
				cap confirm string var `1'
				if _rc == 0{
					gen str `rightLB`rcolsN'' = `1'
				}
				if _rc != 0{
					local f: format `1'
					gen str `rightLB`rcolsN'' = string(`1', "`f'")
					replace `rightLB`rcolsN'' = "" if `rightLB`rcolsN'' == "."
				}
				replace `rightLB`rcolsN'' = "" if (`use' == 3 |  `use' == -2 )
				
				local colName: variable label `1'
				if "`colName'"==""{
					local colName = "`1'"
				}

				// WORK OUT IF TITLE IS BIGGER THAN THE VARIABLE
				// SPREAD OVER UP TO FOUR LINES IF NECESSARY
				local titleln = length("`colName'")
				tempvar tmpln
				gen `tmpln' = length(`rightLB`rcolsN'')
				qui summ `tmpln' if `use' == 1
				local otherln = r(max)
				drop `tmpln'
				// NOW HAVE LENGTH OF TITLE AND MAX LENGTH OF VARIABLE
				local spread = int(`titleln'/`otherln')+1
				if `spread' > 4{
					local spread = 4
				}

				local line = 1
				local end = 0

				gettoken now remain : colName
				while `end' == 0 {
					replace `rightLB`rcolsN'' = `rightLB`rcolsN'' + " " + "`now'" in `line'
					gettoken now remain : remain

					if ("`now'" == "") | (`line' == 4) {
						local end = 1
					}
					if  length("`remain'") > `titleln'/`spread' {
						if `end' == 0 {
							local line = `line' + 1
						}
					}
				}
				if `line' > `maxline'{
					local maxline = `line'
				}
				mac shift
			}
		}

		// now get rid of extra title rows if they weren't used
		if `maxline'==3 {
			drop in 4	
		}
		if `maxline'==2 {
			drop in 3/4 
		}
		if `maxline'==1 {
			drop in 2/4 
		}
		
		sort `id' `groupvar' `use'

		egen `newid' = group(`id')
		replace `id' = `newid'
		drop `newid'

		
		if "`compabs'" != "" {
			replace `id' = _n
		}
				
		if "`compabs'" != "" {
			local borderline = `maxline' + 1.5
		}
		else {
			local borderline = `maxline' + 0.75
		}
		 
		local leftWDtot = 0
		local rightWDtot = 0
		local leftWDtotNoTi = 0

		forvalues i = 1/`lcolsN'{
			getlen `leftLB`i'' `leftWD`i''
			qui summ `leftWD`i'' if `use' != 5 	// DON'T INCLUDE OVERALL STATS AT THIS POINT
			local maxL = r(max)
			local leftWDtotNoTi = `leftWDtotNoTi' + `maxL'
			replace `leftWD`i'' = `maxL'
		}
		tempvar titleLN				// CHECK IF OVERALL LENGTH BIGGER THAN REST OF LCOLS
		getlen `leftLB1' `titleLN'	
		qui summ `titleLN' if `use' == 5
		local leftWDtot = max(`leftWDtotNoTi', r(max))

		forvalues i = 1/`rcolsN'{
			getlen `rightLB`i'' `rightWD`i''
			qui summ `rightWD`i'' if  `use' != 5
			
			replace `rightWD`i'' = r(max)
			local rightWDtot = `rightWDtot' + r(max)
		}
		

		// CHECK IF NOT WIDE ENOUGH (I.E., OVERALL INFO TOO WIDE)
		// LOOK FOR EDGE OF DIAMOND summ `lci' if `use' == ...

		tempvar maxLeft
		getlen `leftLB1' `maxLeft'
		qui count if `use' == 2 | `use' == 5 
		if r(N) > 0 {
			summ `maxLeft' if `use' == 2 | `use' == 5 	// NOT TITLES THOUGH!
			local max = r(max)
			if `max' > `leftWDtotNoTi'{
				// WORK OUT HOW FAR INTO PLOT CAN EXTEND
				// WIDTH OF LEFT COLUMNS AS FRACTION OF WHOLE GRAPH
				local x = `leftWDtot'*(`astext'/100)/(`leftWDtot'+`rightWDtot')
				tempvar y
				// SPACE TO LEFT OF DIAMOND WITHIN PLOT (FRAC OF GRAPH)
				gen `y' = ((100-`astext')/100)*(`lci'-`DXmin') / (`DXmax'-`DXmin') 
				qui summ `y' if `use' == 2 | `use' == 5
				local extend = 1*(r(min)+`x')/`x'
				local leftWDtot = max(`leftWDtot'/`extend',`leftWDtotNoTi') // TRIM TO KEEP ON SAFE SIDE
													// ALSO MAKE SURE NOT LESS THAN BEFORE!
			}

		}
		local LEFT_WD = `leftWDtot'
		local RIGHT_WD = `rightWDtot'
		
		local textWD = (`DXwidth'*(`astext'/(100-`astext'))) /(`leftWDtot' + `rightWDtot')
				
		forvalues i = 1/`lcolsN'{
			gen `left`i'' = `DXmin' - 0.05*(`DXmax' - `DXmin') - `leftWDtot'*`textWD'
			local leftWDtot = `leftWDtot'-`leftWD`i''
		}

		gen `right1' = `DXmax' + 0.05*(`DXmax' - `DXmin')
		forvalues i = 2/`rcolsN'{
			local r2 = `i' - 1
			gen `right`i'' = `right`r2'' + `rightWD`r2''*`textWD'
		}

		local AXmin = `left1'
		*local AXmin = 0.85
		local AXmax = `DXmax' + `rightWDtot'*`textWD'

		// DIAMONDS 
		tempvar DIAMleftX DIAMrightX DIAMbottomX DIAMtopX DIAMleftY1 DIAMrightY1 DIAMleftY2 DIAMrightY2 DIAMbottomY DIAMtopY
		
		//Complete diamond
		gen `DIAMleftX'   = `lci' if `use' == 2 | `use' == 5 
		gen `DIAMleftY1'  = `id' if (`use' == 2 | `use' == 5) 
		gen `DIAMleftY2'  = `id' if (`use' == 2 | `use' == 5) 
		
		gen `DIAMrightX'  = `uci' if (`use' == 2 | `use' == 5)
		gen `DIAMrightY1' = `id' if (`use' == 2 | `use' == 5)
		gen `DIAMrightY2' = `id' if (`use' == 2 | `use' == 5)
		
		gen `DIAMtopX'    = `es' if (`use' == 2 | `use' == 5)
		gen `DIAMbottomY' = `id' - 0.4 if (`use' == 2 | `use' == 5)
		gen `DIAMtopY' 	  = `id' + 0.4 if (`use' == 2 | `use' == 5)

		if `seconddiamond'  {
			tempvar DIAMleftX1 DIAMrightX1
			gen `DIAMleftX1'   = `lci2' if `use' == 2 | `use' == 5 
			gen `DIAMrightX1'  = `uci2' if (`use' == 2 | `use' == 5)
		}
		
		//Incomplete diamonds
		replace `DIAMleftX' = `DXmin' if (`lci' < `DXmin' ) & (`use' == 2 | `use' == 5) //cut the left side to the left limit
		replace `DIAMleftX' = . if (`es' < `DXmin' ) & (`use' == 2 | `use' == 5) //miss it if outside limit
		replace `DIAMleftX' = . if (`df' < 2) & (`use' == 2 | `use' == 5)  //miss it if one study
		
		replace `DIAMleftY1' = `id' + 0.4*(abs((`DXmin' -`lci')/(`es'-`lci'))) if (`lci' < `DXmin' ) & (`use' == 2 | `use' == 5) 
		replace `DIAMleftY1' = . if (`es' < `DXmin' ) & (`use' == 2 | `use' == 5) //miss it if outside limit
	
		replace `DIAMleftY2' = `id' - 0.4*( abs((`DXmin' -`lci')/(`es'-`lci')) ) if (`lci' < `DXmin' ) & (`use' == 2 | `use' == 5) 
		replace `DIAMleftY2' = . if (`es' < `DXmin' ) & (`use' == 2 | `use' == 5) 
		
		//Cutting the right side 
		replace `DIAMrightX' = `DXmax' if (`uci' > `DXmax' ) & (`use' == 2 | `use' == 5) 
		replace `DIAMrightX' = . if (`es' > `DXmax' ) & (`use' == 2 | `use' == 5) 
		
		//If one study, no diamond
		replace `DIAMrightX' = . if (`df' == 1) & (`use' == 2 | `use' == 5) 
	
		replace `DIAMrightY1' = `id' + 0.4*( abs((`uci'-`DXmax' )/(`uci'-`es')) ) if (`uci' > `DXmax' ) & (`use' == 2 | `use' == 5) 
		replace `DIAMrightY1' = . if (`es' > `DXmax' ) & (`use' == 2 | `use' == 5) 

		replace `DIAMrightY2' = `id' - 0.4*( abs((`uci'-`DXmax' )/(`uci'-`es')) ) if (`uci' > `DXmax' ) & (`use' == 2 | `use' == 5) 
		replace `DIAMrightY2' = . if (`es' > `DXmax' ) & (`use' == 2 | `use' == 5) 
			
		replace `DIAMbottomY' = `id' - 0.4*( abs((`uci'-`DXmin' )/(`uci'-`es')) ) if (`es' < `DXmin' ) & (`use' == 2 | `use' == 5) & (abs((`uci'-`DXmin' )/(`uci'-`es')) < 1)
		replace `DIAMbottomY' = `id' - 0.4*( abs((`DXmax' -`lci')/(`es'-`lci')) ) if (`es' > `DXmax' ) & (`use' == 2 | `use' == 5) & abs((`DXmax' -`lci')/(`es'-`lci')) < 1

		replace `DIAMtopY' = `id' + 0.4*( abs((`uci'-`DXmin' )/(`uci'-`es')) ) if (`es' < `DXmin' ) & (`use' == 2 | `use' == 5) & (abs((`uci'-`DXmin' )/(`uci'-`es')) < 1)
		replace `DIAMtopY' = `id' + 0.4*( abs((`DXmax' -`lci')/(`es'-`lci')) ) if (`es' > `DXmax' ) & (`use' == 2 | `use' == 5) & (abs((`DXmax' -`lci')/(`es'-`lci')) < 1)
		
		
		replace `DIAMtopX' = `DXmin'  if (`es' < `DXmin' ) & (`use' == 2 | `use' == 5) 
		replace `DIAMtopX' = `DXmax'  if (`es' > `DXmax' ) & (`use' == 2 | `use' == 5) 
		replace `DIAMtopX' = . if ((`uci' < `DXmin' ) | (`lci' > `DXmax' )) & (`use' == 2 | `use' == 5) //miss it if outside limit
		
		gen `DIAMbottomX' = `DIAMtopX'
	} // END QUI
	
		if "`compabs'" != "" {
			tempvar strgroupvar 
			decode `groupvar', gen(`strgroupvar')
		}
		
	forvalues i = 1/`lcolsN'{
		if "`compabs'" != "" {
			*qui replace `leftLB`i'' = "" if (`expanded')
			if `i'==1 {
				*qui replace `leftLB`i'' = "" if (`order' == 2 & `use'==1) | (`use' == -2 & `id' == `=`startg1'-1')
				*qui replace `leftLB`i'' = "Mean: " + "`groupvar'" + " = " + `strgroupvar' if `use'==2
				qui replace `leftLB`i'' = "Summary: " + "`groupvar'" + " = " + `strgroupvar' if `use'==2
			}
		}
		local lcolCommands`i' "(scatter `id' `left`i'', msymbol(none) mlabel(`leftLB`i'') mlabcolor(black) mlabpos(3) mlabsize(`texts'))"
	}

	forvalues i = 1/`rcolsN' {
		/*if "`compabs'" != "" {
			qui replace `rightLB`i'' = "" if (`expanded') 
		}*/
		local rcolCommands`i' "(scatter `id' `right`i'', msymbol(none) mlabel(`rightLB`i'') mlabcolor(black) mlabpos(3) mlabsize(`texts'))"
	}
	
	/*
	forvalues i = 1/`lcolsN'{
		local lcolCommands`i' "(scatter `id' `left`i'', msymbol(none) mlabel(`leftLB`i'') mlabcolor(black) mlabpos(3) mlabsize(`texts'))"
	}

	forvalues i = 1/`rcolsN' {
		local rcolCommands`i' "(scatter `id' `right`i'', msymbol(none) mlabel(`rightLB`i'') mlabcolor(black) mlabpos(3) mlabsize(`texts'))"
	}
	*/
	
	if `"`diamopts'"' == "" {
		local diamopts "lcolor(red)"
		local diamopts1 "lcolor(black)"
		local diamopts2 "lcolor(orange)"
	}
	else {
		if strpos(`"`diamopts'"',"hor") != 0 | strpos(`"`diamopts'"',"vert") != 0 {
			di as error "Options horizontal/vertical not allowed in diamopts()"
			exit
		}
		if strpos(`"`diamopts'"',"con") != 0{
			di as error "Option connect() not allowed in diamopts()"
			exit
		}
		if strpos(`"`diamopts'"',"lp") != 0{
			di as error "Option lpattern() not allowed in diamopts()"
			exit
		}
		local diamopts `"`diamopts'"'
	}
	
	
	if "`compabs'" != "" {
		/*
		local diamopts0 "lcolor("0 0 0")"
		local diamopts1 "lcolor("255 127 0")"
		*/
		
		local diamopts0 "lcolor("40 40 85")"
		local diamopts1 "lcolor("163 100 249")"
	}
		
	//Box options
	if "`box'" == "" {
		local iw = "[aw = `weights']"
		if `"`boxopts'"' != "" & strpos(`"`boxopts'"',"msy") == 0{
			local boxopts = `"`boxopts' msymbol(square)"' 
		}
		if `"`boxopts'"' != "" & strpos(`"`boxopts'"',"msi") == 0{
			local boxopts = `"`boxopts' msize(0.5)"' 
		}
		if `"`boxopts'"' != "" & strpos(`"`boxopts'"',"mco") == 0{
			local boxopts = `"`boxopts' mcolor("0 0 0")"' 
		}
		if `"`boxopts'"' == "" {
			local boxopts "msymbol(square) msize(.5) mcolor("0 0 0")"
		}
		else{
			local boxopts `"`boxopts'"'
		}
	}
	if ("`box'" != "") {
		local boxopts "msymbol(none)"
		local iw
	}
	
	//Point options
	if "`smooth'" != "" {
		local pointsymbol "msymbol(Oh)"
		local pointcolor "mcolor("128 128 128")"
		local pointsize "msize(small)"
	}
	else {
		local pointsymbol "msymbol(O)"
		local pointcolor "mcolor("0 0 0")"
		local pointsize "msize(vsmall)"
	}
	
	if `"`pointopts'"' != "" & strpos(`"`pointopts'"',"msy") == 0 {
		local pointopts = `"`pointopts' `pointsymbol'"' 
	}
	if `"`pointopts'"' != "" & strpos(`"`pointopts'"',"ms") == 0 {
		local pointopts = `"`pointopts' `pointsize'"' 
	}
	if `"`pointopts'"' != "" & strpos(`"`pointopts'"',"mc") == 0 {
		local pointopts = `"`pointopts' `pointcolor'"' 
	}
	if `"`pointopts'"' == "" {
		local pointopts `"`pointsymbol' `pointsize' `pointcolor'"'
	}
	else{
		local pointopts `"`pointopts'"'
	}
		
	//Smooth Point options

	if `"`smoothpointopts'"' == "" {
		local smoothpointopts "msymbol(D) msize(vsmall) mcolor("0 0 0")"
		if "`compabs'" != "" {
			/*
			local smoothpointopts0 "msymbol(D) msize(vsmall) mcolor("0 0 0")"
			local smoothpointopts1 "msymbol(D) msize(vsmall) mcolor("255 127 0")"
			*/
			local smoothpointopts0 "msymbol(D) msize(vsmall) mcolor("40 40 85")"
			local smoothpointopts1 "msymbol(D) msize(vsmall) mcolor("163 100 249")"
		}
	}
	else {
		local smoothpointopts `"`smoothpointopts'"'
	}
	
	// CI options
	if `"`ciopts'"' == "" {		
		if "`compabs'" != "" {
				/*
				local ciopts0 = `"lcolor("0 0 0")"' 
				local ciopts1 = `" lcolor("255 127 0")"'
				*/
				local ciopts0 = `"lcolor("40 40 85")"' 
				local ciopts1 = `" lcolor("163 100 249")"'
		}
		
		if "`superimpose1'" != "" {
			local ciopts1 = `"lcolor("0 0 0")"' 
			local ciopts2 = `"lwidth(1.5) lcolor(orange)"'
		}
		if "`superimpose2'" != "" {
			local ciopts3 = `"lwidth(0.75) lcolor(cyan)"'
		}
		
		if "`smooth'" != "" {
			local ciopts = `"`ciopts' lwidth(1.25) lcolor(gs13)"'
		}
		else {
			if "`cicolor'" == "" {
				local ciopts = `"`ciopts' lcolor("0 0 0")"' 
			}
			else {
				local cioptsg = `"`ciopts' lcolor("green")"' 
				local cioptsr = `"`ciopts' lcolor("red")"' 
			}			
		}
	}
	else {
		if strpos(`"`ciopts'"',"hor") != 0 | strpos(`"`ciopts'"',"vert") != 0{
			di as error "Options horizontal/vertical not allowed in ciopts()"
			exit
		}
		if strpos(`"`ciopts'"',"con") != 0{
			di as error "Option connect() not allowed in ciopts()"
			exit
		}
		if strpos(`"`ciopts'"',"lp") != 0 {
			di as error "Option lpattern() not allowed in ciopts()"
			exit
		}
	
		local ciopts `"`ciopts'"'
	}
	
	//Smooth ci
	if `"`smoothciopts'"' == "" {
		local smoothciopts "lcolor("0 0 0")"
	}
	else {
		if strpos(`"`smoothciopts'"',"hor") != 0 | strpos(`"`smoothciopts'"',"vert") != 0{
			di as error "Options horizontal/vertical not allowed in ciopts()"
			exit
		}
		if strpos(`"`smoothciopts'"',"con") != 0{
			di as error "Option connect() not allowed in ciopts()"
			exit
		}
		if strpos(`"`smoothciopts'"',"lp") != 0 {
			di as error "Option lpattern() not allowed in ciopts()"
			exit
		}
		if strpos(`"`smoothciopts'"',"lw") == 0 {
				local smoothciopts = `"`smoothciopts' lwidth(.5)"' 
		}

		local smoothciopts `"`smoothciopts'"'
	}
	
	// Arrow options
	if `"`arrowopts'"' == "" {
		if "`smooth'" != "" {
			local arrowopts "mcolor(red) lstyle(none)"
		}
		else {
			local arrowopts "mcolor("0 0 0") lstyle(none)"
		}
	}
	else {
		local forbidden "connect horizontal vertical lpattern lwidth lcolor lsytle"
		foreach option of local forbidden {
			if strpos(`"`arrowopts'"',"`option'")  != 0 {
				di as error "Option `option'() not allowed in arrowopts()"
				exit
			}
		}
		if `"`arrowopts'"' != "" & strpos(`"`arrowopts'"',"mc") == 0{
			local arrowopts = `"`arrowopts' mcolor("0 0 0")"' 
		}
		local arrowopts `"`arrowopts' lstyle(none)"'
	}

	// END GRAPH OPTS

	tempvar tempOv overrallLine ovMin ovMax h0Line
	
	if `"`olineopts'"' == "" {
		local olineopts "lwidth(thin) lcolor(red) lpattern(shortdash)"
	}
	qui summ `id'
	local DYmin = r(min)
	local DYmax = r(max) + 2
	
	qui summ `es' if `use' == 5 
	local overall = r(max)
	if `overall' > `DXmax' | `overall' < `DXmin' | "`ovline'" != "" {	// ditch if not on graph
		local overallCommand ""
	}
	else {
		local overallCommand `" (pci `=`DYmax'-2' `overall' `borderline' `overall', `olineopts') "'
	
	}
	if "`ovline'" != "" {
		local overallCommand ""
	}
	if "`subline'" != "" & "`groupvar'" != "" {
		local sublineCommand ""		
		qui label list `groupvar'
		local nlevels = r(max)
		forvalues l = 1/`nlevels' {
			summ `es' if `use' == 2  & `groupvar' == `l' 
			local tempSub`l' = r(mean)
			qui summ `id' if `use' == 1 & `groupvar' == `l'
			local subMax`l' = r(max) + 1
			local subMin`l' = r(min) - 2
			qui count if `use' == 1 & `groupvar' == `l' 
			if r(N) > 1 {
				local sublineCommand `" `sublineCommand' (pci `subMin`l'' `tempSub`l'' `subMax`l'' `tempSub`l'', `olineopts')"'
			}
		}
	}
	else {
		local sublineCommand ""
	}

	if `"`xline'"' != "" {
		tokenize "`xline'", parse(",")
		if "`logscale'" != "" {
			if "`1'" == "0" {
				local xlineval = ln(`=10^(-`dp')')
			}
			else {
				local xlineval = ln(`1')
			}
		}
		else {
			local xlineval = `1'
		}
		if "`3'" == "" {
			local xlineopts = "`3'"
		}
		else {
			local xlineopts = "lcolor(black)"
		}
		local xlineCommand `" (pci `=`DYmax'-2' `xlineval' `borderline' `xlineval', `xlineopts') "'
	}

	qui {
		//Generate indicator on direction of the off-scale arro
		tempvar rightarrow leftarrow biarrow noarrow rightlimit leftlimit offRhiY offRhiX offRloY offRloX offLloY offLloX offLhiY offLhiX
		gen `rightarrow' = 0
		gen `leftarrow' = 0
		gen `biarrow' = 0
		gen `noarrow' = 0
		
		replace `rightarrow' = 1 if ///
			(round(`uci', 0.001) > round(`DXmax' , 0.001)) & ///
			(round(`lci', 0.001) >= round(`DXmin' , 0.001))  & ///
			(`use' == 1 | `use' == 4) & (`uci' != .) & (`lci' != .)
			
		replace `leftarrow' = 1 if ///
			(round(`lci', 0.001) < round(`DXmin' , 0.001)) & ///
			(round(`uci', 0.001) <= round(`DXmax' , 0.001)) & ///
			(`use' == 1 | `use' == 4) & (`uci' != .) & (`lci' != .)
		
		replace `biarrow' = 1 if ///
			(round(`lci', 0.001) < round(`DXmin' , 0.001)) & ///
			(round(`uci', 0.001) > round(`DXmax' , 0.001)) & ///
			(`use' == 1 | `use' == 4) & (`uci' != .) & (`lci' != .)
			
		replace `noarrow' = 1 if ///
			(`leftarrow' != 1) & (`rightarrow' != 1) & (`biarrow' != 1) & ///
			(`use' == 1 | `use' == 4) & (`uci' != .) & (`lci' != .)	

		replace `lci' = `DXmin'  if (round(`lci', 0.001) < round(`DXmin' , 0.001)) & (`use' == 1 | `use' == 4) 
		replace `uci' = `DXmax'  if (round(`uci', 0.001) > round(`DXmax' , 0.001)) & (`uci' !=.) & (`use' == 1 | `use' == 4) 
		
		//Upper bound not estimable
		replace `rightarrow' = 1 if  (`lci' !=.) & (`uci' ==.) & (`use' == 1 | `use' == 4) 
		replace `uci' = `DXmax'  if  (`lci' !=.) & (`uci' ==.) & (`use' == 1 | `use' == 4) 
		
		replace `lci' = . if (round(`uci', 0.001) < round(`DXmin' , 0.001)) & (`uci' !=. ) & (`use' == 1 | `use' == 4) 
		replace `uci' = . if (round(`lci', 0.001) > round(`DXmax' , 0.001)) & (`lci' !=. ) & (`use' == 1 | `use' == 4)
		replace `es' = . if (round(`es', 0.001) < round(`DXmin' , 0.001)) & (`use' == 1 | `use' == 4) 
		replace `es' = . if (round(`es', 0.001) > round(`DXmax' , 0.001)) & (`use' == 1 | `use' == 4) 
		
		if "`es2'" != "" {
			replace `es2' = . if (round(`es2', 0.001) < round(`DXmin' , 0.001)) & (`use' == 1 | `use' == 4) 
			replace `es2' = . if (round(`es2', 0.001) > round(`DXmax' , 0.001)) & (`use' == 1 | `use' == 4) 
			
			if "`superimpose1'" != "" {
				tempvar rightarrow2 leftarrow2 biarrow2
				gen `rightarrow2' = 0
				gen `leftarrow2' = 0
				gen `biarrow2' = 0
				
				replace `rightarrow2' = 1 if ///
					(round(`uci2', 0.001) > round(`DXmax' , 0.001)) & ///
					(round(`lci2', 0.001) >= round(`DXmin' , 0.001))  & ///
					(`use' == 1 | `use' == 4) & (`uci2' != .) & (`lci2' != .)
					
				replace `leftarrow2' = 1 if ///
					(round(`lci2', 0.001) < round(`DXmin' , 0.001)) & ///
					(round(`uci2', 0.001) <= round(`DXmax' , 0.001)) & ///
					(`use' == 1 | `use' == 4) & (`uci2' != .) & (`lci2' != .)
					
				replace `biarrow2' = 1 if ///
					(round(`lci2', 0.001) < round(`DXmin' , 0.001)) & ///
					(round(`uci2', 0.001) > round(`DXmax' , 0.001)) & ///
					(`use' == 1 | `use' == 4) & (`uci2' != .) & (`lci2' != .)
				
				
				replace `lci2' = `DXmin'  if (round(`lci2', 0.001) < round(`DXmin' , 0.001)) & (`use' == 1 | `use' == 4) 
				replace `uci2' = `DXmax'  if (round(`uci2', 0.001) > round(`DXmax' , 0.001)) & (`uci2' !=.) & (`use' == 1 | `use' == 4) 
				
				replace `lci2' = . if (round(`uci2', 0.001) < round(`DXmin' , 0.001)) & (`uci2' !=. ) & (`use' == 1 | `use' == 4) 
				replace `uci2' = . if (round(`lci2', 0.001) > round(`DXmax' , 0.001)) & (`lci2' !=. ) & (`use' == 1 | `use' == 4)
			}		
			
		}
		
		if "`es3'" != "" {
			replace `es3' = . if (round(`es3', 0.001) < round(`DXmin' , 0.001)) & (`use' == 1 | `use' == 4) 
			replace `es3' = . if (round(`es3', 0.001) > round(`DXmax' , 0.001)) & (`use' == 1 | `use' == 4) 
			
			if "`superimpose2'" != "" {
				tempvar rightarrow3 leftarrow3 biarrow3
				gen `rightarrow3' = 0
				gen `leftarrow3' = 0
				gen `biarrow3' = 0
				
				replace `rightarrow3' = 1 if ///
					(round(`uci3', 0.001) > round(`DXmax' , 0.001)) & ///
					(round(`lci3', 0.001) >= round(`DXmin' , 0.001))  & ///
					(`use' == 1 | `use' == 4) & (`uci3' != .) & (`lci3' != .)
					
				replace `leftarrow3' = 1 if ///
					(round(`lci3', 0.001) < round(`DXmin' , 0.001)) & ///
					(round(`uci3', 0.001) <= round(`DXmax' , 0.001)) & ///
					(`use' == 1 | `use' == 4) & (`uci3' != .) & (`lci3' != .)
					
				replace `biarrow3' = 1 if ///
					(round(`lci3', 0.001) < round(`DXmin' , 0.001)) & ///
					(round(`uci3', 0.001) > round(`DXmax' , 0.001)) & ///
					(`use' == 1 | `use' == 4) & (`uci3' != .) & (`lci3' != .)
				
				
				replace `lci3' = `DXmin'  if (round(`lci3', 0.001) < round(`DXmin' , 0.001)) & (`use' == 1 | `use' == 4) 
				replace `uci3' = `DXmax'  if (round(`uci3', 0.001) > round(`DXmax' , 0.001)) & (`uci3' !=.) & (`use' == 1 | `use' == 4) 
				
				replace `lci3' = . if (round(`uci3', 0.001) < round(`DXmin' , 0.001)) & (`uci3' !=. ) & (`use' == 1 | `use' == 4) 
				replace `uci3' = . if (round(`lci3', 0.001) > round(`DXmax' , 0.001)) & (`lci3' !=. ) & (`use' == 1 | `use' == 4)
			}		
			
		}

		if "`smooth'" != "" {		
			replace `modellci' = `DXmin'  if (round(`modellci', 0.001) < round(`DXmin' , 0.001)) & (`use' == 1) 
			replace `modeluci' = `DXmax'  if (round(`modeluci', 0.001) > round(`DXmax' , 0.001)) & (`modeluci' !=.) & (`use' == 1 ) 
			
			replace `modellci' = . if (round(`modeluci', 0.001) < round(`DXmin' , 0.001)) & (`modeluci' !=. ) & (`use' == 1 ) 
			replace `modeluci' = . if (round(`modellci', 0.001) > round(`DXmax' , 0.001)) & (`modellci' !=. ) & (`use' == 1 )
			replace `modeles' = . if (round(`modeles', 0.001) < round(`DXmin' , 0.001)) & (`use' == 1 ) 
			replace `modeles' = . if (round(`modeles', 0.001) > round(`DXmax' , 0.001)) & (`use' == 1 )
		}
		
		summ `id'
		local xaxislineposition = r(max)

		local xaxis "(pci `xaxislineposition' `DXmin' `xaxislineposition' `DXmax', lwidth(thin) lcolor(black))"
		/*Xaxis 1 title */
		local xaxistitlex `=(`DXmax' + `DXmin')*0.5'
		
		local xaxistitle  (scatteri `=`xaxislineposition' + 2.25' `xaxistitlex' "`sumstat'", msymbol(i) mlabcolor(black) mlabpos(0) mlabsize(`texts'))
		
		/*xticks*/
		local ticksx
		tokenize "`xtick'", parse(",")	
		while "`1'" != "" {
			if "`1'" != "," {
				local ticksx "`ticksx' (pci `xaxislineposition'  `1' 	`=`xaxislineposition'+.25' 	`1' , lwidth(thin) lcolor(black)) "
			}
			macro shift 
		}
		/*labels*/
		local xaxislabels
		tokenize `lblcmd'
		while "`1'" != ""{			
			local xaxislabels "`xaxislabels' (scatteri `=`xaxislineposition'+1' `1' "`2'", msymbol(i) mlabcolor(black) mlabpos(0) mlabsize(`texts'))"
			macro shift 2
		}
		if "`grid'" != "" {
			tempvar gridy gridxmax gridxmin
			
			gen `gridy' = `id' + 0.5
			gen `gridxmax' = `AXmax'
			gen `gridxmin' = `left1'
			local betweengrids "(pcspike `gridy' `gridxmin' `gridy' `gridxmax'  if `use' == 1 , lwidth(vvthin) lcolor(gs12))"	
		}		
	}	// end qui	
	/*===============================================================================================*/
	/*====================================  GRAPH    ================================================*/
	/*===============================================================================================*/
	//Smooth stats
	if "`smooth'" != "" {
		local xboxcenter "`modeles'"
		local smoothcommands1 "(pcspike `id' `modellci' `id' `modeluci' if `use' == 1 , `smoothciopts')"
		local smoothcommands2 "(scatter `id' `modeles' if `use' == 1 , `smoothpointopts')"
		
		if "`compabs'" != "" {
			local smoothcommands10 "(pcspike `id' `modellci' `id' `modeluci' if `use' == 1 & `groupvar'==1  , `smoothciopts0')"
			local smoothcommands11 "(pcspike `id' `modellci' `id' `modeluci' if `use' == 1 & `groupvar'==2 , `smoothciopts1')"
			
			local smoothcommands20 "(scatter `id' `modeles' if `use' == 1 & `groupvar'==1, `smoothpointopts0')"
			local smoothcommands21 "(scatter `id' `modeles' if `use' == 1 & `groupvar'==2, `smoothpointopts1')"
		}
	}
	else {
		local xboxcenter "`es'"
	}
	local bandcolor = "orange green"
	local i = 1
	local area
	while `bands' > 0 {
		local colorb : word `i' of `bandcolor'
		if "`band`i''" != "" {
			tempvar upbound`i' lobound`i'
			tokenize `band`i'', parse(" ")
			
			gen `upbound`i'' = `2'
			gen `lobound`i'' = `1'
			local area "`area' (area  `upbound`i'' `id' if `use' == 1, base(`base') color(`colorb') fintensity(`=10*`i'') horizontal ) "
			local area "`area' (area  `lobound`i'' `id' if `use' == 1, base(`base') color(`colorb') fintensity(`=10*`i'') horizontal ) "			
		}
		local --bands
		local ++i
	}
	
	//Observed CI
	if "`cicolor'" == "" {
		local cicommand "(pcspike `id' `lci' `id' `uci' if `use' == 1 , `ciopts')"
	}
	else {
		local cicommand "(pcspike `id' `lci' `id' `uci' if `use' == 1 & `cicolor' == "green", `cioptsg') (pcspike `id' `lci' `id' `uci' if `use' == 1 & `cicolor' == "red", `cioptsr')"
	}
	
	if "`compabs'" != "" {
			local cicommand1 "(pcspike `id' `lci' `id' `uci' if `use' == 1 & `groupvar'==1, `ciopts0')"
			local cicommand2 "(pcspike `id' `lci' `id' `uci' if `use' == 1 & `groupvar'==2, `ciopts1')"
	} 
	
	if "`superimpose1'" != "" {
			local cicommand1 "(pcspike `id' `lci2' `id' `uci2' if `use' == 1 & mi(`es2') == 0 , `ciopts2')"
			local cicommand2 "(pcspike `id' `lci' `id' `uci' if `use' == 1, `ciopts1')"
	} 
	
	if "`superimpose2'" != "" {
			local cicommand3 "(pcspike `id' `lci3' `id' `uci3' if `use' == 1 & mi(`es3') == 0 , `ciopts3')"
	}
	
	//Diamonds
	local diamondcommand1 "(pcspike `DIAMleftY1' `DIAMleftX' `DIAMtopY' `DIAMtopX' if ( `use' == 2 |`use' == 5) , `diamopts1')"
	local diamondcommand2 "(pcspike `DIAMtopY' `DIAMtopX' `DIAMrightY1' `DIAMrightX' if (`use' == 2 |`use' == 5) , `diamopts1')"
	local diamondcommand3 "(pcspike `DIAMrightY2' `DIAMrightX' `DIAMbottomY' `DIAMbottomX' if (`use' == 2 |`use' == 5) , `diamopts1')"
	local diamondcommand4 "(pcspike `DIAMbottomY' `DIAMbottomX' `DIAMleftY2' `DIAMleftX' if (`use' == 2 |`use' == 5) , `diamopts1')"
	
	if `seconddiamond' {
		local diamondcommand11 "(pcspike `DIAMleftY1' `DIAMleftX1' `DIAMtopY' `DIAMtopX' if ( `use' == 2 |`use' == 5) , `diamopts2')"
		local diamondcommand21 "(pcspike `DIAMtopY' `DIAMtopX' `DIAMrightY1' `DIAMrightX1' if (`use' == 2 |`use' == 5) , `diamopts2')"
		local diamondcommand31 "(pcspike `DIAMrightY2' `DIAMrightX1' `DIAMbottomY' `DIAMbottomX' if (`use' == 2 |`use' == 5) , `diamopts2')"
		local diamondcommand41 "(pcspike `DIAMbottomY' `DIAMbottomX' `DIAMleftY2' `DIAMleftX1' if (`use' == 2 |`use' == 5) , `diamopts2')"
	}

	//legend
		if "`compabs'" != "" {
			local lab0:label `groupvar' 1
			local lab1:label `groupvar' 2
			local legendon "legend(order(1 2) lab(1 "`groupvar' = `lab0'") lab(2 "`groupvar' = `lab1'") bmargin(zero) size(`texts') cols(2) ring(2) position(5))"
		}
		else {
			local legendoff "legend(off)"
		}
	//add text 
	if "`addtextleft'" != "" {
		qui {
			tempvar lasty firstx
			gen `lasty' = `id'[`use'==5] 
			local lasty = _N-1
			local firstx = `left1'
			replace `leftLB1' = "" if `use' == 5
		}
		if "`textopts'" == "" {
			local textopts "place(e) color(red) size(`texts')"
		}
		local textleft "text(`lasty' `firstx' "`addtextleft'",  `textopts')"
	}
	if "`addtextright'" != "" {
		qui {
			tempvar lasty firstx
			gen `lasty' = `id'[`use'==5] 
			local lasty = _N-1
			local firstx = `right1'
			replace `rightLB1' = "" if `use' == 5
		}
		if "`textoptsright'" == "" {
			local textopts "place(e) color(black) size(`texts')"
		}
		local textright "text(`lasty' `firstx' "`addtextright'",  `textopts')"
	}
	
	if "`es2'" != "" {
		if "`superimpose1'" != "" {
			local poppoints "(scatter `id' `es2' if `use' == 1 & mi(`es2') == 0, msymbol(D) msize(medium) mcolor(orange))"
			local arrows2 "(pcarrow `id' `uci2' `id' `lci2' if `leftarrow2' == 1 &  `use' == 1 , `arrowopts')  (pcarrow `id' `lci2' `id' `uci2' if `rightarrow2' == 1 &  `use' == 1, `arrowopts') (pcbarrow `id' `lci2' `id' `uci2' if `biarrow2' == 1 &  `use' == 1, `arrowopts')"
			
			if "`superimpose2'" != "" {
				local poppoints3 "(scatter `id' `es3' if `use' == 1 & mi(`es3') == 0, msymbol(D) msize(medium) mcolor(cyan))"
				local arrows3 "(pcarrow `id' `uci3' `id' `lci3' if `leftarrow3' == 1 &  `use' == 1 , `arrowopts')  (pcarrow `id' `lci3' `id' `uci3' if `rightarrow3' == 1 &  `use' == 1, `arrowopts') (pcbarrow `id' `lci3' `id' `uci3' if `biarrow3' == 1 &  `use' == 1, `arrowopts')"		
			}
		}
		else {
			local es2points "(scatter `id' `es2' if `use' == 1 , msymbol(D) `pointsize' `pointcolor')"
		}
	}
	
	#delimit ;
	twoway
	
	//bands
		`area'
	
	// Draw diamond to make it to construct the legend
		
		`xlineCommand' `xaxis' `xaxistitle' 
		`ticksx' `xaxislabels' 
		
	 /*COLUMN VARIABLES */	 
		`lcolCommands1' `lcolCommands2' `lcolCommands3' `lcolCommands4'  `lcolCommands5'  `lcolCommands6'
		`lcolCommands7' `lcolCommands8' `lcolCommands9' `lcolCommands10' `lcolCommands11' `lcolCommands12'
		`rcolCommands1' `rcolCommands2' `rcolCommands3' `rcolCommands4'  `rcolCommands5'  `rcolCommands6' 
		`rcolCommands7' `rcolCommands8' `rcolCommands9' `rcolCommands10' `rcolCommands11' `rcolCommands12' 
		
	 /*PLOT BOXES AND PUT ALL THE GRAPH OPTIONS IN THERE */ 
		(scatter `id' `xboxcenter' `iw' if `use' == 1, 
			`boxopts'		
			yscale(range(`DYmin' `DYmax') noline reverse)
			ylabel(none) ytitle("")
			xscale(range(`AXmin' `AXmax') noline)
			xlabel(none)
			yline(`borderline', lwidth(thin) lcolor(gs12))
			xtitle("") `legendoff' xtick(""))
			
	 /*HERE ARE GRIDS */
		`betweengrids'			
	 /*HERE ARE THE CONFIDENCE INTERVALS */
	
		`cicommand' `cicommand1' `cicommand2' `cicommand3'
		`smoothcommands1'	`smoothcommands10' `smoothcommands11' `smoothcommands2' `smoothcommands21' `smoothcommands20'	
	 /*ADD ARROWS */
		`arrows2' `arrows3'
		(pcarrow `id' `uci' `id' `lci' if `leftarrow' == 1 &  `use' == 1 , `arrowopts')	
		(pcarrow `id' `lci' `id' `uci' if `rightarrow' == 1 &  `use' == 1, `arrowopts')	
		(pcbarrow `id' `lci' `id' `uci' if `biarrow' == 1 &  `use' == 1, `arrowopts')
		
	 /*DIAMONDS FOR SUMMARY ESTIMATES  */
		`diamondcommand1' 
		`diamondcommand2' 
		`diamondcommand3' 
		`diamondcommand4' 
		
		`diamondcommand11' 
		`diamondcommand21' 
		`diamondcommand31' 
		`diamondcommand41' 
		
	 /*POINTS `poppoints3' `poppoints' */
		`poppoints'
		(scatter `id' `es' if `use' == 1 , `pointopts')
		`es2points' 
				
	//overall & sublines	
		`overallCommand' `sublineCommand'
		
	//Others		
		, `textleft' `textright' `plotopts' legend(off) 
		;
		#delimit cr	
		
		if `"`graphsave'"' != `""' {
			di _n
			noi graph save `graphsave', replace
		}
		restore
end

/*==================================== GETWIDTH  ================================================*/
/*===============================================================================================*/
capture program drop getlen
program define getlen
//From metaprop

qui{

	gen `2' = 0
	count
	local N = r(N)
	forvalues i = 1/`N'{
		local this = `1'[`i']
		local width: _length "`this'"
		replace `2' =  `width' in `i'
	}
} 
end