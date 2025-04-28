cap program drop replicatefit
program define replicatefit

	syntax, truevarcov(name) truemeans(name) seed(string)  
	
	global seed = `seed'	
	global studies =  _N
	
	keep studyid total0 total1 obsevent0 obsevent1
	

		drawnorm logit0i lori, n(`studies') cov(`truevarcov') means(`truemeans') seed(`seed')
		
		gen logit1i = logit0i  + lori
					
		gen p0 = invlogit(logit0i)	
		gen p1 = invlogit(logit1i)

		gen event1 = rbinomial(total1, p1)
		gen event0 = rbinomial(total0, p0)
	

	gen noevent0 = total0 - event0 
	gen noevent1 = total1 - event1
	
	//True parameters
	global mu0_i = `invfn'(`truemeans'[1, 1])
	global or_i = exp(`truemeans'[1, 2])		
	global tausq_i	= `truevarcov'[1,1]	
	global sigmasq_i	= `truevarcov'[2,2]	
	
	//Fit metan
	do "$scriptdir/fitmetanmodels.do"

	//Fit meta
	do "$scriptdir/fitmetamodels.do"

	*Prepare to export data
	drop _*			
	gen sim = $seed
	gen mu0true = $mu0_i
	gen ortrue = $or_i
	gen tausqtrue = $tausq_i
	gen sigmasqtrue = $sigmasq_i
	gen nstudies = $studies
	gen studysize = $studies
	
	*Export data as csv for R
	local simdataname "simdata-`seed'.csv"
	local simdatafile = "$wdir" + "/" + "`simdataname'"
	export delimit using "`simdatafile'",  replace
	
	//Code for R
	local readdata = "data <- read.csv(" +  "'" + "$wdir" + "/" + "`simdataname'" +  "'" + ")"
	local specifyoutfile = "outfile <-" +  "'" + "$outfile" +  "'" 
	local runmodels = "source(" +  "'" + "$scriptdir" + "/" + "fitRmodels.R" + "'" +")"
	local deletedata = "unlink("  +  "'" + "$wdir" + "/" + "`simdataname'" +  "'" + ")"
	
	//write the lines to the R script
	local rfitscriptname "$wdir/runR.R"
	file open rscript using `rfitscriptname',  text write replace //replace the file
	file write rscript  "`readdata'" _n
	file write rscript  "`specifyoutfile'" _n
	file write rscript  "`runmodels'" _n
	file write rscript  "`deletedata'" _n
	file close rscript
	
	*=========Run R script
	rsource using "`rfitscriptname'", noloutput
	
	preserve
	//To long format
	reshape long event total, i(studyid) j(class)

	gen group = "Control" if class==0
	replace group = "Treatment" if class==1	
	
	//Fit binomial and binomial-models in metapreg 
	do "$scriptdir/fitmetapregmodels.do"
	
	//free intercepts models in metapreg
	do "$scriptdir/fitxtrametapregmodels.do"
	
	
	//hexact model in metapreg
	do "$scriptdir/fithexactmetapregmodel.do"
	restore
end	
