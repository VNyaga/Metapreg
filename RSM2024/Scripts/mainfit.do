
global seed = 1
global studies =  _N

//True parameters
global studysize_i = .
global mu0_i = .
global or_i = .
global tausq_i	= .
global sigmasq_i = .
global studysize = .

gen noevent0 = total0 - event0 
gen noevent1 = total1 - event1

gen rid = _n

//Fit metan
do "$scriptdir/fitmetanmodels.do"
*do "$scriptdir/fitmetanormodels.do"
*do "$scriptdir/fitmetanrdmodels.do"

//Fit meta
do "$scriptdir/fitmetamodels.do"
*do "$scriptdir/fitmetaormodels.do"
*do "$scriptdir/fitmetardmodels.do"

*Prepare to export data

drop _*			
gen sim = $seed
gen mu0true = $mu0_i
gen ortrue = $or_i
gen tausqtrue = $tausq_i
gen sigmasqtrue = $sigmasq_i
gen nstudies = $studies
gen studysize = $studysize

*Export data
local simdataname "toR.csv"
local simdatafile = "$wdir" + "/" + "`simdataname'"
export delimit using "`simdatafile'",  replace

//Assempt the lines to write
local readdata = "data <- read.csv(" +  "'" + "$wdir" + "/" + "`simdataname'" +  "'" + ")"
local specifyoutfile = "outfile <-" +  "'" + "$outfile" +  "'" 
*local runormodels = "source(" +  "'" + "$scriptdir" + "/" + "fitRormodels.R" + "'" +")"
*local runrdmodels = "source(" +  "'" + "$scriptdir" + "/" + "fitRrdmodels.R" + "'" +")"
local runmodels = "source(" +  "'" + "$scriptdir" + "/" + "fitRmodels.R" + "'" +")"
local deletedata = "unlink("  +  "'" + "$wdir" + "/" + "`simdataname'" +  "'" + ")"

//write the lines to the rscript
local rfitscriptname "$wdir/runR.R"
file open rscript using `rfitscriptname',  text write replace //replace the file
file write rscript  "`readdata'" _n
file write rscript  "`specifyoutfile'" _n
*file write rscript  "`runormodels'" _n
*file write rscript  "`runrdmodels'" _n
file write rscript  "`runmodels'" _n
file write rscript  "`deletedata'" _n
file close rscript

*=========Run R script
rsource using "`rfitscriptname'"  , noloutput

//To long format
reshape long event total, i(rid) j(code)

gen group = "Control" if code==0
replace group = "Treatment" if code==1	

tostring rid, replace force
cap confirm var studyid

if _rc!=0  {
	gen studyid = "Study " + rid
}

global brun 1
global frun 1
global run 1
global test "yes"

sort studyid group
//Fit metapreg other links
do "$scriptdir/fitmetapregmodels.do"

//Fit metapreg log link
do "$scriptdir/fitmetapregmodels2.do"

//hexact model
do "$scriptdir/fithexactmetapregmodel.do"

//Refit betabin loglog
*do "$scriptdir/refitbetabin.do"

//free intercepts
do "$scriptdir/fitxtrametapregmodels.do"

//free intercept log link
do "$scriptdir/fitxtrametapregmodels2.do"

