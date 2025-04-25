
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

//Use metan commands
do "$scriptdir/fitmetanmodels.do"

//Use meta suite
do "$scriptdir/fitmetamodels.do"

//Prepare to export data to R
drop _*			
gen sim = $seed
gen mu0true = $mu0_i
gen ortrue = $or_i
gen tausqtrue = $tausq_i
gen sigmasqtrue = $sigmasq_i
gen nstudies = $studies
gen studysize = $studysize

*Export data to csv for R to use
local simdataname "toR.csv"
local simdatafile = "$wdir" + "/" + "`simdataname'"
export delimit using "`simdatafile'",  replace

//Code to run in R
local readdata = "data <- read.csv(" +  "'" + "$wdir" + "/" + "`simdataname'" +  "'" + ")"
local specifyoutfile = "outfile <-" +  "'" + "$outfile" +  "'" 
local runmodels = "source(" +  "'" + "$scriptdir" + "/" + "fitRmodels.R" + "'" +")"
local deletedata = "unlink("  +  "'" + "$wdir" + "/" + "`simdataname'" +  "'" + ")"

//write the lines to the  R script
local rfitscriptname "$wdir/runR.R"
file open rscript using `rfitscriptname',  text write replace //replace the file
file write rscript  "`readdata'" _n
file write rscript  "`specifyoutfile'" _n
file write rscript  "`runmodels'" _n
file write rscript  "`deletedata'" _n
file close rscript

*=========Run R script
rsource using "`rfitscriptname'"  , noloutput

//Transform data to long format
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
//Fit binomial and binomial-normal models in metapreg
do "$scriptdir/fitmetapregmodels.do"

//hexact model in metaprge
do "$scriptdir/fithexactmetapregmodel.do"

//free intercepts in metapreg
do "$scriptdir/fitxtrametapregmodels.do"
