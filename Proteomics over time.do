
clear
cd "/Users/larslind/Lind kopior/"
use "PIVUS data 70_75_80 proteomics_long.dta"
sort lpnr time
merge lpnr time using "PIVUS data 70_75_80_basalfil_long.dta"
tab _merge
drop _merge
sort lpnr time
merge lpnr time using "/Users/larslind/Lind kopior/PIVUS data 70_75_80 klin kemi_long.dta"
tab _merge
drop _merge
sort lpnr time
merge lpnr time using "/Users/larslind/Lind kopior/PIVUS data 70_80 fat_DXA_long.dta"
tab _merge
drop _merge



cd "/Users/larslind/Lind kopior/Proteomics/Over time"

*basic table

file open table1 using "basic table_over time.txt", replace write
file write table1 "Variable"  _tab "Mean (SD)" _tab "Mean (SD)" _tab "Mean (SD)"  _n


foreach var of varlist kn bmi blodsocker manuelltsbp bmi ldl hdl GFRcombo hb rkarenu mi sjukhusvrdadstroke diabetes {

	quietly: sum `var' if time==70
	file write table1 "`var'"  _tab (round(r(mean),0.01)) " (" (round(r(sd),0.01)) ")"  _tab
	quietly: sum `var' if time==75
	file write table1  (round(r(mean),0.01)) " (" (round(r(sd),0.01)) ")"  _tab
	quietly: sum `var' if time==80
	file write table1  (round(r(mean),0.01)) " (" (round(r(sd),0.01)) ")"  _n

}
file close table1




**************************förändring över tid bara medicinfria
drop if mediciner01==1


file open table2 using "Proteomics over 10 years_std_ej medicin.txt", replace write
file write table2 "Variable"  _tab "Beta" _tab "SE" _tab "P-value"  _tab "conv"	_n

foreach var of varlist IL_8std-ECPstd {


      
 quietly: xtmixed   `var' time kn lpnr || lpnr: , mle iter(20) // sätt en övre gräns vid 20 iterationer
		local conv = e(converged)
        quietly: test time
        local pval = r(p)
        file write table2 "`var'" _tab (round(_b[time],0.001)) _tab (round(_se[time],0.001))   _tab (`pval') _tab (`conv') _n 
        }

file close table2


*för tabell
use "Proteomics over 10 years_std_ej medicin.dta", clear
sort variable
merge 1:1 variable using "/Users/larslind/Lind kopior/Proteomics/proteinnamn_long_analysis.dta"
tab _merge
sort p3
*sätt in p-value=0 manuellt
gen p4=1/p3
sort p4
keep proteinname beta se pvalue p3
order proteinname beta se pvalue p3





**************************förändring över tid _alla



file open table2 using "Proteomics over 10 years_std_alla.txt", replace write
file write table2 "Variable"  _tab "Beta" _tab "SE" _tab "P-value"  _tab "conv"	_n

foreach var of varlist IL_8std-ECPstd {


      
 quietly: xtmixed   `var' time kn lpnr || lpnr: , mle iter(20) // sätt en övre gräns vid 20 iterationer
		local conv = e(converged)
        quietly: test time
        local pval = r(p)
        file write table2 "`var'" _tab (round(_b[time],0.001)) _tab (round(_se[time],0.001))   _tab (`pval') _tab (`conv') _n 
        }

file close table2


*könsinteraction

gen time_sex=time*kn

file open table2 using "Proteomics over 10 years_std_alla_sexinteraction.txt", replace write
file write table2 "Variable"  _tab "Beta" _tab "SE" _tab "P-value"  _tab "conv"	_n

foreach var of varlist IL_8std-ECPstd {


      
 quietly: xtmixed   `var' time kn time_sex lpnr || lpnr: , mle iter(20) // sätt en övre gräns vid 20 iterationer
		local conv = e(converged)
        quietly: test time_sex
        local pval = r(p)
        file write table2 "`var'" _tab (round(_b[time_sex],0.001)) _tab (round(_se[time_sex],0.001))   _tab (`pval') _tab (`conv') _n 
        }

file close table2

insheet using "/Users/larslind/Lind kopior/Proteomics/Over time/Proteomics over 10 years_std_alla_sexinteraction.txt", clear
sort pvalue
sort variable
merge 1:1 variable using "/Users/larslind/Lind kopior/Proteomics/proteinnamn_long_analysis.dta"
tab _merge
sort pvalue
keep proteinname beta se pvalue   
order proteinname beta se pvalue 


xtmixed  TRAIL_R2std time kn time_sex lpnr || lpnr:
xtmixed  TRAIL_R2std time  lpnr if kn==0 || lpnr:
xtmixed  TRAIL_R2std time  lpnr if kn==1 || lpnr:

xtmixed  GFRcombo time kn time_sex lpnr || lpnr:

xtmixed  hb time kn time_sex lpnr || lpnr:
xtmixed  hb time  lpnr if kn==0 || lpnr:
xtmixed  hb time  lpnr if kn==1 || lpnr:


 




*kolla om beta är lika för alla vs friska
insheet using "Proteomics over 10 years_std_ej medicin.txt", clear
keep variable beta
ren beta beta_frisk
sort variable
outsheet using "beta_friska.txt"

insheet using "Proteomics over 10 years_std_alla.txt", clear
keep variable beta
ren beta beta_alla
sort variable
outsheet using "beta_alla.txt"
save "beta_alla.dta"

insheet using "beta_friska.txt", clear
sort variable
merge 1:1 variable using "beta_alla.dta"
tw scatter beta_alla beta_frisk, xtitle("Beta in the healthy subgroup") ytitle("Beta in the total sample") title("Change over time for 84 proteins")
pwcorr beta_alla beta_frisk



********************************************** Vulcanoplot
*få fram exakta p-värden för vucanoplot
set pformat %5.4e
foreach var of varlist /// 
NT_pro_BNPstd  ///  //4.5e-93 
PTX3std  ///  //2.2e-58
CD40std  ///  //1.8e-42 
CASP_8std  ///  //1.9e-74 
AGRPstd  ///  //1.3e-47 
OPGstd  ///   //2.e-101
CX3CL1std  ///   //1.2e-42
ESM_1std  ///   //2.6e-69 
FGF_23std  ///  //2.6e-62
PAR_1std  ///   //1.4e-84
AMstd  ///   //1.0e-122
IL27_Astd  ///   //8.0e-100 
t_PAstd  ///   //7.9e-89
MMP_12std  ///  //2.7e-56 
SPON1std  ///   //1.5e-53 
U_PARstd  ///  //1.0e-104
GDF_15std  ///  //1.0e-73
VEGF_Dstd  ///  //2.0e-55
TRAILstd  ///  //9.5e-41 
SRCstd  {    //1.6e-41
xtmixed   `var' time kn lpnr || lpnr: , mle iter(20)
}

replace pvalue=4.5e-93 if variable=="4.5e-93"

insheet using "Proteomics over 10 years_std_ej medicin.txt", clear
sort pvalue
replace pvalue=2.475e-37 if pvalue==0
gen p2=1/pvalue
gen p3 =log10(p2)
tw scatter p3 beta
*sedan ditsatta -log_p-värden för hand för 0-orna
save "Proteomics over 10 years_std_ej medicin.dta", replace

*själva plotten
use "Proteomics over 10 years_std_ej medicin.dta", clear
tw scatter p3 beta, xline(0) xtitle("Regression coefficient") ytitle("-10log p-value") yline(3.2) scheme(s2mono) title("Changes in proteins in subject free from medication")


***************************** proteomics vs change in GFR

xtmixed GFRcombo time kn lpnr || lpnr: , mle iter(20)


foreach var of varlist GFRcombo {

sort lpnr time
by lpnr: gen `var'0=`var' if _n==1
replace `var'0=`var'0[_n-1] if `var'0==.
gen `var'change=`var'-`var'0
}


*för tabell
file open table4 using "Proteomics over time vs GFRchange.txt", replace write
file write table4 "Variable" _tab "Beta" _tab "SE" _tab "p-value" 	_n


foreach var of varlist IL_8std-ECPstd  {  
	qui: xtmixed  `var' GFRcombochange GFRcombo0 kn time lpnr || lpnr:, mle iter(20)
	local conv = e(converged)
	qui: test GFRcombochange
	local pval_change = r(p)
	qui: test GFRcombo0
	local pval_0 = r(p)
	qui: lincom GFRcombochange
	local b_change = r(estimate)
	local se_change = r(se)
	qui: lincom GFRcombo0
	local b_0 = r(estimate)
	local se_0 = r(se)
	file write table4 "`var'" _tab (round(`b_change' , 0.001)) _tab (round(`se_change' ,0.001)) _tab (`pval_change')  _n
}
file close table4

insheet using "Proteomics over time vs GFRchange.txt", clear
sort pvalue

set pformat %5.4e
foreach var of varlist  ///
TNF_R2std  ///  //1.8e-61 
CSTBstd ///  //3.8e-65 
FABP4std ///  //7.9e-38
TRAIL_R2std ///  9.3e-62
CD40std ///  //7.8e-61
TNF_R1std ///  //3.3e-73 
FGF_23std ///  //3.2e-50 
RENstd ///  //2.0e-46
PlGFstd {  //2.4e-36
xtmixed  `var' GFRcombochange GFRcombo0 kn time lpnr || lpnr:, mle iter(20)

}


insheet using "Proteomics over time vs GFRchange.txt", clear
gen p2=1/pvalue
gen p3 =log10(p2)
sort pvalue
*sedan ditsatta -log_p-värden för hand för 0-orna
save "Proteomics over time vs GFRchange.dta", replace

*själva plotten
use "Proteomics over time vs GFRchange.dta", clear
tw scatter p3 beta, xline(0) xtitle("Regression coefficient") ytitle("-10log p-value") yline(3.2) scheme(s2mono) title("Changes in proteins vs change in GFR")


***************************** proteomics vs change in Hb

xtmixed hb time kn lpnr || lpnr: , mle iter(20)




foreach var of varlist hb {

sort lpnr time
by lpnr: gen `var'0=`var' if _n==1
replace `var'0=`var'0[_n-1] if `var'0==.
gen `var'change=`var'-`var'0
}

xtmixed GFRcombo hbchange hb0 kn time lpnr || lpnr:, mle iter(20)


*för tabell
file open table4 using "Proteomics over time vs Hbchange.txt", replace write
file write table4 "Variable" _tab "Beta" _tab "SE" _tab "p-value" 	_n


foreach var of varlist IL_8std-ECPstd  {  
	qui: xtmixed  `var' hbchange hb0 kn time lpnr || lpnr:, mle iter(20)
	local conv = e(converged)
	qui: test hbchange
	local pval_change = r(p)
	qui: test hb0
	local pval_0 = r(p)
	qui: lincom hbchange
	local b_change = r(estimate)
	local se_change = r(se)
	qui: lincom hb0
	local b_0 = r(estimate)
	local se_0 = r(se)
	file write table4 "`var'" _tab (round(`b_change' , 0.001)) _tab (round(`se_change' ,0.001)) _tab (`pval_change')  _n
}
file close table4

insheet using "Proteomics over time vs Hbchange.txt", clear
sort pvalue

set pformat %5.4e
foreach var of varlist  ///
TNF_R2std  ///  //1.8e-61 
CSTBstd ///  //3.8e-65 
FABP4std ///  //7.9e-38
TRAIL_R2std ///  9.3e-62
CD40std ///  //7.8e-61
TNF_R1std ///  //3.3e-73 
FGF_23std ///  //3.2e-50 
RENstd ///  //2.0e-46
PlGFstd {  //2.4e-36
xtmixed  `var' GFRcombochange GFRcombo0 kn time lpnr || lpnr:, mle iter(20)

}


insheet using "Proteomics over time vs Hbchange.txt", clear
gen p2=1/pvalue
gen p3 =log10(p2)
sort pvalue
*sedan ditsatta -log_p-värden för hand för 0-orna
save "Proteomics over time vs Hbchange.dta", replace

*själva plotten
use "Proteomics over time vs Hbchange.dta", clear
tw scatter p3 beta, xline(0) xtitle("Regression coefficient") ytitle("-10log p-value") yline(3.2) scheme(s2mono) title("Changes in proteins vs change in hemoglobin")

******************både GFR och Hb

foreach var of varlist hb GFRcombo {

sort lpnr time
by lpnr: gen `var'0=`var' if _n==1
replace `var'0=`var'0[_n-1] if `var'0==.
gen `var'change=`var'-`var'0
}


*för tabell hb
file open table4 using "Proteomics over time vs Hbchange_båda.txt", replace write
file write table4 "Variable" _tab "Beta" _tab "SE" _tab "p-value" 	_n


foreach var of varlist IL_8std-ECPstd  {  
	qui: xtmixed  `var' hbchange GFRcombochange GFRcombo0 hb0  kn time lpnr || lpnr:, mle iter(20)
    local conv = e(converged)
	qui: test hbchange
	local pval_change = r(p)
	qui: test hb0
	local pval_0 = r(p)
	qui: lincom hbchange
	local b_change = r(estimate)
	local se_change = r(se)
	qui: lincom hb0
	local b_0 = r(estimate)
	local se_0 = r(se)
	file write table4 "`var'" _tab (round(`b_change' , 0.001)) _tab (round(`se_change' ,0.001)) _tab (`pval_change')  _n
}
file close table4


*för tabell GFR
file open table4 using "Proteomics over time vs GFRchange_båda.txt", replace write
file write table4 "Variable" _tab "Beta" _tab "SE" _tab "p-value" 	_n


foreach var of varlist IL_8std-ECPstd  {  
	qui: xtmixed  `var' hbchange GFRcombochange GFRcombo0 hb0  kn time lpnr || lpnr:, mle iter(20)
	local conv = e(converged)
	qui: test GFRcombochange
	local pval_change = r(p)
	qui: test GFRcombo0
	local pval_0 = r(p)
	qui: lincom GFRcombochange
	local b_change = r(estimate)
	local se_change = r(se)
	qui: lincom GFRcombo0
	local b_0 = r(estimate)
	local se_0 = r(se)
	file write table4 "`var'" _tab (round(`b_change' , 0.001)) _tab (round(`se_change' ,0.001)) _tab (`pval_change')  _n
}
file close table4

*kombinerad tabell
insheet using "Proteomics over time vs Hbchange_båda.txt", clear
ren beta beta_hb
ren se se_hb
ren pvalue pvalue_hb
sort variable
save "Proteomics over time vs Hbchange_båda.dta", replace
insheet using "Proteomics over time vs GFRchange_båda.txt", clear
sort variable
merge variable using "Proteomics over time vs Hbchange_båda.dta"
tab _merge
drop _merge
sort variable
merge 1:1 variable using "/Users/larslind/Lind kopior/Proteomics/proteinnamn_long_analysis.dta"
tab _merge
keep proteinname beta- pvalue_hb
order proteinname beta- pvalue_hb

*för tabell Over time justerad för GFR och Hb


foreach var of varlist hb GFRcombo {

sort lpnr time
by lpnr: gen `var'0=`var' if _n==1
replace `var'0=`var'0[_n-1] if `var'0==.
gen `var'change=`var'-`var'0
}


file open table4 using "Proteomics over time just GFRchange_hb.txt", replace write
file write table4 "Variable" _tab "Beta" _tab "SE" _tab "p-value" 	_n


foreach var of varlist IL_8std-ECPstd  {  
	

	 quietly: xtmixed   `var' time kn lpnr hbchange GFRcombochange GFRcombo0 hb0 || lpnr: , mle iter(20) // sätt en övre gräns vid 20 iterationer
		local conv = e(converged)
        quietly: test time
        local pval = r(p)
        file write table4 "`var'" _tab (round(_b[time],0.001)) _tab (round(_se[time],0.001))   _tab (`pval')  _n 
        }

file close table4

insheet using "Proteomics over time just GFRchange_hb.txt", clear
sort pvalue
 


*könsskillnader hela samplet
file open table3 using "Proteomics over time interactions_sex.txt", replace write
file write table3 "Variable" _tab "Beta" _tab "SE" _tab "P-value"  	_n

foreach var of varlist IL_8std-ECPstd  { 

		sort lpnr time
		tempvar baseline // baseline är en temporär variabel som aldrig sparas i data
		quietly: bysort lpnr: gen `baseline'=`var'[1] // genererar baselinevariabel
      
        quietly: xtmixed  `var' c.time##c.kn `baseline' lpnr || lpnr: , mle iter(20)
		local conv = e(converged)
        quietly: test c.time#c.kn
        local pval = r(p)
		quietly: lincom c.time#c.kn // detta blir samma som koefficienten i detta fall
		local beta = r(estimate)
		local se = r(se)
        file write table3 "`var'" _tab (round(`beta' , 0.001))   _tab (round(`se' ,0.001))  _tab (`pval')  _n
        }

file close table3

insheet using "Proteomics over time interactions_sex.txt", clear
sort pvalue


xtmixed  MMP_3std c.time `baseline' lpnr if kn==1 || lpnr: , mle iter(20)
xtmixed  MMP_3std c.time `baseline' lpnr if kn==0 || lpnr: , mle iter(20)


********************************* Change in BMI vs proteomet
sort lpnr
by lpnr: gen motion=.
by lpnr: replace motion=1 if lugnmotion<2 & lugnmotion!=.
by lpnr: replace motion=2 if lugnmotion>=2 & lugnmotion!=.
by lpnr: replace motion=3 if hrdmotion>0 & hrdmotion!=.
by lpnr: replace motion=4 if hrdmotion>2 & hrdmotion!=.

sort time
by time: sum motion

sort lpnr time

*skapa change over time-variabler

foreach var of varlist bmi whratio {

sort lpnr time
by lpnr: gen `var'0=`var' if _n==1
replace `var'0=`var'0[_n-1] if `var'0==.
gen `var'change=`var'-`var'0
}


*för tabell
file open table4 using "Proteomics over time bmi_confounders.txt", replace write
file write table4 "Variable" _tab "bmichange_b" _tab "bmichange_se" _tab "bmichange_p"  _tab "bmi0_b" _tab "bmi0_se" _tab "bmi0_p" _tab "conv"	_n


foreach var of varlist IL_8std-ECPstd {  
	qui: xtmixed  `var' bmichange bmi0 kn time betablock diuretika statiner insulin oralaantidiabetika rkarenu motion GFRcombo lpnr || lpnr:, mle iter(20)
	local conv = e(converged)
	qui: test bmichange
	local pval_change = r(p)
	qui: test bmi0
	local pval_0 = r(p)
	qui: lincom bmichange
	local b_change = r(estimate)
	local se_change = r(se)
	qui: lincom bmi0
	local b_0 = r(estimate)
	local se_0 = r(se)
	file write table4 "`var'" _tab (round(`b_change' , 0.001)) _tab (round(`se_change' ,0.001)) _tab (`pval_change') _tab (round(`b_0' , 0.001)) _tab (round(`se_0' ,0.001))  _tab (`pval_0') _tab (`conv') _n
}
file close table4

cd "/Users/larslind/Lind kopior/Proteomics/Over time"
insheet using "Proteomics over time bmi_confounders.txt", clear
sort bmichange_p
sort variable
merge 1:1 variable using "/Users/larslind/Lind kopior/Proteomics/proteinnamn_long_analysis.dta"
tab _merge
sort bmichange_p
keep proteinname bmichange_b bmichange_se bmichange_p   
order proteinname bmichange_b bmichange_se bmichange_p  

save "Proteomics over time bmi_confounders.dta", replace

	xtmixed LEPstd bmichange bmi0 kn time betablock diuretika statiner insulin oralaantidiabetika rkarenu motion GFRcombo lpnr || lpnr:, mle iter(20)
	
*kolla könsskillnader	

sort lpnr
gen bmich_kn=bmichange*kn

file open table4 using "Proteomics over time bmi_confounders_sexinteraction.txt", replace write
file write table4 "Variable" _tab "bmichange_b" _tab "bmichange_se" _tab "bmichange_p"  _tab "bmi0_b" _tab "bmi0_se" _tab "bmi0_p" _tab "conv"	_n


foreach var of varlist IL_8std-ECPstd {  
	qui: xtmixed  `var' bmich_kn bmichange bmi0 kn time betablock diuretika statiner insulin oralaantidiabetika rkarenu motion GFRcombo lpnr || lpnr:, mle iter(20)
	local conv = e(converged)
	qui: test bmich_kn
	local pval_change = r(p)
	qui: test bmi0
	local pval_0 = r(p)
	qui: lincom bmich_kn
	local b_change = r(estimate)
	local se_change = r(se)
	qui: lincom bmi0
	local b_0 = r(estimate)
	local se_0 = r(se)
	file write table4 "`var'" _tab (round(`b_change' , 0.001)) _tab (round(`se_change' ,0.001)) _tab (`pval_change') _tab (round(`b_0' , 0.001)) _tab (round(`se_0' ,0.001))  _tab (`pval_0') _tab (`conv') _n
}
file close table4

insheet using "Proteomics over time bmi_confounders_sexinteraction.txt", clear

	
	
	
********************************* Change in WHR vs proteomet

sort lpnr
by lpnr: gen motion=.
by lpnr: replace motion=1 if lugnmotion<2 & lugnmotion!=.
by lpnr: replace motion=2 if lugnmotion>=2 & lugnmotion!=.
by lpnr: replace motion=3 if hrdmotion>0 & hrdmotion!=.
by lpnr: replace motion=4 if hrdmotion>2 & hrdmotion!=.

sort time
by time: sum motion



sort lpnr time

*skapa change over time-variabler

foreach var of varlist bmi whratio {

sort lpnr time
by lpnr: gen `var'0=`var' if _n==1
replace `var'0=`var'0[_n-1] if `var'0==.
gen `var'change=`var'-`var'0
}


*för tabell
file open table4 using "Proteomics over time WHR_confounders.txt", replace write
file write table4 "Variable" _tab "bmichange_b" _tab "bmichange_se" _tab "bmichange_p"  _tab "bmi0_b" _tab "bmi0_se" _tab "bmi0_p" _tab "conv"	_n


foreach var of varlist IL_8std-ECPstd {  
	qui: xtmixed  `var' whratiochange whratio0 kn time betablock diuretika statiner insulin oralaantidiabetika rkarenu motion GFRcombo lpnr || lpnr:, mle iter(20)
	local conv = e(converged)
	qui: test whratiochange
	local pval_change = r(p)
	qui: test whratio0
	local pval_0 = r(p)
	qui: lincom whratiochange
	local b_change = r(estimate)
	local se_change = r(se)
	qui: lincom whratio0
	local b_0 = r(estimate)
	local se_0 = r(se)
	file write table4 "`var'" _tab (round(`b_change' , 0.001)) _tab (round(`se_change' ,0.001)) _tab (`pval_change') _tab (round(`b_0' , 0.001)) _tab (round(`se_0' ,0.001))  _tab (`pval_0') _tab (`conv') _n
}
file close table4

cd "/Users/larslind/Lind kopior/Proteomics/Over time"
insheet using "Proteomics over time WHR_confounders.txt", clear
sort bmichange_p
sort variable
merge 1:1 variable using "/Users/larslind/Lind kopior/Proteomics/proteinnamn_long_analysis.dta"
tab _merge
sort bmichange_p
keep proteinname bmichange_b bmichange_se bmichange_p   
order proteinname bmichange_b bmichange_se bmichange_p 

ren bmichange_b whrchange_b
ren bmichange_se whrchange_se
ren bmichange_p whrchange_p

sort whrchange_p
save "Proteomics over time WHR_confounders.dta", replace
sort 

*välja ut till figurer
use "Proteomics over time bmi_confounders.dta", clear
sort proteinname
merge 1:1 proteinname using "Proteomics over time WHR_confounders.dta"
sort bmichange_p
sort whrchange_p


*kolla könsinteractions
sort lpnr
gen whrch_kn=whratiochange*kn

file open table4 using "Proteomics over time whr_confounders_sexinteraction.txt", replace write
file write table4 "Variable" _tab "bmichange_b" _tab "bmichange_se" _tab "bmichange_p"  _tab "bmi0_b" _tab "bmi0_se" _tab "bmi0_p" _tab "conv"	_n


foreach var of varlist IL_8std-ECPstd {  
	qui: xtmixed  `var' whrch_kn whratiochange whratio0 kn time betablock diuretika statiner insulin oralaantidiabetika rkarenu motion GFRcombo lpnr || lpnr:, mle iter(20)
	local conv = e(converged)
	qui: test whrch_kn
	local pval_change = r(p)
	qui: test whratio0
	local pval_0 = r(p)
	qui: lincom whrch_kn
	local b_change = r(estimate)
	local se_change = r(se)
	qui: lincom whratio0
	local b_0 = r(estimate)
	local se_0 = r(se)
	file write table4 "`var'" _tab (round(`b_change' , 0.001)) _tab (round(`se_change' ,0.001)) _tab (`pval_change') _tab (round(`b_0' , 0.001)) _tab (round(`se_0' ,0.001))  _tab (`pval_0') _tab (`conv') _n
}
file close table4

insheet using "Proteomics over time whr_confounders_sexinteraction.txt", clear

variable
LEPstd
PRLstd
MMP_3std

xtmixed  LEPstd  whratiochange whratio0  time betablock diuretika statiner insulin oralaantidiabetika rkarenu motion GFRcombo lpnr if kn==0 || lpnr:, mle iter(20)
xtmixed  LEPstd  whratiochange whratio0  time betablock diuretika statiner insulin oralaantidiabetika rkarenu motion GFRcombo lpnr if kn==1 || lpnr:, mle iter(20)

xtmixed  PRLstd  whratiochange whratio0  time betablock diuretika statiner insulin oralaantidiabetika rkarenu motion GFRcombo lpnr if kn==0 || lpnr:, mle iter(20)
xtmixed  PRLstd  whratiochange whratio0  time betablock diuretika statiner insulin oralaantidiabetika rkarenu motion GFRcombo lpnr if kn==1 || lpnr:, mle iter(20)

xtmixed  MMP_3std  whratiochange whratio0  time betablock diuretika statiner insulin oralaantidiabetika rkarenu motion GFRcombo lpnr if kn==0 || lpnr:, mle iter(20)
xtmixed  MMP_3std  whratiochange whratio0  time betablock diuretika statiner insulin oralaantidiabetika rkarenu motion GFRcombo lpnr if kn==1 || lpnr:, mle iter(20)


********figurer för delta BMI vs delta protein


foreach var of varlist bmi whratio {

sort lpnr time
by lpnr: gen `var'0=`var' if _n==1
replace `var'0=`var'0[_n-1] if `var'0==.
gen `var'change=`var'-`var'0
}

******** skapar pred-fil.
** bmi0 sätts till medianen och förändringen till tisspecifika p25, p50 och p75. Se help(collapse) för andra alternativ
** Skapa kn och sätt till 0 eller 1, det spelar ingen roll men måste finnas för att prediktionerna ska kunna göras
preserve
collapse (median) bmi0 (p25) bmichange1 = bmichange (p50) bmichange2 = bmichange (p75) bmichange3 = bmichange, by(time)
reshape long bmichange, i(time) j (change) // variabeln change innehåller 1-3 som motsvarar p25, p50 och p75
gen kn = 1
save pred_bmi, replace
restore

foreach var of varlist IL_8std-ECPstd  {  
qui: xtmixed  `var' bmichange bmi0 kn time   || lpnr:

	estimates store est
	preserve
	use pred_bmi, clear
	estimates restore est
	predict xb_`var', xb
	save pred_bmi, replace
	restore
}


use "pred_bmi.dta", clear

tw line xb_ESM_1std time if change == 1, sort  || ///
line xb_ESM_1std time if change == 2,  || ///
line xb_ESM_1std time if change == 3, ///
scheme(s2mono)  ///
xlabel(70(5)80) title("Change in ESM-1 vs change in BMI") ///
legend(lab(1 25th percentile) lab(2 50th percentile) lab(3 75th percentile)) xtitle("Age") ytitle("ESM-1") name(ESM1)

tw line xb_ST2std time if change == 1, sort  || ///
line xb_ST2std time if change == 2,  || ///
line xb_ST2std time if change == 3, ///
scheme(s2mono)  ///
xlabel(70(5)80) title("Change in ST2 vs change in BMI") ///
legend(lab(1 25th percentile) lab(2 50th percentile) lab(3 75th percentile)) xtitle("Age") ytitle("ST2") name(ST2)

tw line xb_LEPstd time if change == 1, sort  || ///
line xb_LEPstd time if change == 2,  || ///
line xb_LEPstd time if change == 3, ///
scheme(s2mono)  ///
xlabel(70(5)80) title("Change in leptin vs change in BMI") ///
legend(lab(1 25th percentile) lab(2 50th percentile) lab(3 75th percentile)) xtitle("Age") ytitle("Leptin") name(BMI_Lep)

tw line xb_FABP4std time if change == 1, sort  || ///
line xb_FABP4std time if change == 2,  || ///
line xb_FABP4std time if change == 3, ///
scheme(s2mono)  ///
xlabel(70(5)80) title("Change in FABP-4 vs change in BMI") ///
legend(lab(1 25th percentile) lab(2 50th percentile) lab(3 75th percentile)) xtitle("Age") ytitle("FABP-4") name(BMI_FABP4)

graph combine ESM1 ST2, altshrink



********figurer för delta WHR vs delta protein


foreach var of varlist bmi whratio {

sort lpnr time
by lpnr: gen `var'0=`var' if _n==1
replace `var'0=`var'0[_n-1] if `var'0==.
gen `var'change=`var'-`var'0
}

******** skapar pred-fil.
** bmi0 sätts till medianen och förändringen till tisspecifika p25, p50 och p75. Se help(collapse) för andra alternativ
** Skapa kn och sätt till 0 eller 1, det spelar ingen roll men måste finnas för att prediktionerna ska kunna göras
preserve
collapse (median) whratio0 (p25) whratiochange1 = whratiochange (p50) whratiochange2 = whratiochange (p75) whratiochange3 = whratiochange, by(time)
reshape long whratiochange, i(time) j (change) // variabeln change innehåller 1-3 som motsvarar p25, p50 och p75
gen kn = 1
save pred_whratio, replace
restore

foreach var of varlist IL_8std-ECPstd  {  
qui: xtmixed  `var' whratiochange whratio0 kn time   || lpnr:

	estimates store est
	preserve
	use pred_whratio, clear
	estimates restore est
	predict xb_`var', xb
	save pred_whratio, replace
	restore
}


use "pred_whratio.dta", clear

tw line xb_CASP_8std time if change == 1, sort  || ///
line xb_CASP_8std time if change == 2,  || ///
line xb_CASP_8std time if change == 3, ///
scheme(s2mono)  ///
xlabel(70(5)80) title("Change in CASP-8 vs change in WHR") ///
legend(lab(1 25th percentile) lab(2 50th percentile) lab(3 75th percentile)) xtitle("Age") ytitle("CASP-8") name(CASP8)

tw line xb_CTSL1std time if change == 1, sort  || ///
line xb_CTSL1std time if change == 2,  || ///
line xb_CTSL1std time if change == 3, ///
scheme(s2mono)  ///
xlabel(70(5)80) title("Change in CTSL-1 vs change in WHR") ///
legend(lab(1 25th percentile) lab(2 50th percentile) lab(3 75th percentile)) xtitle("Age") ytitle("CTSL-1") name(CTSL12)

tw line xb_LEPstd time if change == 1, sort  || ///
line xb_LEPstd time if change == 2,  || ///
line xb_LEPstd time if change == 3, ///
scheme(s2mono)  ///
xlabel(70(5)80) title("Change in leptin vs change in WHR") ///
legend(lab(1 25th percentile) lab(2 50th percentile) lab(3 75th percentile)) xtitle("Age") ytitle("Leptin") name(WHR_Lep)

tw line xb_FABP4std time if change == 1, sort  || ///
line xb_FABP4std time if change == 2,  || ///
line xb_FABP4std time if change == 3, ///
scheme(s2mono)  ///
xlabel(70(5)80) title("Change in FABP-4 vs change in WHR") ///
legend(lab(1 25th percentile) lab(2 50th percentile) lab(3 75th percentile)) xtitle("Age") ytitle("FABP-5") name(WHR_FABP)

graph combine CASP8 CTSL12, altshrink

graph combine BMI_Lep WHR_Lep BMI_FABP4 WHR_FABP, altshrink


foreach var of varlist imtnmedian-ldsd{
ren `var' `var'sin
}



****** Figurer på BMI och WHR over time i R

library(haven)
PIVUS_data_70_75_80_risk_factors_long <- read_dta("Lind kopior/PIVUS data/Long format/PIVUS data 70_75_80 risk factors_long.dta")

library(lattice)
bwtheme <- standard.theme("pdf", color=FALSE)


 densityplot(~ bmi, groups = time, data = PIVUS_data_70_75_80_risk_factors_long,
             plot.points = FALSE, auto.key = list(space = "right", columns = 1), 
             par.settings=bwtheme, xlab = "BMI")
			 
			 
			
			


*************** DXA**************************************************

*change in fat mass



sort lpnr
by lpnr: gen motion=.
by lpnr: replace motion=1 if lugnmotion<2 & lugnmotion!=.
by lpnr: replace motion=2 if lugnmotion>=2 & lugnmotion!=.
by lpnr: replace motion=3 if hrdmotion>0 & hrdmotion!=.
by lpnr: replace motion=4 if hrdmotion>2 & hrdmotion!=.

ren fett_total fett_total_1000
sort lpnr
by lpnr: gen fett_total=fett_total_1000/1000

ren fettrunkleg fettrunkleg_oneleg
gen fettrunkleg=fett_trunk/(2*fett_leg)



sort lpnr time

*skapa change over time-variabler

foreach var of varlist fett_total fett_arm fett_leg fett_trunk fettrunkleg {

sort lpnr time
by lpnr: gen `var'0=`var' if _n==1
replace `var'0=`var'0[_n-1] if `var'0==.
gen `var'change=`var'-`var'0
}


*för tabell
file open table4 using "Proteomics over time fatmass_confounders.txt", replace write
file write table4 "Variable" _tab "bmichange_b" _tab "bmichange_se" _tab "bmichange_p"  _tab "bmi0_b" _tab "bmi0_se" _tab "bmi0_p" _tab "conv"	_n


foreach var of varlist IL_8std-ECPstd {  
	qui: xtmixed  `var' fett_totalchange fett_total0 kn time betablock diuretika statiner insulin oralaantidiabetika rkarenu motion GFRcombo lpnr || lpnr:, mle iter(20)
	local conv = e(converged)
	qui: test fett_totalchange
	local pval_change = r(p)
	qui: test fett_total0
	local pval_0 = r(p)
	qui: lincom fett_totalchange
	local b_change = r(estimate)
	local se_change = r(se)
	qui: lincom fett_total0
	local b_0 = r(estimate)
	local se_0 = r(se)
	file write table4 "`var'" _tab (round(`b_change' , 0.001)) _tab (round(`se_change' ,0.001)) _tab (`pval_change') _tab (round(`b_0' , 0.001)) _tab (round(`se_0' ,0.001))  _tab (`pval_0') _tab (`conv') _n
}
file close table4

cd "/Users/larslind/Lind kopior/Proteomics/Over time"
insheet using "Proteomics over time fatmass_confounders.txt", clear
sort bmichange_p
sort variable
merge 1:1 variable using "/Users/larslind/Lind kopior/Proteomics/proteinnamn_long_analysis.dta"
tab _merge
sort bmichange_p
keep proteinname bmichange_b bmichange_se bmichange_p   
order proteinname bmichange_b bmichange_se bmichange_p 


*change in fat mass trunk/legratio



sort lpnr
by lpnr: gen motion=.
by lpnr: replace motion=1 if lugnmotion<2 & lugnmotion!=.
by lpnr: replace motion=2 if lugnmotion>=2 & lugnmotion!=.
by lpnr: replace motion=3 if hrdmotion>0 & hrdmotion!=.
by lpnr: replace motion=4 if hrdmotion>2 & hrdmotion!=.





sort lpnr time

*skapa change over time-variabler

foreach var of varlist fett_total fett_arm fett_leg fett_trunk fettrunkleg {

sort lpnr time
by lpnr: gen `var'0=`var' if _n==1
replace `var'0=`var'0[_n-1] if `var'0==.
gen `var'change=`var'-`var'0
}


*för tabell
file open table4 using "Proteomics over time fettrunkleg_confounders_twoleg.txt", replace write
file write table4 "Variable" _tab "bmichange_b" _tab "bmichange_se" _tab "bmichange_p"  _tab "bmi0_b" _tab "bmi0_se" _tab "bmi0_p" _tab "conv"	_n


foreach var of varlist IL_8std-ECPstd {  
	qui: xtmixed  `var' fettrunklegchange fettrunkleg0 kn time betablock diuretika statiner insulin oralaantidiabetika rkarenu motion GFRcombo lpnr || lpnr:, mle iter(20)
	local conv = e(converged)
	qui: test fettrunklegchange
	local pval_change = r(p)
	qui: test fettrunkleg0
	local pval_0 = r(p)
	qui: lincom fettrunklegchange
	local b_change = r(estimate)
	local se_change = r(se)
	qui: lincom fettrunkleg0
	local b_0 = r(estimate)
	local se_0 = r(se)
	file write table4 "`var'" _tab (round(`b_change' , 0.001)) _tab (round(`se_change' ,0.001)) _tab (`pval_change') _tab (round(`b_0' , 0.001)) _tab (round(`se_0' ,0.001))  _tab (`pval_0') _tab (`conv') _n
}
file close table4

cd "/Users/larslind/Lind kopior/Proteomics/Over time"
insheet using "Proteomics over time fettrunkleg_confounders_twoleg.txt", clear
sort bmichange_p
sort variable
merge 1:1 variable using "/Users/larslind/Lind kopior/Proteomics/proteinnamn_long_analysis.dta"
tab _merge
sort bmichange_p
keep proteinname bmichange_b bmichange_se bmichange_p   
order proteinname bmichange_b bmichange_se bmichange_p 
