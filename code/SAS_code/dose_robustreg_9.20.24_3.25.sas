/*******************************************
This file uses robust regression
to test dose-response associations.

Robust reg does not have the option to list
multiple dependent variables, so i had to trick 
SAS into running many at once by transposing the CpG
vars from wide to long, and then running analyses by 
the long CpG variable

robust reg gives diff results with every run
must save to hard data file and then use that hard
data file for analysis going forward
*************************************************/

data addrecdrg; *bringing in recurrent drug dx;
	set canolder.cannabiscombined_2021_02_11; 
	keep snum recdrg2645_2;
	proc sort; by snum;

data addfh;
	set canolder.cannabiscogn_rev18Aug21;
	keep snum prsub1;
	proc sort; by snum;

data maindata;
	retain snum cg08923376 cg01472075  cg10616121 cg13337507 cg16609878 cg24339197; 
	set methyl.mastersept24;
	proc sort; by snum;
	
data one;
	merge addrecdrg addfh maindata; by snum;
	
	*To make sure models 1-3 have same N;
	if dailypotgp1845 ne . and pc1_age45--pc32_age45 ne . and sex ne . and ntrphlsp45 ne . and lymphcytsp45 ne . and mncytsp45 ne . and esinphlsp45 ne . and bsphlsp45 and SESchildhd ne . and lscuw311 ne . and prsub1 ne . and smkdxgp1845 ne . and alcdxgrp1845 ne . and recdrg2645_2 ne . then nm_candose =1;

	proc means data=one; var ntrphlsp45 lymphcytsp45  mncytsp45  esinphlsp45  bsphlsp45; *for supp table of measures; run;
	proc means n mean data=one; var cg08923376--cg27564939; where nm_candose = 1; run; *model Ns;
	proc freq; tables dailypotgp1845; where cg05575921 ne . ; run; *N for cannabis groups;
	proc freq; tables smkdxgp1845; where cg05575921 ne . ; run; *N for tobacco groups;
	proc freq; tables alcdxgrp1845; where cg05575921 ne . ; run; *N for alcohol groups;
	proc freq; tables  recdrg2645_2;where cg05575921 ne . ; run; *N for pers drug groups reported in supp table of measures;



*trick to run many ROBUST regs at once: 
transpose wide to long, and sort by cpg so i can run robust regs by cpg
https://communities.sas.com/t5/SAS-Communities-Library/How-do-I-write-a-macro-to-run-multiple-regressions/ta-p/223663;
proc transpose data=one out=long prefix=cg;
by snum;
var cg08923376--cg27564939;
run;

data three;
	set long;
	proc sort; by snum;

*getting a count by snum var of the cpg, so i can use this as my "by" variable in robust reg;
data four;
	set three;
	by snum;
	if first.snum then cpg_index = 1;
	else cpg_index+1;
	run;
 

*dropping cgs from original data set to merge this with the transposed data and not be confusing that we have long and wide cgs;
data five;
	set one; drop cg08923376--cg27564939;
	proc sort; by snum;

*merging transposed data with original data with untransposed (i.e., wide) covariates;
data six;
	merge four five; by snum; run;

	proc sort; by cpg_index;
	*analytic Ns for dose-response;
	proc freq; tables dailypotgp1845; by cpg_index; where nm_candose=1;run;

/*******************************
	Can Model 1
*******************************/

	/*******************************************************************
	OLS diagnostics 
	***********************************************************************/


		ods exclude all; *this speeds analysis by not showing output, but still dumps output into a dataset;

		*ods exclude none; *to show output, which slows down analyses;
		*ods graphics on; *turn on to get plots;
		*ods html path='c:\temp' (url=none) file='regdiagnostic.html';  *here i output the results and diagnostic plots into one file saved in submission 2 folder;

		proc reg data=six plots=diagnostics; model cg1= dailypotgp1845 pc1_age45--pc32_age45 sex ntrphlsp45 lymphcytsp45 mncytsp45 esinphlsp45 bsphlsp45; where nm_candose = 1; by cpg_index;
		output out=t student=res cookd=cookd h=lev;
		run;
		quit;

		ods exclude none;

	*sgplots show some high leverage and high residuals, which is why i use robust regression;
	data diag;
		set t; 
		resid_sq=res*res;

	proc sgplot data=diag;
	scatter y=lev x=resid_sq; by cpg_Index; run; *high leverage obs;

	
	proc sort data=t; *this shows which obs are downweighted;
	by wgt;
	proc print data=t;
	var snum cg1 res cookd lev wgt; run;

	proc corr data=t3; var wgt; with cg1 res cookd lev; run; *this shows weight is correlated with cooksd.;

	proc print data=t;
	where cookd>4/819; *conventional cutoff for large cooksd = 4/n;
	var snum cg1 _NAME_ cookd; by cpg_index; run;

	
	
		*robust reg diagnostics for CpGs - plot of residuals. ;
		proc robustreg data=six method=mm (asympcov=H1) plots=qqplot; model cg1= dailypotgp1845 pc1_age45--pc32_age45 sex ntrphlsp45 lymphcytsp45 mncytsp45 esinphlsp45 bsphlsp45; where nm_candose = 1;  by cpg_index;
		output out=robust_diagnostics residual=residuals predicted=fitted_values; run;

		*this tells us outliers or poor fit;
		proc sgplot data=robust_diagnostics;
	   scatter x=fitted_values y=residuals;
	   refline 0 / axis=y lineattrs=(color=red);
		run;

		*this tells us if residuals are normally distributed;
		proc univariate data=robust_diagnostics;
		 var residuals;
		 histogram / normal;
		qqplot / normal;
		run;

	*/

	*robust reg results;

	ods exclude all; *must do this to exclude printing output, bc printing output prolongs the run;
	
	ods output ParameterEstimates=parest_canregrob;
		proc robustreg data=six method=mm (asympcov=H1); model cg1= dailypotgp1845 pc1_age45--pc32_age45 sex ntrphlsp45 lymphcytsp45 mncytsp45 esinphlsp45 bsphlsp45; where nm_candose = 1; by cpg_index;run;
		quit; 

*robust reg output data, limiting results to can reg predictor (i.e., excluding covariate estimates);
data seven;
	set parest_canregrob;
	if Parameter="DailyPotGp1845"; run;

*need to get name of cpg, so keeping name for one participant who has all 238 (because long file means only need one) to merge back with results;
data eight; 
	set four;
	keep _name_ cpg_index;
	if snum=[removed to avoid publishing SNUMs;
	proc sort; by cpg_index;

data nine; *can select this for merging can m1 results with results of other models;
	merge eight seven; by cpg_index; 
	
	rename estimate=estimate_canm1 lowercl= lowercl_canm1 uppercl= uppercl_canm1  probchisq=canm1_probt;


		proc sort; by cpg_index; run;

data ten; *for selecting and sorting by can model 1 p value;
	set nine;
	keep _name_ canm1_probt cpg_index;
	proc sort; by cpg_index;
	run;

data tena;
	merge six ten; by cpg_index; 

	if canm1_probt ge .05 then cg1 = -99; *setting cpgs to missing for model 2 analysis if model 1 p GE .05. The CpGs with p >.05 are not tested. Setting them to -99 helps for data entry;

/*************************
	Can model 2
**************************/
		ods output ParameterEstimates=parest_canregrob_m2;
		proc robustreg data=tena method=mm (asympcov=H1); model cg1= dailypotgp1845 pc1_age45--pc32_age45 sex ntrphlsp45 lymphcytsp45 mncytsp45 esinphlsp45 bsphlsp45 SESchildhd lscuw311 prsub1; where nm_candose = 1; by cpg_index;run;
		quit; 

data eleven;
	set parest_canregrob_m2;
	if Parameter="DailyPotGp1845"; run;

data twelve; *Select this for merging can model 2 results with other models;
	merge eleven eight ; by cpg_index;

	rename estimate=estimate_canm2 lowercl= lowercl_canm2 uppercl= uppercl_canm2  probchisq=canm2_probt;

	
		proc sort; by cpg_index; run;

data twelvea; *for selecting and sorting by can model 2 p values;
	set twelve;
	keep cpg_index canm2_probt; 
	proc sort; by cpg_index;

data twelveb;
	merge six twelvea; by cpg_index;
	if  canm2_probt = . or canm2_probt ge .05 then cg1=-99; *setting cpgs to missing for model 3 analysis if model 2 p value was GE .05;

/*************************
	Can model 3
**************************/
		ods output ParameterEstimates=parest_canregrob_m3;
		proc robustreg data=twelveb method=mm (asympcov=H1); model cg1= dailypotgp1845 pc1_age45--pc32_age45 sex ntrphlsp45 lymphcytsp45 mncytsp45 esinphlsp45 bsphlsp45 SESchildhd lscuw311 prsub1 smkdxgp1845 alcdxgrp1845 recdrg2645_2 ; where nm_candose = 1; by cpg_index;run;
		quit; 

data thirteen;
	set parest_canregrob_m3;
	if Parameter="DailyPotGp1845";
	run;

data fourteen; *can model 3 results, for merging with other results;
	merge thirteen eight ten; by cpg_index;

	rename estimate=estimate_canm3 lowercl= lowercl_canm3 uppercl= uppercl_canm3  probchisq=canm3_probt;
	


	proc sort; by cpg_index; run;

/******************************
	Tob model 1
*********************************/
		ods output ParameterEstimates=parest_tobregrob_m1;
		proc robustreg data=six method=mm (asympcov=H1); model cg1= smkdxgp1845 pc1_age45--pc32_age45 sex ntrphlsp45 lymphcytsp45 mncytsp45 esinphlsp45 bsphlsp45; where nm_candose = 1; by cpg_index;run;
		quit;

data fifteen;
	set parest_tobregrob_m1;
	if Parameter="SmkDxgp1845";
	run;

data sixteen; *select to merge tob m1 results with other models;
	merge fifteen eight ten; by cpg_index;

	rename estimate=estimate_tobm1 lowercl= lowercl_tobm1 uppercl= uppercl_tobm1  probchisq=tobm1_probt;


		proc sort; by cpg_index; run;

data sixteena;
	set sixteen; keep tobm1_probt cpg_index;
	proc sort; by cpg_index;

data sixteenb;
	merge six sixteena; by cpg_index;
	if tobm1_probt ge .05 then cg1 = -99; *doing this so I don't run a test of any CPG that was non-sig in model 1;

/******************************
	Tob model 2
*********************************/
		ods output ParameterEstimates=parest_tobregrob_m2;
		proc robustreg data=sixteenb method=mm (asympcov=H1); model cg1= smkdxgp1845 pc1_age45--pc32_age45 sex ntrphlsp45 lymphcytsp45 mncytsp45 esinphlsp45 bsphlsp45 SESchildhd lscuw311 prsub1; where nm_candose = 1; by cpg_index;run;
		quit;


data seventeen;
	set parest_tobregrob_m2;
	if Parameter="SmkDxgp1845";
	run;

data eighteen; *tob m2 results for merging with other results;
	merge seventeen eight ; by cpg_index;

	rename estimate=estimate_tobm2 lowercl= lowercl_tobm2 uppercl= uppercl_tobm2  probchisq=tobm2_probt;
	
		proc sort; by cpg_index; run;

data eighteena;
	set eighteen;
	keep tobm2_probt cpg_index; *for selecting and sorting by tob m2 p value;

data eighteenb;
	merge eighteena six; by cpg_index;
	if tobm2_probt= . or tobm2_probt ge .05 then cg1 = -99;


/******************************
	Tob model 3
*********************************/
		ods output ParameterEstimates=parest_tobregrob_m3;
		proc robustreg data=eighteenb method=mm (asympcov=H1); model cg1= smkdxgp1845 pc1_age45--pc32_age45 sex ntrphlsp45 lymphcytsp45 mncytsp45 esinphlsp45 bsphlsp45 SESchildhd lscuw311 prsub1 dailypotgp1845 alcdxgrp1845 recdrg2645_2; where nm_candose = 1; by cpg_index;run;
		quit;

data nineteen;
	set parest_tobregrob_m3;
	if Parameter="SmkDxgp1845";
	run;

data twenty; *tob model 3 results for mering with other data sets;
	merge nineteen eight ; by cpg_index;

	rename estimate=estimate_tobm3 lowercl= lowercl_tobm3 uppercl= uppercl_tobm3  probchisq=tobm3_probt;


	proc sort; by cpg_index; run;

/*Hard data set with results. 
data methyl.dose_mmrobust_Sept24; *merging model cannabis models 1-3 with model 1,2,3 tobacco;
	merge twenty eighteen sixteen fourteen twelve nine; by cpg_index;
	
	
	drop parameter df chisq stderr;

	proc contents data=methyl.dose_mmrobust_Sept24;run;

	run;

*/


data ins;
/*merge twenty eighteen sixteen fourteen twelve nine; by cpg_index; */
 set methyl.dose_mmrobust_Sept24; *select if running the hard dataset;


	

	/***********************************
	Table 2 of can dose-respon models 1-3, sorted by can model 1 p
	Sept 20 2024
	****************************************************/
	ods exclude none;
	*tabled results;
	proc sort; by canm1_probt;

	proc print; var _name_ estimate_canm1 lowercl_canm1 uppercl_canm1 canm1_probt;where canm1_probt <.05; run;
	proc print; var _name_ estimate_canm2 lowercl_canm2 uppercl_canm2 canm2_probt ; where canm1_probt <.05; run;
	proc print; var _name_ estimate_canm3 lowercl_canm3 uppercl_canm3 canm3_probt; where canm1_probt <.05;

	run;

	*N that are sig in each model;
	proc print; var cpg_index _name_; where canm1_probt <.05; *52;
		proc print; var cpg_index _name_; where canm1_probt<.000203; run;
	proc print; var cpg_index _name_; where canm1_probt <.05 and canm2_probt <.05;  run;
		proc print; var cpg_index _name_; where canm1_probt<.05 and canm2_probt <.00096; run;
	proc print; var cpg_index _name_; where canm1_probt <.05 and canm2_probt<.05 and canm3_probt <.05; 
		proc print; var cpg_index _name_; where canm1_probt <.05 and canm2_probt<.05 and canm3_probt <.00111; 
	run;


	/**********************************************************************
	Table 3 of tob dose-response models 1-3, sorted by tob model 1 p
	**********************************************************************/
	proc sort; by tobm1_probt;

	proc print; var _name_ estimate_tobm1 lowercl_tobm1 uppercl_tobm1 tobm1_probt;where tobm1_probt <.05; run;
	proc print; var _name_ estimate_tobm2 lowercl_tobm2 uppercl_tobm2 tobm2_probt ; where tobm1_probt <.05; run;
	proc print; var _name_ estimate_tobm3 lowercl_tobm3 uppercl_tobm3 tobm3_probt; where tobm1_probt <.05; run;

	run;
		*N that are sig in each model;
	proc print; var cpg_index _name_; where tobm1_probt <.05;  run;
		proc print; var cpg_index; where tobm1_probt <.000203;  run;
	proc print; var cpg_index; where tobm1_probt <.05 and tobm2_probt <.05;  run;
			proc print; var cpg_index; where tobm1_probt <.05 and tobm2_probt <.000877;  run;
	proc print; var cpg_index; where tobm1_probt <.05 and tobm2_probt<.05 and tobm3_probt <.05;  run;
	proc print; var cpg_index _name_; where tobm1_probt <.05 and tobm2_probt<.05 and tobm3_probt <.00106; run;
	run;

	run;

	

	
