/*************************************************************
Creating one hard dataset from all Karen's data sets
to test age-38 long-term cannabis users vs. non-users
(Supplemental Table 6)
************************************************************/



/*Bringing in methylation beta values data at age 38*/

data one;
	set methyl.p38canbeta; *this data file has NOT been updated to add CpGs from Fang and Garrett;

	/*proc contents; run; */*This shows which are character;



	proc sort; by snum;

data oneins;
	set methyl.age38cg230; *bringing in age 38 cg23079012, which is from fang/garret updated lit search and sig at age 45;
	proc sort; by snum;

data oneb;
	merge one oneins; by snum; 

	*standardizing all age-38 cpgs;
	proc standard m=0 std=1 out=stdcpg38; var cg01035812--cg23079012;
	proc univariate; var cg01035812--cg23079012;run;
	proc means data=stdcpg38; var cg01035812--cg23079012;
	proc sort data=stdcpg38; by snum; run;

data two; *bringing in age 38 principal components;

	set methyl.p38pcs; 	
	proc contents; 
	proc sort; by snum;
	run;



data three;/*Bringing in exposures and covariates*/
	set canolder.cannabiscogn_rev18Aug21;
	keep snum comp0_ltuser2021 comp1_noSUD2021 Comp2_LTtob2021_REV Comp3_LTalc2021REV comp4_recreatn2021 comp5_quit2021 potdxgrp1845 DailyPotGp1845
	 alcdxgrp1845  prsub1 TotTob1845lmt TotAlc1845lmt TotDrg2645lmt TotMar1845lmt marfreq45 smokes45 weeklyalc45 sex SESchildhd lscuw311  seennot45;

	proc sort; by snum;

data four; *bringing in vars to create an age-38 tobacco quit group;
	set methyl.methylproj_phenotypes_dec2023;
	keep snum dxtob18 dxtob21 dxtob26 dxtob32 dxTob38Dsm4 dxTob45Dsm4 cigsday18 cigsday21 CigsDay26 CigsDay32 CigsDay38 CigsDay45;



	proc sort; by snum;


data five; *bringing in the correct smoking variable;
		set canolder.cannabiscombined_2021_02_11; 
		keep snum smkdxgp1845;

		proc sort; by snum;

data six; *brining in white blood cells;
	set methyl.age38wbc;
	keep snum neutrophils38 lymphocytes38 monocytes38 eosinophils38 basophils38 ntrphlsp38 lymphcytsp38 mncytsp38 esinphlsp38  bsphlsp38;

ntrphlsp38=neutrophils38*1;
lymphcytsp38 = lymphocytes38*1;
mncytsp38 =  monocytes38*1;
esinphlsp38 = eosinophils38*1;
 bsphlsp38 = basophils38*1;

	proc freq; tables ntrphlsp38 lymphcytsp38 mncytsp38 esinphlsp38 bsphlsp38; run;

	proc contents; run;
	proc sort; by snum;

data seven; *merging age-38 std CpGs, exposures, tobacco quit group, smoking var, white blood cells;
	merge stdcpg38 two three four five six ; by snum; 


	*long-term users vs squeaky clean;
		if comp0_ltuser2021 = 1 then Lt_clean = 0;
		if comp1_noSUD2021 = 1 then lt_clean = 1;

	
		*long-term users vs tob;
		if comp0_ltuser2021 = 1 then lt_tob=0;
		if Comp2_LTtob2021_REV = 1 then lt_tob = 1;

		*long-term users vs alc ;
		if comp0_ltuser2021 = 1 then lt_alc=0;
		if Comp3_LTalc2021REV = 1 then lt_alc = 1;

		*long-term users vs. rec;
		if comp0_ltuser2021 = 1 then lt_rec = 0;
		if comp4_recreatn2021 = 1 then lt_rec = 1;

		*long-term users vs. quit;
		if comp0_ltuser2021 = 1 then lt_quit = 0;
		if comp5_quit2021 = 1 then lt_quit = 1;

		*Making a grouping var for three independent groups: long-term cannabis, nonusers, tobacco, quitters;
		if  comp0_ltuser2021 = 1 then group = 1; *long-term cannabis;
		if  comp1_noSUD2021 = 1 then group = 2; *nonusers;
		if  Comp2_LTtob2021_REV = 1 then group = 3; *long-term tobacco;
		if comp5_quit2021 then group = 4; *cannabis quitters; 

		
	*creating a tobacco quitters group;
	nmisstob1845 = nmiss (of dxtob18, dxtob21, dxtob26, dxtob32, dxTob38Dsm4, dxTob45Dsm4);
	if nmisstob1845 ne . and nmisstob1845 le 2 then counttob1845 = sum (of dxtob18, dxtob21, dxtob26, dxtob32, dxTob38Dsm4, dxTob45Dsm4);
	if counttob1845 ge 1 and cigsday45 = 0 then tobquit = 1;
	if cigsday45 = . then tobquit = .;



data age38;
	set seven; run;
		proc sort; by snum;

data age38a; *bringing in age 38 lt can group;
	set methyl.canadd;
	proc sort; by snum;

data age38b;
 	merge age38 age38a;
	
	*lt can vs non - age 38;
	if ltcan_38 = 1 then d1=1;
	if group = 2 then d1=0;

	*table s6 means;
	proc sort data=age38b; by d1; proc means data=age38b; var cg05575921 cg21566642 cg03636183 cg21161138 cg01940273 cg23079012; by d1; where pc1_age38 ne . and sex ne . and ntrphlsp38 ne . and lymphcytsp38 ne . and mncytsp38 ne . and esinphlsp38 ne . and bsphlsp38 ne .;  
	run;

	proc freq; tables d1; where  cg05575921 ne . and ntrphlsp38 ne . ; run; *shows analytic N (n = 76) lt can at age 38;

	*Table s6 - testing the 9 cpgs from age 45 that were robust to cov adjustment. only 6 are present at age 38;
	proc robustreg data=age38b method=mm (asympcov=H1) plots=qqplot; model cg05575921= d1  pc1_age38--pc32_age38 sex ntrphlsp38 lymphcytsp38 mncytsp38 esinphlsp38 bsphlsp38;  
		output out=residuals p=predicted r=residual;
		run;
		quit; 

	proc robustreg data=age38b method=mm (asympcov=H1) plots=qqplot; model cg21566642= d1  pc1_age38--pc32_age38 sex ntrphlsp38 lymphcytsp38 mncytsp38 esinphlsp38 bsphlsp38;  
		output out=residuals p=predicted r=residual;
		run;
		quit; 

	proc robustreg data=age38b method=mm (asympcov=H1) plots=qqplot; model cg03636183= d1  pc1_age38--pc32_age38 sex ntrphlsp38 lymphcytsp38 mncytsp38 esinphlsp38 bsphlsp38;  
		output out=residuals p=predicted r=residual;
		run;
		quit; 

	proc robustreg data=age38b method=mm (asympcov=H1) plots=qqplot; model cg21161138= d1   pc1_age38--pc32_age38 sex ntrphlsp38 lymphcytsp38 mncytsp38 esinphlsp38 bsphlsp38;  
		output out=residuals p=predicted r=residual;
		run;
		quit; 


	proc robustreg data=age38b method=mm (asympcov=H1) plots=qqplot; model cg01940273= d1   pc1_age38--pc32_age38 sex ntrphlsp38 lymphcytsp38 mncytsp38 esinphlsp38 bsphlsp38;  
		output out=residuals p=predicted r=residual;
		run;
		quit; 


	proc robustreg data=age38b method=mm (asympcov=H1) plots=qqplot; model cg17739917= d1   pc1_age38--pc32_age38 sex ntrphlsp38 lymphcytsp38 mncytsp38 esinphlsp38 bsphlsp38;  
		output out=residuals p=predicted r=residual;
		run;
		quit; 

	proc robustreg data=age38b method=mm (asympcov=H1) plots=qqplot; model cg05086879= d1   pc1_age38--pc32_age38 sex ntrphlsp38 lymphcytsp38 mncytsp38 esinphlsp38 bsphlsp38;  
		output out=residuals p=predicted r=residual;
		run;
		quit; 

	proc robustreg data=age38b method=mm (asympcov=H1) plots=qqplot; model cg02978227= d1   pc1_age38--pc32_age38 sex ntrphlsp38 lymphcytsp38 mncytsp38 esinphlsp38 bsphlsp38;  
		output out=residuals p=predicted r=residual;
		run;
		quit; 

	proc robustreg data=age38b method=mm (asympcov=H1) plots=qqplot; model cg23079012= d1   pc1_age38--pc32_age38 sex ntrphlsp38 lymphcytsp38 mncytsp38 esinphlsp38 bsphlsp38;  
		output out=residuals p=predicted r=residual;
		run;
		quit; 
