/*******************************************
Qualitative group comparisons

This file uses robust ANOVA with MM est
to test group comparisons.

THIS COMPARES LT CAN (group 1), NON (group 2), LT TOB (group 3)

Robust reg does not have the option to list
multiple dependent variables, so i had to trick 
SAS into running many at once by transposing the CpG
vars from wide to long, and then running analyses by 
the long CpG variable


CRUDE MEANS ARE SHOWN AT BEGINNGING OF THIS FILE
RESULTS OF STATS TESTS ARE SHOWN AT THE END OF THE FILE

THE RESULTS FROM MM ESTIMATION VARY WITH EACH RUN.
THEREFORE, USE THE OUTPUTTED DATA SET 
QUAL_COMBINED TO GO BACK AND LOOK AT RESULTS
*************************************************/

data adddx45;*bringing in can, alc, and drug dep at age 45;
	set canolder.cannabisprojfile_2020_07_31; 
	keep snum  dxmar45 dxal45d4 dxdrg45m;
	proc sort; by snum;

data addearlypot; *bringing in can use before age 18;
	set dundata.cannabis_IQ_final_July2012 ;
	keep snum pot1315;
	proc sort; by snum;

data addfh; *bringing in family history of substance dependence;
	set canolder.cannabiscogn_rev18Aug21;
	keep snum prsub1;
	proc sort; by snum;

data maindata; *bringing in main data;
	retain snum cg08923376 cg01472075  cg10616121 cg13337507 cg16609878 cg24339197; 
	set methyl.mastersept24;
	proc sort; by snum;

data age38lt; *bringing in age 38 lt can group;
	set methyl.canadd;
	proc sort; by snum;

	ods graphics off;

data mjlife; *Bringing in data on mj frequency from 18-45;
	set canolder.cannabis_firstlook_2020_07_27;
	keep snum marfrq18t marfrq21t marfrq26t marfrq32t marfrq38  MarFreq45; 
	proc sort; by snum;


data addrecdrg; *bringing in recurrent drug dx;
	set canolder.cannabiscombined_2021_02_11; 
	keep snum recdrg2645_2;
	proc sort; by snum;

data one;
	merge adddx45 addearlypot addfh maindata age38lt mjlife addrecdrg; by snum; 
	
	*creating a variable indicating snums with more than 2/3 of the Cpgs used here. n=818;
	if nmiss (of cg08923376 --cg27564939) ge 164 then havecpgdata = 0;
	else havecpgdata=1;
	
	*creating a var indicating not missing sex, PCs, and white blood cell counts, so crude means match N in robustreg;
	if pc1_age45--pc32_age45 ne . and sex ne . and ntrphlsp45 ne . and lymphcytsp45 ne . and mncytsp45 ne . and esinphlsp45 ne . and bsphlsp45 ne .
	then nm_qual = 1;

	*creating a var for lifetime history of sub dx;
	if TotMar1845lmt ge 1 then lifetime_candx = 1; else if totmar1845lmt ne . then lifetime_candx = 0;
	if TotAlc1845lmt ge 1 then lifetime_alcdx = 1; else if totalc1845lmt ne . then lifetime_alcdx = 0;
	if TotTob1845lmt ge 1 then lifetime_tobdx = 1; else if tottob1845lmt ne . then lifetime_tobdx = 0;
	if TotDrg2645lmt ge 1 then lifetime_drgdx = 1; else if totdrg2645lmt ne . then lifetime_drgdx = 0;

	*creating a var for reg can use at age 45;
	if marfreq45 ge 208 then regmar_45 = 1;
	if marfreq45 ne . and marfreq45 lt 208 then regmar_45 = 0;

	*Weekly cannabis at age 45;
	if marfreq45 ge 52 then weeklycannabis =1;
	if marfreq45 ne . and marfreq45 lt 52 then  weeklycannabis = 0;
	
	*to cross check with results of transpose approach.  obs=311;
	if group = 1 then d1 = 1;
	if group IN (2,3) then d1=0; *d1 =long term cannabis vs non and tob;

	if group = 3 then d2 = 1; 
	if group IN (1,2) then d2=0; *d2 = tob vs long term can and non;

	*when did quitter quit - reviewer 3;

		array higgs (6) marfrq18t marfrq21t marfrq26t marfrq32t marfrq38  MarFreq45;
		array higgs2 (6) any18 any21 any26 any32 any38 any45;
		do i = 1 to 6;
		if higgs(i) gt 0 then higgs2(i) = 1;
		if higgs(i) =0 then higgs2(i) = 0;
		end;

	
			if group=4 and any38 = 1 then quitby45 = 1;*n=23 quit by age 45;
			if group=4 and any38 = 0 and any32 = 1 then quitby38 = 1; *n=14 quit by age 38;
			if group=4 and any38 = 0 and any32 = 0 and any26 = 1 then quitby32 =1; *n=12 quit by age 32;
			if group=4 and any38 = 0 and any32 = 0 and any26 = 0 and any21=1 then quitby26 = 1;;  *n=6 quit by age 26;
			if group=4 and any38 = 0 and any32 = 0 and any26 = 0 and any21=0 and any18 =1 then quitby21=1; *n=4 quit by age 21;
			if group=4 and quitby45 ne 1 and quitby38 ne 1 and quitby32 ne 1 and quitby26 ne 1 and  quitby21 ne 1 then leftover=1; 

			if group = 1 then quitgroup = 0; *has not quit;
			if quitby45 = 1 then quitgroup = 1; *quit by 45;
			if quitby38=1 then quitgroup = 2;*quit by 38;
			if quitby32 = 1 then quitgroup = 3; *quit by 32;
			if quitby26 = 1 or quitby21 = 1 then quitgroup = 4; *quit by 21/26;
		
			proc print; var marfrq18t any18 marfrq21t any21 marfrq26t any26 marfrq32t any32 marfrq38 any38 marfreq45 any45; run;
			proc freq; tables quitgroup; run;
			proc freq; tables quitby45 quitby38 quitby32 quitby26 quitby21 leftover; run;
			proc print; var snum any18 any21 any26 any32 any38 any45; where leftover = 1; run;

			proc sort; by quitgroup;
			proc means; var marfrq18t marfrq21t marfrq26t marfrq32t marfrq38  MarFreq45; by quitgroup; run;

			proc print; var snum any18 any21 any26 any32 any38 any45; where group=4 and any38 = 1; run; *n=23 quit by age 45;
			proc print; var snum any18 any21 any26 any32 any38 any45; where group=4 and any38 = 0 and any32 = 1; run; *n=14 quit by age 38;
			proc print; var snum any18 any21 any26 any32 any38 any45; where group=4 and any38 = 0 and any32 = 0 and any26 = 1; run; *n=12 quit by age 32;
			proc print; var snum any18 any21 any26 any32 any38 any45; where group=4 and any38 = 0 and any32 = 0 and any26 = 0 and any21=1; run; *n=6 quit by age 26;
			proc print; var snum any18 any21 any26 any32 any38 any45; where group=4 and any38 = 0 and any32 = 0 and any26 = 0 and any21=0 and any18 =1; run; *n=4 quit by age 21;

			proc glm; class quitgroup; model cg05575921 cg21566642 cg03636183 cg21161138 cg01940273 cg17739917 cg05086879 cg02978227 cg23079012 = quitgroup; means quitgroup; run;

	/*	*Table S9 - Reviewer 3: testing quit length as a predictor of 9 robust CpGs;
	
			ods output ParameterEstimates=quitlength1;
		proc robustreg data=one method=mm (asympcov=H1); model cg05575921= quitgroup pc1_age45--pc32_age45 sex ntrphlsp45 lymphcytsp45 mncytsp45 esinphlsp45 bsphlsp45 SESchildhd lscuw311 prsub1 smkdxgp1845 alcdxgrp1845 ;  run;
		quit; 

					ods output ParameterEstimates=quitlength2;
		proc robustreg data=one method=mm (asympcov=H1); model cg21566642= quitgroup pc1_age45--pc32_age45 sex ntrphlsp45 lymphcytsp45 mncytsp45 esinphlsp45 bsphlsp45 SESchildhd lscuw311 prsub1 smkdxgp1845 alcdxgrp1845 recdrg2645_2 ;  run;
		quit; 

					ods output ParameterEstimates=quitlength3;
		proc robustreg data=one method=mm (asympcov=H1); model cg03636183= quitgroup pc1_age45--pc32_age45 sex ntrphlsp45 lymphcytsp45 mncytsp45 esinphlsp45 bsphlsp45 SESchildhd lscuw311 prsub1 smkdxgp1845 alcdxgrp1845 recdrg2645_2;  run;
		quit; 

					ods output ParameterEstimates=quitlength4;
		proc robustreg data=one method=mm (asympcov=H1); model cg21161138= quitgroup pc1_age45--pc32_age45 sex ntrphlsp45 lymphcytsp45 mncytsp45 esinphlsp45 bsphlsp45 SESchildhd lscuw311 prsub1 smkdxgp1845 alcdxgrp1845 recdrg2645_2;  run;
		quit; 

					ods output ParameterEstimates=quitlength5;
		proc robustreg data=one method=mm (asympcov=H1); model cg01940273= quitgroup pc1_age45--pc32_age45 sex ntrphlsp45 lymphcytsp45 mncytsp45 esinphlsp45 bsphlsp45 SESchildhd lscuw311 prsub1 smkdxgp1845 alcdxgrp1845 recdrg2645_2 ;  run;
		quit; 

					ods output ParameterEstimates=quitlength6;
		proc robustreg data=one method=mm (asympcov=H1); model cg17739917 = quitgroup pc1_age45--pc32_age45 sex ntrphlsp45 lymphcytsp45 mncytsp45 esinphlsp45 bsphlsp45 SESchildhd lscuw311 prsub1 smkdxgp1845 alcdxgrp1845 recdrg2645_2;  run;
		quit; 

					ods output ParameterEstimates=quitlength7;
		proc robustreg data=one method=mm (asympcov=H1); model cg05086879 = quitgroup pc1_age45--pc32_age45 sex ntrphlsp45 lymphcytsp45 mncytsp45 esinphlsp45 bsphlsp45 SESchildhd lscuw311 prsub1 smkdxgp1845 alcdxgrp1845 recdrg2645_2;  run;
		quit; 

					ods output ParameterEstimates=quitlength8;
		proc robustreg data=one method=mm (asympcov=H1); model cg02978227 = quitgroup pc1_age45--pc32_age45 sex ntrphlsp45 lymphcytsp45 mncytsp45 esinphlsp45 bsphlsp45 SESchildhd lscuw311 prsub1 smkdxgp1845 alcdxgrp1845 recdrg2645_2 ;  run;
		quit; 

					ods output ParameterEstimates=quitlength9;
		proc robustreg data=one method=mm (asympcov=H1); model cg23079012= quitgroup pc1_age45--pc32_age45 sex ntrphlsp45 lymphcytsp45 mncytsp45 esinphlsp45 bsphlsp45 SESchildhd lscuw311 prsub1 smkdxgp1845 alcdxgrp1845 recdrg2645_2;  run;
		quit; 
*/


	proc freq; tables seennot45*havecpgdata; run; *of the 938 seen at 45, 818 have methylation data;

	*agreement between age 38 and age 45 lt cannabis groups;
	proc freq; tables ltcan_38*ltcan_45/agree; where havecpgdata=1;run; *reviewer 1;


	*Supp Fig 2 on groups and overlap;
	proc freq; tables group tobquit; where havecpgdata=1; 
	proc freq; tables group*tobquit; where havecpgdata=1;run;
	
	
	*methods section on groups;

		*LT users;
		proc freq; table sex  pot1315 dailypotgp1845 potdxgrp1845 regmar_45 smkdxgp1845 smokes45;  where  group = 1 and havecpgdata=1;
		proc means n mean median; var marfreq45 dailypotgp1845; where group=1 and havecpgdata=1;
		run;
	

		*LT tobacco users;
		proc freq; tables sex marfreq45 weeklycannabis potdxgrp1845; where group=3 and havecpgdata=1; run;

		*Nonusers;
		proc freq; tables sex smkdxgp1845 potdxgrp1845 totalc1845lmt totdrg2645lmt; where group=2 and havecpgdata=1; run;

		*cannabis quitters;
		proc freq; tables sex marfreq45  smkdxgp1845 smokes45; where group=4 and havecpgdata=1; 
		proc freq; tables potdxgrp1845*dailypotgp1845; where group=4 and havecpgdata=1; run;

		*tobacco quitters;
		proc freq; tables sex cigsday45 smkdxgp1845 potdxgrp1845 regmar_45; where tobquit=1 and havecpgdata=1; run;


	*Table s4 on cohort and groups ;
	*Full cohort;
	proc means; var SESchildhd lscuw311 prsub1 marfreq45 ; where havecpgdata=1;
	proc freq; tables sex regmar_45 dxmar45  smokes45 dxtob45dsm4  weeklyalc45 dxal45d4 dxdrg45m dailypotgp1845 lifetime_candx lifetime_tobdx lifetime_alcdx lifetime_drgdx; where havecpgdata=1;
	run;

	*Long-term cannabis users;
	proc means; var SESchildhd lscuw311 prsub1 marfreq45 ; where havecpgdata=1 and group=1;
	proc freq; tables  sex regmar_45 dxmar45  smokes45 dxtob45dsm4  weeklyalc45 dxal45d4 dxdrg45m dailypotgp1845 lifetime_candx lifetime_tobdx lifetime_alcdx lifetime_drgdx; where havecpgdata=1 and group=1;
	run;

	*Nonusers;
	proc means; var SESchildhd lscuw311 prsub1 marfreq45 ; where havecpgdata=1 and group=2;
	proc freq; tables sex regmar_45 dxmar45  smokes45 dxtob45dsm4  weeklyalc45 dxal45d4 dxdrg45m dailypotgp1845 lifetime_candx lifetime_tobdx lifetime_alcdx lifetime_drgdx; where havecpgdata=1 and group=2;
	run;

	*Long-term tobacco users;
	proc means; var SESchildhd lscuw311 prsub1 marfreq45 ; where havecpgdata=1 and group=3;
	proc freq; tables sex regmar_45 dxmar45  smokes45 dxtob45dsm4  weeklyalc45 dxal45d4 dxdrg45m dailypotgp1845 lifetime_candx lifetime_tobdx lifetime_alcdx lifetime_drgdx; where havecpgdata=1 and group=3;
	run;

	*Cannabis quitters;
	proc means; var SESchildhd lscuw311 prsub1 marfreq45 ; where havecpgdata=1 and group=4;
	proc freq; tables sex regmar_45 dxmar45  smokes45 dxtob45dsm4  weeklyalc45 dxal45d4 dxdrg45m dailypotgp1845 lifetime_candx lifetime_tobdx lifetime_alcdx lifetime_drgdx; where havecpgdata=1 and group=4;
	run;

	*Tobacco quitters;
	proc means; var SESchildhd lscuw311 prsub1 marfreq45 ; where havecpgdata=1 and tobquit=1;
	proc freq; tables sex  regmar_45 dxmar45  smokes45 dxtob45dsm4  weeklyalc45 dxal45d4 dxdrg45m dailypotgp1845 lifetime_candx lifetime_tobdx lifetime_alcdx lifetime_drgdx; where havecpgdata=1 and tobquit=1;
	run;

	*Crude means for table 1. I got the order of the CPGs from running the regressions below.
	They are sorted by p value of can vs non and then tob vs non.;
	proc sort; by group; 
	proc means n mean std clm; var  
cg05575921
cg21566642
cg01940273
cg17739917
cg21161138
cg03636183
cg25189904
cg18110140
cg05086879
cg23079012
cg02978227
cg18387338
cg09935388
cg23916896
cg05009104
cg26286198
cg15088912
cg19089201
cg17350345
cg26764244
cg22124000
cg14219071
cg12510044
cg21852554
cg16276224
cg21491555
cg15802887
cg13966609
cg18880190
cg22572071
cg04864586
cg05922265
cg03093806
cg27055782
cg08153621
cg00813162
cg27209861
cg11175241
cg08688629
cg17464820
cg15204119
cg08839808
cg03802952
cg02235741
cg10328583
cg17250160
cg12028375
cg10620881
cg09501516
cg21960184
cg09825346
cg25069772
cg21998512
cg00250546
cg20958467
;
    by group; where nm_qual = 1;
	run;


*Test of variance diff (ANOVA assumption of equal variances) for LT vs. non;
proc glm data=one outstat= hovtestrev1;
 class group;
 model cg08923376--cg27564939 = group ; means group  / hovtest = levene; where group IN (1,2); run;

	*count of LT vs non CpGs with different variances at alpha adjusted for 246 tests;
	data hov_ltnon;
	set hovtestrev1; if _SOURCE_="group" and _type_ = "SS3";
	if prob <.00021; *CpGs with p <.00021, the threshold for multiple testing, indicating sig diff in variance of CpGs for LT vs non;
	run;

*Test of variance diff for LT vs Tob;
proc glm data=one outstat= hovtestrev1a;
 class group;
 model cg08923376--cg27564939 = group ; means group  / hovtest = levene; where group IN (1,3); run;

 	*count of LT vs tob CpGs with different variances at alpha adjusted for 246 tests;
	data hov_ttob;
	set hovtestrev1a; if _SOURCE_="group" and _type_="SS3";
	if prob <.00021; 
	run;

	*crude means for figure 3/table s7 - can quit non;
	proc means data=one n mean clm; var
cg05575921
cg21566642
cg03636183
cg21161138
cg01940273
cg17739917
cg05086879
cg02978227
cg23079012
; by group; where nm_qual = 1; run;

	
	*Crude means for figure 3/table s8 - tob quit v non;

	proc sort data=one; by group;
	proc means data=one n mean clm; var
cg05575921
cg21566642
cg03636183
cg21161138
cg01940273
cg17739917
cg05086879
cg02978227
cg23079012
cg18110140
cg09935388
cg25189904
cg05009104
cg23916896
cg18387338
cg15088912
cg19089201
; by group; where nm_qual = 1; run;

 proc means data=one n mean clm; var
cg05575921
cg21566642
cg03636183
cg21161138
cg01940273
cg17739917
cg05086879
cg02978227
cg23079012
cg18110140
cg09935388
cg25189904
cg05009104
cg23916896
cg18387338
cg15088912
cg19089201; where tobquit = 1 and nm_qual = 1;
 run;




	proc sort; by snum;

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
 
*need covariates from original data. dropping cgs from original data set to merge this with the transposed data and not be confusing that we have long and wide cpgs;
data five;
	set one; drop cg08923376--cg27564939;
	proc sort; by snum;

*merging transposed data with original data with untransposed (i.e., wide) covariates;

	/********************************
	LT v NON
	*********************************/
data six;
	merge four five; by snum; 

	/*********************************
	dummy coding group:
	D1 - LT CAN v NON
	D2 - LT TOB V NON

	***************************************/
	if group = 1 then d1 = 1;
	if group IN (2,3) then d1=0;

	if group = 3 then d2 = 1;
	if group IN (1,2) then d2=0;


	proc sort; by cpg_index;

	/*Table 1 - see results in methyl.qualcombined_sept24*/
	*ods exclude none; 
	ods exclude all; *this code speeds analyses. ran this as exclude none to confirm analytic Ns. N=311 (sum of lt, non, tob), which matches Ns from table of crude means;
	ods output ParameterEstimates=canvall;
		proc robustreg data=six method=mm (asympcov=H1) plots=qqplot; model cg1= d1 d2  pc1_age45--pc32_age45 sex ntrphlsp45 lymphcytsp45 mncytsp45 esinphlsp45 bsphlsp45;  by cpg_index;
		output out=residuals p=predicted r=residual;
		run;
		quit; 

		proc univariate data=one; var cg21566642; histogram; run;

	*Diagnostics;
		proc sort data=residuals; by cpg_index; 
	proc sgplot data=residuals;
   scatter x=predicted y=residual;
   refline 0 / axis=y lineattrs=(color=red); where cg1 ne .; by cpg_index; run;

data seven; *Selecting only can v non;
	set canvall; if Parameter="d1";

	proc sort; by cpg_index;

	run;

*need to get name of cpg, so keeping name for one participant (because long file means only need one) to merge back with results;
data eight; 
	set four;
	keep _name_ cpg_index;
	if snum=[removed to avoid publishing SNUMs];
	proc sort; by cpg_index;

data nine; *this is  data set of can vs non robust reg tests;
	merge seven eight; by cpg_index; 

	run;
	
	/***************************************
	LT tob vs non- table 1

	*****************************************/
data ten;
	set canvall;
	if Parameter="d2";

	proc sort; by cpg_index;

	run;

data eleven;
	merge ten eight; by cpg_index; *this is  data set of tob vs non robust reg tests;

	

/**********************************************

		LT Can v LT Tob
************************************************/
data twelve;
	set six;

	*Dummy code comparing can and non to tob. select d1_a for can vs tob;

	if group = 1 then d1_a = 1;
	if group IN (2,3) then d1_a=0;
	
	if group = 2 then d2_a = 1;
	if group IN (1,3) then d2_a = 0;

	ods exclude all;
	ods output ParameterEstimates=canvtob;
		proc robustreg data=twelve method=mm (asympcov=H1); model cg1= d1_a d2_a  pc1_age45--pc32_age45 sex ntrphlsp45 lymphcytsp45 mncytsp45 esinphlsp45 bsphlsp45; by cpg_index;
	run;
		quit; 

	proc sort; by cpg_index;

data thirteen;
	set work.canvtob;
	if Parameter="d1_a"; *can v tob;

	proc sort; by cpg_index;

data fourteen;
	merge thirteen eight; by cpg_index; *this is dataset of can v tob robust reg tests;



	/*****************************************
	LT Can vs Quit - 
	Table S7
	*******************************************/

data fifteen;
	set six; 


	*creating dummies to compare Quit to LT and Non.
	So quit is the reference group.
	canquit_d1 = LT vs. Quit
	Canquit_d2 = Non vs. Quit;

	if group =1 then canquit_d1 = 1;
	if group IN (2,4) then canquit_d1 = 0;

	if group =2 then canquit_d2 = 1;
	if group IN (1,4) then canquit_d2 =0;


	proc sort; by cpg_index;


	ods exclude all;
	ods output ParameterEstimates=canquit;
		proc robustreg data=fifteen method=mm (asympcov=H1); model cg1= canquit_d1 canquit_d2  pc1_age45--pc32_age45 sex ntrphlsp45 lymphcytsp45 mncytsp45 esinphlsp45 bsphlsp45;  by cpg_index;
	run;
		quit; 

data sixteen; *bringing in the cpg name;
	merge canquit eight; by cpg_index;


data seventeen;
	set sixteen; *this is dataset of robust reg results of can quit vs lt can;
	if Parameter="canquit_d1"; *selecting LT can vs Quit;

	
	*Reverse scoring param est so it is quit minus lt;
	rec_estimate=estimate*-1;
	rec_uppercl= lowercl*-1;
	rec_lowercl=uppercl*-1;

	


data eighteen;
	set sixteen; *this is data set of non vs. can quit;
	if Parameter="canquit_d2"; *selecting NOn vs Quit;

	
	*Reverse scoring param est so it is quit minus non;
	rec_estimate=estimate*-1;
	rec_uppercl= lowercl*-1;
	rec_lowercl=uppercl*-1;



/***********************************
	Tob v Quit
 Table s8
**************************/
data nineteen;
	set six;

	
	*tob quit vs. lt tob and non. 
	tob quit is ref group.
	d1 = lt tob vs quit
	d2=non vs quit;
	if group=3 then tobquit_d1 = 1;
	if group =2 or tobquit = 1 then tobquit_d1 = 0;

	if group=2 then tobquit_d2 = 1;
	if group=3 or tobquit=1 then tobquit_d2 = 0;

	proc sort; by cpg_index;

	ods exclude all;
	ods output ParameterEstimates=tobquit;
		proc robustreg data=nineteen method=mm (asympcov=H1); model cg1= tobquit_d1 tobquit_d2 pc1_age45--pc32_age45 sex ntrphlsp45 lymphcytsp45 mncytsp45 esinphlsp45 bsphlsp45;  by cpg_index;
	run;
		quit; 

data twenty;
	merge tobquit eight;
	by cpg_index;
	
data twentyone;
	set twenty; *this is data set of robust reg results of tob quitters vs long-term tob;
	if Parameter="tobquit_d1"; *this is lt tob minus quitters;

	*Reverse scoring param est so it is quit minus lt tob;
	rec_estimate=estimate*-1;
	rec_uppercl= lowercl*-1;
	rec_lowercl=uppercl*-1;


data twentytwo;
	set twenty; *this is dataset of robust reg results of tob quitters vs non;
	if Parameter="tobquit_d2"; *this is non minus tob quitters;

	*Reverse scoring param est so it is quit minus non;
	rec_estimate=estimate*-1;
	rec_uppercl= lowercl*-1;
	rec_lowercl=uppercl*-1;



*Merging all results of qual comparisons. To do that, I must rename the parameter estimates and p values;
data qual_cannon;
	set nine;

	Can_minusnon_est=estimate;
	can_minusnon_lcl = lowercl;
	can_minusnon_ucl = uppercl;
	can_minusnon_p = probchisq;

drop parameter df estimate stderr lowercl uppercl chisq probchisq ;

proc sort; by cpg_index;
run;

data qual_tobnon;
	set eleven;

	tob_minusnon_est=estimate;
	tob_minusnon_lcl = lowercl;
	tob_minusnon_ucl = uppercl;
	tob_minusnon_p = probchisq;

drop parameter df estimate stderr lowercl uppercl chisq probchisq ;

proc sort; by cpg_index;
run;

data qual_cantob;
	set fourteen;

	can_minustob_est=estimate;
	can_minustob_lcl = lowercl;
	can_minustob_ucl = uppercl;
	can_minustob_p = probchisq;

drop parameter df estimate stderr lowercl uppercl chisq probchisq ;

proc sort; by cpg_index;
run;

data qual_canquitcan;
	set seventeen;

	canquit_minuscan_est=rec_estimate;
	canquit_minuscan_lcl=rec_lowercl;
	canquit_minuscan_ucl=rec_uppercl;
	canquit_minuscan_p = probchisq; 

drop parameter df estimate stderr lowercl uppercl chisq probchisq  rec_estimate rec_lowercl rec_uppercl;

proc sort; by cpg_index;
run;

data qual_canquitnon;
	set eighteen;

	canquit_minusnon_est=rec_estimate;
	canquit_minusnon_lcl=rec_lowercl;
	canquit_minusnon_ucl=rec_uppercl;
	canquit_minusnon_p = probchisq; 

drop parameter df estimate stderr lowercl uppercl chisq probchisq rec_estimate rec_lowercl rec_uppercl;

proc sort; by cpg_index;
run;

data qual_tobquittob;
	set twentyone;

	tobquit_minustob_est=rec_estimate;
	tobquit_minustob_lcl=rec_lowercl;
	tobquit_minustob_ucl=rec_uppercl;
	tobquit_minustob_p = probchisq; 

drop parameter df estimate stderr lowercl uppercl chisq probchisq  rec_estimate rec_lowercl rec_uppercl; 

proc sort; by cpg_index;
run;

data qual_tobquitnon;
	set twentytwo;

	tobquit_minusnon_est=rec_estimate;
	tobquit_minusnon_lcl=rec_lowercl;
	tobquit_minusnon_ucl=rec_uppercl;
	tobquit_minusnon_p = probchisq; 

drop parameter df estimate stderr lowercl uppercl chisq probchisq  rec_estimate rec_lowercl rec_uppercl;

proc sort; by cpg_index;
run;

*Creating one file of all qual findings ;
/*
data methyl.qual_combined_sept24; *THIS file is super important because robust reg results change with every run. So, all results from robust reg analyses must be based off this hard
							data set of results;
	merge qual_cannon qual_canquitcan qual_canquitnon qual_cantob qual_tobnon qual_tobquitnon qual_tobquittob; by cpg_index; run;
*/

data ins;
	set methyl.qual_combined_sept24;

	*Results paragraph 1 - summary of sig <.05 and p<.00020 findings;
	*Table S5;
	ods exclude none;
	proc sort; by can_minusnon_p; 
	proc print; var _name_ can_minusnon_est can_minusnon_lcl can_minusnon_ucl can_minusnon_p; where can_minusnon_p <.05; run; 
	proc print; var _name_ can_minusnon_est can_minusnon_lcl can_minusnon_ucl can_minusnon_p; where can_minusnon_p <.00020; run; 

	proc sort; by can_minusnon_p tob_minusnon_p;
	proc print; var _name_ tob_minusnon_est tob_minusnon_lcl tob_minusnon_ucl tob_minusnon_p; where tob_minusnon_p <.05; run; 
	proc print; var _name_ tob_minusnon_est tob_minusnon_lcl tob_minusnon_ucl tob_minusnon_p; where tob_minusnon_p <.00020; run; 

	proc sort; by can_minusnon_p;
	proc print; var _name_ can_minustob_est can_minustob_lcl can_minustob_ucl can_minustob_p; where can_minustob_p <.05; run; 
	proc print; var _name_ can_minustob_est can_minustob_lcl can_minustob_ucl can_minustob_p; where can_minustob_p <.00020; run; 
	run;

	ods exclude none;

	*Table 1 results for lt can, lt tob, nonuser comparisons, adjusted mean diffs, for the p<.05 sig can vs non, sorted by can vs non p;
	proc sort; by can_minusnon_p;
	proc print; var _name_ can_minusnon_est can_minusnon_lcl can_minusnon_ucl can_minusnon_p tob_minusnon_est 
				tob_minusnon_lcl tob_minusnon_ucl tob_minusnon_p can_minustob_est can_minustob_lcl can_minustob_ucl can_minustob_p;  where can_minusnon_p < .05; 

	run;

	*table 1 results for the 1 additional sig tob vs non;
	proc sort; by tob_minusnon_p;
		proc print; var _name_ tob_minusnon_est 
				tob_minusnon_lcl tob_minusnon_ucl tob_minusnon_p ;  where can_minusnon_p ge .05 and tob_minusnon_p <.00020;


	run;

	*Figure 3 tests, supp tables 7 and 8  9.25.24;
	/*This is the order - the first 9 are for cannabis (sig qual and quant) and the second 17 are for tob.


cannabis
1	cg05575921
2	cg21566642
3	cg03636183
4	cg21161138
5	cg01940273
6	cg17739917
7	cg05086879
8	cg02978227
9	cg23079012

tobacco
cg05575921
cg21566642
cg01940273
cg03636183
cg17739917
cg21161138
cg18110140
cg09935388
cg25189904
cg05086879
cg02978227
cg23079012
cg05009104
cg23916896
cg18387338
cg15088912
cg19089201

	***********************************/
	*Can quitters;
	proc print data=ins; var _name_ canquit_minuscan_est canquit_minuscan_lcl canquit_minuscan_ucl canquit_minuscan_p canquit_minusnon_est canquit_minusnon_lcl 
					canquit_minusnon_ucl canquit_minusnon_p; where _name_="cg05575921"; run;
	proc print data=ins; var  _name_ canquit_minuscan_est canquit_minuscan_lcl canquit_minuscan_ucl canquit_minuscan_p canquit_minusnon_est canquit_minusnon_lcl 
					canquit_minusnon_ucl canquit_minusnon_p; where _name_="cg21566642"; run;
	proc print data=ins; var _name_ canquit_minuscan_est canquit_minuscan_lcl canquit_minuscan_ucl canquit_minuscan_p canquit_minusnon_est canquit_minusnon_lcl 
					canquit_minusnon_ucl canquit_minusnon_p; where _name_="cg03636183"; run;
	proc print data=ins; var _name_ canquit_minuscan_est canquit_minuscan_lcl canquit_minuscan_ucl canquit_minuscan_p canquit_minusnon_est canquit_minusnon_lcl 
					canquit_minusnon_ucl canquit_minusnon_p; where _name_="cg21161138"; run;
	proc print data=ins; var _name_ canquit_minuscan_est canquit_minuscan_lcl canquit_minuscan_ucl canquit_minuscan_p canquit_minusnon_est canquit_minusnon_lcl 
					canquit_minusnon_ucl canquit_minusnon_p; where _name_="cg01940273"; run;
	proc print data=ins; var _name_ canquit_minuscan_est canquit_minuscan_lcl canquit_minuscan_ucl canquit_minuscan_p canquit_minusnon_est canquit_minusnon_lcl 
					canquit_minusnon_ucl canquit_minusnon_p; where _name_="cg17739917"; run;
	proc print data=ins; var _name_ canquit_minuscan_est canquit_minuscan_lcl canquit_minuscan_ucl canquit_minuscan_p canquit_minusnon_est canquit_minusnon_lcl 
					canquit_minusnon_ucl canquit_minusnon_p; where _name_="cg05086879"; run;
	proc print data=ins; var  _name_ canquit_minuscan_est canquit_minuscan_lcl canquit_minuscan_ucl canquit_minuscan_p canquit_minusnon_est canquit_minusnon_lcl 
					canquit_minusnon_ucl canquit_minusnon_p; where _name_="cg02978227"; run;
	proc print data=ins; var  _name_ canquit_minuscan_est canquit_minuscan_lcl canquit_minuscan_ucl canquit_minuscan_p canquit_minusnon_est canquit_minusnon_lcl 
					canquit_minusnon_ucl canquit_minusnon_p; where _name_="cg23079012"; run;

	*Tob quitters;
	proc print data=ins; var _name_ tobquit_minustob_est tobquit_minustob_lcl tobquit_minustob_ucl tobquit_minustob_p tobquit_minusnon_est tobquit_minusnon_lcl 
					tobquit_minusnon_ucl tobquit_minusnon_p; where _name_="cg05575921"; run;
	proc print data=ins; var  _name_ tobquit_minustob_est tobquit_minustob_lcl tobquit_minustob_ucl tobquit_minustob_p tobquit_minusnon_est tobquit_minusnon_lcl 
					tobquit_minusnon_ucl tobquit_minusnon_p; where _name_="cg21566642"; run;
	proc print data=ins; var _name_ tobquit_minustob_est tobquit_minustob_lcl tobquit_minustob_ucl tobquit_minustob_p tobquit_minusnon_est tobquit_minusnon_lcl 
					tobquit_minusnon_ucl tobquit_minusnon_p; where _name_="cg03636183"; run;
	proc print data=ins; var _name_ tobquit_minustob_est tobquit_minustob_lcl tobquit_minustob_ucl tobquit_minustob_p tobquit_minusnon_est tobquit_minusnon_lcl 
					tobquit_minusnon_ucl tobquit_minusnon_p; where _name_="cg21161138"; run;
	proc print data=ins; var _name_ tobquit_minustob_est tobquit_minustob_lcl tobquit_minustob_ucl tobquit_minustob_p tobquit_minusnon_est tobquit_minusnon_lcl 
					tobquit_minusnon_ucl tobquit_minusnon_p; where _name_="cg01940273"; run;
	proc print data=ins; var _name_ tobquit_minustob_est tobquit_minustob_lcl tobquit_minustob_ucl tobquit_minustob_p tobquit_minusnon_est tobquit_minusnon_lcl 
					tobquit_minusnon_ucl tobquit_minusnon_p; where _name_="cg17739917"; run;
	proc print data=ins; var _name_ tobquit_minustob_est tobquit_minustob_lcl tobquit_minustob_ucl tobquit_minustob_p tobquit_minusnon_est tobquit_minusnon_lcl 
					tobquit_minusnon_ucl tobquit_minusnon_p; where _name_="cg05086879"; run;
	proc print data=ins; var  _name_ tobquit_minustob_est tobquit_minustob_lcl tobquit_minustob_ucl tobquit_minustob_p tobquit_minusnon_est tobquit_minusnon_lcl 
					tobquit_minusnon_ucl tobquit_minusnon_p; where _name_="cg02978227"; run;
	proc print data=ins; var  _name_ tobquit_minustob_est tobquit_minustob_lcl tobquit_minustob_ucl tobquit_minustob_p tobquit_minusnon_est tobquit_minusnon_lcl 
					tobquit_minusnon_ucl tobquit_minusnon_p; where _name_="cg23079012"; run;
	proc print data=ins; var  _name_ tobquit_minustob_est tobquit_minustob_lcl tobquit_minustob_ucl tobquit_minustob_p tobquit_minusnon_est tobquit_minusnon_lcl 
					tobquit_minusnon_ucl tobquit_minusnon_p; where _name_="cg18110140"; run;
	proc print data=ins; var  _name_ tobquit_minustob_est tobquit_minustob_lcl tobquit_minustob_ucl tobquit_minustob_p tobquit_minusnon_est tobquit_minusnon_lcl 
					tobquit_minusnon_ucl tobquit_minusnon_p; where _name_="cg09935388"; run;
	proc print data=ins; var  _name_ tobquit_minustob_est tobquit_minustob_lcl tobquit_minustob_ucl tobquit_minustob_p tobquit_minusnon_est tobquit_minusnon_lcl 
					tobquit_minusnon_ucl tobquit_minusnon_p; where _name_="cg25189904"; run;
	proc print data=ins; var  _name_ tobquit_minustob_est tobquit_minustob_lcl tobquit_minustob_ucl tobquit_minustob_p tobquit_minusnon_est tobquit_minusnon_lcl 
					tobquit_minusnon_ucl tobquit_minusnon_p; where _name_="cg05009104"; run;
	proc print data=ins; var  _name_ tobquit_minustob_est tobquit_minustob_lcl tobquit_minustob_ucl tobquit_minustob_p tobquit_minusnon_est tobquit_minusnon_lcl 
					tobquit_minusnon_ucl tobquit_minusnon_p; where _name_="cg23916896"; run;
	proc print data=ins; var  _name_ tobquit_minustob_est tobquit_minustob_lcl tobquit_minustob_ucl tobquit_minustob_p tobquit_minusnon_est tobquit_minusnon_lcl 
					tobquit_minusnon_ucl tobquit_minusnon_p; where _name_="cg18387338"; run;
	proc print data=ins; var  _name_ tobquit_minustob_est tobquit_minustob_lcl tobquit_minustob_ucl tobquit_minustob_p tobquit_minusnon_est tobquit_minusnon_lcl 
					tobquit_minusnon_ucl tobquit_minusnon_p; where _name_="cg15088912"; run;
	proc print data=ins; var  _name_ tobquit_minustob_est tobquit_minustob_lcl tobquit_minustob_ucl tobquit_minustob_p tobquit_minusnon_est tobquit_minusnon_lcl 
					tobquit_minusnon_ucl tobquit_minusnon_p; where _name_="cg19089201"; run;





/*Figure 3 - A, B, C*/

data ins2;
	set ins;
*Volcano plots;
	ods exclude none;
	
	/*can vs non*/
			 logp_canvnon = log10(can_minusnon_p)*-1;
			 if logp_canvnon >100 then logp_canvnon=100; *recoding one cpg p value (cg055..) to 100 bc it was outlier and fig looked odd;
			

		*pvalue groups -- bonferroni sig or not, used in sgplot to show different colors for each group;
		if can_minusnon_p <.00020 then pgroup_canvnon =1;
		if can_minusnon_p <.05 and can_minusnon_p ge .00020 then pgroup_canvnon =2; *sig at p<.05 but not at p<.0002;
		if can_minusnon_p ge .05 then pgroup_canvnon = 3;
		
							

		*making new var names to refer to log p for sgplot to run;
		if pgroup_canvnon = 1 then bonp_canvnon =logp;
		if pgroup_canvnon = 2 then sigp_canvnon = logp;
		if pgroup_canvnon = 3 then nonsigp_canvnon = logp;

		label bonp_canvnon ='-log10 (P-value)' sigp_canvnon='-log10 (P-value)' nonsigp_canvnon ='-log10 (P-value)';;

		
		*checking to make sure no issue with log p values;
		proc means n range min max data=ins2; var logp_canvnon; run;

		*volcano plot example: https://www.lexjansen.com/phuse-us/2020/dv/DV14.pdf 
						example https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4391864/;
		*volcano plot code: https://blogs.sas.com/content/graphicallyspeaking/2016/05/23/ctspedia-clinical-graphs-volcano-plot/;
		*used volcano plot code: https://communities.sas.com/t5/Graphics-Programming/Group-data-on-a-plot-by-multiple-grouping-variables/td-p/175541;

	ODS GRAPHICS/WIDTH=4IN HEIGHT=3IN;
	

	*Fig 2a - Can vs non ;
	proc sgplot; *correct number of obs;
	styleattrs datacontrastcolors=(blue darkred gray );
	  label logp_canvnon='-log10 (P-value)';
	  label can_minusnon_est='Mean Difference (Long-term Cannabis Users vs Non-Users)';
	  styleattrs datasymbols=(circlefilled);
	  scatter x=can_minusnon_est y=logp_canvnon /group=pgroup_canvnon /*datalabel=label */; *use this code to get datalabels;
	  refline 0 / axis=x lineattrs=(pattern=shortdash);
	  	  yaxis values = (0,10, 20, 30, 40, 50, 60,70, 80, 90, 100);
	  xaxis values=(-2.8, -2.6, -2.4, -2.2, -2.0, -1.8, -1.6, -1.4, -1.2, -1.0, -0.8, -0.6, -0.4, -0.2, 0, .2, .4, .6, 0.8, 1.0) ;
	  keylegend 'a' / across=1 position=right valueattrs=(size=6);
	run;
	


data ins3;
	set ins;
/*tob v non*/
	 logp_tobvnon = log10(tob_minusnon_p)*-1;
			 if logp_tobvnon >100 then logp_tobvnon=100; *recoding cpg p values to 100 bc outlier ;
			/*proc print; var _name_ tob_minusnon_p logp_tobvnon; run;*/
			if tob_minusnon_p <.00021 then label =_name_;

		*pvalue groups -- bonferroni sig or not, used in sgplot to show different colors for each group;
		if tob_minusnon_p <.00020 then pgroup_tobvnon =1;
		if tob_minusnon_p <.05 and tob_minusnon_p ge .00020 then pgroup_tobvnon =2; *sig at p<.05 but not at p<.0002;
		if tob_minusnon_p ge .05 then pgroup_tobvnon = 3;
		
							

		*making new var names to refer to log p for sgplot to run;
		if pgroup_tobvnon = 1 then bonp_tobvnon =logp_tobvnon;
		if pgroup_tobvnon = 2 then sigp_tobvnon = logp_tobvnon;
		if pgroup_tobvnon = 3 then nonsigp_tobvnon = logp_tobvnon;

		label bonp_tobvnon ='-log10 (P-value)' sigp_tobvnon='-log10 (P-value)' nonsigp_tobvnon ='-log10 (P-value)';;

		
		*checking to make sure no issue with log p values;
		proc means n range min max data=ins3; var logp_tobvnon; run;

	ODS GRAPHICS/WIDTH=4IN HEIGHT=3IN;
	

	*Fig 2b - Tob vs non ;
	proc sgplot ; 
	styleattrs datacontrastcolors=(blue darkred gray);
	  label logp_tobvnon='-log10 (P-value)';
	  label tob_minusnon_est='Mean Difference (Long-term Tobacco Users vs Non-Users)';
	  styleattrs datasymbols=(circlefilled);
	  scatter x=tob_minusnon_est y=logp_tobvnon /group=pgroup_tobvnon /*datalabel=label */; *use this code to get datalabels;
	  refline 0 / axis=x lineattrs=(pattern=shortdash);
	  	  yaxis values = (0,10, 20, 30, 40, 50, 60,70, 80, 90, 100);
	  xaxis values=(-2.8, -2.6, -2.4, -2.2, -2.0, -1.8, -1.6, -1.4, -1.2, -1.0, -0.8, -0.6, -0.4, -0.2, 0, .2, .4, .6, 0.8, 1.0) ;
	  keylegend 'a' / across=1 position=right valueattrs=(size=6);
	run;

data ins4;
	set ins;

/*can v tob*/

			 logp_canminustob = log10(can_minustob_p)*-1;
			 if logp_canminustob >100 then logp_canminustob=100; *recoding one logp_canminustob to 100 bc it was outlier and fig looked odd;
			/*proc print; var _name_ can_minustob_p logp_canminustob; run;*/
			if can_minustob_p <.00021 then label =_name_;

		*pvalue groups -- bonferroni sig or not, used in sgplot to show different colors for each group;
		if can_minustob_p <.00021 then pgroup_canminustob =1;
		if can_minustob_p <.05 and can_minustob_p ge .00021 then pgroup_canminustob =2; *sig at p<.05 but not at p<.0002;
		if can_minustob_p ge .05 then pgroup_canminustob = 3;
		
							

		*making new var names to refer to log p for sgplot to run;
		if pgroup_canminustob = 1 then bonp_canminustob =logp_canminustob;
		if pgroup_canminustob = 2 then sigp_canminustob = logp_canminustob;
		if pgroup_canminustob = 3 then nonsigp_canminustob = logp_canminustob;

		label bonp_canminustob ='-log10 (P-value)' sigp_canminustob='-log10 (P-value)' nonsigp_canminustob ='-log10 (P-value)';;


		*volcano plot example: https://www.lexjansen.com/phuse-us/2020/dv/DV14.pdf 
						example https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4391864/;
		*volcano plot code: https://blogs.sas.com/content/graphicallyspeaking/2016/05/23/ctspedia-clinical-graphs-volcano-plot/;
		*used volcano plot code: https://communities.sas.com/t5/Graphics-Programming/Group-data-on-a-plot-by-multiple-grouping-variables/td-p/175541;

		*checking to make sure no issue with log p values;
		proc means n range min max data=ins4; var logp_canminustob; run;

	ODS GRAPHICS/WIDTH=4IN HEIGHT=3IN;
	
	ods exclude none;

	*Figure 2c;
	proc sgplot ; 
	styleattrs datacontrastcolors=(blue gray darkred );
	  label logp_canminustob='-log10 (P-value)';
	  label can_minustob_est='Mean Difference (Long-term Cannabis Users vs Long-Term Tobacco Users)';
	  styleattrs datasymbols=(circlefilled);
	  scatter x=can_minustob_est y=logp_canminustob /group=pgroup_canminustob /*datalabel=label */; *use this code to get datalabels;
	  refline 0 / axis=x lineattrs=(pattern=shortdash);
	  	  yaxis values = (0,10, 20, 30, 40, 50, 60,70, 80, 90, 100);
	  xaxis values=(-2.8, -2.6, -2.4, -2.2, -2.0, -1.8, -1.6, -1.4, -1.2, -1.0, -0.8, -0.6, -0.4, -0.2, 0, .2, .4, .6, 0.8, 1.0) ;
	  keylegend 'a' / across=1 position=right valueattrs=(size=6);
	run;
