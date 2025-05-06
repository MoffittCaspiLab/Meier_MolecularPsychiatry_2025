*Supplemental Table 2 on Probes -
 list of cpgs, chromosome, nearest Gene, LOcation

 Correlation between dose-response probe p value and reliability;

data one;
	set methyl.cpg_genomicnotes_aug24; 
	*this file is updated to include CpGs from Fang and Garrett. it has 272 probes because there are duplicates ;


	*creating a variable for chromosome number;



	chr=substr(Genomic_Location__hg19_,4,2);*getting just the chromosome number 1-23;
	chromosomei = compress(chr, ':'); *removing colon;
	chrom = chromosomei;
	if chrom = '' then chromosomei = 'NA';
	*if chromosomei = "NA" then chromo = '';  *Removing NAs;
	*can't make numeric due to presence of x chromosome;

	*creating a variable for start position;
	startpointi=scan(Genomic_Location__hg19_,-1,':'); *taking the first block after the colon;
	*startpointii = scan(startpointi, 1,'-'); *taking the first block before the dash;
	*if startpointii = "NA" then startpointii = '';*Removing NAs;
	*startpnt = compress(translate(startpointii,"",',')); *removing commas;

	*Shortening the reliability name of 450-epic;
	rel450 = Reliability__450K_EPIC_;
	drop Reliability__450K_EPIC_;
	
		new_probeid = compress (probeid,'cg'); *removing cg from probeid so i can sort by probe;


	proc sort; by new_probeid;

	proc print; var new_probeid manuscript beadchip; by new_probeid; run; *this shows me duplicate probes are within study because the study had more than one exposure. the exception is tobacco cpg 05575921, for which 3 studies found this probe was associated with cannabis (osbrone, nannin, garrett);

	proc contents data=one; run; *probeid is character with length 10 and format $10. new_probeid is character with length 10. this is important for merging with data on regression results;

		proc print data=one; var probeid; where manuscript="Nannini"; run;*n=201 for nannini but not unique -- bc 4 predictors;

data two;
	set one; proc sort nodupkey; by new_probeid ; *removing duplicates, leaves me with 251.  which I will need to add back in if i want to show notes section on exposures and I will need to add back that multiple studies had cg05575921;

	proc print; var probeid; where manuscript="Nannini"; run; *actually 182 when count cg05575921;
	proc print; var probeid; where manuscript="Nannini" and  mean_45_beta_value ne .; run;*this shows 179 but there are actually 180 bc 1 is 05575921 which is dup removed from nannini;
	*Table S2 - limit to the 246 on Dunedin age-45 beadchip;
	proc print; var probeid manuscript; where mean_45_beta_value ne .; run;
	proc print; var probeid beadchip Genomic_Location__hg19_ nearest_gene Location_relative_to_nearest_gen  Mean_45_Beta_value  rel450;  where mean_45_beta_value ne .; run; *Supp table 2 -- I printed this, put it in excel and then added Osborne and Nannini and garrett on cg05565921;
	
	proc means n mean std median min max; var rel450; where mean_45_beta_value ne .; run; *I only have reliability data for 156 ;

	proc sort; by new_probeid;

	proc means n mean std median min max; var rel450;  where mean_45_beta_value ne .; run;*Average reliability of probes;


*Correlation between probe p-value and reliability. 
	Step 1 = pull in dose-response results;
data three; 
	set methyl.dose_mmrobust_Sept24; *pulling in regression output from persistent reg can use, Model 1 from all 246 probes;
	

	probeid = _NAME_;
	new_probeid = compress (probeid,'cg'); *removing cg from probeid so i can sort by probe;

		proc sort; by new_probeid;
		run;

data four;
	merge two three; by new_probeid;

	*taking the absolute value of effect sizes to correlate with reliability data;
	abs_canm1 = abs(estimate_canm1);
	abs_tobm1 = abs(estimate_tobm1);


	*reported in discussion. only report for model 1, since models 2 and 3 drop nonsig CPGs (for p>.05);
	proc corr; var rel450; with canm1_probt abs_canm1  tobm1_probt abs_tobm1; run;

	proc corr; var canm1_probt abs_canm1; run;
