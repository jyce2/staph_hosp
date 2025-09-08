*********************************************************************
*  Title:    	  Risk Factors of Hospitalization among S. aureus-infected 
*				  individuals in Fulton County, Georgia in 2017        
*                                                                    
*  Description:   Project - Optional Step 4 Predictive Modeling - 
*									1) Impute missing values 
*									2) Variable Screening
*									3) Survey-weighted Logistic LASSO Regression 
*									using training set.
*									4) Evaluate lowest average squared error (ASE)
*									from validation set.
*									5) Compare LASSO model to logistic regression 
*      							 	using AUC-ROC.
*
*  Name:          Joyce C
*
*  Date:          9/2025                          
*------------------------------------------------------------------- 
*				       
*  Language:      SAS, VERSION 9.4  
*
*  Input:         Original CSV dataset
*
*  Output:        PDF 
*                                                                    
********************************************************************;
%let homefolder = BIOS992;
options nodate;
ods pdf file = "~/&homefolder/Project_Modeling_Step4_Output.pdf";

/* Import original CSV file */
proc import datafile="~/&homefolder/original_staph.csv"
    out=original_staph
    dbms=CSV
    replace;
    getnames=yes;
run;

/* ===================================================== */
/* Imputing Missing Values                 			     */
/* ===================================================== */

/* Identify variables with missing values */
proc means data=original_staph nmiss nway noprint;
	var _numeric_;
	output out=count_missing  nmiss= / autoname;
run;

/* Identify rows with missing values*/
data original_staph_mi(drop=i);
   set original_staph;
   /* create missing indicator variables */
   array mi{*} mi_newethnic mi_READMIT 
   			   mi_BOIL mi_CIRR
  			   mi_CPD11 mi_RENAL
   			   mi_CSBREAK9 mi_HEART 
   			   mi_CTD11 mi_SMOKER
   			   mi_CVA mi_CYSTIC9
   			   mi_DULCER7 mi_DEMENT9
   			   mi_DIABETES mi_HEMAP9
   			   mi_HIV mi_FLU
   			   mi_IVDU mi_TUMOR9
   			   mi_MI11 mi_OBESITY
   			   mi_DRUG7 mi_PEPTIC9
   			   mi_PVD mi_SSABC mi_SSAIW mi_SSBI mi_SSCEL mi_SSCHR mi_SSHER
			   mi_SSINF mi_SSMAS mi_SSMYO mi_SSNF mi_SSPUS mi_CLINDA mi_OX mi_SXT mi_VANCO;
   /* select variables with missing values */
   array x{*}  newethnic READMIT 
   			   BOIL CIRR
  			   CPD11 RENAL
   			   CSBREAK9 HEART 
   			   CTD11 SMOKER
   			   CVA CYSTIC9
   			   DULCER7 DEMENT9
   			   DIABETES HEMAP9
   			   HIV FLU
   			   IVDU TUMOR9
   			   MI11 OBESITY
   			   DRUG7 PEPTIC9
   			   PVD SSABC SSAIW SSBI SSCEL SSCHR SSHER
			   SSINF SSMAS SSMYO SSNF SSPUS CLINDA OX SXT VANCO;
   do i=1 to dim(mi);
      mi{i}=(x{i}=.);
      nummiss+mi{i};
   end;
run;

/* Impute missing values with median */
proc stdize data=original_staph_mi method=median
            reponly out=original_staph_imputed;
   			var newethnic READMIT 
   			   BOIL CIRR
  			   CPD11 RENAL
   			   CSBREAK9 HEART 
   			   CTD11 SMOKER
   			   CVA CYSTIC9
   			   DULCER7 DEMENT9
   			   DIABETES HEMAP9
   			   HIV FLU
   			   IVDU TUMOR9
   			   MI11 OBESITY
   			   DRUG7 PEPTIC9
   			   PVD SSABC SSAIW SSBI SSCEL SSCHR SSHER
			   SSINF SSMAS SSMYO SSNF SSPUS CLINDA OX SXT VANCO ;
run;
options nolabel;

/* Check imputed values */
title 'Verify imputed values, check non-missing';
proc means data=original_staph_imputed median nmiss maxdec=1;
   var newethnic READMIT 
   			   BOIL CIRR
  			   CPD11 RENAL
   			   CSBREAK9 HEART 
   			   CTD11 SMOKER
   			   CVA CYSTIC9
   			   DULCER7 DEMENT9
   			   DIABETES HEMAP9
   			   HIV FLU
   			   IVDU TUMOR9
   			   MI11 OBESITY
   			   DRUG7 PEPTIC9
   			   PVD SSABC SSAIW SSBI SSCEL SSCHR SSHER
			   SSINF SSMAS SSMYO SSNF SSPUS CLINDA OX SXT VANCO;
run;
options label;
title;

/* ===================================================== */
/* Reducing Redundancy by Clustering Variables */
/* ===================================================== */

/* Identify variables by name and type */
proc contents data=original_staph_imputed(drop=nummiss) out=varinfo(keep=name type) noprint;
run;

/* Save numeric variables as macro list */
proc sql noprint;
	select name into :var_list separated by ' '
	from varinfo
	where type = 1 AND name NOT LIKE "mi_%"; /* Numeric */
quit;

/* Check numeric variables */
%put &var_list;

/* Prepare variable clusters based on eigenvalue */
ods select none;
ods output clusterquality=work.summary
           rsquare=work.clusters;
proc varclus data=original_staph
             hi maxeigen=0.70;
   			 var &var_list;
run;

ods select all;

data _null_;
   set work.summary;
   call symput('nvar',compress(NumberOfClusters));
run;

/* Use cluster correlation metric to remove redundant variables*/
title1 "Variables by Cluster";
proc print data=work.clusters noobs label split='*';
   where NumberOfClusters=&nvar;
   var Cluster Variable RSquareRatio;
   label RSquareRatio="1 - RSquare*Ratio"; /*High correlation with own cluster / Low correlation with other clusters -> Lower Ratio preferred to select representative variable of a cluster*/
run;

title1 "Variation Explained by Clusters";
proc print data=work.summary label;
run;
title1 ;


/* ========================================================== */
/* Variable Screening using Rank Correlation                  */
/* Check monotonic associations with target response variable */
/* ========================================================== */

/* Save remaining variables into a variable list */
%let ex_reduced= cacase hosp_onset newwhite bsi LTCYR causal11 
				 iab LTACYR11 newage renal smoker othsite priorinvasive
				 bji MI11 UTIT HEMAP9 ICU16 newethnic PLEURAL CSF SEX
				 BOIL mrsafinal CSBREAK9 SST CTD11 ethnicity_real
				 OTHERrace HOMELESS9 INCERC9 OBESITY HIV OX CNS TUMOR9
				 BONE UND CIRR BODYSITE CPD11 SXT SURGYR9
				 IVDU PVD CVI DISLTACH10 PNE PEPTIC9 HEMMALIG CVA
				 CVC9 DIABETES CLINDA DEMENT9 HEART DISLTC8 
				 VANCO UNKRACE OTHPOS DULCER7;
	


/* Compute Spearman (monotonic) & Hoeffding (nonlinear) correlations between
predictor variables and hospitazliation */
ods select none;
ods output spearmancorr=work.spearman
           hoeffdingcorr=work.hoeffding;
proc corr data=original_staph_imputed
          spearman hoeffding;
  		  var HOSPITAL;
 		  with &ex_reduced;
run;

ods select all;

proc sort data=work.spearman;
    by variable;
run;

proc sort data=work.hoeffding;
    by variable;
run;

data work.correlations;
   attrib variable length=$32;
   merge work.spearman(rename=
         (hospital=scorr phospital=spvalue))
         work.hoeffding
         (rename=(hospital=hcorr phospital=hpvalue));
   by variable;
   scorr_abs=abs(scorr);
   hcorr_abs=abs(hcorr);
run;

/* Order correlations by rank */
proc rank data=correlations 
          out=correlations1 descending;
    var scorr_abs hcorr_abs;
    ranks ranksp rankho;
run;

proc sort data=correlations1;
   by ranksp;
run;

title1 "Rank of Spearman Correlations and Hoeffding Correlations";
proc print data=correlations1 label split='*';
   var variable ranksp rankho scorr spvalue hcorr hpvalue;
   label ranksp='Spearman rank*of variables'
         scorr='Spearman Correlation'
         spvalue='Spearman p-value'
         rankho='Hoeffding rank*of variables'
         hcorr='Hoeffding Correlation'
         hpvalue='Hoeffding p-value';
run;

%global vref href;
proc sql noprint;
   select min(ranksp) into :vref 
   from (select ranksp 
         from work.correlations1 
         having spvalue > .5);
   select min(rankho) into :href 
   from (select rankho
         from work.correlations1
         having hpvalue > .5);
quit;

/* Screen variables that are low rank in both Spearman and Hoeffding */
title1 "Scatter Plot of the Ranks of Spearman vs. Hoeffding";
proc sgplot data=correlations1;
   refline &vref / axis=y;
   refline &href / axis=x;
   scatter y=ranksp x=rankho / datalabel=variable;
   yaxis label="Rank of Spearman";
   xaxis label="Rank of Hoeffding";
run;


/* Save screened variables based on low s_rank and h_rank*/
proc sql number;
	select variable into :screened separated by ' '
	from correlations1
	where ranksp <= 10 and rankho < 14;
quit;

%global screened;
%put &screened;

/* ===================================================== */
/* LASSO Logistic regression model 					    */
/* Training & Validation datasets			             */
/* ===================================================== */

/* Sort data by strata variables */
proc sort data=original_staph_imputed out=sort_staph;
   by invasive cacase;
run;

/* Stratified sampling w/o replacement by strata (balanced);
split data into two: training set 70% + validation set 30% */
proc surveyselect noprint data=sort_staph
                  samprate=0.7 out=sample seed=512
                  outall stratumseed=restore;
   				  strata invasive cacase;
run;



/* Fit logistic regression model w/ lasso regularization on training set */
proc sort data=sample;
	by cacase;
run;
ods output LassoSelectionDetails=summary;
title "Model Selection: LASSO logistic regression";
proc hpgenselect data=sample technique=nrridg;	  		 
	by cacase;
	partition role=selected(train='1' validate='0');	 /*partition dataset 70% training + 30% validation*/
    class und(ref='0') 
		  disltc8(ref='0') 
		  diabetes(ref='0')
		  renal(ref='0')
		  sst(ref='0')
		  mrsafinal(ref='0')
		  bsi(ref='0') / param=glm;						  /* class parameterization = glm (Reference-cell coding, less than full rank) */ 
    model hospital(event='1') = und disltc8 diabetes renal sst mrsafinal newage bsi/  
    	  distribution=binary
    	  link=logit
    	  cl;
    selection method=lasso(choose=validate) details=all; /* choose model with lowest average standard error (ASE) in validation model */
    performance details;
    weight weight;										   /* apply survey-weights */
    id uniqueid cacase hospital weight und disltc8 diabetes renal sst mrsafinal bsi newage ;
    output out=pred p=phat role=role  ; 				 /* output predicted values, separated by role: train=1, validate=2 */
run;

/* Rename variables for series plot */
proc sql;
    create table summary2 as
    select cacase, step, descr, lambda, ase, validase, 
    	case 
    		when lambda = max(lambda) then round(max(validase), 0.001)
    		when lambda = min(lambda) then round(min(validase), 0.001)
    		else .
    	end as vase_min_max,
        case
            when lowcase(descr) like "%und%" then 'Undetermined site of Infection'
            when lowcase(descr) like "%disltc8%" then 'Long-term facility'
            when lowcase(descr) like "%diabetes%" then 'Diabetes'
            when lowcase(descr) like "%renal%" then 'Renal disease'
            when lowcase(descr) like "%sst%" then 'Skin-Soft Tissue'
            when lowcase(descr) like "%mrsafinal%" then 'MRSA Status'
            when lowcase(descr) like "%newage%" then 'Patient Age'
            when lowcase(descr) like "%bsi%" then 'Septicemia'
            else descr
        end as newdescr
	from summary
quit;


/* Split summary output by community-onset vs. healthcare exposure */
data co(where=(cacase=1)) ha(where=(cacase=0));
	set summary2;
    if ^missing(descr) or step >= 18 then do;
        ase_range = round(ase, 0.001);
        vase_range = round(validase, 0.001);
        output;
    end;
    else do;
        ase_range = .;
        vase_range = .;
        lambda = . ;
        output;
    end;
run;


/* Plot training vs. validation ASE */
ods graphics on / width = 8in;
title 'Predictors of Community-onset Infection: Training vs. Validation ASE from LASSO model';
proc sgplot data=co;
	series x=lambda y=ase / markers datalabel=newdescr legendlabel="Training Average Squared Error (ASE)";
	series x=lambda y=validase /markers datalabel=vase_range legendlabel="Validation Average Squared Error (ASE)";
	xaxis type=log reverse;
run;

title 'Predictors of Prior Healthcare Infection: Training vs. Validation ASE from LASSO model';
proc sgplot data=ha;
	series x=lambda y=ase / markers datalabel=newdescr legendlabel="Training Average Squared Error (ASE)";
	series x=lambda y=validase /markers datalabel=vase_range legendlabel="Validation Average Squared Error (ASE)";
	xaxis type=log reverse;
run;
title;

ods graphics on/ reset;
/* ROC comparison */
/* Baseline logistic vs. LASSO logistic models */
ods select ROCOverlay ROCAssociation;
proc logistic data=pred;
	weight weight;
	Baseline_Logistic_Model:model hospital(event='1')= und disltc8 diabetes renal sst mrsafinal newage bsi;
	roc 'LASSO Logistic Model' pred=phat ;
run; 

/*  Split the data into training and validation sets */
data training(drop=selected SelectionProb SamplingWeight)
     validate(drop=selected SelectionProb SamplingWeight);
   set sample;
   if selected then output training;
   else output validate;
run;

/* Stratify by healthcare exposure vs. community-associated to generate separate models */
proc sort data=training; by cacase; run;
proc sort data=validate; by cacase; run;

/* ROC comparison */
/* Standard logistic: Training vs. Validation Set */
title 'Baseline Logistic Model: Training vs. Validation Set';
ods output OddsRatios=or_estimates;
proc logistic data=training;
	weight weight;
	class und(ref='0') 
		  disltc8(ref='0') 
		  diabetes(ref='0')
		  renal(ref='0')
		  sst(ref='0')
		  mrsafinal(ref='0')
		  bsi(ref='0') / param=glm;	
	Training: model hospital(event='1')= und disltc8 diabetes renal sst mrsafinal newage bsi / technique=newton outroc=troc;
	score data=validate out=valpred outroc=vroc;
run; 

/* Rename variables for forest plot */
proc sql;
    create table editor as
    select effect, oddsratioest, lowercl, uppercl, 
    cat(put(round(oddsratioest, .01), 5.2),' (', put(round(lowercl, .01),5.2), ', ', put(round(uppercl,.01),5.2), ')') as cl,
        case
            when lowcase(effect) like "%und%" then 'Undetermined site of Infection'
            when lowcase(effect) like "%disltc8%" then 'Long-term facility'
            when lowcase(effect) like "%diabetes%" then 'Diabetes'
            when lowcase(effect) like "%renal%" then 'Renal disease'
            when lowcase(effect) like "%sst%" then 'Skin-Soft Tissue'
            when lowcase(effect) like "%mrsafinal%" then 'MRSA Status'
            when lowcase(effect) like "%newage%" then 'Patient Age'
            when lowcase(effect) like "%bsi%" then 'Septicemia'
            else effect
        end as neweffect
	from or_estimates;
quit;

proc sort data=editor out=editor;
	by oddsratioest;
run;

* Forest plot of Standard Logistic Model w/ predictors in order;
proc sgplot data=editor noautolegend;
    scatter y=neweffect x=oddsratioest / markerattrs=(size=3 symbol=CircleFilled color=blue);
    highlow y=neweffect low=lowercl high=uppercl / lineattrs=(color=gray thickness=1) type=line clipcapshape=serif highcap=serif lowcap=serif;
    xaxis type=log logbase=10 logstyle=linear values=(0.1 0.7 1 5 10 15 25) label="Adjusted Odds Ratio";
    yaxis discreteorder=data grid label="Risk Factor" type=discrete;
    yaxistable cl / label="aOR 95% CI" valueattrs=(color=black size=7)
       labelattrs=(color=black size=7);
    title1 "Weighted Logistic Model: Risk Factors of Hospitalization from S. aureus Infection";
    title2 "Odds Ratios with 95% Confidence Intervals";
    refline 1 / axis=x lineattrs=(thickness=1 color=darkgrey pattern=dash);  
run;


ods pdf close;