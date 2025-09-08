*********************************************************************
*  Title:    	  Risk Factors of Hospitalization among S. aureus-infected 
*				  individuals in Fulton County, Georgia in 2017        
*                                                                    
*  Description:   Project - Step 3 (Random-effect Logistic Model 
									using relevant variables from medical literature)
*
*  Name:          Joyce C
*
*  Date:          4/30/2025                               
*------------------------------------------------------------------- 
*				       
*  Language:      SAS, VERSION 9.4  
*
*  Input:         SAS dataset
*
*  Output:        PDF 
*                                                                    
********************************************************************;
options nodate;
ods pdf file = "~/BIOS992/Project_Modeling_Step3_Output.pdf";
libname s "~/BIOS992";


/* Change response variable to indicator */
data staph;
	set s.staph;
	if strip(lowcase(hosp))='yes' then hosp_num=1;
	else hosp_num = 0;
run;

/* Stratify dataset */
data comm hop;
	set s.staph;
	if co="Community-onset" then output comm;
	else output hop;
run;

/* Check dataset based on altered variables*/
proc surveyfreq data=staph;
	table hosp*hosp_num ;
	table co;
	weight weight;
run;

proc freq data=staph;
	table hosp*hosp_num / missing list;
	table co / missing list;
	weight weight;
run;

proc freq data=comm;
	table co / missing list;
	weight weight;
run;

proc freq data=hop;
	table co / missing list;
	weight weight;
run;


/* Adjusted Analyses */
* Weighted logistic regression, unstratified ;
title1 'SAS Modeling Output';
ods output OddsRatios=or_tot;
title2 'Baseline Weighted Logistic Regression Model';
proc surveylogistic data=staph varmethod=taylor;
	strata invasive /list ;
	class wound(ref='No') mrsafinal(ref='MSSA') kidney(ref='No') diabetes(ref='No') smoker(ref='No') bsi(ref='No') / param=ref;
	model  hosp(event='Yes') = smoker mrsafinal kidney diabetes bsi wound newage/  technique=newton link=logit rsquare corrb;
	weight weight;
	output out=pred_b p=phat_b ; 
	store model_b;
run;
title;

/* Empirical logit plots - check linearity by age */
proc rank data=pred_b out=ranks groups=40;
  var newage;
  ranks bin;
run;

/* create bins */
proc means data=ranks noprint nway;
  class bin;
  var hosp_num newage;
  output out=binned_means sum(hosp_num)=hosp_num mean(newage)=meanage;
run;

/* elogit calculates the estimated probability based on smoother addition */
data binned_means;
   set binned_means;
   elogit=log((hosp_num+(sqrt(_FREQ_ )/2))/
          ( _FREQ_ - hosp_num+(sqrt(_FREQ_ )/2)));
run;

/* Plot - Check linear relationship */
title1 "Empirical Logit against Age";
proc sgplot data=binned_means;
   reg y=elogit x=meanage /
       curvelabel="Linear Relationship?"
       curvelabelloc=outside
       lineattrs=(color=ligr);
   series y=elogit x=meanage;
   label newage = 'Age (years)';
run;
title1;


/* Plot - Check linear relationship using bins*/
title1 "Empirical Logit against Binned Age";
proc sgplot data=binned_means;
   reg y=elogit x=bin /
       curvelabel="Linear Relationship?"
       curvelabelloc=outside
       lineattrs=(color=ligr);
   series y=elogit x=bin;
   label meanage = 'Age (years)';
run;
title1;

/* Test linear association between logit and age */
title 'Linear Association between logit and age';
proc corr data=binned_means;
    var elogit meanage; 
run;

title 'Linear Regression between logit and age';
proc reg data=binned_means;
    model elogit = meanage;
run;
title;

/* Test linear association between logit and bin */
title 'Linear Association between logit and age';
proc corr data=binned_means;
    var elogit bin; 
run;

title 'Linear Regression between logit and age';
proc reg data=binned_means;
    model elogit = bin;
run;
title;


***************************************************************************************************************;


/* rename variables */
proc sql;
    create table editor as
    select effect, oddsratioest, lowercl, uppercl, 
    cat(put(round(oddsratioest, .01), 5.2),' (', put(round(lowercl, .01),5.2), ', ', put(round(uppercl,.01),5.2), ')') as cl,
        case
            when lowcase(effect) like "%smoke%" then 'Smoker'
            when lowcase(effect) like "%mrsa%" then 'MRSA'
            when lowcase(effect) like "%kidney%" then 'Dialysis'
            when lowcase(effect) like "%diabetes%" then 'Diabetes'
            when lowcase(effect) like "%bsi%" then 'Septicemia'
            when lowcase(effect) like "%wound%" then 'Wound'
            when lowcase(effect) like "%age%" then 'Patient age'
            else effect
        end as neweffect
	from or_tot;
quit;

* Forest plot;
proc sgplot data=editor noautolegend;
    scatter y=neweffect x=oddsratioest / markerattrs=(size=3 symbol=CircleFilled color=blue);
    highlow y=neweffect low=lowercl high=uppercl / lineattrs=(color=gray thickness=1) type=line clipcapshape=serif highcap=serif lowcap=serif;
    xaxis type=log logbase=10 logstyle=linear label="Adjusted Odds Ratio";
    yaxis discreteorder=data grid label="Risk Factor" type=discrete;
    yaxistable cl / label="aOR 95% CI" valueattrs=(color=black size=7)
       labelattrs=(color=black size=7);
    title1 "Baseline Weighted Logistic Model: Risk Factors of Hospitalization";
    title2 "Odds Ratios with 95% Confidence Intervals";
    refline 1 / axis=x lineattrs=(thickness=1 color=darkred pattern=dash);  
run;


/* Stratified analysis*/
* Weighted logistic regression, subset to comm;
ods output OddsRatios=or_comm;
title 'Weighted Logistic Stratified: Community-onset model';
proc surveylogistic data=comm varmethod=taylor;
	strata invasive /list ;
	class mrsafinal(ref='MSSA') kidney(ref='No') diabetes(ref='No') smoker(ref='No') bsi(ref='No') wound(ref='No') / param=ref;
	model  hosp(event='Yes') = smoker mrsafinal kidney diabetes bsi wound newage /  technique=newton link=logit rsquare;
	weight weight;
	output out=pred_c p=phat_c; 
	store model_c;
run;

* Weighted logistic regression, subset to hop;
ods output OddsRatios=or_hop;
title 'Weighted Logistic Stratified: Hospital-associated model';
proc surveylogistic data=hop varmethod=taylor;
	strata invasive /list ;
	class mrsafinal(ref='MSSA') kidney(ref='No') diabetes(ref='No') smoker(ref='No') bsi(ref='No') wound(ref='No') / param=ref;
	model  hosp(event='Yes') = smoker mrsafinal kidney diabetes bsi wound newage /  technique=newton link=logit rsquare;
	weight weight;
	output out=pred_d p=phat_d; 
	store model_d;
run;

* Rename variables;
proc sql;
    create table editor_merge as
    select effect, oddsratioest, lowercl, uppercl, 'Community-onset' as ds,
    cat(put(round(oddsratioest, .01),5.2),' (', put(round(lowercl, .01),5.2), ', ', put(round(uppercl,.01),5.2), ')') as cl,
        case
            when lowcase(effect) like "%smoke%" then 'Smoker'
            when lowcase(effect) like "%mrsa%" then 'MRSA'
            when lowcase(effect) like "%kidney%" then 'Dialysis'
            when lowcase(effect) like "%diabetes%" then 'Diabetes'
            when lowcase(effect) like "%bsi%" then 'Septicemia'
            when lowcase(effect) like "%wound%" then 'Wound'
            when lowcase(effect) like "%age%" then 'Patient age'
            else effect
        end as neweffect
	from or_comm
	union all
	select effect, oddsratioest, lowercl, uppercl, 'Healthcare-associated' as ds,
    cat(put(round(oddsratioest, .01),5.2),' (', put(round(lowercl, .01),5.2), ', ', put(round(uppercl,.01),5.2), ')') as cl,
        case
            when lowcase(effect) like "%smoke%" then 'Smoker'
            when lowcase(effect) like "%mrsa%" then 'MRSA'
            when lowcase(effect) like "%kidney%" then 'Dialysis'
            when lowcase(effect) like "%diabetes%" then 'Diabetes'
            when lowcase(effect) like "%bsi%" then 'Septicemia'
            when lowcase(effect) like "%wound%" then 'Wound'
            when lowcase(effect) like "%age%" then 'Patient age'
            else effect
        end as neweffect
	from or_hop
	order by ds, oddsratioest desc;
quit;

ods graphics on / height = 5in width = 8in;
* Odds ratio plot;
proc sgplot data=editor_merge noautolegend;
    scatter y=neweffect x=oddsratioest / group=ds groupdisplay=cluster grouporder=descending markerattrs=(size=8 symbol=CircleFilled);
    highlow y=neweffect low=lowercl high=uppercl / group=ds groupdisplay=cluster grouporder=descending lineattrs=(color=gray thickness=1) type=line clipcapshape=serif highcap=serif lowcap=serif;
    xaxis type=log logstyle=logexpand label="Adjusted Odds Ratio";
    yaxis discreteorder=data grid label="Risk Factor" type=discrete fitpolicy=split;
    yaxistable cl / class=ds classdisplay=stack colorgroup=ds labeljustify=left label="aOR 95% CI" valueattrs=(color=black)
       labelattrs=(color=black);
    title1 "Stratified Weighted Logistic Model: Risk Factors of Hospitalization";
    title2 "Odds Ratios with 95% Confidence Intervals";
    refline 1 / axis=x lineattrs=(thickness=1 color=darkgray pattern=dash);  
run;


/* Modeling Correlated Binary Outcomes Through Hierarchical Logistic Regression Models */
* Random-Intercept Model = simplest correlated model of GLMM
* Response: hospitalization;
* Fixed effects: age sex smoker mrsafinal kidney diabetes bsi wound;
* Model has hospital ID as the intercept
 - addresses the conditional distribution of the probability of outcome given the random effects 
 - models the conditional mean of a binary response given hospital ID
 - subject-specific model, E[yij | gj], j= 1 to hospitalID, i=1 to patient
 - Observations within a cluster are assumed independent given the random cluster effect (conditional independence);


* Non-stratified Generalized Linear Mixed Model;
/*Weighted logistic GLMM */
ods output OddsRatios=or_glm;
title1 'Weighted Logistic Random-Intercept Model';
title2 'Conditional on Hospital Cluster';
proc glimmix data=s.staph method=quad empirical=mbn;
   *effect newage_3 = spline(newage);
   class txhosp mrsafinal(ref='MSSA') kidney(ref='No') diabetes(ref='No') smoker(ref='No') wound(ref='No') bsi(ref='No');
   model hosp(event='Yes') = smoker mrsafinal kidney diabetes wound bsi newage / dist=binary corrb chisq
   									  link=logit
   									  obsweight=weight 
   									  ddfm=betwithin
   									  solution
   									  or;
   random intercept  / subject=txhosp
   					   solution 
   					   type=vc
   					   cl;
   covtest / wald;
   ods output SolutionR=random;
   output out=glimmix_out pred(noilink)=pred pred(ilink)=phat_g pred(noblup)=pred_marg resid=resid resid(noblup)=res_marg student=student student(noblup)=student_marg;
run;
title;


/* transpose ds for */
proc transpose data=or_glm out=or_glm_wide(rename=(_NAME_=effect));
	var smoker mrsafinal kidney diabetes wound bsi newage;
	copy estimate lower upper;
run;

/* rename variables */
proc sql;
    create table or_glm_wide_name as
    select effect, estimate, lower, upper, 
    cat(put(round(estimate, .01),5.2), ' (', put(round(lower, .01),5.2), ', ', put(round(upper,.01), 5.2), ')') as cl,
        case
            when lowcase(effect) like "%smoke%" then 'Smoker'
            when lowcase(effect) like "%mrsa%" then 'MRSA'
            when lowcase(effect) like "%kidney%" then 'Dialysis'
            when lowcase(effect) like "%diabetes%" then 'Diabetes'
            when lowcase(effect) like "%bsi%" then 'Septicemia'
            when lowcase(effect) like "%wound%" then 'Wound'
            when lowcase(effect) like "%age%" then 'Patient age'
            else effect
        end as neweffect
	from or_glm_wide
	order by estimate desc;
quit;


/* OR plot random-effects logistic model */
proc sgplot data=or_glm_wide_name noautolegend;
    scatter y=neweffect x=estimate / markerattrs=(size=8 symbol=CircleFilled);
    highlow y=neweffect low=lower high=upper / lineattrs=(color=gray thickness=1) type=line clipcapshape=serif highcap=serif lowcap=serif;
    xaxis type=log logstyle=linear label="Adjusted Odds Ratio";
    yaxis discreteorder=data grid label="Risk Factor" type=discrete fitpolicy=split;
    yaxistable cl / label="aOR 95% CI" valueattrs=(color=black)
       labelattrs=(color=black);
    title1 "Weighted Logistic Random-Intercept Model: Risk Factors of Hospitalization";
    title2 "Adjusted Odds Ratio with 95% Confidence Intervals";
    refline 1 / axis=x lineattrs=(thickness=1 color=grey pattern=dash);  
    format lower 5.2 upper 5.2;
run;


***********************************************************************************************;

/* GLMM */ 
/* Subset by CO vs. HA */
ods output OddsRatios=or_glm_comm;
title 'Weighted Logistic Random-Intercept Stratified Model: Community-Onset';
proc glimmix data=comm method=quad empirical=mbn;
   class txhosp mrsafinal(ref='MSSA') kidney(ref='No') diabetes(ref='No') smoker(ref='No') wound(ref='No') bsi(ref='No');
   *effect newage_3 = spline(newage);
   model hosp(event='Yes') = smoker mrsafinal kidney diabetes wound bsi newage/ dist=binary 
   									  link=logit
   									  obsweight=weight 
   									  ddfm=betwithin
   									  solution
   									  or;
   random intercept  / subject=txhosp
   					   solution 
   					   type=vc
   					   cl;
   covtest / wald;
   output out=glimmix_out_h pred(noilink)=pred pred(ilink)=phat_h pred(noblup)=pred_marg resid=resid resid(noblup)=res_marg student=student student(noblup)=student_marg;
run;
title;



ods output OddsRatios=or_glm_hop;
title 'Weighted Logistic Random-Intercept Stratified Model: Hospital-Associated';
proc glimmix data=hop method=quad empirical=mbn;
   *effect newage_3 = spline(newage);
   class txhosp mrsafinal(ref='MSSA') kidney(ref='No') diabetes(ref='No') smoker(ref='No') wound(ref='No') bsi(ref='No');
   model hosp(event='Yes') = smoker mrsafinal kidney diabetes wound bsi newage/ dist=binary
   									  link=logit
   									  obsweight=weight 
   									  ddfm=betwithin
   									  solution
   									  or;
   random intercept  / subject=txhosp
   					   solution 
   					   type=vc
   					   cl;
   covtest / wald;
   output out=glimmix_out_i pred(noilink)=pred pred(ilink)=phat_i pred(noblup)=pred_marg resid=resid resid(noblup)=res_marg student=student student(noblup)=student_marg;
run;
title;

/* transpose 2 ds*/
proc transpose data=or_glm_comm out=or_glm_comm_w(rename=(_NAME_=effect));
	var smoker mrsafinal kidney diabetes wound bsi newage;
	copy estimate lower upper;
run;

proc transpose data=or_glm_hop out=or_glm_hop_w(rename=(_NAME_=effect));
	var smoker mrsafinal kidney diabetes wound bsi newage;
	copy estimate lower upper;
run;


* Rename variables;
proc sql;
    create table or_glm_merge as
    select effect, estimate format=5.2, lower format=5.2, upper format=5.2, 'Community-onset' as ds,
    cat(put(round(estimate, .01), 5.2),' (', put(round(lower, .01),5.2), ', ', put(round(upper,.01),5.2), ')') as cl,
        case
            when lowcase(effect) like "%smoke%" then 'Smoker'
            when lowcase(effect) like "%mrsa%" then 'MRSA'
            when lowcase(effect) like "%kidney%" then 'Dialysis'
            when lowcase(effect) like "%diabetes%" then 'Diabetes'
            when lowcase(effect) like "%bsi%" then 'Septicemia'
            when lowcase(effect) like "%wound%" then 'Wound'
            when lowcase(effect) like "%age%" then 'Patient age'
            else effect
        end as neweffect
	from or_glm_comm_w
	union all 
	select effect, estimate format=5.2, lower format=5.2, upper format=5.2, 'Healthcare-associated' as ds,
    cat(put(round(estimate, .01), 5.2),' (', put(round(lower, .01), 5.2), ', ', put(round(upper,.01),5.2), ')') as cl,
        case
            when lowcase(effect) like "%smoke%" then 'Smoker'
            when lowcase(effect) like "%mrsa%" then 'MRSA'
            when lowcase(effect) like "%kidney%" then 'Dialysis'
            when lowcase(effect) like "%diabetes%" then 'Diabetes'
            when lowcase(effect) like "%bsi%" then 'Septicemia'
            when lowcase(effect) like "%wound%" then 'Wound'
            when lowcase(effect) like "%age%" then 'Patient age'
            else effect
        end as neweffect
	from or_glm_hop_w
	order by ds, estimate desc;
quit;


* Forest plot;
ods graphics on / height = 5in width = 8in;
proc sgplot data=or_glm_merge;
    scatter y=neweffect x=estimate / group=ds groupdisplay=cluster grouporder=descending markerattrs=(size=8 symbol=CircleFilled);
    highlow y=neweffect low=lower high=upper / group=ds groupdisplay=cluster grouporder=descending lineattrs=(color=gray thickness=1) type=line clipcapshape=serif highcap=serif lowcap=serif;
    xaxis type=log logstyle=logexpand label="Adjusted Odds Ratios"; 
    yaxis discreteorder=data label="Risk Factor" type=discrete;
    yaxistable cl / class=ds classdisplay=stack colorgroup=ds labeljustify=left label="Adj. OR, 95% CI" valueattrs=(color=black size=10)
       labelattrs=(color=black size=10);
    title1 "Stratified Random-Intercept Logistic Model: Risk Factors of Hospitalization ";
    title2 "Odds Ratios with 95% Confidence Intervals";
    refline 1 / axis=x lineattrs=(thickness=1 color=darkgray pattern=dash);  
    format estimate 5.2;
run;
title;
********************************************************************************************************************************************;

/* Studentized Binned Residual analysis */
/* Binning the predicted probabilities into frequency bins */
proc rank data=glimmix_out out=binned_data groups=50;
  var pred;
  ranks bin;
run;

/* Calculating average residuals and predicted probabilities per bin */
proc means data=binned_data noprint nway;
  class bin;
  var pred student;
  output out=binned_means mean(pred)=avg_pred_prob mean(student)=avg_s_resid stderr=stderr_raw_resid;
run;

ods graphics / height = 10in;
/* Plotting the binned residuals */
proc sgplot data=binned_means;
  scatter x=avg_pred_prob y=avg_s_resid;
  refline 0 / axis=y;
  refline 2 / axis=y;
  refline -2 / axis=y;
  yaxis values=(-4 -2 0 2 4) label="Average Studentized Residual";
  title "Binned Residuals vs. Predicted Probability: Weighted Logistic GLMM (Full)";
  xaxis label="Average Predicted Probability";
run;


/* Studentized Binned Residual analysis by Strata */
/* Binning the predicted probabilities into frequency bins */
proc rank data=glimmix_out_h out=binned_data_h groups=50;
  var pred;
  ranks bin;
run;

/* Calculating average residuals and predicted probabilities per bin */
proc means data=binned_data_h noprint nway;
  class bin;
  var pred student;
  output out=binned_means_h mean(pred)=avg_pred_prob mean(student)=avg_s_resid stderr=stderr_raw_resid;
run;

ods graphics / reset;
/* Plotting the binned residuals */
proc sgplot data=binned_means_h;
  scatter x=avg_pred_prob y=avg_s_resid;
  refline 0 / axis=y;
  refline 2 / axis=y;
  refline -2 / axis=y;
  yaxis values=(-4 -2 0 2 4) label="Average Studentized Residual";
  title "Binned Residuals vs. Predicted Probability: Community-Onset Stratum";
  xaxis label="Average Predicted Probability";
run;

/* Studentized Binned Residual analysis by Strata */
/* Binning the predicted probabilities into frequency bins */
proc rank data=glimmix_out_i out=binned_data_i groups=50;
  var pred;
  ranks bin;
run;

/* Calculating average residuals and predicted probabilities per bin */
proc means data=binned_data_i noprint nway;
  class bin;
  var pred student;
  output out=binned_means_i mean(pred)=avg_pred_prob mean(student)=avg_s_resid stderr=stderr_raw_resid;
run;

ods graphics / reset;
/* Plotting the binned residuals */
proc sgplot data=binned_means_i;
  scatter x=avg_pred_prob y=avg_s_resid;
  refline 0 / axis=y;
  refline 2 / axis=y;
  refline -2 / axis=y;
  yaxis values=(-4 -2 0 2 4) label="Average Studentized Residual";
  title "Binned Residuals vs. Predicted Probability: Hospital-Associated Straum";
  xaxis label="Average Predicted Probability";
run;



***********************************************************************************;

/* Collect predicted probabilities from models into one dataset */
proc sql;
	create table pred as 
	select a.hosp, 
		   a.newage, 
		   a.weight, 
		   a.mrsafinal,
		   a.co,
		   a.bsi,
		   a.kidney,
		   a.diabetes,
		   a.wound,
		   a.smoker,
		   b.phat_b, c.phat_c, d.phat_d, g.phat_g, h.phat_h, i.phat_i,
		CASE WHEN c.phat_c IS NOT NULL THEN c.phat_c
			WHEN d.phat_d IS NOT NULL THEN d.phat_d
			ELSE .  
			END AS phat_e,
		CASE WHEN h.phat_h IS NOT NULL THEN h.phat_h
			WHEN i.phat_i IS NOT NULL THEN i.phat_i
			ELSE .
			END AS phat_j
	from s.staph as a
		left join pred_b as b
			on a.uniqueid = b.uniqueid
		left join pred_c as c
			on a.uniqueid = c.uniqueid
		left join pred_d as d
			on a.uniqueid = d.uniqueid
		left join glimmix_out as g
			on a.uniqueid = g.uniqueid
		left join glimmix_out_h as h
			on a.uniqueid = h.uniqueid
		left join glimmix_out_i as i
			on a.uniqueid = i.uniqueid
	order by a.uniqueid;
quit;

/* ROC model comparison */
ods select ROCOverlay ROCAssociation;
proc logistic data=pred;
	weight weight;
	Baseline_Weighted_Logistic_Model:model hosp(event='yes')=phat_b;
	roc 'Stratified Weighted Logistic Model' pred=phat_e;
	roc 'Random-effects Weighted Logistic Model' pred =phat_g;
	roc 'Stratified, Random-effects Weighted Logistic Model' pred = phat_j;
run; 


/* Residual analysis of random-effect logistic model*/
proc rank data=glimmix_out out=bin_pred groups=500;
  var pred;
  ranks bin;
run;

/* Calculating average residuals and predicted probabilities per bin */
proc means data=bin_pred noprint nway;
  class bin;
  var pred student;
  output out=binned_means mean(pred)=avg_pred_prob mean(student)=avg_s_resid stderr=stderr_raw_resid;
run;

ods graphics / reset;
/* Plot binned residuals by the predicted probabilities*/
proc sgplot data=binned_means;
  scatter x=avg_pred_prob y=avg_s_resid;
  refline 0 / axis=y;
  refline 2 / axis=y;
  refline -2 / axis=y;
  yaxis values=(-4 to 4 by 2) label="Average Studentized Residual";;
  title "Binned residual plot of Random-effects Weighted logistic Model (Bins=500)";
  xaxis label="Average Predicted Probability";
run;

/* Plot binned residuals on a QQplot to check normality */
proc univariate data=binned_means normaltest;
    var avg_s_resid;
    qqplot avg_s_resid / normal(mu=0 sigma= 0.65) square;
    histogram avg_s_resid / normal;
    id bin;
run;


/* Locate influential points, or studentized residuals +/-2 above mean*/
title 'Influential points';
proc sql;
	select g.uniqueid,
		   g.hospid label="Hospital ID",
		   g.factype label="Facility type",
		   g.txhosp label="Treated Hospital",
		   g.hosp label="Hospitalized", 
		   g.newage label="Age (years)", 
		   g.invasive label="Invasive status",
		   g.mrsafinal label="Strain",
		   g.co label="Setting",
		   g.bsi label="Bloodstream infection",
		   g.kidney label ="Dialysis status",
		   g.diabetes label= "Diabetes status",
		   g.wound label="Wound status",
		   g.smoker label= "Smoker status",
		   b.*
	from bin_pred as g
	left join binned_means as b
	on g.bin = b.bin
	where (b.avg_s_resid > 2) or (b.avg_s_resid < -2);
run;

***********************************************************************************;
/* Prepare log-odds plot by hospital */
data odds;
	set random;
	odds = exp(estimate);
	l_odds = exp(lower);
	u_odds = exp(upper);
	odds_95 = cat("(", put(round(l_odds, .01),5.2), ", ", put(round(u_odds, .01),5.2),")");
	txhosp_name = scan(subject, 2);
run;

proc sql;
	CREATE TABLE odds2 AS
		SELECT * 
		FROM odds
		LEFT JOIN (SELECT txhosp, count(txhosp) AS count_txhosp
					FROM s.staph
					GROUP BY txhosp) AS sub
		ON odds.txhosp_name = sub.txhosp;	
quit;

proc sort data=odds2 out=sort_odds;
	by descending odds;
run;

/* Odds from Odds(Intercept) Random-effects by Hospital Cluster*/
ods graphics / height=8in;
proc sgplot data=sort_odds noautolegend;
    scatter y=txhosp_name x=odds / markerattrs=(size=3 symbol=CircleFilled color=blue);
    highlow y=txhosp_name low=l_odds high=u_odds / lineattrs=(color=gray thickness=1) type=line clipcapshape=serif highcap=serif lowcap=serif;
     xaxis type=log logstyle=linear label="Odds of Hospitalization";
    yaxis discreteorder=data grid label="Hospital ID" type=discrete fitpolicy=none;
    title "Odds of Hospitalization from Hospital Cluster";
    yaxistable odds/ label="Odds"  valueattrs=(color=black)
     labelattrs=(color=black);
    yaxistable odds_95/ label="95% CI"  valueattrs=(color=black)
     labelattrs=(color=black);
    yaxistable count_txhosp / label="n" valueattrs=(color=black)
     labelattrs=(color=black);
     format odds 5.2;
    refline 1 / axis=x lineattrs=(thickness=1 color=darkred pattern=dash);  
run;

ods pdf close;
