
libname s "~/folder";

ods graphics on/ reset;

/* Examine dataset*/
/*ods exclude EngineHost;
proc contents data=s.staph ;
run;

/* Stratify dataset */
data comm hop;
	set s.staph;
	if co="Community-onset" then output comm;
	else output hop;
run;


/* Adjusted Analyses */
* Weighted logistic regression, unstratified ;
title1 'SAS Modeling Output';
ods output OddsRatios=or_tot;
title2 'Baseline Weighted Logistic Regression Model';
proc surveylogistic data=s.staph varmethod=taylor;
	strata invasive /list ;
	class mrsafinal(ref='MSSA') kidney(ref='No') diabetes(ref='No') smoker(ref='No') bsi(ref='No') wound(ref='No') / param=ref;
	model  hosp(event='yes') = smoker mrsafinal kidney diabetes bsi wound newage /  technique=newton link=logit rsquare corrb;
	weight weight;
	output out=pred_b p=phat_b; 
	store model_b;
run;
title;


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
	model  hosp(event='yes') = smoker mrsafinal kidney diabetes bsi wound newage /  technique=newton link=logit rsquare;
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
	model  hosp(event='yes') = smoker mrsafinal kidney diabetes bsi wound newage /  technique=newton link=logit rsquare;
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

* Odds ratio plot;
proc sgplot data=editor_merge noautolegend;
    scatter y=neweffect x=oddsratioest / group=ds groupdisplay=cluster grouporder=descending markerattrs=(size=8 symbol=CircleFilled);
    highlow y=neweffect low=lowercl high=uppercl / group=ds groupdisplay=cluster grouporder=descending lineattrs=(color=gray thickness=1) type=line clipcapshape=serif highcap=serif lowcap=serif;
    xaxis type=log logstyle=linear label="Adjusted Odds Ratio";
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
proc glimmix data=s.staph method=quad empirical=mbn plots=residualpanel (conditional marginal);
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
   output out=glimmix_out pred(noilink)=pred pred(ilink)=phat_g;
run;
title;

/* transpose ds*/
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
proc glimmix data=comm method=quad empirical=mbn plots=residualpanel (conditional marginal);
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
    output out=glimmix_out_h pred(noilink)=pred pred(ilink)=phat_h;
run;
title;

ods output OddsRatios=or_glm_hop;
title 'Weighted Logistic Random-Intercept Stratified Model: Hospital-Associated';
proc glimmix data=hop method=quad empirical=mbn plots=residualpanel (conditional marginal);
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
    output out=glimmix_out_i pred(noilink)=pred pred(ilink)=phat_i;
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
ods graphics on/ width = 8in;
proc sgplot data=or_glm_merge;
    scatter y=neweffect x=estimate / group=ds groupdisplay=cluster grouporder=descending markerattrs=(size=8 symbol=CircleFilled);
    highlow y=neweffect low=lower high=upper / group=ds groupdisplay=cluster grouporder=descending lineattrs=(color=gray thickness=1) type=line clipcapshape=serif highcap=serif lowcap=serif;
    xaxis type=log logstyle=linear label="Adjusted Odds Ratios"; 
    yaxis discreteorder=data label="Risk Factor" type=discrete;
    yaxistable cl / class=ds classdisplay=stack colorgroup=ds labeljustify=left label="Adj. OR, 95% CI" valueattrs=(color=black size=10)
       labelattrs=(color=black size=10);
    title1 "Stratified Random-Intercept Logistic Model: Risk Factors of Hospitalization ";
    title2 "Odds Ratios with 95% Confidence Intervals";
    refline 1 / axis=x lineattrs=(thickness=1 color=darkgray pattern=dash);  
    format estimate 5.2;
run;
title;

***********************************************************************************;

/* Collect predicted probabilities from models into one dataset */
proc sql;
	create table pred as 
	select a.hosp, a.weight, b.phat_b, c.phat_c, d.phat_d, g.phat_g, h.phat_h, i.phat_i,
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


/* extra 
proc sgplot data=glimmix_out ;
    hbox p / category=txhosp;
    yaxis label="Predicted Probability of Hospitalization";
    title "Predicted Risk by Hospital";
run;
*/

/* Training and test dataset */ 

/* Scoring model on a new dataset */
/*proc plm restore model_b;
	score data=.  out=scored_. predicted/ ilink;
run; 
