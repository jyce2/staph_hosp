libname s "~/folder";

/* Exploratory Data Analysis */
/* Histograms */
ods graphics on/reset;
ods graphics on/ height = 5in;
title1 'SAS Exploratory Data Analysis (EDA) Output';
title2 'Age Distribution by Hospitalization';
proc sgplot data=s.staph;
	histogram newage / group=hosp weight=weight transparency=0.8 scale=count;
run;

proc sgplot data=s.staph noautolegend;
	scatter x=hosp y=newage / jitter freq=weight;
	vbox newage / nofill nomean category = hosp;
run;
title;


title 'Age Distribution by Sex';
proc sgplot data=s.staph;
	histogram newage / group=sex weight=weight transparency=0.8 scale=count;
run;

proc sgplot data=s.staph noautolegend;
	scatter x=sex y=newage / jitter freq=weight;
	vbox newage / nofill nomean category = sex weight=weight;
run;
title;

title 'Age Distribution by BSI';
proc sgplot data=s.staph;
	histogram newage / group=bsi weight=weight transparency=0.8 scale=count;
run;

proc sgplot data=s.staph noautolegend;
	scatter x=bsi y=newage / jitter freq=weight;
	vbox newage / nofill nomean category = bsi weight=weight;
run;
title;

title 'Age Distribution by Dialysis';
proc sgplot data=s.staph;
	histogram newage / group=kidney weight=weight transparency=0.8 scale=count;
run;

proc sgplot data=s.staph noautolegend;
	scatter x=kidney y=newage / jitter freq=weight;
	vbox newage / nofill nomean category=kidney weight=weight;
run;

title 'Age Distribution by Wound';
proc sgplot data=s.staph;
	histogram newage / group=wound weight=weight transparency=0.8 scale=count;
run;

proc sgplot data=s.staph;
	scatter x=wound y=newage / jitter freq=weight;
	vbox newage / nofill nomean category=wound weight=weight;
run;
title;

title 'Age Distribution by Diabetes';
proc sgplot data=s.staph noautolegend;
	histogram newage / group=diabetes weight=weight transparency=0.8 scale=count;
run;

proc sgplot data=s.staph;
	scatter x=diabetes y=newage / jitter freq=weight;
	vbox newage / nofill nomean category=diabetes weight=weight;
run;
title;

title 'Age Distribution by MRSA';
proc sgplot data=s.staph;
	histogram newage / group=mrsafinal weight=weight transparency=0.8 scale=count;
run;

proc sgplot data=s.staph noautolegend;
	scatter x=mrsafinal y=newage / jitter freq=weight;
	vbox newage / nofill nomean category=mrsafinal weight=weight;
run;
title;

title 'Age Distribution by Smoker';
proc sgplot data=s.staph;
	histogram newage / group=smoker weight=weight transparency=0.8 scale=count;
run;

proc sgplot data=s.staph noautolegend;
	scatter x=smoker y=newage / jitter freq=weight;
	vbox newage / nofill nomean category=smoker weight=weight;
run;
title;


title 'Age Distribution by Setting';
proc sgplot data=s.staph;
	histogram newage / group=co weight=weight transparency=0.8 scale=count;
	label co = 'Setting';
run;

proc sgplot data=s.staph noautolegend;
	scatter x=co y=newage / jitter freq=weight;
	vbox newage / nofill nomean category=co weight=weight;
run;
title;


* sort &topic by decreasing median;
proc sql;
create table sorted as 
select *, median(newage) as med_age
	from s.staph
    group by txhosp
    order by med_age;
quit;

title 'Weighted Summary statistics: Overall Age';
proc means data=sorted n mean median stddev qrange stddev maxdec=1 nway;
	var newage / weight=weight;
run;

title 'Summary statistics: Age by Hospital ID';
proc means data=sorted n mean median stddev qrange stddev maxdec=1 nway;
	class txhosp / order=data;
	var newage / weight=weight;
run;

ods trace on;
title 'Summary statistics: Age by Hospital ID';
proc surveymeans data=s.staph ;
	domain hosp;
	var newage;
	weight weight;
run;


ods graphics on / reset;
ods graphics / height = 10in width=8in;
title 'Median Age across Hospitals';
proc sgplot data=sorted;
	scatter y=txhosp x=newage / jitter 
	filledoutlinedmarkers 
      	markerfillattrs=(color=blue) 
      	markeroutlineattrs=(color=black thickness=1)
      	markerattrs=(symbol=circlefilled size=5)
      	transparency=0.5;
	hbox newage / nomean nofill 
		displaystats=(n median) category=txhosp weight=weight 
		medianattrs=(color=blue);
	yaxis discreteorder=data label ='Hospital ID' fitpolicy=none 
		labelattrs=(size=10pt color=black) 
		valueattrs=(size=7pt color=blue) ;
	xaxis values=(0 10 20 30 40 50 60 70 80 90 105) 
		label = 'Age (years)';
	refline 45 / axis=x labelloc=inside labelpos=auto 
		label='Overall Median' labelattrs=(size=8 color=darkred) 
		lineattrs=(color=darkred thickness=1pt pattern=dash);
run;
title;
/*************************************************************************************************************************************************************************/

ods graphics / reset;
%macro desc(dom=, label=);
title "Weighted 2-sample t-test: Age by &label.";
proc ttest data=s.staph sided=2 test=diff;
	class &dom;
	var newage;
	weight weight;
	label &dom = "&label.";
title;
%mend desc; 

%desc(dom=co, label=Setting);
%desc(dom=hosp, label=Hospitalized);
%desc(dom=sex, label=Patient sex);
%desc(dom=bsi, label=Bloodstream Infection);
%desc(dom=kidney, label=Dialysis);
%desc(dom=diabetes, label=Diabetes);
%desc(dom=wound, label=Wound);
%desc(dom=mrsafinal, label=MRSA status);
%desc(dom=smoker, label=Smoker);


ods proctitle;
/* weighted chi-square test w/ unadjusted odds ratios */
title 'Weighted Chi-square test by hosp w/ unadjusted OR';
proc surveyfreq data=s.staph;
	table (sex mrsafinal kidney diabetes smoker bsi wound)*hosp / nopercent lrchisq or;
	weight weight;
run;
title;

title 'Weighted Chi-square test by setting w/ unadjusted OR';
proc surveyfreq data=s.staph;
	table (sex mrsafinal kidney diabetes smoker bsi wound)*co / nopercent lrchisq or;
	weight weight;
run;
title;


* Adjust for confounding = co indicator;
proc sort data=s.staph out=sort_staph;
	by descending hosp mrsafinal;
run;

* 3x2 table;
title 'MRSA*Hospitalization adjusted RR (of CO/HA infections)';
ods graphics on;
proc freq data=sort_staph order=data;
	table co*mrsafinal*hosp / nopercent measures cmh plots=oddsratioplot ;
	weight weight;
run;

* 2x2 table;
title 'MRSA*Hospitalization unadjusted RR';
proc freq data=sort_staph order=data;
	table mrsafinal*hosp / nopercent nocol measures oddsratio plots=all;
	weight weight;
run;


/*Variables of interest: 
HOSP 
mrsafinal kidney diabetes smoker bsi wound (binary predictors)
newage (continous predictors)
invasive weight (strata, weight)
hospid txhosp studyid (cluster, pt)*/

/* Unadjusted 1-way Analyses */

title 'Table 1 Weighted Categorical Statistics';
proc surveyfreq data = s.staph varmethod=taylor; *Table 2;
	weight weight;
	tables (sex mrsafinal kidney diabetes smoker bsi wound)*hosp / chisq or ;
run;

ods exclude SummaryPanel;
title 'Table 2 Weighted Interval Statistics';
proc surveymeans data = s.staph varmethod=taylor median mean var nmiss lclm uclm; 
	weight weight;
	class hosp;
	var newage ;
run;

ods exclude SummaryPanel;
ods select TestsForNormality Plots;
title 'Check Age distribution for normality';
proc univariate data=s.staph plots normal normaltest;
	var newage;
	*weight weight;
run;


