/* fitting beta-binomial and GEE models */

data sim;
    infile "SimData.txt" firstobs=2;
    input ID ClusterSize NResp Trt $ Dataset;
run;
    
proc sort data=sim;
    by Dataset descending Trt;
run;    
    
/* GEE */
data simlong;
    set sim;
    do  j=1 to ClusterSize;
        Y = (j <= NResp); output;
    end;
run;    
        
ods listing close;    
proc gee data=simlong descending ;
    by Dataset;
    ods output GEEEmpPEst = geeest;
    class Trt ID / order=data;
    model Y = Trt/dist=binomial link=log;
    repeated subject=id / type=cs;
run;    
ods listing;
    
proc export data=geeest(keep = Dataset Parm Level1 Estimate) file="GEEest.txt" DBMS=TAB replace; run;    
proc print data=geeest(obs=20); run;    
    
    
/* Beta-binomial */    
ods listing close;  
proc fmm data=sim order=data  tech=nrridg;
    by Dataset;
    ods output ParameterEstimates = est FitStatistics=fit;
    class Trt;
    model NResp / ClusterSize = Trt/dist=betabinomial link=log parms(-0.7 -2 -1 -0.5 02);
run;    
ods listing;
  
   
proc export data=est(keep = Dataset Effect Trt Estimate) file="BBest.txt" DBMS=TAB replace; run;       
proc print data=est(obs=20); run;    
    
