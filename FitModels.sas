/* fitting beta-binomial and GEE models */

data sim;
    infile "SimData.txt" firstobs=2;
    input ID ClusterSize NResp Trt $;
run;
    
proc sort data=sim;
    by descending Trt;
run;    
    
ods listing close;  
proc fmm data=sim order=data;
    ods output ParameterEstimates = est FitStatistics=fit;
    class Trt;
    model NResp / ClusterSize = Trt/dist=betabinomial link=log;
run;    
ods listing;
    
    
proc print data=fit; run;
    
proc print data=est; run;    
    

data simlong;
    set sim;
    do  j=1 to ClusterSize;
        Y = (j <= NResp); output;
    end;
run;    
        
ods listing close;    
proc gee data=simlong descending;
    ods output GEEEmpPEst = geeest;
    class Trt ID / order=data;
    model Y = Trt/dist=binomial link=log;
    repeated subject=id / type=cs;
run;    
ods listing;
    
proc print data=geeest; run;    
    
