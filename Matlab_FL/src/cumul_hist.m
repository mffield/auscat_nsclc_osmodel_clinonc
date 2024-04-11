function [CumHist,CumHistLevels,NumSamples,minmax] = cumul_hist(var, NrPtPerBin)

NumSamples = length(var);
CumHistLevels=NrPtPerBin/NumSamples;
CumHistLevels=0:CumHistLevels:1;
CumHistLevels=prctile(var,100*CumHistLevels);
CumHist=arrayfun(@(x) sum(var<x),CumHistLevels);
CumHist(1)=0;
CumHist(end)=NumSamples;
CumHist=CumHist/NumSamples;
minmax = [min(var) max(var)];
    
end