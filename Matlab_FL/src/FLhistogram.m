function histstruct = FLhistogram(p, NrPtPerBin)

NumSamples = length(p);
CumHistLevels=NrPtPerBin/NumSamples;
CumHistLevels=0:CumHistLevels:1;
CumHistLevels=prctile(p,100*CumHistLevels);
CumHist=arrayfun(@(x) sum(p<x),CumHistLevels);
CumHist(1)=0;
CumHist(end)=NumSamples;
histstruct.CumHist=CumHist/NumSamples;
histstruct.CumHistLevels = CumHistLevels;
histstruct.NumSamples=NumSamples;
histstruct.minmax = [min(p) max(p)];

end