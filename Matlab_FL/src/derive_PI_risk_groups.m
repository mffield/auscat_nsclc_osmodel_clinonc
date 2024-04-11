function model = derive_PI_risk_groups(model, pi_str)


model.PI.hist=[];
for j=1:length(model.client(1).validation_result)
    model.PI.rawhist = [];
    for i=1:length(model.client)
        NS = model.client(i).validation_result(j).(pi_str).NumSamples;
        CH = model.client(i).validation_result(j).(pi_str).CumHist;
        CHL = model.client(i).validation_result(j).(pi_str).CumHistLevels;
        model.PI.rawhist = [model.PI.rawhist;[[0;diff(CH*NS)],CHL]];
    end
    [chl, ord]=sort(model.PI.rawhist(:,2));
    model.PI.rawhist = model.PI.rawhist(ord,:);
    cumhist = cumsum(model.PI.rawhist(:,1))/sum(model.PI.rawhist(:,1));
    model.PI.histlevel = 0:0.005:1; %0.005
    model.PI.hist = [model.PI.hist; interp1(model.PI.rawhist(:,2),cumhist,model.PI.histlevel)];
end
% [m,i1]=min(abs(mean(model.PI.hist,'omitnan')-0.25)); [m,i2]=min(abs(mean(model.PI.hist,'omitnan')-0.75)); 
% model.settings.PI_risk_group = [0 model.PI.histlevel(i1);... 
%                                 model.PI.histlevel(i1) model.PI.histlevel(i2);...
%                                 model.PI.histlevel(i2) 1];
                            
% % [m,i1]=min(abs(mean(model.PI.hist)-0.5));
% model.settings.PI_risk_group = [0 model.PI.histlevel(i1);... 
%                                 model.PI.histlevel(i2) 1];
%                             
                                                       
[m,i1]=min(abs(median(model.PI.hist,'omitnan')-0.25)); [m,i2]=min(abs(median(model.PI.hist,'omitnan')-0.75)); 
model.settings.PI_risk_group = [0 model.PI.histlevel(i1);... 
                                model.PI.histlevel(i1) model.PI.histlevel(i2);...
                                model.PI.histlevel(i2) 1];

model.settings.PI_median = median(model.PI.hist,'omitnan');

[~,ii]=min(abs(median(model.PI.hist,'omitnan')-0.5));
model.settings.PI_median = model.PI.histlevel(ii);
  