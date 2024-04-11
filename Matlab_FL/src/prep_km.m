function [KM_b, KM_t, Nm, Nriskgroup] = prep_km(model, treat_cohort)


Nm = length(model.client(1).validation_km);
Nriskgroup = length(model.client(1).validation_km(1,1).(treat_cohort));
KM_b=struct; %bootstraps
KM_t=struct; %total
Nclient = length(model.client);
for i=1:Nclient
    for g = 1:Nriskgroup
        KM_b(i,g).D = zeros(size(model.client(1).validation_km(1).(treat_cohort)(1).D,1),Nm);
        KM_b(i,g).N = zeros(size(model.client(1).validation_km(1).(treat_cohort)(1).D,1),Nm);
        KM_b(i,g).num_pts = 0;
        for mn = 1:Nm
            KM_b(i,g).D(:,mn) = model.client(i).validation_km(mn).(treat_cohort)(g).D;
            KM_b(i,g).N(:,mn) = model.client(i).validation_km(mn).(treat_cohort)(g).N;
            KM_b(i,g).num_pts(:,mn) = model.client(i).validation_km(mn).(treat_cohort)(g).num_pts;
        end
    end
end
for g = 1:Nriskgroup
    KM_b(Nclient+1,g).D = zeros(size(model.client(1).validation_km(1).(treat_cohort)(1).D,1),Nm);
    KM_b(Nclient+1,g).N = zeros(size(model.client(1).validation_km(1).(treat_cohort)(1).D,1),Nm);
    KM_b(Nclient+1,g).num_pts = 0;
    for i=1:Nclient
        KM_b(Nclient+1,g).D = KM_b(Nclient+1,g).D + KM_b(i,g).D;
        KM_b(Nclient+1,g).N = KM_b(Nclient+1,g).N + KM_b(i,g).N;
        KM_b(Nclient+1,g).num_pts = KM_b(Nclient+1,g).num_pts + KM_b(i,g).num_pts;
    end
end

% [median, 95%CI -> [0.025, 0.975]] [2.5,50,97.5]
prcintervals = [2.5,50,97.5];%[2.5,50,97.5]
for i=1:(Nclient+1) 
    for j=1:Nriskgroup
        KM_t(i,j).D = prctile(KM_b(i,j).D',prcintervals)';
        KM_t(i,j).N = prctile(KM_b(i,j).N',prcintervals)';
        KM_t(i,j).num_pts = prctile(KM_b(i,j).num_pts',prcintervals);
    end
end


end
