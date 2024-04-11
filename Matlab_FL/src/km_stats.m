function [KM] = km_stats(data,g,num_group,model)

    KM=struct;
    for k=1:num_group

        risk_group_surv = data.event_time(ismember(g,k))/12;
        risk_group_vital = data.censor(ismember(g,k));
        if ~isempty(risk_group_surv)
            [temp_km_D, temp_km_N] = emp_km_curve(risk_group_surv,risk_group_vital,model.settings.km_periods,model.settings.km_limit);
        else
            temp_km_D = zeros(model.settings.km_periods*model.settings.km_limit,1);
            temp_km_N = zeros(model.settings.km_periods*model.settings.km_limit,1);
        end
        KM(k).D = temp_km_D;
        KM(k).N = temp_km_N;
        KM(k).num_pts = length(risk_group_surv);

    end
    
end