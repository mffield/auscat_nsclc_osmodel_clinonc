function [central_stats]=central_stats_collection(model,stat_string)

stats_collect = {'N','pts_curative','pts_palliative','pts_ruled_out','pts_select_cur_no_imaging','pts_select_pall_no_imaging','lost_followup',...
    'pts_select_cur','pts_select_pall','train_outcomes_curative','valid_outcomes_curative','train_outcomes_palliative','valid_outcomes_palliative'};
for i = 1: length(stats_collect)
    if isfield(model.client(1).(stat_string),stats_collect{i})
        stats.(stats_collect{i}) = [];
        for j=1:length(model.client)
            stats.(stats_collect{i}) = [stats.(stats_collect{i}) model.client(j).(stat_string).(stats_collect{i})];
        end
    end
end

if isfield(model.client(1).(stat_string),'localN')
    stats.sum = []; stats.localN = [];
    for j=1:length(model.client)
        stats.sum = [stats.sum model.client(j).(stat_string).localSum];
        stats.localN = [stats.localN model.client(j).(stat_string).localN];
    end
    stats.Mu = sum(stats.sum,2)./sum(stats.localN,2);
    
    if isfield(model.client(1).(stat_string),'localStdSum')
        stats.stdsum = [];
        for j=1:length(model.client)
            stats.stdsum = [stats.stdsum model.client(j).(stat_string).localStdSum];
        end        
        stats.Std = sqrt(sum(stats.stdsum,2)./((sum(stats.localN,2))-1));
    end

end


vars = {'fev','dlco','smoking_pack_yr','age_start_rt','tumour_volume','presc_dose_ttl','eqd2','survival','time_to_treat','time_to_cardiac','event_time'};%
for i=1:length(vars)
    if isfield(model.client(1).(stat_string),vars{i})
        stats.(vars{i}).sum = 0; stats.(vars{i}).N = 0; stats.(vars{i}).Nvec = [];
        for j=1:length(model.client)
            stats.(vars{i}).sum = stats.(vars{i}).sum + model.client(j).(stat_string).(vars{i}).sum;
            stats.(vars{i}).N = stats.(vars{i}).N + model.client(j).(stat_string).(vars{i}).N;
            stats.(vars{i}).Nvec = [stats.(vars{i}).Nvec model.client(j).(stat_string).(vars{i}).N];
        end
        stats.(vars{i}).mean = stats.(vars{i}).sum/stats.(vars{i}).N;
        if isfield(model.client(1).(stat_string).(vars{i}),'stdsum')
            stats.(vars{i}).stdsum = 0;
            for j=1:length(model.client)
                stats.(vars{i}).stdsum = stats.(vars{i}).stdsum + model.client(j).(stat_string).(vars{i}).stdsum;
            end
            stats.(vars{i}).std = real(sqrt(stats.(vars{i}).stdsum/stats.(vars{i}).N-1));
        end
    end
    
end

vars = {'overall_stage','histology','ecog','gender','tumour_loc','tumour_grade','smoking','weightloss',...
    'marital_status','alcoholism','laterality','technique',...
    'surgery','prev_body_rt','cardiac_comorbidity','t_stage','n_stage','rt_intent'};
for i=1:length(vars)
    if isfield(model.client(1).(stat_string),vars{i})
        stats.(vars{i}).num = zeros(length(model.client(1).(stat_string).(vars{i}).num),length(model.client));
        stats.(vars{i}).category = model.client(1).(stat_string).(vars{i}).category;
        for j=1:length(model.client)
            stats.(vars{i}).num(:,j) = model.client(j).(stat_string).(vars{i}).num;
        end
    end
end
central_stats = stats;
