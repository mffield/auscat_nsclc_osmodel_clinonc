function [Data, columnHeaders]=loadLungData(config, model)

Data=[]; columnHeaders=[];

if strcmp(config.dataType,'csv')
    
    % File is expected to contain:
    % Comma separated variables - at least columns survival time and vital status  
    Data = readtable(config.csv_location);
    columnHeaders = Data.Properties.VariableNames;
    date_vars = {'death_date','diagnosis_date','first_RT_date','last_RT_date','censor_date','death_date_clinic','death_date_rbdm',...
        'censor_date_clinic','censor_date_rbdm','cardiac_event_date'};
    for j=1:length(date_vars)
        if ismember(date_vars{j},Data.Properties.VariableNames)
            if ~isa(Data.(date_vars{j})(1),'datetime')
                tmp1 = cellfun(@(x) strsplit(x,'T'),Data.(date_vars{j}),'UniformOutput',false);
                Data.(date_vars{j}) = datetime(cellfun(@(x) x{:,1},tmp1,'UniformOutput',false));
            end
        end
    end
    
elseif strcmp(config.dataType,'sparql')
    
    %Config expected to contain file for sparql query and endpoint detail
    %or alternatively the query is sent through a message.
    if isfield(model, 'sparql_query')
        sparql_query = model.sparql_query;
    end
    [columnHeaders, Data, extra] = sparql_execute(config.endpoint_location, sparql_query, 1);

end


missing = false(size(Data));
for i=1:size(Data,2)
    if isa(Data.(Data.Properties.VariableNames{i}),'cell')
        if ~contains(Data.Properties.VariableNames{i},'vital')
            if contains(Data.Properties.VariableNames{i},'date')
                isadatetime = (cell2mat(cellfun(@(x) isa(x,'datetime'),Data.(Data.Properties.VariableNames{i}),'UniformOutput',false)));
                if any(isadatetime)
                    for j=1:length(isadatetime)
                        if ~isadatetime(j)
                            Data.(Data.Properties.VariableNames{i})(j)={NaT};
                        end
                    end
                    Data.(Data.Properties.VariableNames{i}) = [Data.(Data.Properties.VariableNames{i}){:}]';
                end   
%                 Data.(Data.Properties.VariableNames{i}) = datetime(Data.(Data.Properties.VariableNames{i}),'InputFormat','yyyy-MM-dd''T''HH:mm:SS');
            else
                missing(:,i) = cell2mat(cellfun(@(x) any(isnan(x)),Data.(Data.Properties.VariableNames{i}),'UniformOutput',false));
                Data.(Data.Properties.VariableNames{i})(missing(:,i))=repmat({''},sum(missing(:,i)),1);
                % impute with mode is missing here otherwise we have hard constraint of no missing data at all
            end
        end
    elseif isa(Data.(Data.Properties.VariableNames{i}),'double')
        if ~contains(Data.Properties.VariableNames{i},'survival')
            missing(:,i) = isnan(Data.(Data.Properties.VariableNames{i}));
            % impute with mean if missing here otherwise we have hard constraint of no missing data at all
        end
    end
end
% remove_data = any(missing,2);
% remove_data = isnan(Data.tumour_volume);
% Data(remove_data,:)=[];

Data.overall_stage = cell(size(Data,1),1); Data.overall_stage_num = zeros(size(Data,1),1);
for i=1:size(Data,1)

    [stage_string, num]=overallstage(Data.t_stage{i}, Data.n_stage{i}, Data.m_stage{i}, Data.stage_group{i});
    Data.overall_stage{i} = stage_string;
    Data.overall_stage_num(i) = num;
    if ismember(stage_string,{'Stage I','Stage IA','Stage IB'})
        Data.overall_stage_num(i) = 1;
    elseif ismember(stage_string,{'Stage II','Stage IIA','Stage IIB'})
        Data.overall_stage_num(i) = 2;
    elseif ismember(stage_string,{'Stage III','Stage IIIA','Stage IIIB'})
        Data.overall_stage_num(i) = 3;
    elseif ismember(stage_string,{'Stage IV','Stage IVA','Stage IVB'})
        Data.overall_stage_num(i) = 4;
    end
end

%%% prepare ecog status
datanotisnan = cell2mat(cellfun(@(x) ~any(isnan(x)),Data.ecog,'UniformOutput',false));
Data.ecog(datanotisnan) = cellfun(@(x) str2double(x), regexp(Data.ecog(datanotisnan),'\d+(\.)?(\d+)?','match') ,'UniformOutput',false);
Data.ecog(cell2mat(cellfun(@(x) isempty(x), Data.ecog,'UniformOutput',false)))={NaN};
Data.ecog = cell2mat(Data.ecog);

% prepare histology
Data.histology(strcmp(Data.histology,'')) = {'?'};

for i=1:size(Data,1)
    if (isnat(Data.death_date(i)) & (year(Data.censor_date(i))<2018))
        Data.censor_date(i) = datetime(2017,12,31);
    end
end
% variables: censor, event_time
if ~isfield(model.settings,'death_date')
    model.settings.death_date = 'death_date';
    model.settings.censor_date = 'censor_date';
end

Data.vital_status = isnat(Data.(model.settings.death_date));
for i=1:size(Data,1)

    if strcmp(model.target,'cardiac_event')
        Data.censor = isnat(Data.cardiac_event_date)|isnat(Data.(model.settings.death_date));
        if ~isnat(Data.cardiac_event_date(i))
            Data.event_time(i) = ((hours(Data.cardiac_event_date(i) - Data.first_RT_date(i))/24)/365.25)*12;
        elseif Data.censor(i)
            Data.event_time(i) = ((hours(Data.(model.settings.censor_date)(i) - Data.first_RT_date(i))/24)/365.25)*12;
        else
            Data.event_time(i) = ((hours(Data.(model.settings.death_date)(i) - Data.first_RT_date(i))/24)/365.25)*12;
        end
    else
        % survival and propensity model (event time not used)
        Data.censor = isnat(Data.(model.settings.death_date));
        if ~isnat(Data.(model.settings.death_date)(i))
            Data.event_time(i) = ((hours(Data.(model.settings.death_date)(i) - Data.first_RT_date(i))/24)/365.25)*12;
        else
            Data.event_time(i) = ((hours(Data.(model.settings.censor_date)(i) - Data.first_RT_date(i))/24)/365.25)*12;
        end
    end
    if Data.vital_status(i)
        Data.survival(i) = ((hours(Data.(model.settings.censor_date)(i) - Data.first_RT_date(i))/24)/365.25)*12;
    else
        Data.survival(i) = ((hours(Data.(model.settings.death_date)(i) - Data.first_RT_date(i))/24)/365.25)*12;
    end
end

% Time difference between diagnosis and RT date as a feature
Data.time_to_treat = hours(Data.first_RT_date - Data.diagnosis_date)/24;
Data.time_to_treat((abs(Data.time_to_treat)>2000) | (Data.time_to_treat<-30)) = nan;
Data.time_to_treat((Data.time_to_treat<0)) = -1;


if ismember('multiple_centre_flag',Data.Properties.VariableNames)
    Data.time_to_cardiac = hours(Data.cardiac_event_date - Data.first_RT_date)/24;
    Data.cardiac_comorbidity(cell2mat(cellfun(@(x) isempty(x),Data.multiple_centre_flag,'UniformOutput',false)))={''};
end

Data.presc_dose_tx = Data.presc_dose_ttl./Data.presc_fractions;
alphabetaratio=10;
Data.eqd2 = Data.presc_dose_ttl.*((Data.presc_dose_tx + alphabetaratio)/(2 + alphabetaratio));
Compare = [Data.eqd2,Data.presc_dose_ttl,Data.deliv_dose_ttl,Data.presc_dose_ttl-Data.deliv_dose_ttl];
%figure; hist(Compare(Compare(:,4)~=0,4),20)
% curative_dose_threshold=44.99;
% survival_threshold = 24; % in months



Data.t_stage(ismember(Data.t_stage,'T0 Stage Finding'))={'TX Stage Finding'};

Data.rt_intent = repmat({'Curative'},size(Data,1),1);
Data.rt_intent(Data.(model.settings.dose_feature) < model.settings.curative_dose_threshold) = repmat({'Palliative'},sum(Data.(model.settings.dose_feature) < model.settings.curative_dose_threshold),1);
Data.rt_intent(Data.(model.settings.dose_feature) < model.settings.curative_dose_threshold) = repmat({'Palliative'},sum(Data.(model.settings.dose_feature) < model.settings.curative_dose_threshold),1);

category_var = {'histology','overall_stage','gender', 'smoking', 'n_stage', 't_stage','ecog','weightloss','laterality',...
    'tumour_loc','tumour_grade','marital_status','alcoholism','cardiac_comorbidity','technique','rt_intent','cause_death'};



%%%%% function on selection criteria analysis - returns numbers of
%%%%% patients and relevant data sets (curative and palliative)

[data_curative, data_palliative]=apply_lung_study_selection_criteria(Data, model, true);

%%%% run again but for constrained time period patients between 2011-2015(or 2016)
%         timespan_filter = year(Data.first_RT_date)>2010 & year(Data.first_RT_date)<2016;
%         Data5yr = Data(timespan_filter,:);
%         [data_curative5yr, data_palliative5yr]=apply_lung_study_selection_criteria(Data5yr, curative_dose_threshold, survival_threshold, true);
%
%         timespan_filter = year(Data.first_RT_date)>2010 & year(Data.first_RT_date)<2020;
%         Data6yr = Data(timespan_filter,:);
%         [data_curative6yr, data_palliative6yr]=apply_lung_study_selection_criteria(Data6yr, curative_dose_threshold, survival_threshold, true);

D = Data;
Dcp = [data_curative; data_palliative];
Data=[];
Data.table = Dcp;
Data.originaltable = D;

Data.data = prepLungData(Data.table);



end
