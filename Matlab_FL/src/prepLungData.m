function newdata = prepLungData(data)

survival_threshold=24;

columns = data.Properties.VariableNames;

newdata = data;
newdata.fev_missing = zeros(size(data,1),1);
newdata.fev_missing = double(isnan(data.fev));

newdata.tumour_volume = log(data.tumour_volume);
newdata.tumour_volume_missing = zeros(size(data,1),1);
newdata.tumour_volume_missing = double(isnan(data.tumour_volume));


% newdata.tumour_volume = -newdata.tumour_volume;

if ismember('rt_intent',columns)
    newdata.rt_intent = double(contains(data.rt_intent,'Curative'));
end

if ismember('gender',columns)
    newdata.gender = double(contains(data.gender,'Male'));
end
if ismember('laterality',columns) % introduce a missing
    newdata.lateral = ones(size(data,1),1);
    newdata.lateral(contains(data.laterality,'Left')) = 1;
    newdata.lateral(contains(data.laterality,'Right')) = 0;

    newdata.lateral(~contains(data.laterality,'Right') & ~contains(data.laterality,'Left')) = NaN;
end

% tumour lobe
if ismember('tumour_loc',columns) % introduce a missing
    newdata.tumourloc_middle = zeros(size(data,1),1);
    newdata.tumourloc_middle(contains(data.tumour_loc,'Middle')) = 1;
    newdata.tumourloc_lower = zeros(size(data,1),1);
    newdata.tumourloc_lower(contains(data.tumour_loc,'Lower')) = 1;
    newdata.tumourloc_bronchus = zeros(size(data,1),1);
    newdata.tumourloc_bronchus(contains(data.tumour_loc,'bronchus')) = 1;
    newdata.tumourloc_missing = double(strcmp(data.tumour_loc,''));

    newdata.tumourloc_middle(strcmp(data.tumour_loc,'')) = NaN;
    newdata.tumourloc_lower(strcmp(data.tumour_loc,'')) = NaN;
    newdata.tumourloc_bronchus(strcmp(data.tumour_loc,'')) = NaN;
end


if ismember('cardiac_comorbidity',columns)
    newdata.cardiac_comorbidity = zeros(size(data,1),1);
    newdata.cardiac_comorbidity = double(contains(data.cardiac_comorbidity,'NO'));
    
    newdata.cardiac_comorbidity(strcmp(data.cardiac_comorbidity,''))=NaN;
end
if ismember('histology',columns)
    newdata.hist_adeno = zeros(size(data,1),1);
    newdata.hist_adeno(contains(data.histology,'Adenocar')) = 1;
    newdata.hist_squamous = zeros(size(data,1),1);
    newdata.hist_squamous(contains(data.histology,'Squamous')) = 1;
    newdata.hist_largecell = zeros(size(data,1),1);
    newdata.hist_largecell(contains(data.histology,'Large ')) = 1;
end
if ismember('t_stage',columns)
    newdata.t2 = zeros(size(data,1),1); newdata.t2 = double(contains(data.t_stage,'T2')); newdata.t2(strcmp(data.t_stage,''))=NaN;
    newdata.t3 = zeros(size(data,1),1); newdata.t3 = double(contains(data.t_stage,'T3')); newdata.t3(strcmp(data.t_stage,''))=NaN;
    newdata.t4 = zeros(size(data,1),1); newdata.t4 = double(contains(data.t_stage,'T4')); newdata.t4(strcmp(data.t_stage,''))=NaN;
    newdata.tx = zeros(size(data,1),1); newdata.tx = double(contains(data.t_stage,'TX')); newdata.tx(strcmp(data.t_stage,''))=NaN;
end
if ismember('n_stage',columns)
    newdata.n1 = zeros(size(data,1),1); newdata.n1 = double(contains(data.n_stage,'N1')); newdata.n1(strcmp(data.n_stage,''))=NaN;
    newdata.n2 = zeros(size(data,1),1); newdata.n2 = double(contains(data.n_stage,'N2')); newdata.n2(strcmp(data.n_stage,''))=NaN;
    newdata.n3 = zeros(size(data,1),1); newdata.n3 = double(contains(data.n_stage,'N3')); newdata.n3(strcmp(data.n_stage,''))=NaN;
    newdata.nx = zeros(size(data,1),1); newdata.nx = double(contains(data.n_stage,'NX')); newdata.nx(strcmp(data.n_stage,''))=NaN;
end
if ismember('overall_stage',columns)
    newdata.sg2 = zeros(size(data,1),1); newdata.sg2 = double(strcmp(data.overall_stage,'Stage IIA')|strcmp(data.overall_stage,'Stage IIB'));
    newdata.sg3a = zeros(size(data,1),1); newdata.sg3a = double(strcmp(data.overall_stage,'Stage IIIA'));
    newdata.sg3b = zeros(size(data,1),1); newdata.sg3b = double(strcmp(data.overall_stage,'Stage IIIB'));
end

if ismember('ecog',columns) % reference point is ecog0
    newdata.ecog1 = double(data.ecog==1); newdata.ecog1(isnan(data.ecog))=NaN;
    newdata.ecog2 = double(data.ecog>=2); newdata.ecog2(isnan(data.ecog))=NaN;
%     newdata.ecog3 = double(data.ecog>=3);
    newdata.ecog_missing = double(isnan(data.ecog));
end
if ismember('tumour_grade',columns)
    newdata.tg_missing = strcmp(data.tumour_grade,'')|contains(data.tumour_grade,'X');
    newdata.tg2 = zeros(size(data,1),1); newdata.tg2 = double(contains(data.tumour_grade,'2'));
    newdata.tg2(newdata.tg_missing)=NaN;
    newdata.tg3 = zeros(size(data,1),1); newdata.tg3 = double(contains(data.tumour_grade,'3')); %%%%%%%%%%%%%%%%%%%%%%%********
    newdata.tg3(newdata.tg_missing)=NaN;
    newdata.tg4 = zeros(size(data,1),1); newdata.tg4 = double(contains(data.tumour_grade,'4'));
    newdata.tg4(newdata.tg_missing)=NaN;
%     newdata.tgx = zeros(size(data,1),1); newdata.tgx = double(contains(data.tumour_grade,'X'));
    
end
if ismember('smoking',columns)
    newdata.smoking_missing = double(strcmp(data.smoking,''));
    newdata.smoking_never = zeros(size(data,1),1);
    newdata.smoking_never = double(contains(data.smoking,'Never')); newdata.smoking_never(strcmp(data.smoking,''))=NaN;
    newdata.smoking_former = zeros(size(data,1),1);
    newdata.smoking_former = double(contains(data.smoking,'Former')); newdata.smoking_former(strcmp(data.smoking,''))=NaN;
end
if ismember('marital_status',columns)
%     newdata.marital_status = zeros(size(data,1),1);
%     newdata.marital_status = double(contains(data.marital_status,'Separated') | contains(data.marital_status,'Widow') | contains(data.marital_status,'Divorce') | contains(data.marital_status,'Never'));

    newdata.marital_never = zeros(size(data,1),1);
    newdata.marital_never = double(contains(data.marital_status,'Never')); newdata.marital_never(strcmp(data.marital_status,''))=NaN;
    newdata.marital_div = zeros(size(data,1),1);
    newdata.marital_div = double(contains(data.marital_status,'Divorce')); newdata.marital_div(strcmp(data.marital_status,''))=NaN;
    newdata.marital_sep = zeros(size(data,1),1);
    newdata.marital_sep = double(contains(data.marital_status,'Separated')); newdata.marital_sep(strcmp(data.marital_status,''))=NaN;
    newdata.marital_widow = zeros(size(data,1),1);
    newdata.marital_widow = double(contains(data.marital_status,'Widow')); newdata.marital_widow(strcmp(data.marital_status,''))=NaN;
end
if ismember('weightloss',columns)
    if isa(data.weightloss,'double')
        newdata.weightloss5=NaN(size(data,1),1);
        newdata.weightloss5_10=NaN(size(data,1),1);
        newdata.weightloss10=NaN(size(data,1),1);
    else
        newdata.weightloss5 = zeros(size(data,1),1);
        newdata.weightloss5 = double(contains(data.weightloss,'Less than 5%')); newdata.weightloss5(strcmp(data.weightloss,''))=NaN;
        newdata.weightloss5_10 = zeros(size(data,1),1);
        newdata.weightloss5_10 = double(contains(data.weightloss,'5 - 10%')); newdata.weightloss5_10(strcmp(data.weightloss,''))=NaN;
        newdata.weightloss10 = zeros(size(data,1),1);
        newdata.weightloss10 = double(contains(data.weightloss,'Greater than or equal to 10%')); newdata.weightloss10(strcmp(data.weightloss,''))=NaN;
    end
end

if ismember('alcoholism',columns)
    newdata.alcoholism = nan(size(data,1),1);
    if iscell(data.alcoholism(1))
        newdata.alcoholism = zeros(size(data,1),1);
        newdata.alcoholism = double(contains(data.alcoholism,'True'));
        newdata.alcoholism(strcmp(data.alcoholism,''))=NaN;
    end
end




if ~any(ismember(data.Properties.VariableNames,'death_date_rbdm'))
    newdata.vital_status = isnat(newdata.death_date);
    for i=1:length(newdata.vital_status)
        if newdata.vital_status(i)
            newdata.survival(i) = ((hours(newdata.censor_date(i) - newdata.first_RT_date(i))/24)/365.25)*12;
        else
            newdata.survival(i) = ((hours(newdata.death_date(i) - newdata.first_RT_date(i))/24)/365.25)*12;
        end
    end
    newdata.censored = newdata.vital_status & (((hours(newdata.censor_date - newdata.first_RT_date)/24)/365.25)*12 < survival_threshold);
   
    
else

    newdata.vital_status = isnat(newdata.death_date);
    newdata.vital_status_clinic = isnat(newdata.death_date_clinic);
    newdata.vital_status_rbdm = isnat(newdata.death_date_rbdm);
    %survival time in months
    for i=1:length(newdata.vital_status)
        if newdata.vital_status(i)
            newdata.survival(i) = ((hours(newdata.censor_date(i) - newdata.first_RT_date(i))/24)/365.25)*12;
        else
            newdata.survival(i) = ((hours(newdata.death_date(i) - newdata.first_RT_date(i))/24)/365.25)*12;
        end
        if newdata.vital_status_clinic(i)
            newdata.survival_clinic(i) = ((hours(newdata.censor_date_clinic(i) - newdata.first_RT_date(i))/24)/365.25)*12;
        else
            newdata.survival_clinic(i) = ((hours(newdata.death_date_clinic(i) - newdata.first_RT_date(i))/24)/365.25)*12;
        end
        if newdata.vital_status_rbdm(i)
            newdata.survival_rbdm(i) = ((hours(newdata.censor_date_rbdm(i) - newdata.first_RT_date(i))/24)/365.25)*12;
        else
            newdata.survival_rbdm(i) = ((hours(newdata.death_date_rbdm(i) - newdata.first_RT_date(i))/24)/365.25)*12;
        end
    end
    newdata.censored = newdata.vital_status & (((hours(newdata.censor_date - newdata.first_RT_date)/24)/365.25)*12 < survival_threshold);
end
% Time difference between diagnosis and RT date as a feature
newdata.RTdiagnosis_time_difference = ((hours(newdata.first_RT_date - newdata.diagnosis_date)/24)/365.25)*12;

end