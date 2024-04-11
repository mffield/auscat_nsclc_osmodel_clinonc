function [stats] = cohort_stats_lung(Data, model, central_stats, normalise)

% descriptive stats for table 1 of article (2003-2017) & (2011-2015) & (2016-2017)
% tumour_volume, fev1, ecog, age, gender, dose, histology, stages
% (TNM and overall), smoking, surgery, 2yr survival (yes, no,
% censored, lost to followup), previous body RT,

verbose=true;
disp('****STATS***')
[data_curative, data_palliative, dataset_total, stats]=apply_lung_study_selection_criteria(Data, model, verbose);


continuous_var = {'fev', 'dlco', 'smoking_pack_yr', 'tumour_volume', 'age_start_rt',...
    'presc_dose_ttl', 'eqd2', 'survival', 'time_to_treat', 'time_to_cardiac', 'event_time'}; %
for i=1:length(continuous_var)
    if ismember(continuous_var{i},Data.Properties.VariableNames)
        stats.(continuous_var{i}).mean = mean(Data.(continuous_var{i})(~isnan(Data.(continuous_var{i}))));
        stats.(continuous_var{i}).std = std(Data.(continuous_var{i})(~isnan(Data.(continuous_var{i}))));
        stats.(continuous_var{i}).sum = sum(Data.(continuous_var{i})(~isnan(Data.(continuous_var{i}))));
        stats.(continuous_var{i}).N = sum(~isnan(Data.(continuous_var{i})));
        stats.(continuous_var{i}).median = median(Data.(continuous_var{i})(~isnan(Data.(continuous_var{i}))));
        if isfield(model,'central_stats')
            stats.(continuous_var{i}).stdsum = sum((Data.(continuous_var{i})(~isnan(Data.(continuous_var{i}))) - model.central_stats.(continuous_var{i}).mean).^2);
        end
    end
end

stats.event_time.event_thresh = sum(Data.event_time>=model.settings.survival_threshold);
stats.event_time.event_thresh_cens = sum(isnat(Data.death_date) & Data.event_time<model.settings.survival_threshold);  

stats.survival.surv2yr = sum(Data.survival>=model.settings.survival_threshold);
stats.survival.surv2yr_cens = sum(isnat(Data.death_date) & Data.survival<model.settings.survival_threshold);  



NrPtPerBin=10;
for i=1:length(continuous_var)
    if any(ismember(Data.Properties.VariableNames,continuous_var{i}))
        NumSamples = sum(~isnan(Data.(continuous_var{i})));
        CumHistLevels=NrPtPerBin/NumSamples;
        CumHistLevels=0:CumHistLevels:1;
        CumHistLevels=prctile(Data.(continuous_var{i})(~isnan(Data.(continuous_var{i}))),100*CumHistLevels);
        CumHist=arrayfun(@(x) sum(Data.(continuous_var{i})(~isnan(Data.(continuous_var{i})))<x),CumHistLevels);
        CumHist(1)=0;
        CumHist(end)=NumSamples;
        stats.(continuous_var{i}).CumHist=CumHist/NumSamples;
        stats.(continuous_var{i}).CumHistLevels = CumHistLevels;
    end
end



stats.presc_dose_ttl.num_category = [sum(Data.presc_dose_ttl<40); sum(Data.presc_dose_ttl>=40 & Data.presc_dose_ttl<48);
    sum(Data.presc_dose_ttl>=48 & Data.presc_dose_ttl<56); sum(Data.presc_dose_ttl>=57 & Data.presc_dose_ttl<64);
    sum(Data.presc_dose_ttl>=64 & Data.presc_dose_ttl<68);sum(Data.presc_dose_ttl>=68)];
stats.presc_dose_ttl.category = {'<40Gy','>=40Gy and <48Gy','>=48Gy and <57Gy','>=57Gy and <64Gy','>=64Gy and <68Gy','>=68Gy'};

stats.eqd2.num_category = [sum(Data.eqd2<40); sum(Data.eqd2>=40 & Data.eqd2<48);
    sum(Data.eqd2>=48 & Data.eqd2<57); sum(Data.eqd2>=57 & Data.eqd2<64);
    sum(Data.eqd2>=64 & Data.eqd2<68);sum(Data.eqd2>=68)];
stats.eqd2.category = {'<40Gy','>=40Gy and <48Gy','>=48Gy and <57Gy','>=57Gy and <64Gy','>=64Gy and <68Gy','>=68Gy'};



date_ranges = 1995:2020;
date_var = {'first_RT_date'};
for i=1:length(date_var)
    if any(ismember(Data.Properties.VariableNames, date_var{i}))
        stats.(date_var{i}).date_ranges = date_ranges;
        stats.(date_var{i}).date_hist = histc(year(Data.(date_var{i})),date_ranges);
        for j=1:length(date_ranges)
            stats.(date_var{i}).survival(j,1) = sum(Data.survival(year(Data.(date_var{i}))==date_ranges(j))>=model.settings.survival_threshold);
            stats.(date_var{i}).survival(j,2) = sum(year(Data.(date_var{i}))==date_ranges(j));
        end
    end
end

date_ranges = 1995:2020;
date_var = {'death_date_clinic', 'death_date_rbdm', 'death_date'};
for i=1:length(date_var)
    if any(ismember(Data.Properties.VariableNames, date_var{i}))
        stats.(date_var{i}).date_ranges = date_ranges;
        stats.(date_var{i}).date_hist = histc(year(Data.(date_var{i})),date_ranges);
    end
end

LinkedData = Data(cell2mat(cellfun(@(x) ~isempty(x),Data.multiple_centre_flag,'UniformOutput',false)),:);
stats.num_unq_clinic = sum(~isnat(LinkedData(isnat(LinkedData.death_date_rbdm),:).death_date_clinic) & year(LinkedData(isnat(LinkedData.death_date_rbdm),:).death_date_clinic)<2018);
stats.num_unq_rbdm = sum(~isnat(LinkedData(isnat(LinkedData.death_date_clinic),:).death_date_rbdm));

% category = {'YES','NO'};
% for i=1:length(category)
%     stats.multiple_centre_flag.num(i) = sum(ismember(Data.multiple_centre_flag,category{i}));
% end
% stats.multiple_centre_flag.category = category;


%categorical variables
category = {'Adenocarcinoma','Non-Small Cell Lung Carcinoma','Large Cell Carcinoma','Squamous Cell Carcinoma','Small Cell Carcinoma'};
for i=1:length(category)
    stats.histology.num(i) = sum(ismember(Data.histology,category{i}));
end
stats.histology.category = category;

Data.overall_stage(ismember(Data.overall_stage,'Stage IA'))={'Stage I'};
Data.overall_stage(ismember(Data.overall_stage,'Stage IB'))={'Stage I'};
Data.overall_stage(ismember(Data.overall_stage,'Stage IIA'))={'Stage II'};
Data.overall_stage(ismember(Data.overall_stage,'Stage IIB'))={'Stage II'};
% Data.overall_stage(ismember(Data.overall_stage,'Stage IIIA'))={'Stage III'};
% Data.overall_stage(ismember(Data.overall_stage,'Stage IIIB'))={'Stage III'};
Data.overall_stage(ismember(Data.overall_stage,'Stage IIIC'))={'Stage IIIB'};
Data.overall_stage(ismember(Data.overall_stage,'Stage IVA'))={'Stage IV'};
Data.overall_stage(ismember(Data.overall_stage,'Stage IVB'))={'Stage IV'};
Data.overall_stage(ismember(Data.overall_stage,'Stage IVC'))={'Stage IV'};

category = {'Stage I','Stage II','Stage IIIA','Stage IIIB','Stage IV'};
for i=1:length(category)
    stats.overall_stage.num(i) = sum(ismember(Data.overall_stage,category{i}));
end
stats.overall_stage.category = category;


category = {'N0 Stage Finding','N1 Stage Finding','N2 Stage Finding','N3 Stage Finding','NX Stage Finding'};
for i=1:length(category)
    stats.n_stage.num(i) = sum(ismember(Data.n_stage,category{i}));
end
stats.n_stage.category = category;


Data.t_stage(ismember(Data.t_stage,'T1a Stage Finding'))={'T1 Stage Finding'};
Data.t_stage(ismember(Data.t_stage,'T1b Stage Finding'))={'T1 Stage Finding'};
Data.t_stage(ismember(Data.t_stage,'T1c Stage Finding'))={'T1 Stage Finding'};
Data.t_stage(ismember(Data.t_stage,'T2a Stage Finding'))={'T2 Stage Finding'};
Data.t_stage(ismember(Data.t_stage,'T2b Stage Finding'))={'T2 Stage Finding'};
Data.t_stage(ismember(Data.t_stage,'T2c Stage Finding'))={'T2 Stage Finding'};
category = {'T0 Stage Finding','T1 Stage Finding','T2 Stage Finding','T3 Stage Finding','T4 Stage Finding','TX Stage Finding'};
for i=1:length(category)
    stats.t_stage.num(i) = sum(ismember(Data.t_stage,category{i}));
end
stats.t_stage.category = category;


category = [0:5];
for i=1:length(category)
    stats.ecog.num(i) = sum(ismember(Data.ecog,category(i)));
end
stats.ecog.category = {'ECOG 0','ECOG 1','ECOG 2','ECOG 3','ECOG 4','ECOG 5'};

% 
if isa(Data.gender(1),'double')
    category = {1;0};
else
    category = {'Male','Female'};
end
for i=1:length(category)
    stats.gender.num(i) = sum(ismember(Data.gender,category{i}));
end
stats.gender.category = {'Male','Female'};


category = {'Left','Right'};
for i=1:length(category)
    stats.laterality.num(i) = sum(ismember(Data.laterality,category{i}));
end
stats.laterality.category = category;

category = {'Upper lobe of lung','Middle lobe of right lung','Main bronchus','Lower lobe of lung'};
for i=1:length(category)
    stats.tumour_loc.num(i) = sum(ismember(Data.tumour_loc,category{i}));
end
stats.tumour_loc.category = category;


category = {'Grade 1','Grade 2','Grade 3','Grade 4','Grade X'};
for i=1:length(category)
    stats.tumour_grade.num(i) = sum(ismember(Data.tumour_grade,category{i}));
end
stats.tumour_grade.category = category;


category = {'None','Less than 5%','5 - 10%','Greater than or equal to 10%'};
for i=1:length(category)
    stats.weightloss.num(i) = sum(ismember(Data.weightloss,category{i}));
end
stats.weightloss.category = category;


category = {'Married','Divorced','Never Married','Separated','Widowed'};
for i=1:length(category)
    stats.marital_status.num(i) = sum(ismember(Data.marital_status,category{i}));
end
stats.marital_status.category = category;

category = {'Never Smoker','Current Smoker','Former Smoker'};
for i=1:length(category)
    stats.smoking.num(i) = sum(ismember(Data.smoking,category{i}));
end
stats.smoking.category = category;

category = {'Endoscopic Resection','Lobectomy of Lung','Pneumonectomy','Segmental resection of lung','Wedge Resection of Lung'};
for i=1:length(category)
    stats.surgery.num(i) = sum(ismember(Data.lung_surgery,category{i}));
end
stats.surgery.category = category;

category = {'YES','NO'};
for i=1:length(category)
    stats.prev_body_rt.num(i) = sum(ismember(Data.has_prev_body_rt,category{i}));
end
stats.prev_body_rt.category = category;

% category = {'YES','NO'};
% for i=1:length(category)
%     stats.multiple_centre_flag.num(i) = sum(ismember(Data.multiple_centre_flag,category{i}));
% end
% stats.multiple_centre_flag.category = category;

category = {'YES','NO'};
for i=1:length(category)
    stats.cardiac_comorbidity.num(i) = sum(ismember(Data.cardiac_comorbidity,category{i}));
end
stats.cardiac_comorbidity.category = category;

category = {'True','False'};
for i=1:length(category)
    stats.alcoholism.num(i) = sum(ismember(Data.alcoholism,category{i}));
end
stats.alcoholism.category = category;

category = {'Parallel-Opposed Fields Radiation Therapy','3-Dimensional Conformal Radiation Therapy','Intensity-Modulated Radiation Therapy','Volume Modulated Arc Therapy'};
for i=1:length(category)
    stats.technique.num(i) = sum(ismember(Data.technique,category{i}));
end
stats.technique.category = category;

category = {'Curative','Palliative'};
for i=1:length(category)
    stats.rt_intent.num(i) = sum(ismember(Data.rt_intent,category{i}));
end
stats.rt_intent.category = {'Curative','Palliative'};
stats.rt_intent.date_ranges = 1995:2020; 
for i=1:length(stats.rt_intent.date_ranges)
    stats.rt_intent.intent_hist(i,:)=[sum(ismember(Data.rt_intent(year(Data.first_RT_date)==stats.rt_intent.date_ranges(i)),'Curative')),...
        sum(ismember(Data.rt_intent(year(Data.first_RT_date)==stats.rt_intent.date_ranges(i)),'Palliative'))];
end

if ismember({'cause_of_death'},Data.Properties.VariableNames)

    Data.cause_of_death(ismember(Data.cause_of_death,{'Malignant neoplasm of bronchus and lung','Malignant neoplasm of trachea, bronchus, and lung'}))={'Lung cancer'};
    Data.cause_of_death(ismember(Data.cause_of_death,{'Malignant melanoma of skin','Malignant neoplasm of anus and anal canal',...
        'Malignant neoplasm of bladder','Malignant neoplasm of brain','Malignant neoplasm of breast','Malignant neoplasm of cervix uteri',...
        'Malignant neoplasm of colon','Malignant neoplasm of kidney, except renal pelvis','Malignant neoplasm of larynx',...
        'Malignant neoplasm of liver and intrahepatic bile ducts','Malignant neoplasm of meninges','Malignant neoplasm of oesophagus',...
        'Malignant neoplasm of other and ill-defined digestive organs','Malignant neoplasm of other and ill-defined sites in the lip, oral cavity and pharynx',...
        'Malignant neoplasm of other and unspecified parts of mouth','Malignant neoplasm of other and unspecified parts of tongue',...
        'Malignant neoplasm of pancreas','Malignant neoplasm of prostate','Malignant neoplasm of rectum','Malignant neoplasm of stomach',...
        'Malignant neoplasm of tonsil','Malignant neoplasm without specification of site','Malignant neoplasms of independent (primary) multiple sites',...
        'Mesothelioma','Multiple myeloma and malignant plasma cell neoplasms','Other malignant neoplasms of skin'}))={'Other cancer'};
    Data.cause_of_death(ismember(Data.cause_of_death,{'Other acute ischaemic heart diseases','Acute myocardial infarction',...
        'Chronic ischaemic heart disease','Heart failure'}))={'Cardiac event'};
    Data.cause_of_death(ismember(Data.cause_of_death,{'Pneumonitis due to solids and liquids','Other interstitial pulmonary diseases',...
        'Other chronic obstructive pulmonary disease','Other and unspecified types of non-Hodgkin''s lymphoma','Emphysema'}))={'Other'};
    
    Data.cause_of_death(isnat(Data.death_date)) = {'NA'};
    Data.cause_of_death(cell2mat(cellfun(@(x) isempty(x),Data.cause_of_death,'UniformOutput',false)) & ~isnat(Data.death_date)) = {'Unknown'};
    
    category = {'Lung cancer','Other cancer','Cardiac event','Other','Unknown','NA'};
    for i=1:length(category)
        stats.cause_of_death.num(i) = sum(ismember(Data.cause_of_death,category{i}));
    end
    stats.cause_of_death.category = category;
    

end


end