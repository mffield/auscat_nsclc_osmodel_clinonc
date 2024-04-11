function [Data_curative, Data_palliative, DataSet, stats]=apply_lung_study_selection_criteria(Data, model, verbose)


Data.censored = Data.vital_status & (((hours(Data.censor_date - Data.first_RT_date)/24)/365.25)*12 < model.settings.survival_threshold);
Data.lost_followup = ((hours(Data.censor_date - Data.first_RT_date)/24)/365.25 < 2) & isnat(Data.death_date) & (hours(datetime() - Data.first_RT_date)/24)/365.25 > 2.5;

DA = Data(:,{'histology','overall_stage','eqd2', 'eqd2', 'censored','tumour_volume'});

%SBRT if prescribed fractional dose is above 6Gy and total dose above 25Gy

selectDA = zeros(size(DA));
selectDA(ismember(Data.histology',{'Adenocarcinoma','Non-Small Cell Lung Carcinoma','Large Cell Carcinoma','Squamous Cell Carcinoma'}),1)=1;
selectDA(ismember(Data.histology,{'Small Cell Carcinoma'}),1)=2;
selectDA(ismember(Data.overall_stage,{'Stage I','Stage IA','Stage IB','Stage II','Stage IIA','Stage IIB','Stage III','Stage IIIA','Stage IIIB','Stage Not IV'}),2)=1;
selectDA(ismember(Data.overall_stage,{'Stage IV','Stage IVA','Stage IVB'}),2)=2;
selectDA(Data.(model.settings.dose_feature)>=model.settings.curative_dose_threshold,3)=1;
selectDA(Data.(model.settings.dose_feature)<model.settings.curative_dose_threshold,3)=2;
selectDA(ismember(Data.lung_surgery,{''}),4)=1;
selectDA(~ismember(Data.lung_surgery,{''}),4)=2;
selectDA(ismember(Data.has_prev_body_rt,{'NO'}),5)=1;
selectDA(ismember(Data.has_prev_body_rt,{'YES'}),5)=2;
selectDA(~((Data.presc_dose_ttl./Data.presc_fractions >= 6) & (Data.presc_dose_ttl >= 25)),6) = 1; %SBRT
selectDA((Data.presc_dose_ttl./Data.presc_fractions >= 6) & (Data.presc_dose_ttl >= 25),6) = 2; %SBRT
selectDA(Data.censored==0,7)=1;
selectDA(Data.censored==1,7)=2;
selectDA(~isnan(Data.tumour_volume),8)=1;

pts_curative = sum(selectDA(:,3)==1);
pts_palliative = sum(selectDA(:,3)==2);
pts_ruled_out = sum(any(selectDA(:,[1 2 4 5 6])==2,2));

pts_select_cur_no_imaging = sum(all(selectDA(:,1:6)==1,2));
Data_curative_no_image = Data(isnan(Data.tumour_volume) & all(selectDA(:,1:6)==1,2),:);

TNMHist_missing = sum((selectDA(:,1)==0) | (selectDA(:,2)==0));

Data_curative = Data(all(selectDA(:,[1 2 3 4 5 6 ])==1,2),:); % [1 2 3 4 5 6 8]

yrfilt = (year(Data.first_RT_date)>2010) & (year(Data.first_RT_date)<2020);

disp('Exclusion')
disp(['Histology:' num2str(sum((selectDA(yrfilt,1)==0) | (selectDA(yrfilt,1)==2)))]) %exclude histology
disp(['Staging:' num2str(sum((selectDA(yrfilt,2)==0) | (selectDA(yrfilt,2)==2)))]) %exclude staging
disp(['Surgery:' num2str(sum((selectDA(yrfilt,4)==0) | (selectDA(yrfilt,4)==2)))]) %exclude surgery
disp(['Prev RT:' num2str(sum((selectDA(yrfilt,5)==0) | (selectDA(yrfilt,5)==2)))]) %exclude prev RT
disp(['SBRT:' num2str(sum((selectDA(yrfilt,6)==0) | (selectDA(yrfilt,6)==2)))]) %exclude SBRT

numpatinperiod = sum((year(Data.first_RT_date)>2010) & (year(Data.first_RT_date)<2020));
disp(['Num patients in period 2011-2019:' num2str(numpatinperiod)])


% change to select palliative cohort
selectDA(Data.(model.settings.dose_feature)<model.settings.curative_dose_threshold,3)=1;
selectDA(Data.(model.settings.dose_feature)>=model.settings.curative_dose_threshold,3)=2;
Data_palliative = Data(all(selectDA(:,[1 2 3 4 5 6 ])==1,2),:); %[1 2 3 4 5 6 8]
pts_select_pall_no_imaging = sum(all(selectDA(:,1:6)==1,2));

stats.pts_inelig_stage = sum(any(selectDA(:,2)~=1,2));
stats.pts_inelig_hist = sum(any(selectDA(:,1)~=1,2));
stats.pts_inelig_surgery = sum(any(selectDA(:,4)~=1,2));
stats.pts_inelig_prevRT = sum(any(selectDA(:,5)~=1,2));
stats.pts_inelig_sbrt = sum(any(selectDA(:,6)~=1,2));
stats.pts_noimage = sum(any(selectDA(:,8)==0,2));



selectDA(Data.(model.settings.dose_feature)>=model.settings.curative_dose_threshold,3)=1;
selectDA(Data.(model.settings.dose_feature)<model.settings.curative_dose_threshold,3)=2;

N = size(Data,1);

lost_followup = ((hours(Data.censor_date - Data.first_RT_date)/24)/365.25 < 2) & isnat(Data.death_date) & (hours(datetime() - Data.first_RT_date)/24)/365.25 > 2.5;

if verbose
disp(['Total number of patients: ' num2str(N)])
disp(['Total number of curative patients: ' num2str(pts_curative)])
disp(['Total number of palliative patients: ' num2str(pts_palliative)])

% how many lost to follow up? not dead and last contact date < 2yrs from RT and RT was more than 2yrs ago.

disp(['Number of patients lost to follow up: ' num2str(sum(lost_followup))])

disp(['Number of patients that definitively did not fit selection criteria: ' num2str(pts_ruled_out)])
disp(['Number of patients that may fit selection criteria: ' num2str(N-pts_ruled_out)])
disp(['Number of patients that we are not sure fit selection criteria based on missing data: ' num2str(N - pts_ruled_out - pts_select_cur_no_imaging - pts_select_pall_no_imaging)])

disp(['Number of patients that fit the selection criteria of model (curative) without imaging: ' num2str(pts_select_cur_no_imaging)])
disp(['Number of patients that fit the selection criteria (palliative) without imaging: ' num2str(pts_select_pall_no_imaging)])

disp(['Number of patients that fit the selection criteria of model (curative): ' num2str(size(Data_curative,1))])
disp(['Number of patients that fit the selection criteria (palliative): ' num2str(size(Data_palliative,1))])
end

stats.N = N;
stats.pts_curative = pts_curative;
stats.pts_palliative = pts_palliative;
stats.lost_followup = sum(lost_followup);
stats.pts_ruled_out = pts_ruled_out;
stats.pts_select_cur_no_imaging = pts_select_cur_no_imaging;
stats.pts_select_pall_no_imaging = pts_select_pall_no_imaging;
stats.pts_select_cur = size(Data_curative,1);
stats.pts_select_pall = size(Data_palliative,1);

stats.TNMHist_missing = TNMHist_missing;


DataSet = [Data_curative; Data_palliative];

end