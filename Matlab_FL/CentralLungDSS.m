% Central algorithm

% This script is intended to run first, it will wait for clients to
% connect to the analysis over the network. At each data client
% LocalLungModel.m will be run and it will establish connection and respond
% to requests from this script.


model = struct;
model.features = {'age_start_rt','gender','tumour_volume','time_to_treat','fev','dlco'...
        ,'cardiac_comorbidity'...
    ,'ecog1','ecog2'...
    ,'eqd2'...
    ,'t2','t3','t4','tx'...
    ,'n1','n2','n3','nx'...
    ,'sg2','sg3a','sg3b'...
    ,'lateral'...
    ,'hist_adeno','hist_squamous','hist_largecell'...
    ,'tg2','tg3','tg4'...
    ,'tumourloc_middle','tumourloc_lower','tumourloc_bronchus'...
    ,'smoking_never','smoking_former'...
    ,'smoking_pack_yr'...
    ,'marital_never','marital_div','marital_sep','marital_widow'...
    ,'alcoholism'...
    ,'weightloss5','weightloss5_10','weightloss10'...
    };

%%%% Missing <50%
model.features = {'age_start_rt','gender','tumour_volume','time_to_treat'...
        ,'cardiac_comorbidity'...
    ,'ecog1','ecog2'...
    ,'eqd2'...
    ,'t2','t3','t4','tx'...
    ,'n1','n2','n3','nx'...
    ,'sg2','sg3a','sg3b'...
    ,'lateral'...
    ,'hist_adeno','hist_squamous','hist_largecell'...
    ,'tumourloc_middle','tumourloc_lower','tumourloc_bronchus'...
    ,'marital_never','marital_div','marital_sep','marital_widow'...
    };

%%%% Missing <50% & P-value < 0.4 -------- AUC 0.68
model.features = {'age_start_rt','gender','tumour_volume','ecog2','eqd2'...
            ,'tumourloc_lower','tumourloc_bronchus','hist_adeno','n1','marital_sep','marital_div'...
     };
% % % %

% Propensity model
model.features = {'age_start_rt','gender','tumour_volume','ecog2'...
            ,'tumourloc_lower','tumourloc_bronchus','hist_adeno','n1','marital_sep','marital_div'...
     };
 

% Overall survival outcome model analysis: 
%    (1) RBDM outcome time span of data should be 2006-2015 with 2-year follow-up to Dec 2017
%    (2) Clinical outcome time span of data should be 2006-2015 with 2-year follow-up to Dec 2017
%    (3) Final time span of data should be 2010-2017 and 2018-2020 for validation.

model.target='treatment_intent'; %treatment_intent/survival
clients_names = {'client1','client2','client3','client4'};
clients_ids = {'Clinic A','Clinic B','Clinic C','Clinic D'}; % for deidentifying clients in output
model.sparql_query = readTextFile('.\sparql_queries\lung_dss.sparql'); %Query for lung model across 4 centres
model.algorithm = 'logreg';


model.settings.curative_dose_threshold=44.99;
model.settings.survival_threshold = 24; % in months
model.settings.calibration_bins = 0:0.05:1;  
model.settings.randseed = 0+[1:200]; 
Nboot = length(model.settings.randseed);
model.settings.predict_thresh=linspace(0,1,200);
model.settings.km_periods = 26;    
model.settings.km_limit = 5;
model.settings.death_date = 'death_date'; %death_date_clinic, death_date_rbdm, death_date
model.settings.censor_date = 'censor_date'; %censor_date_clinic, censor_date_rbdm, censor_date
model.settings.dose_feature = 'eqd2';

model.settings.timespan = [2010 2017 2020]; % training t(1)->t(2) and testing t(2)inclusive->t(3)
%Mode1 1 -> [2010 2017 2020]
%Model 2 (comparing RBDM/clinic) [2005 2016]

model.settings.rrss_prop=[0.6 0.4];
model.settings.SMOTE_balance=0;
model.settings.upsample_balance=0;
%Bootstrap models Nboot=100 in earlier time period and test in later period
model.settings.validation='time_bootstrap'; %'time_rrss';'time_bootstrap';bootstrap
% model.settings.varselect = ff2n(size(model.features,2));
% model.settings.varselect = model.settings.varselect(2:end,:);
% model.settings.varselect = [eye(size(model.features,2));ones(1,size(model.features,2))];
model.settings.varselect = ones(1,size(model.features,2));
model.settings.max_itr=10;
model.settings.lambda = 0.01; 
model.settings.km_conf=false;

% Port must align with local MPI service
mpi = ServerConnectionSetup('MPIClient','8084',strjoin(clients_names,','),true);

model.w = zeros(Nboot*size(model.settings.varselect,1),size(model.features,2));

model=initialise_network_data(mpi, model);

pcorr = pcorr_central(model);

DescriptiveStats;

model.itr=1;
iter=1;
maxiter = 15;
rng(200); 
model.w = 0.1*[randn(Nboot*size(model.settings.varselect,1),1),...
                randn(Nboot*size(model.settings.varselect,1),size(model.features,2))];
model.u=0.1*randn(size(model.w)); 
model.g=0.1*randn(size(model.w));
model.kappa=zeros(size(model.w,1),1);
model.convergence = false(Nboot*size(model.settings.varselect,1),1);
model.tol = 1e-8;
L=[]; E=[]; VL=[]; VE=[]; W=[model.w];

% 
while ~all(model.convergence) && iter < maxiter

    model.state = 'updateModel2';
    central_message = constructJSONMessage(model,[{'algorithm'},{'state'},{'w'},{'u'},{'Mu'},{'Std'}]);
    model_ = model;
    [model, model.itr] = Broadcast_Receive_from_MPI(mpi.svc, central_message, mpi.num_clients, model.itr, model_);
    [model, lik, error, vlik, verr] = update_logreg_model(model);
    W(:,:,iter)=model.w;
    L = [L, lik]; E = [E, error]; VL = [VL, vlik]; VE = [VE, verr];
    iter = iter+1;
end

model.state = 'testModelAUC2';
central_message = constructJSONMessage(model,[{'algorithm'},{'state'},{'w'},{'settings'},{'Mu'},{'Std'}]);
model_ = model;
[model, model.itr] = Broadcast_Receive_from_MPI(mpi.svc, central_message, mpi.num_clients, model.itr, model_);

featurelabels = model.features;
figure('color','w'); 
boxplot(model.w,'labels',[{'constant'},featurelabels],'PlotStyle','compact','orientation', 'horizontal');grid on
title('Model coefficients')

[AUC, TPR, FPR]=ROCfromConfMat3(model);
plotROC(AUC, TPR, FPR, [clients_ids, 'Overall'])

model.propensity.w = model.w;
model.propensity.features = model.features;
model.propensity.Mu = model.Mu;
model.propensity.Std = model.Std;

model = rmfield(model,{'Mu', 'Std'});
model.target='survival'; %treatment_intent/survival/cardiac_event
model.features = {'age_start_rt','gender','tumour_volume','ecog2','eqd2'...
            ,'tumourloc_lower','tumourloc_bronchus','hist_adeno','n1','marital_sep','marital_div'...
     };
model.settings.varselect = ones(1,size(model.features,2));
model.w = zeros(Nboot*size(model.settings.varselect,1),size(model.features,2));
model=initialise_network_data(mpi, model);

model.itr=1;
iter=1;
maxiter = 15;
model.w = 0.1*[randn(Nboot*size(model.settings.varselect,1),1),...
                randn(Nboot*size(model.settings.varselect,1),size(model.features,2))];
model.u=0.1*randn(size(model.w)); 
model.g=0.1*randn(size(model.w));
model.kappa=zeros(size(model.w,1),1);
model.convergence = false(Nboot*size(model.settings.varselect,1),1);
model.tol = 1e-8;
L=[]; E=[]; VL=[]; VE=[]; W=[model.w];

while ~all(model.convergence) && iter < maxiter

    model.state = 'updateModel2';
    central_message = constructJSONMessage(model,[{'algorithm'},{'state'},{'w'},{'u'},{'Mu'},{'Std'}]);
    model_ = model;
    [model, model.itr] = Broadcast_Receive_from_MPI(mpi.svc, central_message, mpi.num_clients, model.itr, model_);
    [model, lik, error, vlik, verr] = update_logreg_model(model);
    W(:,:,iter)=model.w;
    L = [L, lik]; E = [E, error]; VL = [VL, vlik]; VE = [VE, verr];
    iter = iter+1;
end




model.state = 'testModelAUC2';
central_message = constructJSONMessage(model,[{'algorithm'},{'state'},{'w'},{'settings'},{'Mu'},{'Std'},{'propensity'}]);
model_ = model;
[model, model.itr] = Broadcast_Receive_from_MPI(mpi.svc, central_message, mpi.num_clients, model.itr, model_);

featurelabels = {'Age at start of RT','Gender','Tumour volume','ECOG PS>=2','EQD2','Location: Lower Lobe',...
    'Location: Main bronchus','Adenocarcinoma','N stage=1','Marital - separated','Marital - divorced'};
% featurelabels = model.features;
figure('color','w'); 
boxplot(model.w,'labels',[{'constant'},featurelabels],'PlotStyle','compact','orientation', 'horizontal');grid on
title('Model coefficients')


figure('color','w'); 
boxplot(wl,'labels',[{'constant'},featurelabels([10 11 9:-1:1])],'PlotStyle','compact','orientation', 'horizontal');grid on
title('Model coefficients')


% figure; boxplot(model.w,'labels',[{'const'},model.features],'orientation', 'horizontal');grid on
 xlabel('Coefficients')
 
% figure('color','w'); boxplot(model.w,'labels',[{'constant','Age at start RT','Gender','Tumour volume','Cardiac comorbidity','ECOG PS >= 2'}],'PlotStyle','compact','orientation', 'horizontal');grid on
% figure('color','w'); boxplot(model.w,'labels',[{'constant','Age at start RT','Gender','Tumour volume','Cardiac comorbidity','Marital widow','Alcoholism', 'Tumour in lower lobe','Weight loss 5-10%'}],'PlotStyle','compact','orientation', 'horizontal');grid on
% xlabel('Coefficients')
% custom_save_fig(gcf, '.\figures\LogReg_weights1', {'png'}, {'900'})

%Compute BIC
valloglik = zeros(length(model.client(1).info),length(model.client(1).info(1).valloglik)); localNval = zeros(size(valloglik));
for j=1:length(model.client(1).info)
    for i=1:length(model.client)
         valloglik(j,:) = valloglik(j,:) + model.client(i).info(j).valloglik';
         localNval(j,:) = localNval(j,:) + model.client(i).info(j).localN';
    end
end
ms_vl = [mean(valloglik,2) std(valloglik,[],2) sum(model.settings.varselect,2) mean(localNval,2)];
[~,ordcomp]=sort(ms_vl(:,3)); ms_vl = ms_vl(ordcomp,:);
[mx, mxi]=max(ms_vl(:,1)); mx-ms_vl(mxi,2)
for i=1:max(ms_vl(:,3)); disp(max(ms_vl(ms_vl(:,3)==i,1))); end
BIC = (ms_vl(:,3)+1).*log(ms_vl(:,4)) - 2*ms_vl(:,1)
%figure; plot(BIC)

%figure('color','w'); boxplot(valloglik(ordcomp,:)','PlotStyle','compact','orientation', 'horizontal')
% figure; boxplot(model.info.P,'labels',[{'const'},model.features],'PlotStyle','compact','orientation', 'horizontal');grid on;

%Plot ROC and compute AUC
[AUC, TPR, FPR]=ROCfromConfMat3(model);
plotROC(AUC, TPR, FPR, [clients_ids, 'Overall'])
% custom_save_fig(gcf, '.\figures\ROCcurve1', {'png'}, {'900'})
charrel = [];
for i=1:length(model.client)
    disp(['AUC: ' model.client(i).client_name ': Mean: ' num2str(mean(AUC(i,:))) ', SD: ' num2str(std(AUC(i,:)))...
        '. Median: ' num2str(median(AUC(i,:))) ' [' num2str(prctile(AUC(i,:),2.5)) '-' num2str(prctile(AUC(i,:),97.5)) ']'])
    disp(['C-index: ' model.client(i).client_name ': Mean ' num2str(mean([model.client(i).validation_result.CHarrell_prob])) ', SD ' num2str(std([model.client(i).validation_result.CHarrell_prob]))...
        '. Median: ' num2str(median([model.client(i).validation_result.CHarrell_prob])) ' [' num2str(prctile([model.client(i).validation_result.CHarrell_prob],2.5)) '-' num2str(prctile([model.client(i).validation_result.CHarrell_prob],97.5)) ']'])
    charrel = [charrel, [model.client(i).validation_result.CHarrell_prob]];
end
disp(['Overall: Mean ' num2str(mean(AUC(end,:))) ', SD ' num2str(std(AUC(end,:)))...
    '. Median: ' num2str(median(AUC(end,:))) ' [' num2str(prctile(AUC(end,:),2.5)) '-' num2str(prctile(AUC(end,:),97.5)) ']'])
[median(charrel), prctile(charrel,2.5), prctile(charrel,97.5)]
% figure('Color','w'); boxplot(AUC','labels',[clients_ids {'Overall'}],'orientation', 'horizontal');grid on
% title('AUC');xlabel('AUC')

model = derive_PI_risk_groups(model, 'pi_cur');
% figure;plot(model.PI.histlevel,mean(model.PI.hist))
% figure;bar(model.PI.histlevel(2:end),diff(nanmean(model.PI.hist)))

FLHist_cur = prepFLhist(model, 'prop_cur');
FLHist_pal = prepFLhist(model, 'prop_pal');

matchPal = zeros(size(FLHist_pal));matchCur = zeros(size(FLHist_cur));

for m = 1:Nboot
    % match propensity
    matching_hist = min([sum(FLHist_pal(:,:,m),2) sum(FLHist_cur(:,:,m),2)],[],2);
    total_matching_samples = sum(matching_hist);
    for h=1:length(matching_hist)
        if matching_hist(h)==sum(FLHist_pal(h,:,m),2)
            matchPal(h,:,m) = FLHist_pal(h,:,m);
        elseif matching_hist(h) < sum(FLHist_pal(h,:,m),2)
            matchPal(h,:,m) = round(matching_hist(h)*FLHist_pal(h,:,m)/sum(FLHist_pal(h,:,m),2));
            if sum(matchPal(h,:,m)) > matching_hist(h)
                [mx,mxi]=max(matchPal(h,:,m));
                matchPal(h,mxi,m)=mx-1;
            elseif sum(matchPal(h,:,m)) < matching_hist(h)
                [mx,mxi]=max(matchPal(h,:,m));
                matchPal(h,mxi,m)=mx+1;
            end
        end
        if matching_hist(h)==sum(FLHist_cur(h,:,m),2)
            matchCur(h,:,m) = FLHist_cur(h,:,m);
        elseif matching_hist(h) < sum(FLHist_cur(h,:,m),2)
            matchCur(h,:,m) = round(matching_hist(h)*FLHist_cur(h,:,m)/sum(FLHist_cur(h,:,m),2));

            if sum(matchCur(h,:,m)) > matching_hist(h)
                [mx,mxi]=max(matchCur(h,:,m));
                matchCur(h,mxi,m)=mx-1;
            elseif sum(matchCur(h,:,m)) < matching_hist(h)
                [mx,mxi]=max(matchCur(h,:,m));
                matchCur(h,mxi,m)=mx+1;
            end
        end
    end
    
end

%[matching_hist sum(matchPal(:,:,m),2) sum(matchCur(:,:,m),2)] 
[matching_hist sum(FLHist_pal(:,:,m),2) sum(FLHist_cur(:,:,m),2)] 

model.prop_match.palliative = matchPal;
model.prop_match.curative = matchCur;

figure;
XD=[];
XData = [0.05:0.05:1]-0.025; 
YData = sum(FLHist_pal(:,:,m),2);
for i=1:length(XData)
    XD = [XD; XData(i)*ones(YData(i),1)];
end
h1 = histfit(XD,20,'kernel');
XD1 = h1(2).XData; YD1 = h1(2).YData;

XD=[];YData = sum(FLHist_cur(:,:,m),2);
for i=1:length(XData)
    XD = [XD; XData(i)*ones(YData(i),1)];
end
h1 = histfit(XD,20,'kernel');
XD2 = h1(2).XData; YD2 = h1(2).YData;

XD=[];YData = sum(matchPal(:,:,m),2);
for i=1:length(XData)
    XD = [XD; XData(i)*ones(YData(i),1)];
end
h1 = histfit(XD,20,'kernel');
XD3 = h1(2).XData; YD3 = h1(2).YData;
close(gcf)

bw = 0.3; lw = 3;
figure('Color','w');bar([0.05:0.05:1]-0.015,sum(FLHist_pal(:,:,m),2),bw,'facecolor','r', 'edgecolor','none','handlevisibility','off'); hold on
bar([0.05:0.05:1]-0.025,sum(FLHist_cur(:,:,m),2),bw,'facecolor','g', 'edgecolor','none','handlevisibility','off')
%bar([0.05:0.05:1]-0.035,sum(matchPal(:,:,m),2),bw,'facecolor','k', 'edgecolor','none')
line(XD1, YD1, 'linewidth', lw, 'color', 'r'); hold on; box off;
line(XD2, YD2, 'linewidth', lw, 'color', 'g');
line(XD3, YD3, 'linewidth', 5, 'color', 'k');
xlabel('Propensity score'); ylabel('Number of patients')
legend({'Palliative RT','Curative RT','Matched Cohort'})
set(gca,'xlim',[0 1])




model.state = 'km_curves';
central_message = constructJSONMessage(model,[{'algorithm'},{'state'},{'w'},{'prop_match'},{'settings'}]);
model_ = model;
[model, model.itr] = Broadcast_Receive_from_MPI(mpi.svc, central_message, mpi.num_clients, model.itr, model_);


labels = {'high','medium','low'};
% labels = {'high','low'};
KM=struct; bootKM=struct;
for m = 1:length(model.client(1).validation_km) % The three external validations

    for g = 1:length(model.client(1).validation_km(m).KM)
        KM(g).D = zeros(size(model.client(1).validation_km(m).KM(1).D));
        KM(g).N = zeros(size(model.client(1).validation_km(m).KM(1).N));
        KM(g).num_pts = 0;
        for i=1:length(model.client)
            KM(g).D = KM(g).D+model.client(i).validation_km(m).KM(g).D;
            KM(g).N = KM(g).N+model.client(i).validation_km(m).KM(g).N;
            KM(g).num_pts = KM(g).num_pts + model.client(i).validation_km(m).KM(g).num_pts;
        end
        bootKM.KM(g).D(:,m) = KM(g).D;
        bootKM.KM(g).N(:,m) = KM(g).N;
        bootKM.KM(g).num_pts(m) = KM(g).num_pts;
    end
%     [kmstat,p_value]=km_curves(KM, labels, model)
end
for g=1:length(KM)
    KM(g).D = median(bootKM.KM(g).D,2); %mean
    KM(g).N = median(bootKM.KM(g).N,2);
    KM(g).num_pts = round(median(bootKM.KM(g).num_pts));
end
% [kmstat,p_value]=km_curves(KM, labels, model);
% title('Kaplan Meier curves for prognostic groups')
% grid(gca,'on');
% custom_save_fig(gcf, '.\figures\KMcurve1', {'png'}, {'600'})


% %KM curves based on Stage
% KM=struct; labels = {'Stage I','Stage II','Stage IIIA','Stage IIIB','Stage IV','Stage missing'};
% for g = 1:length(model.client(1).stats.KM_stage_group)
%     KM(g).D = zeros(size(model.client(1).stats.KM_stage_group(1).D));
%     KM(g).N = zeros(size(model.client(1).stats.KM_stage_group(1).D));
%     KM(g).num_pts = 0;
%     for i=1:length(model.client)
%         KM(g).D = KM(g).D+model.client(i).stats.KM_stage_group(g).D;
%         KM(g).N = KM(g).N+model.client(i).stats.KM_stage_group(g).N;
%         KM(g).num_pts = KM(g).num_pts + model.client(i).stats.KM_stage_group(g).num_pts;
%     end
% end
% [kmstat,p_value]=km_curves(KM, labels, model);

%P-values from Hessian matrix
if isfield(model,'info'); model = rmfield(model,'info'); end
model.info = model.client(1).info(1);
for j=1:length(model.client(1).info)
    varnum = sum(model.settings.varselect(j,:))+1;
    model.info(j)=model.client(1).info(j);
    modelfields = fieldnames(model.client(1).info(1));
    for i=1:length(model.client)
        for k=1:length(modelfields)
            model.info(j).(modelfields{k}) = model.info(j).(modelfields{k}) + model.client(i).info(j).(modelfields{k});
        end 
    end
end
for j=1:length(model.client(1).info)
    varnum = sum(model.settings.varselect(j,:))+1;
    for i=1:Nboot
        se = sqrt(diag(inv(-(reshape(model.info(j).hessian(i,:),varnum,varnum)))))';
        zz = model.w(((j-1)*Nboot)+i,1:varnum)./se;
        model.info(j).P(i,:) = erfc(-zz/sqrt(2))/2;
        model.info(j).P(i,:) = 2*tcdf(-abs(zz),model.info(j).localN(i)-1); 
    end
end
% figure; boxplot(model.info.P,'labels',[{'const'},model.features],'PlotStyle','compact','orientation', 'horizontal');grid on;
% title('P-value significance of features from Hessian')
model.features([median(model.info.P(:,2:end))<=0.4])
model.features([median(model.info.P(:,2:end))>0.2])

for j=1:length(model.info)
    Np = model.info(j).Npos;
    NN = model.info(j).localN;
    
    oll = Np.*log(Np./NN)  +  (NN-Np).*log(ones(length(Np),1) - (Np./NN));
    fit = model.info(j).loglik;
    chisq=2*(fit-oll);
    ((oll-fit)./oll)';
    pvalue = chi2pdf(chisq,repmat(size(model.info(j).gradient,2)-1,length(chisq),1));
end



save(['.\models\Model_' date() char(datetime(datetime,'Format','_ddMMyy_HHmmSS'))])



% labels={'good prognosis, curative','poor prognosis, curative','good prognosis, palliative','poor prognosis, palliative'}; 
treatment_groups = {'curative', 'palliative'};

treatment_cohorts = {'KM','KM_pal'};
[KM_b_cur, KM_t_cur]=prep_km(model,treatment_cohorts{1});
[KM_b_pal, KM_t_pal, Nm, Nriskgroup]=prep_km(model,treatment_cohorts{2});

[overall_pal, overall_t_pal, ~, ~]=prep_km(model,'overall_pal');
[overall_cur, overall_t_cur, ~, ~]=prep_km(model,'overall_cur');
[pal_dss, pal_t_dss, ~, ~]=prep_km(model,'overall_pal_dss');
[cur_dss, cur_t_dss, ~, ~]=prep_km(model,'overall_cur_dss');

labels={'curative','palliative', 'curative DSS', 'palliative DSS'}; 
[kmstat,p_value]=km_curves_bootstrap([overall_cur(5,:) overall_pal(5,:) cur_dss(5,:) pal_dss(5,:)], labels, model);

[prop_pal_dss, prop_t_pal_dss, ~, ~]=prep_km(model,'prop_pal_dss');
[prop_cur_dss, prop_t_cur_dss, ~, ~]=prep_km(model,'prop_cur_dss');
[prop_pal, prop_t_pal, ~, ~]=prep_km(model,'prop_pal');
[prop_cur, prop_t_cur, ~, ~]=prep_km(model,'prop_cur');
[overall_prop_pal, ~, ~, ~]=prep_km(model,'overall_prop_pal');
[overall_prop_cur, ~, ~, ~]=prep_km(model,'overall_prop_cur');
[cur_low, ~, ~, ~]=prep_km(model,'cur_low');
[cur_high, ~, ~, ~]=prep_km(model,'cur_high');


% KM curves for above and below median prognostic index
model.settings.km_conf=true;
labels={'curative-high', 'curative-low'}; 
[kmstat,p_value]=km_curves_bootstrap([cur_low(5,:) cur_high(5,:)], labels, model);

% KM curves of all risk groups and treatments
model.settings.km_conf=false;
labels={'curative-high', 'curative-moderate', 'curative-low', 'palliative-high ', 'palliative-moderate', 'palliative-low'}; 
[kmstat,p_value]=km_curves_bootstrap([prop_cur(5,:) prop_pal(5,:)], labels, model);

model.settings.km_conf=true;
n=5;
labels={'curative-low', 'curative-moderate', 'curative-high', 'palliative-low ', 'palliative-moderate', 'palliative-high'}; 
[kmstat,p_value, pv]=km_curves_bootstrap([KM_b_cur(n,3) KM_b_cur(n,2) KM_b_cur(n,1) KM_b_pal(n,3) KM_b_pal(n,2) KM_b_pal(n,1)], labels, model);

% KM curves of all risk groups and treatments with CI
model.settings.km_conf=true;
labels={'curative-low', 'curative-moderate', 'curative-high', 'palliative-low ', 'palliative-moderate', 'palliative-high'}; 
[~,p_value,pv]=km_curves_bootstrap([prop_cur(5,3) prop_cur(5,2) prop_cur(5,1) prop_pal(5,3) prop_pal(5,2) prop_pal(5,1)], labels, model);


prop_cur_sc = prop_cur(5,:);
prop_cur_ov = prop_cur_sc(1);
prop_pal_sc = prop_pal(5,:);
prop_pal_ov = prop_pal_sc(1);

for m=1:Nboot
    sc1 = (prop_pal(5,3).num_pts(m)+prop_cur(5,3).num_pts(m))/prop_cur(5,3).num_pts(m);
    prop_cur_sc(:,3).num_pts(m) = round(prop_cur(5,3).num_pts(m)*sc1);
    prop_cur_sc(:,3).N(:,m) = round(prop_cur(5,3).N(:,m)*sc1);
    prop_cur_sc(:,3).D(:,m) = round(prop_cur(5,3).D(:,m)*sc1);
    
    sc1 = (prop_pal(5,2).num_pts(m)+prop_cur(5,2).num_pts(m))/prop_cur(5,2).num_pts(m);
    prop_cur_sc(:,2).num_pts(m) = round(prop_cur(5,2).num_pts(m)*sc1);
    prop_cur_sc(:,2).N(:,m) = round(prop_cur(5,2).N(:,m)*sc1);
    prop_cur_sc(:,2).D(:,m) = round(prop_cur(5,2).D(:,m)*sc1);
    
    prop_cur_ov.num_pts(m) = prop_cur_sc(1).num_pts(m)+prop_cur_sc(2).num_pts(m)+prop_cur_sc(3).num_pts(m);
    prop_cur_ov.N(:,m) = prop_cur_sc(1).N(:,m) + prop_cur_sc(2).N(:,m)+prop_cur_sc(3).N(:,m);
    prop_cur_ov.D(:,m) = prop_cur_sc(1).D(:,m) + prop_cur_sc(2).D(:,m)+prop_cur_sc(3).D(:,m);
     
%     prop_cur_ov.num_pts(m) = prop_cur_sc(2).num_pts(m)+prop_cur_sc(3).num_pts(m);
%     prop_cur_ov.N(:,m) = prop_cur_sc(2).N(:,m)+prop_cur_sc(3).N(:,m);
%     prop_cur_ov.D(:,m) = prop_cur_sc(2).D(:,m)+prop_cur_sc(3).D(:,m);
% 
%     sc1 = (prop_cur(5,1).num_pts(m)+prop_pal(5,1).num_pts(m))/prop_pal(5,1).num_pts(m);
%     prop_pal_ov.num_pts(m) = round(prop_pal(5,1).num_pts(m)*sc1);
%     prop_pal_ov.N(:,m) = round(prop_pal(5,1).N(:,m)*sc1);
%     prop_pal_ov.D(:,m) = round(prop_pal(5,1).D(:,m)*sc1);

end

%KM curves comparing overall curative/palliative changes with reclassifying risk groups
[~,p_value]=km_curves_bootstrap([overall_prop_cur(5,:) overall_prop_pal(5,:) prop_cur_ov prop_pal_ov], {'curative', 'palliative', 'curative DSS', 'palliative DSS'}, model);

%KM curves comparing only overall curative changes with reclassifying risk groups
[~,p_value]=km_curves_bootstrap([overall_prop_cur(5,:) prop_cur_ov], {'curative', 'curative DSS'}, model);


% model.settings.km_conf=true;
% labels={'high risk','moderate risk','low risk'}; 
% % labels={'high risk','low risk'}; 
% [~,p_value]=km_curves_bootstrap(KM_b_cur(5,:), labels, model);
% set(gcf,'Position',[750 218 488 347]);
% [kmstat,p_value]=km_curves_bootstrap(KM_b_pal(5,:), labels, model);
% set(gcf,'Position',[750 218 488 347]);



num_curve = length(KM_t_cur(5,:));
logrank_stat=zeros(num_curve);
p_value=zeros(num_curve);
for j=1:num_curve
    for k=1:num_curve
        [logrank_stat(j,k), p_value(j,k)] = KM_logrank(KM_t_cur(5,j).N(:,2),KM_t_cur(5,j).D(:,2),...
                                                       KM_t_pal(5,k).N(:,2),KM_t_pal(5,k).D(:,2),0.05,...
                                                       [num2str(j) ' ' num2str(k)], [num2str(j) ' ' num2str(k)]);
    end
end



A = p_value; t = 'Log-rank';
for i=1:length(A); A(i,i)=1; end
figure('Color','w');
image(A'*128)
colormap(autumn(128))
% set(gca,'XTick',[],'YTick',[])
labels={'curative-low', 'curative-moderate', 'curative-high', 'palliative-low ', 'palliative-moderate', 'palliative-high'}; 
set(gca, 'YTicklabel', labels, 'XTicklabel',labels)
xtickangle(30)
[x,y] = meshgrid(1:size(pv,1),1:size(pv,1));
text(x(:),y(:),num2str(A(:),2),'HorizontalAlignment','center')
title(t)

% figure('Color','w');
% plotCumulHist(model,'validation_result','pi_pal', 'r');
% plotCumulHist(model,'validation_result','pi_cur', 'g');
% legend({'Palliative RT', 'Curative RT'})
% xlabel('Prognostic Index'); ylabel('Cumulative probability');
% grid on; box off;
% 
% figure('Color','w');
% plotCumulHist(model,'validation_hist','prop_pal', 'r');
% plotCumulHist(model,'validation_hist','prop_cur', 'g');
% legend({'Palliative RT', 'Curative RT'})
% xlabel('Prognostic Index'); ylabel('Cumulative probability');
% grid on; box off;
% 
% figure('Color','w');
% plotCumulHist(model,'validation_hist','prop_pal_rg_low', 'g');
% plotCumulHist(model,'validation_hist','prop_pal_rg_med', 'b');
% plotCumulHist(model,'validation_hist','prop_pal_rg_high', 'r');
% box off; grid on;
% title('Palliative RT')
% xlabel('Prognostic Index'); ylabel('Cumulative probability');
% legend({'Low Risk', 'Medium Risk', 'High Risk'})
% 
% figure('Color','w');
% plotCumulHist(model,'validation_hist','prop_cur_rg_low', 'g');
% plotCumulHist(model,'validation_hist','prop_cur_rg_med', 'b');
% plotCumulHist(model,'validation_hist','prop_cur_rg_high', 'r');
% box off; grid on;
% xlabel('Prognostic Index'); ylabel('Cumulative probability');
% title('Curative RT')
% legend({'Low Risk', 'Medium Risk', 'High Risk'})



model.state = 'calibration';
central_message = constructJSONMessage(model,[{'algorithm'},{'state'},{'w'},{'settings'}]);
model_ = model;
[model, model.itr] = Broadcast_Receive_from_MPI(mpi.svc, central_message, mpi.num_clients, model.itr, model_);



% Number of patients per risk group
numrisk=zeros(length(model.client),length(model.client(1).validation_km(1).KM),length(model.client(1).validation_km));
propsurvrisk=zeros(length(model.client),length(model.client(1).validation_km(1).KM),length(model.client(1).validation_km));
for i = 1:length(model.client)
    for m = 1:length(model.client(i).validation_km)
        for g = 1:length(model.client(i).validation_km(1).KM)
            numrisk(i,g,m) = numrisk(i,g,m) + model.client(i).validation_km(m).KM(g).num_pts;
            propsurvrisk(i,g,m) = propsurvrisk(i,g,m) + model.client(i).validation_km(m).KM(g).N(26);
        end
    end
    
end

nummodel = length(model.client(1).validation_result);
calibration_data = zeros(length(model.client(1).validation_result(1).actual_outcome_dist),2,length(model.client),nummodel);

for j=1:length(model.client(1).validation_result)
    for i=1:length(model.client)
        calibration_data(:,1,i,j) = model.client(i).validation_result(j).actual_outcome_dist;
        calibration_data(:,2,i,j) = model.client(i).validation_result(j).predicted_outcome_dist; 
    end
end


overall_calib = sum(calibration_data,4)/nummodel;
overall_calib = cumsum(sum(sum(calibration_data,4)/nummodel,3));

calib_mean = squeeze(sum(calibration_data,3));
calib_mean = calib_mean(1:1:end,:,:);
[fc_all, b] = plotModelCalibration(calib_mean, model.settings.calibration_bins(1:1:end-1)+0.05);

size(calib_mean)
size(calibration_data)

for k=1:size(model.settings.PI_risk_group,1)
    cp(k) = mean(model.settings.PI_risk_group(k,:));
end

for j=1:length(model.client(1).validation_result)
    for i=1:length(model.client)
        calibration_rg(1,:,j,i) = model.client(i).validation_result(j).calib_rg_act;
        calibration_rg(2,:,j,i) = model.client(i).validation_result(j).calib_rg_pred;
        calibration_rg(3,:,j,i) = model.client(i).validation_result(j).calib_x;
    end
end

% calibration_rg = squeeze(calibration_rg(:,:,:,4));%for one clinic
calibration_rg = squeeze(sum(calibration_rg,4))/4;%for overall

obs_prob_ci = squeeze(prctile(calibration_rg(1,:,:)./calibration_rg(2,:,:),[2.5,50,97.5],3));
calibration_rg = prctile(calibration_rg,[2.5 50 97.5],3);
yy = obs_prob_ci';
xx = [squeeze(calibration_rg(3,:,:))]';

x = xx(2,:); y = yy(2,:);
errorbar(xx(2,:),yy(2,:),yy(2,:)-yy(1,:),yy(3,:)-yy(2,:),xx(2,:)-xx(1,:),xx(3,:)-xx(2,:),'s','MarkerSize',10,...
    'MarkerEdgeColor',[0.4660 0.6740 0.1880],'MarkerFaceColor',[0.4660 0.6740 0.1880],'color','k','linewidth',2)
plot([0 1],[0 1],'-k')


[b,bint,r,rint,stats] = regress(y',[ones(length(x),1) x']);
disp(b)
brier_score = sum((x'-y').^2)/length(x);
disp(brier_score)
plot([0:0.01:1],b(1)+b(2)*[0:0.01:1],'--k','linewidth',2)
set(gca,'xlim',[0 1], 'ylim',[0 1])
annotation('textbox',[.12 .67 .3 .3],'String',{['Slope: ' num2str(b(2),'%.3f')],['Intercept: ' num2str(b(1),'%.3f')],['R^2: ' num2str(stats(1),'%.3f')], ['Brier: ' num2str(brier_score,'%.3f')]},'FitBoxToText','on');

%
surv_years = zeros(26,2);
for i=1:4
    surv_years = surv_years + model.client(i).stats_cur.first_RT_date.survival;
end
figure; bar(model.client(i).stats_cur.first_RT_date.date_ranges,surv_years(:,1)./surv_years(:,2))



model.state='terminate_algorithm';
central_message = constructJSONMessage(model,[{'algorithm'},{'state'}]);
model_ = model;
[model, model.itr] = Broadcast_Receive_from_MPI(mpi.svc, central_message, mpi.num_clients, model.itr, model_);


