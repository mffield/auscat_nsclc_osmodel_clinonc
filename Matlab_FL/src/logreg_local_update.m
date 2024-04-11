function [local_message,model,mpi,itr]=logreg_local_update(mpi,itr,central_message,model,data)
% Logistic regression - implemented through conjugate gradient ascent
%
% beta: step size
% u: search direction
% g: gradient of the weights
% H: Hessian
%
local_message=[];

message_tokens = strsplit(central_message,';');
model.state = message_tokens(ismember(message_tokens,'state')+1);

[model, received_parameters] = readJSONMessage('local',model,central_message,model.alg_client_index);

curative_thresh = model.settings.curative_dose_threshold;
timespan = model.settings.timespan;
surv_thresh = model.settings.survival_threshold;
dose_feature = model.settings.dose_feature;

if strcmp(model.state,'initialise')
         
    if isfield(model,'info'); model = rmfield(model,'info'); end
    model = init_logreg_model('local',model,model.features,[],[],[]);
    [model, ~] = readJSONMessage('local',model,central_message,model.alg_client_index);
    % %%%%%% Imputation 
    % First approach is a direct imputation using all local data not
    % separate into train/validation - do this when data set is loaded into memory
    
    %%%%%%%%%%%%%%%%%%%%%%%%
    
    % Second approach is to impute on training data locally. 
    % Third approach is a network wide imputation.
        
    Data = data.data(:,model.features);
    Data = table2array(Data);
    
    [model.stats_original]=cohort_stats_lung(data.originaltable, model, [], false);
    [model.stats]=cohort_stats_lung(data.table, model, [], false);
    
    [itr,model,central_message,mpi]=sendMessage(mpi,model,itr, {'client_name','stats','stats_original'});
    [model, ~] = readJSONMessage('local',model,central_message,model.alg_client_index);
    
    censorfilter = ~(data.table.vital_status & (data.table.survival<surv_thresh));
    yearfilter = (year(data.table.first_RT_date)>timespan(1)) & (year(data.table.first_RT_date)<timespan(2));
    yearfilter2 = (year(data.originaltable.first_RT_date)>timespan(1)) & (year(data.originaltable.first_RT_date)<timespan(2));
    yearfilter3 = (year(data.table.first_RT_date)>timespan(1)) & (year(data.table.first_RT_date)<timespan(3));
    
    %%%%%
    yearfilter = yearfilter3;
    %%%%%
    
    cur_filter=(data.data.(dose_feature) > curative_thresh);
    
    pcorrcoeff = corr(Data(cur_filter,:),(data.data.event_time(cur_filter)>=surv_thresh)); 
    for i = 1:size(Data,2)
        mm = ~isnan(Data(:,i)) & censorfilter & cur_filter & yearfilter;
        y = double(((data.data.event_time(mm)>=surv_thresh)));
        x = Data(mm,i);
        model.pcorr.n(i) = sum(mm);
        model.pcorr.xs(i) = sum(x); model.pcorr.xs2(i) = sum(x.^2); 
        model.pcorr.ys(i) = sum(y); model.pcorr.ys2(i) = sum(y.^2);
        model.pcorr.xys(i) = sum(x.*y);
    end
    disp('Pearson correlation:')
    disp(pcorrcoeff')
       
    if isfield(model,'Mu')
        model = rmfield(model,{'Mu','Std'});
    end
    
    
    [model.stats_original]=cohort_stats_lung(data.originaltable(yearfilter2,:), model, 'central_stats_original', false);
    [model.stats]=cohort_stats_lung(data.table(yearfilter,:), model, 'central_stats', true);
    [model.stats_cur]=cohort_stats_lung(data.table(yearfilter & (data.table.(dose_feature)>curative_thresh)...
        ,:), model,[], false);
    model.stats_cur=normalise_features(model,data.data,model.stats_cur);
%     [model.stats_pal]=cohort_stats_lung(data.data((ff_yr>timespan(1)) & (ff_yr<timespan(2)) & (data.data.(dose_feature)<=curative_thresh)...
%         ,:), model, 'central_stats', true);
    
    [itr,model,central_message,mpi]=sendMessage(mpi,model,itr, {'client_name','stats','stats_original','stats_cur'});
    [model, ~] = readJSONMessage('local',model,central_message,model.alg_client_index);
    
    model.alg_client_index = find(ismember(model.settings.client_names,model.client_name));
    
    if isfield(model.stats_cur,'localN')
        model.localN = model.stats_cur.localN;
        model.localSum = model.stats_cur.localSum;
    end
    
    [model.stats_original]=cohort_stats_lung(data.originaltable(yearfilter2,:), model, 'central_stats_original', false);
    [model.stats]=cohort_stats_lung(data.table(yearfilter,:), model, 'central_stats', true);
    [model.stats_cur]=cohort_stats_lung(data.table(yearfilter & cur_filter,:), model, 'central_stats_cur', true);
    model.stats_cur=normalise_features(model,data.data,model.stats_cur);
    [model.stats_pal]=cohort_stats_lung(data.table(yearfilter & ~cur_filter,:), model, 'central_stats', true);
    
%     testfilter = data.test.index;
%     y_train_c = table2array(data.data(yearfilter & cur_filter & ~testfilter,'event_time'))>=surv_thresh;
%     y_valid_c = table2array(data.data(yearfilter & cur_filter,'event_time'))>=surv_thresh;
%     y_train_p = table2array(data.data(yearfilter & ~cur_filter & ~testfilter,'event_time'))>=surv_thresh;
%     y_valid_p = table2array(data.data(yearfilter & ~cur_filter,'event_time'))>=surv_thresh;
%     
%     y_test_c = table2array(data.data(yearfilter & cur_filter & testfilter,'event_time'))>=surv_thresh;
%     y_test_p = table2array(data.data(yearfilter & ~cur_filter & testfilter,'event_time'))>=surv_thresh;
    
    
    %Balance of data sets
%     model.stats.train_outcomes_curative = [sum(y_train_c==1) sum(y_train_c==0)];
%     model.stats.valid_outcomes_curative = [sum(y_valid_c==1) sum(y_valid_c==0)];
%     model.stats.train_outcomes_palliative = [sum(y_train_p==1) sum(y_train_p==0)];
%     model.stats.valid_outcomes_palliative = [sum(y_valid_p==1) sum(y_valid_p==0)];
%     
%     model.stats.test_outcomes_curative = [sum(y_test_c==1) sum(y_test_c==0)];
%     model.stats.test_outcomes_palliative = [sum(y_test_p==1) sum(y_test_p==0)];

    
    stagegroup = zeros(length(data.originaltable.overall_stage),1);
    stagegroup(ismember(data.originaltable.overall_stage,'Stage I'))=1;
    stagegroup(ismember(data.originaltable.overall_stage,'Stage IA'))=1;
    stagegroup(ismember(data.originaltable.overall_stage,'Stage IB'))=1;
    stagegroup(ismember(data.originaltable.overall_stage,'Stage IIA'))=2;
    stagegroup(ismember(data.originaltable.overall_stage,'Stage IIB'))=2;
    stagegroup(ismember(data.originaltable.overall_stage,'Stage IIIA'))=3;
    stagegroup(ismember(data.originaltable.overall_stage,'Stage IIIB'))=4;
    stagegroup(ismember(data.originaltable.overall_stage,'Stage IIIC'))=4;
    stagegroup(ismember(data.originaltable.overall_stage,'Stage IV'))=5;
    stagegroup(ismember(data.originaltable.overall_stage,'Stage IVA'))=5;
    stagegroup(ismember(data.originaltable.overall_stage,'Stage IVB'))=5;
    stagegroup(ismember(data.originaltable.overall_stage,'Stage IVC'))=5;
    stagegroup(ismember(data.originaltable.overall_stage,'missing'))=6; 
    stagegroup(ismember(data.originaltable.overall_stage,'Stage Not IV'))=6; 
    [model.stats.KM_stage_group] = km_stats(data.originaltable(yearfilter2,:),stagegroup(yearfilter2),6,model);

    
    local_message = {'client_name','stats','stats_original','stats_cur','stats_pal','pcorr'};

    

elseif strcmp(model.state,'updateModel2') 
    
    
    curative_thresh = model.settings.curative_dose_threshold;
    timespan = model.settings.timespan;
    surv_thresh = model.settings.survival_threshold;
    
    if min(size(model.settings.varselect))==1; model.settings.varselect=model.settings.varselect(:)'; end
    data.normdata = data.data;
    Nboot = length(model.settings.randseed);
    numModel = size(model.settings.varselect,1);
%     model.loglik = zeros(numModel*Nboot,1);
%     model.valloglik = zeros(numModel*Nboot,1);
%     model.gradient = zeros(numModel,size(model.w,2));
%     model.hessian = zeros(numModel,(size(model.w,2))^2);
%     model.localN = zeros(numModel*Nboot,1);
%     model.Npos = zeros(numModel*Nboot,1);
%     model.localN_val = zeros(numModel*Nboot,1);
%     model.Npos_val = zeros(numModel*Nboot,1);
    
    x = data.data(:,model.features);

    %%%%% Impute with local mean 
    x(:,model.features) = fillmissing(x(:,model.features),'constant',model.Mu(:)');
    x_table = x;
    x = (table2array(x) - repmat(model.Mu(:)',size(x,1),1))./repmat(model.Std(:)',size(x,1),1);

    N = size(x,1);
    X = [ones(N,1),x]';
    
    if (strcmp(model.target,'survival') || strcmp(model.target,'cardiac_event')) 
        y = double(data.data.event_time>=surv_thresh);
    elseif strcmp(model.target,'treatment_intent')
        y = data.data.rt_intent; 
    end
    y(y==0)=-1; y=y(:)';
    
    
    
    for j=1:numModel
        w=model.w(((j-1)*Nboot+1):j*Nboot,:);
        if min(size(model.settings.varselect))==1
            varselect = logical(model.settings.varselect);
        else
            varselect = logical(model.settings.varselect(j,:));
        end
        
        w = w(:,[true true(1,sum(varselect))]);
        Xt = X([true varselect(:)'],:);
        
        tt = find(all(Xt'==0)); %| (std(Xt')<0.01)
        if ~isempty(tt)
            for i=1:length(tt)
                Xt(tt(i),:)=Xt(tt(i),:) + 0.05*randn(1,size(Xt,2));
            end
        end
        
        
        
        for i=1:Nboot
            % curative/palliative filter?
            %& (data.data.(dose_feature) > curative_thresh)
            
            train_index = false(size(data.data,1),1); train_index(data.train(i).index)=true; 
            val_index = false(size(data.data,1),1); val_index(data.valid(i).index)=true; 
            if (strcmp(model.target,'survival') || strcmp(model.target,'cardiac_event'))
            	filter = ~(data.data.censor & (data.data.event_time<surv_thresh)) & (data.data.(dose_feature) > curative_thresh);
                train_index = train_index & filter;            
                val_index = val_index & filter;
            end
            
            % Address data imbalance in training - check y(train_index)
            % balance and upsample smaller group

%             data.data(train_index,:)
            
            if model.settings.SMOTE_balance
                [Xb, yb] = SMOTE(Xt(:,train_index)', y(train_index)');
                yb(yb==0)=-1; yb=yb(:)'; Xb=Xb';
            elseif model.settings.upsample_balance
                nneg = sum(y(train_index)==-1); npos = sum(y(train_index)==1);
                percent_balance = 100*[nneg npos]/length(y);
                if percent_balance(1)<(percent_balance(2)-10)
                    dataind = find(y(train_index)==-1);
                    dataind_kcopies = repmat(dataind,20,1); %increase up to max 20 times
                    dataind_kcopies = dataind_kcopies(randperm(numel(dataind_kcopies)));
                    train_index = [dataind_kcopies(1:npos) find(y(train_index)==1)];
                elseif percent_balance(2)<(percent_balance(1)-10)
                    dataind = find(y(train_index)==1);
                    dataind_kcopies = repmat(dataind,20,1); %increase up to max 20 times
                    dataind_kcopies = dataind_kcopies(randperm(numel(dataind_kcopies)));
                    train_index = [dataind_kcopies(1:nneg) find(y(train_index)==-1)];
                end
                Xb = Xt(:,train_index);
                yb = y(train_index);yb(yb==0)=-1;yb=yb(:)'; 
            else
                Xb = Xt(:,train_index);
                yb = y(train_index);
            end
            
            model.info(j).balance = sum(yb==1)/length(yb);
%             
%             t = templateTree('Reproducible',true); % For reproducibiliy of random predictor selections
%             Mdl = fitcensemble(Xb',yb','Method','Bag','Learners',t)
%             
%             nTrees=300;
%             B = TreeBagger(nTrees,Xb',yb', 'Method', 'classification'); 
%             predChar1 = B.predict(Xt(:,val_index)');  % Predictions is a char though. We want it to be a number.
%             c = str2double(predChar1);
%             consistency=sum(c==y(val_index)')/length(y(val_index));


            [gradient, hessian, LL] = logreg_evaluate_grad(w(i,:),Xb,yb, model.settings.lambda);
            
            [model.info(j).loglik(i), model.info(j).trainSumError(i), prob_tr] = test_logreg(w(i,:), Xb,yb);
            model.info(j).gradient(i,:) = gradient(:)';
            model.info(j).hessian(i,:) = hessian(:)';
            model.info(j).localN(i) = sum(train_index);
            model.info(j).localN_val(i) = sum(val_index);
            model.info(j).Npos(i) = sum(y(train_index)>0);
            model.info(j).Npos_val(i) = sum(y(val_index)>0);
            
%             [model.loglik(((j-1)*Nboot)+i), model.trainSumError(((j-1)*Nboot)+i), prob_tr] = test_logreg(w(i,:), Xt(:,train_index),y(train_index));
%             model.gradient(((j-1)*Nboot)+i,:) = [gradient(:)' zeros(1,size(model.gradient,2)-length(gradient(:)'))];
%             model.hessian(((j-1)*Nboot)+i,:) = [hessian(:)' zeros(1,size(model.hessian,2)-length(hessian(:)'))];
%             
%             model.localN(((j-1)*Nboot)+i) = sum(train_index);
%             model.localN_val(((j-1)*Nboot)+i) = sum(val_index);
%             model.Npos(((j-1)*Nboot)+i) = sum(y(train_index)>0);
%             model.Npos_val(((j-1)*Nboot)+i) = sum(y(val_index)>0);
            
            if sum(val_index)>0
%                 [model.valloglik(((j-1)*Nboot)+i), model.valSumError(((j-1)*Nboot)+i),  prob_val] = test_logreg(w(i,:), Xt(:,val_index),y(val_index));
                [model.info(j).valloglik(i), ~,  prob_val] = test_logreg(w(i,:), Xt(:,val_index),y(val_index));
                %model.info(j).valSumError(i)
            else
%                 model.valloglik(((j-1)*Nboot)+i) = NaN;
                model.info(j).valloglik(i) = NaN;
            end
            
        end
        
        
    end    
    
 local_message = {'client_name','info'}; 
   
 
 
elseif strcmp(model.state,'testModelAUC2')    
    
    
    if min(size(model.settings.varselect))==1; model.settings.varselect=model.settings.varselect(:)'; end
    data.normdata = data.data;
    Nboot = length(model.settings.randseed);
    numModel = size(model.settings.varselect,1);
    model.valloglik = zeros(numModel*Nboot,1);
    
    curfilter = (data.data.(dose_feature) > curative_thresh);
    palfilter = (data.data.(dose_feature) <= curative_thresh);
    filter = ~(data.data.censor & (data.data.event_time<surv_thresh));
    
    temp = data.data.(dose_feature)(palfilter);
    data.data.(dose_feature)(palfilter) = 60;

    x = data.data(:,model.features);
    x(:,model.features) = fillmissing(x(:,model.features),'constant',model.Mu(:)');
    x = (table2array(x) - repmat(model.Mu(:)',size(x,1),1))./repmat(model.Std(:)',size(x,1),1);
       
    N = size(x,1);
    X = [ones(N,1),x]';
    
    if isfield(model,'propensity')
        % data for propensity model
        px = data.data(:,model.propensity.features);
        px(:,model.propensity.features) = fillmissing(px(:,model.propensity.features),'constant',model.propensity.Mu(:)');
        px = (table2array(px) - repmat(model.propensity.Mu(:)',size(px,1),1))./repmat(model.propensity.Std(:)',size(px,1),1);
        PX = [ones(size(px,1),1),px]';
        
    end
    
    
    if (strcmp(model.target,'survival') || strcmp(model.target,'cardiac_event')) 
        y = double(data.data.event_time>=surv_thresh);
    elseif strcmp(model.target,'treatment_intent')
        y = data.data.rt_intent; 
    end
    y(y==0)=-1; y=y(:)';
    
    
    
    for j=1:numModel
        w=model.w(((j-1)*Nboot+1):j*Nboot,:);
        if min(size(model.settings.varselect))==1
            varselect = logical(model.settings.varselect);
        else
            varselect = logical(model.settings.varselect(j,:));
        end

        w = w(:,[true true(1,sum(varselect))]);
        Xt = X([true varselect(:)'],:);
        for i=1:Nboot
            modelind =((j-1)*Nboot)+i;

            model.training_result(i,j).confusionMatrix=zeros(2,2,length(model.settings.predict_thresh));
            model.training_result(i,j).predicted_outcome_dist=zeros(length(model.settings.calibration_bins)-1,size(model.w,2));
            model.training_result(i,j).actual_outcome_dist=zeros(length(model.settings.calibration_bins)-1,size(model.w,2));
            model.validation_result(i,j).confusionMatrix=zeros(2,2,length(model.settings.predict_thresh));
            model.validation_result(i,j).predicted_outcome_dist=zeros(length(model.settings.calibration_bins)-1,size(model.w,2));
            model.validation_result(i,j).actual_outcome_dist=zeros(length(model.settings.calibration_bins)-1,size(model.w,2));
           
            val_ind = false(size(data.data,1),1); 
            if strcmp(model.settings.validation,'time_bootstrap') || strcmp(model.settings.validation,'time')
                val_ind(data.test.index)=true;
            else
                val_ind(data.valid(i).index)=true;
            end 
%             train_index = false(size(data.data,1),1);
%             train_index(data.train(i).index)=true;
%             train_index(data.valid(i).index)=true;
%             val_ind = train_index;
            
            val_index = val_ind;
            val_index2 = val_ind;
            if (strcmp(model.target,'survival') || strcmp(model.target,'cardiac_event')) 

                val_index = val_ind & filter & curfilter;
                val_index2 = val_ind & curfilter; % no censor

            end

            if sum(val_index)>0
                [~, ~,  model.pi(i,j).pi] = test_logreg(w(i,:), Xt,y);
                [model.validation_result(i,j).valloglik, model.validation_result(i,j).valSumError,  prob_val] = test_logreg(w(i,:), Xt(:,val_index),y(val_index));
                [~, hessian, ~] = logreg_evaluate_grad(w(i,:),Xt(:,val_index),y(val_index),model.settings.lambda);
                model.validation_result(i,j).hessian = hessian;
            else
                model.validation_result(i,j).valloglik = NaN;
            end
            
            for k=1:length(model.settings.predict_thresh)
                y_outcome = ones(size(prob_val)); % assign the class labels based on a chosen threshold (default to zero)
                y_outcome(prob_val<model.settings.predict_thresh(k)) = -1;
                y_outcome(prob_val>=model.settings.predict_thresh(k)) = 1;
                acc = sum(y_outcome==y(val_index))/length(y(val_index)); % compute the prediction accuracy
                [ConfMat, ~] = confusionmat(y(val_index),y_outcome,'order',[1 -1]); % compute the confusion matrix
                model.validation_result(i,j).confusionMatrix(:,:,k) = ConfMat;
            end
%             [~,~,~,model.validation_result(modelind).AUC_train]=perfcurve(y(train_index),prob_val,1);
            [~,~,~,model.validation_result(i,j).AUC_val]=perfcurve(y(val_index),prob_val,1);
            model.validation_result(i,j).CHarrell_prob=concordance_index(data.data.event_time(val_index2),~data.data.censor(val_index2),model.pi(i,j).pi(val_index2));
            model.pi(i,j).pi_cur=prob_val;
            model.pi(i,j).pi_pal=model.pi(i,j).pi(val_ind & filter & palfilter);
            model.validation_result(modelind).predicted_outcome_dist = histcounts(prob_val,model.settings.calibration_bins);
            model.validation_result(modelind).actual_outcome_dist = histcounts(prob_val(y(val_index)==1),model.settings.calibration_bins);
            
            
            NrPtPerBin=5;
            model.validation_result(modelind).pi_cur = FLhistogram(prob_val, NrPtPerBin);
            model.validation_result(modelind).pi_pal = FLhistogram(model.pi(i,j).pi_pal, NrPtPerBin);
            model.validation_result(modelind).pi = FLhistogram(model.pi(i,j).pi, NrPtPerBin);

            [model.validation_result(modelind).pi_cur.CumHist,...
            model.validation_result(modelind).pi_cur.CumHistLevels,...
            model.validation_result(modelind).pi_cur.NumSamples,...
            model.validation_result(modelind).pi_cur.minmax] = cumul_hist(prob_val, NrPtPerBin);
            
            [model.validation_result(modelind).pi_pal.CumHist,...
            model.validation_result(modelind).pi_pal.CumHistLevels,...
            model.validation_result(modelind).pi_pal.NumSamples,...
            model.validation_result(modelind).pi_pal.minmax] = cumul_hist(model.pi(i,j).pi_pal, NrPtPerBin);
        
            [model.validation_result(modelind).pi.CumHist,...
            model.validation_result(modelind).pi.CumHistLevels,...
            model.validation_result(modelind).pi.NumSamples,...
            model.validation_result(modelind).pi.minmax] = cumul_hist(model.pi(i,j).pi, NrPtPerBin);
        
            if isfield(model,'propensity')
                [~, ~,  prop_prob] = test_logreg(model.propensity.w(i,:), PX,zeros(1,size(PX,2)));
                model.pi(i,j).prop_pi=prop_prob;
                [model.validation_result(modelind).prop_cur.CumHist,...
                    model.validation_result(modelind).prop_cur.CumHistLevels,...
                    model.validation_result(modelind).prop_cur.NumSamples,...
                    model.validation_result(modelind).prop_cur.minmax] = cumul_hist(prop_prob(val_ind & curfilter), NrPtPerBin);
        
                [model.validation_result(modelind).prop_pal.CumHist,...
                    model.validation_result(modelind).prop_pal.CumHistLevels,...
                    model.validation_result(modelind).prop_pal.NumSamples,...
                    model.validation_result(modelind).prop_pal.minmax] = cumul_hist(prop_prob(val_ind & palfilter), NrPtPerBin);
            
                binwidth = 0.05;
                [model.validation_result(modelind).prop_cur.PHist,model.validation_result(modelind).prop_cur.PHistLevels]=histcounts(prop_prob(val_ind & curfilter),0:binwidth:1); 
                [model.validation_result(modelind).prop_pal.PHist,model.validation_result(modelind).prop_pal.PHistLevels]=histcounts(prop_prob(val_ind & palfilter),0:binwidth:1); 
            
            end
        
        
        end
    end
    data.data.(dose_feature)(palfilter) = temp;        
    local_message = {'client_name','validation_result'};
        
       


elseif strcmp(model.state,'km_curves')   

    % assumes that the state preceding was 'testModelAUC2'
    Nboot = length(model.settings.randseed);
    
    if min(size(model.settings.varselect))==1; model.settings.varselect=model.settings.varselect(:)'; end
    numModel = size(model.settings.varselect,1);
    
    for j=1:numModel
        for i=1:Nboot
            modelind =((j-1)*Nboot)+i;
    
            val_index = false(size(data.data,1),1); 
            if strcmp(model.settings.validation,'time_bootstrap') || strcmp(model.settings.validation,'time')
                val_index(data.test.index)=true;
            else
                val_index(data.valid(i).index)=true;
            end
            %             test_index = false(size(data.data,1),1); test_index(data.test.index)=true;

            tr_index = false(size(data.data,1),1); 
            tr_index(data.train(i).index)=true;
            
            curfilt = (data.table.(dose_feature) > curative_thresh);
            palfilt = (data.table.(dose_feature) <= curative_thresh);
 
            model.pi(modelind).g = zeros(size(model.pi(modelind).pi));
            rg = model.settings.PI_risk_group;
            for k=1:size(rg,1)
                model.pi(modelind).g((model.pi(modelind).pi >= rg(k,1)) & (model.pi(modelind).pi<=rg(k,2))) = k;
            end
            
            [model.validation_km(modelind).KM] = km_stats(data.table(val_index & curfilt,:),model.pi(modelind).g(val_index & curfilt),size(rg,1),model);
            [model.validation_km(modelind).KM_pal] = km_stats(data.table(val_index & palfilt,:),model.pi(modelind).g(val_index & palfilt),size(rg,1),model);
            
            [model.validation_km(modelind).overall_cur] = km_stats(data.table(val_index & curfilt,:),ones(size(data.table(val_index & curfilt,:),1),1),1,model);
            [model.validation_km(modelind).overall_pal] = km_stats(data.table(val_index & palfilt,:),ones(size(data.table(val_index & palfilt,:),1),1),1,model);
            
%             f = find(model.pi(modelind).g'==2);
%             f = f(randperm(length(f),round(0.11*length(f))));
%             g = zeros(length(model.pi(modelind).g),1); g(f)=2; %|(g==2)
            val_cur = data.table(val_index & curfilt & ~((model.pi(modelind).g'==1)),{'censor','event_time'});
            val_cur2 = data.table(val_index & curfilt & ~(model.pi(modelind).g'==1),{'censor','event_time'});
            val_pal = data.table(val_index & palfilt & ~(((model.pi(modelind).g'==size(rg,1))|(model.pi(modelind).g'==2))),{'censor','event_time'});
            
%             N_pl = size(data.table(val_index & palfilt & ((model.pi(modelind).g'==size(rg,1))|(model.pi(modelind).g'==2)),:),1);
%             N_ch = size(data.table(val_index & curfilt & ((model.pi(modelind).g'==1)),:),1);
%             
%             
%             highrisk = model.pi(modelind).g(:)==1;
%             lowrisk = model.pi(modelind).g(:)==size(rg,1);
%             cur_low = data.data(~data.data.censor &~val_index & curfilt & lowrisk,:); 
%             pal_high = data.data(~data.data.censor & ~val_index & palfilt & highrisk,:); 
%             
%             gmm_cur_low = train_gmm1(cur_low, model);
%             if N_pl>0
%                 val_cur = [val_cur; sample_gmm(gmm_cur_low, N_pl)];
%             end
%             gmm_pal_high = train_gmm1(pal_high, model);
%             if N_ch>0
%                 val_pal = [val_pal; sample_gmm(gmm_pal_high, N_ch)];
%             end
            
            [model.validation_km(modelind).overall_cur_dss] = km_stats(val_cur,ones(size(val_cur,1),1),1,model);
            [model.validation_km(modelind).overall_pal_dss] = km_stats(val_pal,ones(size(val_pal,1),1),1,model);

            prop_cur = data.table(val_index & curfilt,{'censor','event_time'});
            prop_pal = data.table(val_index & palfilt,{'censor','event_time'});
            
            pmc = model.prop_match.curative(:,model.alg_client_index,i);
            pmp = model.prop_match.palliative(:,model.alg_client_index,i);
            
            binwidth = 0.05;        
            [cur_N,cur_lvl,cur_bin]=histcounts(model.pi(i,j).prop_pi(val_index & curfilt),0:binwidth:1); 
            [pal_N,pal_lvl,pal_bin]=histcounts(model.pi(i,j).prop_pi(val_index & palfilt),0:binwidth:1); 
            
            pal_match = zeros(1,sum(val_index & palfilt));
            cur_match = zeros(1,sum(val_index & curfilt));
            for m=1:length(pmp)
                x = find(pal_bin==m);
                if x
                    x = x(1:min([length(x) pmp(m)]));
                    pal_match(x)=1;
                end
                x = find(cur_bin==m);
                if x                
                    x = x(1:min([length(x) pmc(m)]));
                    cur_match(x)=1;
                end
            end
            %[sum(pal_match) sum(cur_match)]
            
            prop_val_cur = prop_cur(logical(cur_match),:);
            prop_val_pal = prop_pal(logical(pal_match),:);
            
            RG_cur = model.pi(modelind).g(val_index & curfilt);RG_cur = RG_cur(logical(cur_match));
            RG_pal = model.pi(modelind).g(val_index & palfilt);RG_pal = RG_pal(logical(pal_match));
            PI_cur = model.pi(modelind).pi(val_index & curfilt);PI_cur = PI_cur(logical(cur_match));
            PI_pal = model.pi(modelind).pi(val_index & palfilt);PI_pal = PI_pal(logical(pal_match));
            
            model.validation_hist(modelind).prop_cur = FLhistogram(PI_cur, 5);
            model.validation_hist(modelind).prop_pal = FLhistogram(PI_pal, 5);
            
            
            model.validation_hist(modelind).prop_cur_rg_high = FLhistogram(PI_cur(RG_cur==1), 2);
            model.validation_hist(modelind).prop_cur_rg_med = FLhistogram(PI_cur(RG_cur==2), 2);
            model.validation_hist(modelind).prop_cur_rg_low = FLhistogram(PI_cur(RG_cur==3), 2);
            model.validation_hist(modelind).prop_pal_rg_high = FLhistogram(PI_pal(RG_pal==1), 2);
            model.validation_hist(modelind).prop_pal_rg_med = FLhistogram(PI_pal(RG_pal==2), 2);
            model.validation_hist(modelind).prop_pal_rg_low = FLhistogram(PI_pal(RG_pal==3), 2);
            
            PI_median = model.settings.PI_median;
            model.validation_km(modelind).cur_high = km_stats(data.table(val_index & curfilt & (model.pi(modelind).pi>=PI_median)',:),ones(sum(val_index & curfilt & (model.pi(modelind).pi>=PI_median)'),1),1,model);
            model.validation_km(modelind).cur_low = km_stats(data.table(val_index & curfilt & (model.pi(modelind).pi<PI_median)',:),ones(sum(val_index & curfilt & (model.pi(modelind).pi<PI_median)'),1),1,model);
          
            
            %model.validation_km(modelind).prop_cur_dss = km_stats(prop_val_cur,ones(size(prop_val_cur,1),1),1,model);
            %model.validation_km(modelind).prop_pal_dss = km_stats(prop_val_pal,ones(size(prop_val_pal,1),1),1,model);

            model.validation_km(modelind).prop_cur = km_stats(prop_val_cur,RG_cur,size(rg,1),model);
            model.validation_km(modelind).prop_pal = km_stats(prop_val_pal,RG_pal,size(rg,1),model);
            
            model.validation_km(modelind).overall_prop_cur = km_stats(prop_val_cur,ones(size(prop_val_cur,1),1),1,model);
            model.validation_km(modelind).overall_prop_pal = km_stats(prop_val_pal,ones(size(prop_val_pal,1),1),1,model);

            
            prop_val_cur = prop_val_cur(~(RG_cur==1),:);
            prop_val_pal = prop_val_pal((RG_pal==1),:);
            
            model.validation_km(modelind).prop_cur_dss = km_stats(prop_val_cur,ones(size(prop_val_cur,1),1),1,model);
            model.validation_km(modelind).prop_pal_dss = km_stats(prop_val_pal,ones(size(prop_val_pal,1),1),1,model);

            
        end
    end
            
    
%     temp = g;
%     temp(ismember(g,model_binary_grouping{j}))=1;
%     temp(~ismember(g,model_binary_grouping{j}))=2;
%     g = temp;
%     [model.validation_result(i+1,j).KM_palliative] = km_stats(data.table(yearfilter & ~curfilter,:),g(yearfilter & ~curfilter),2,model);
    

    local_message = {'client_name','validation_km','validation_hist'};
    
    
elseif strcmp(model.state,'calibration')    
    
    
    if min(size(model.settings.varselect))==1; model.settings.varselect=model.settings.varselect(:)'; end
    Nboot = length(model.settings.randseed);
    numModel = size(model.settings.varselect,1);
    
    y = double(data.data.event_time>=surv_thresh);
    y(y==0)=-1; y=y(:)';
    
    curfilter = (data.data.(dose_feature) > curative_thresh);
    palfilter = (data.data.(dose_feature) <= curative_thresh);
    filter = ~(data.data.censor & (data.data.event_time<surv_thresh));
    
    for j=1:numModel
        w=model.w(((j-1)*Nboot+1):j*Nboot,:);
        if min(size(model.settings.varselect))==1
            varselect = logical(model.settings.varselect);
        else
            varselect = logical(model.settings.varselect(j,:));
        end

        for i=1:Nboot
            modelind =((j-1)*Nboot)+i;

            model.validation_result(i,j).confusionMatrix=zeros(2,2,length(model.settings.predict_thresh));
            model.validation_result(i,j).predicted_outcome_dist=zeros(length(model.settings.calibration_bins)-1,size(model.w,2));
            model.validation_result(i,j).actual_outcome_dist=zeros(length(model.settings.calibration_bins)-1,size(model.w,2));
             
            val_ind = false(size(data.data,1),1); 
            if strcmp(model.settings.validation,'time_bootstrap') || strcmp(model.settings.validation,'time')
                val_ind(data.test.index)=true;
            else
                val_ind(data.valid(i).index)=true;
            end 
            val_index = val_ind & filter & curfilter;
            val_index2 = val_ind & curfilter; % no censor
%             train_index = false(size(data.data,1),1);
%             train_index(data.train(i).index)=true;
%             train_index = train_index & filter;
%             val_index = train_index & filter;

            model.validation_result(modelind).predicted_outcome_dist = histcounts(model.pi(i,j).pi_cur,model.settings.calibration_bins);
            model.validation_result(modelind).actual_outcome_dist = histcounts(model.pi(i,j).pi_cur(y(val_index)==1),model.settings.calibration_bins);
            %risk group calibration - rows is risk groups and the columns
            %are the two values: predicted and actual.
            
            rg = unique(model.settings.PI_risk_group(:));
            model.validation_result(modelind).calib_rg_pred = histcounts(model.pi(i,j).pi_cur, rg);
            model.validation_result(modelind).calib_rg_act = histcounts(model.pi(i,j).pi_cur(y(val_index)==1), rg);

            for v = 1:length(rg)-1
                model.validation_result(modelind).calib_x(v) = median(model.pi(i,j).pi_cur((model.pi(i,j).pi_cur>rg(v)) & (model.pi(i,j).pi_cur<=rg(v+1))));
            end
        end
    end
            
    local_message = {'client_name','validation_result'};
       
    
    
    
elseif strcmp(model.state,'externvalidate')   

    
    Nboot = length(model.settings.randseed);
    y = double(data.data.survival>=surv_thresh);
    y(y==0)=-1; y=y(:)';
    filter = ~(data.data.vital_status & (data.data.survival<surv_thresh));
    yearfilter = (year(data.data.first_RT_date)>timespan(1)) & (year(data.data.first_RT_date)<timespan(2));
    filter = filter & yearfilter;
    curfilter=(data.data.(dose_feature) > curative_thresh);
    
    model_str = {'maastro2009','LCPI','stage'}; model_scaling=[1,28,3];
    model_cat = [4,4,3];
    model_binary_grouping = {[1,2],[1],[1,2]};
    
    for j = 1:3
        for i = 1:Nboot
            model.validate_compare(i,j).confusionMatrix=zeros(2,2,length(model.settings.predict_thresh));
            model.validate_compare(i,j).predicted_outcome_dist=zeros(length(model.settings.calibration_bins)-1,1);
            model.validate_compare(i,j).actual_outcome_dist=zeros(length(model.settings.calibration_bins)-1,1);
            
            val_index = false(size(data.data,1),1); val_index(data.valid(i).index)=true;
            val_index = val_index & filter;
            
            p=data.table.([model_str{j} '_score'])(val_index)/model_scaling(j);
            g=data.table.([model_str{j} '_group'])(val_index);
            
            [CM, acc, predicted_outcome_dist, actual_outcome_dist] = model_scoring(p, y(val_index), model.settings.predict_thresh,model.settings.calibration_bins);
                        
            model.validate_compare(i,j).confusionMatrix = CM;
            model.validate_compare(i,j).acc = acc;
            model.validate_compare(i,j).predicted_outcome_dist = predicted_outcome_dist;
            model.validate_compare(i,j).actual_outcome_dist = actual_outcome_dist;
            if length(unique(y(val_index)))>1
                [~,~,~,model.validate_compare(i,j).AUC]=perfcurve(y(val_index),p,1);
            else
                model.validate_compare(i,j).AUC = -99; % sentinel value to discount from analysis
            end
            model.validate_compare(i,j).CHarrell_prob=CHarrellIndex(data.table.survival(val_index),~data.table.vital_status(val_index),p);
            model.validate_compare(i,j).CHarrell_group=CHarrellIndex(data.table.survival(val_index),~data.table.vital_status(val_index),g);
%             [model.validate_compare(i,j).KM] = km_stats(data.table(val_index,:),g,model_cat(j),model);
        end
        model.validate_compare(i+1,j).AUC = -99;
        p=data.table.([model_str{j} '_score'])/model_scaling(j);
        g=data.table.([model_str{j} '_group']);
        [CM, acc, predicted_outcome_dist, actual_outcome_dist] = model_scoring(p(yearfilter), y(yearfilter), model.settings.predict_thresh,model.settings.calibration_bins);
        model.validate_compare(i+1,j).confusionMatrix = CM;
        model.validate_compare(i+1,j).acc = acc;
        model.validate_compare(i+1,j).predicted_outcome_dist = predicted_outcome_dist;
        model.validate_compare(i+1,j).actual_outcome_dist = actual_outcome_dist;
        if length(unique(y(yearfilter)))>1
            [~,~,~,model.validate_compare(i+1,j).AUC]=perfcurve(y(yearfilter),p(yearfilter),1);
        end

      
        model.validate_compare(i+1,j).CStat_prob = concordance_index(data.table.survival(yearfilter),~data.table.vital_status(yearfilter),p(yearfilter));
        model.validate_compare(i+1,j).CStat_group = concordance_index(data.table.survival(yearfilter),~data.table.vital_status(yearfilter),(model_cat(j)+1)-g(yearfilter));
        model.validate_compare(i+1,j).CStat_prob_c=concordance_index(data.table.survival(yearfilter& curfilter),~data.table.vital_status(yearfilter& curfilter),p(yearfilter& curfilter));
        model.validate_compare(i+1,j).CStat_group_c=concordance_index(data.table.survival(yearfilter& curfilter),~data.table.vital_status(yearfilter& curfilter),(model_cat(j)+1)-g(yearfilter & curfilter));
        model.validate_compare(i+1,j).CStat_prob_p=concordance_index(data.table.survival(yearfilter& ~curfilter),~data.table.vital_status(yearfilter&~curfilter),p(yearfilter&~curfilter));
        model.validate_compare(i+1,j).CStat_group_p=concordance_index(data.table.survival(yearfilter& ~curfilter),~data.table.vital_status(yearfilter&~curfilter),(model_cat(j)+1)-g(yearfilter & ~curfilter));

 
        [model.validate_compare(i+1,j).KM] = km_stats(data.table(yearfilter,:),g(yearfilter),model_cat(j),model);
        
        temp = g;
        temp(ismember(g,model_binary_grouping{j}))=1;
        temp(~ismember(g,model_binary_grouping{j}))=2;
        g = temp;
        [model.validate_compare(i+1,j).KM_curative] = km_stats(data.table(yearfilter & curfilter,:),g(yearfilter & curfilter),2,model);
        [model.validate_compare(i+1,j).KM_palliative] = km_stats(data.table(yearfilter & ~curfilter,:),g(yearfilter & ~curfilter),2,model);
        
        
        
    end
    local_message = {'client_name','validate_compare'};
    
    
   
  

else
    disp('Model state not found.')
end
