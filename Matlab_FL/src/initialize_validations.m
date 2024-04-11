function [model, data] = initialize_validations(model, data, type, stratifying_index)

[N, ~] = size(data.table);

if strcmp(type,'crossval')
    
    data.val_index = zeros(N,1);
    if isempty(stratifying_index)
        stratifying_index = ones(N,1);
    end
    if min(stratifying_index)==0; stratifying_index=stratifying_index+1; end
    index1=[];index2=[];
        % partition data while stratifying for propensity of outcome
    for j=1:length(unique(stratifying_index))
        num_pos = length(find(y(stratifying_index==j)==1));
        num_neg = length(find(y(stratifying_index==j)==-1));

        temp = [repmat(1:model.cross_val,[1 floor(num_pos/model.cross_val)])'; [1:(num_pos-floor(num_pos/model.cross_val)*model.cross_val)]'];
        if length(temp)<num_pos; temp = [temp;ones(num_pos-length(temp),1)]; end
        index1 = [index1; find(stratifying_index==j & y==1), temp];
        
        temp = [repmat(1:model.cross_val,[1 floor(num_neg/model.cross_val)])'; [1:(num_neg-floor(num_neg/model.cross_val)*model.cross_val)]'];
        if length(temp)<num_neg; temp = [temp;ones(num_neg-length(temp),1)]; end
        index2 = [index2; find(stratifying_index==j & y==-1), temp];
    end

    index1(:,2) = index1(randperm(size(index1,1)),2);
    index2(:,2) = index2(randperm(size(index2,1)),2);
    data.val_index(index1(:,1))=index1(:,2);
    data.val_index(index2(:,1))=index2(:,2);
    
%     ind1 = [find(y==1), [repmat(1:model.cross_validation,[1 floor(num_pos/model.cross_validation)])'; [1:(num_pos-floor(num_pos/model.cross_validation)*model.cross_validation)]']];
%     ind2 = [find(y==-1), [repmat(1:model.cross_validation,[1 floor(num_neg/model.cross_validation)])'; [1:(num_neg-floor(num_neg/model.cross_validation)*model.cross_validation)]']];
%     %randomize the entries
%     ind1(:,2) = ind1(randperm(size(ind1,1)),2);
%     ind2(:,2) = ind2(randperm(size(ind2,1)),2);
%     data.val_index(ind1(:,1))=ind1(:,2);
%     data.val_index(ind2(:,1))=ind2(:,2);
    
elseif strcmp(type,'rrss')
       
    
    rng(model.randseed(1));
    [N, ~] = size(data.table);
    y = double(data.table.survival>=model.survival_threshold);
    y(y==0)=-1;
    data.val_index = zeros(N,1); % partition data while stratifying for propensity of outcome
    
    if isempty(stratifying_index)
        stratifying_index = ones(N,1);
    end
    if min(stratifying_index)==0; stratifying_index=stratifying_index+1; end
    index1=[];index2=[];
    for j=1:length(unique(stratifying_index))
        
        num_pos = length(find(y(stratifying_index==j)==1));
        num_neg = length(find(y(stratifying_index==j)==-1));
        
        temp = [ones(floor(num_pos*model.rrss_prop(1)),1); 2*ones(floor(num_pos*model.rrss_prop(2)),1)];
        if length(temp)<num_pos; temp = [temp;ones(num_pos-length(temp),1)]; end
        index1 = [index1; find(stratifying_index==j & y==1), temp];
        
        temp = [ones(floor(num_neg*model.rrss_prop(1)),1); 2*ones(floor(num_neg*model.rrss_prop(2)),1)];
        if length(temp)<num_neg; temp = [temp;ones(num_neg-length(temp),1)]; end
        index2 = [index2; find(stratifying_index==j & y==-1), temp];
    end
    
    index1(:,2) = index1(randperm(size(index1,1)),2);
    index2(:,2) = index2(randperm(size(index2,1)),2);
    data.val_index(index1(:,1))=index1(:,2);
    data.val_index(index2(:,1))=index2(:,2);

    data.train.index=data.val_index==1;
    data.valid.index=data.val_index==2;

elseif strcmp(type,'time_rrss')
    
    rng(model.randseed(1));
    data.data = data.data((year(data.data.first_RT_date)>model.settings.timespan(1)) & (year(data.data.first_RT_date)<model.settings.timespan(2)),:);

    stratifying_index = stratifying_index((year(data.data.first_RT_date)>model.settings.timespan(1)) & (year(data.data.first_RT_date)<model.settings.timespan(2)));
    [N, ~] = size(data.data);
    y = double(data.data.survival>=model.settings.survival_threshold);
    y(y==0)=-1;
    data.val_index = zeros(N,1); % partition data while stratifying for propensity of outcome
    
    if isempty(stratifying_index)
        stratifying_index = ones(N,1);
    end
    if min(stratifying_index)==0; stratifying_index=stratifying_index+1; end
    index1=[];index2=[];
    for j=1:length(unique(stratifying_index))
        
        num_pos = length(find(y(stratifying_index==j)==1));
        num_neg = length(find(y(stratifying_index==j)==-1));
        
        temp = [ones(floor(num_pos*model.settings.rrss_prop(1)),1); 2*ones(floor(num_pos*model.settings.rrss_prop(2)),1)];
        if length(temp)<num_pos; temp = [temp;ones(num_pos-length(temp),1)]; end
        index1 = [index1; find(stratifying_index==j & y==1), temp];
        
        temp = [ones(floor(num_neg*model.settings.rrss_prop(1)),1); 2*ones(floor(num_neg*model.settings.rrss_prop(2)),1)];
        if length(temp)<num_neg; temp = [temp;ones(num_neg-length(temp),1)]; end
        index2 = [index2; find(stratifying_index==j & y==-1), temp];
    end
    
    index1(:,2) = index1(randperm(size(index1,1)),2);
    index2(:,2) = index2(randperm(size(index2,1)),2);
    data.val_index(index1(:,1))=index1(:,2);
    data.val_index(index2(:,1))=index2(:,2);

    data.train.index=data.val_index==1;
    data.valid.index=data.val_index==2;
    data.test.index=(year(data.data.first_RT_date)>=model.settings.timespan(2)) & (year(data.data.first_RT_date)<model.settings.timespan(3));

    
elseif strcmp(type,'time_bootstrap')
     
    D = data.table((year(data.table.first_RT_date)>model.settings.timespan(1)) & (year(data.table.first_RT_date)<model.settings.timespan(2)),:);
    D_index = find((year(data.table.first_RT_date)>model.settings.timespan(1)) & (year(data.table.first_RT_date)<model.settings.timespan(2)));
    [N, ~] = size(D);
    Nboot = length(model.settings.randseed);
    
    for i=1:Nboot
        
        rng(model.settings.randseed(i)); % + model.localseed
        if model.settings.randseed(i)==0
            train_index=(1:N);
            valid_index=(1:N);
        else
            valid_index=[];
            while isempty(valid_index)
                train_index=unique(ceil(N*rand(1,N)));
                valid_index=setdiff((1:N),train_index);
            end
        end
        data.train(i).index = false(size(data.table,1),1);
        data.valid(i).index = false(size(data.table,1),1);
        data.train(i).index(D_index(train_index)) = true;
        data.valid(i).index(D_index(valid_index)) = true;
    end    
    data.test.index=(year(data.table.first_RT_date)>=model.settings.timespan(2)) & (year(data.table.first_RT_date)<model.settings.timespan(3));

    
elseif strcmp(type, 'time')
    
    data.val_index = zeros(N,1); % partition data on time (data may not be stratified, but may use propensity matching).
    % timespan has three entries (1) -> (2) is training set. (2) inclusive -> (3) is testing.
    
    data.train.index=(year(data.table.first_RT_date)>model.settings.timespan(1)) & (year(data.table.first_RT_date)<model.settings.timespan(2));
    data.valid.index=data.train.index;
%     data.valid.index=(year(data.table.first_RT_date)>=model.settings.timespan(2)) & (year(data.table.first_RT_date)<model.timespan(3));
    data.test.index=(year(data.table.first_RT_date)>=model.settings.timespan(2)) & (year(data.table.first_RT_date)<model.settings.timespan(3));
    
elseif strcmp(type, 'bootstrap')
    
    D = data.table((year(data.table.first_RT_date)>model.settings.timespan(1)) & (year(data.table.first_RT_date)<model.settings.timespan(2)),:);
    D_index = find((year(data.table.first_RT_date)>model.settings.timespan(1)) & (year(data.table.first_RT_date)<model.settings.timespan(2)));
    [N, ~] = size(D);
    Nboot = length(model.settings.randseed);
    
    for i=1:Nboot
        
        rng(model.settings.randseed(i));
        if model.settings.randseed(i)==0
            train_index=(1:N);
            valid_index=(1:N);
        else
            valid_index=[];
            while isempty(valid_index)
                train_index=(ceil(N*rand(1,N)));
                valid_index=setdiff((1:N),train_index);
            end
        end
        data.train(i).index = false(size(data.table,1),1);
        data.valid(i).index = false(size(data.table,1),1);
        data.train(i).index(D_index(train_index)) = true;
        data.valid(i).index(D_index(valid_index)) = true;

    end
    
elseif strcmp(type, 'bootstrap_strat')
    
    Nboot = length(model.settings.randseed);
    yearfilter = (year(data.table.first_RT_date)>model.settings.timespan(1)) & (year(data.table.first_RT_date)<model.settings.timespan(2));

    D = data.table(yearfilter & (data.table.(model.settings.dose_feature)>model.settings.curative_dose_threshold),:);
    D_index = find(yearfilter & (data.table.(model.settings.dose_feature)>model.settings.curative_dose_threshold));
    [N, ~] = size(D);
    
    for i=1:Nboot
        
        rng(model.settings.randseed(i));
        valid_index=[];
        while isempty(valid_index)
            train_index=ceil(N*rand(1,N));
            valid_index=setdiff((1:N),train_index);
        end
        data.train(i).index=D_index(train_index);
        data.valid(i).index=D_index(valid_index);
    end
    

end