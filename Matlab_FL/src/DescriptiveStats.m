

%%% Report on data descriptive stats

format compact
stats_str = 'stats_cur';

for i=1:length(clients_names)
    indnew = find(ismember({model.client.client_name},clients_names{i}));
    model.client_new(i) = model.client(indnew);
end
model.client = model.client_new;

Yo = []; Ys = [];
for i = 1:length(model.client)
    Yo = [Yo model.client(i).stats_original.first_RT_date.date_hist];
    Ys = [Ys model.client(i).(stats_str).first_RT_date.date_hist];
end
figure('Color','w'); b = bar(model.client(1).(stats_str).first_RT_date.date_ranges(16:end-1),Ys(16:end-1,:),'stacked');
colors = distinguishable_colors(10); colors(4,:)=[];colors(4,:)=[];
for i=1:4;  b(i).FaceColor=colors(i,:);  end
xlabel('Years'); ylabel('Number of patients per clinic'); box off;
legend(clients_ids,'Location','northwest')

% Cur_Ys = Ys;
% Cur_pt_per_yr = pt_per_yr;

% Pal_Ys = Ys;
% Pal_pt_per_yr = pt_per_yr;

% 
% figure('Color','w'); b = bar(model.client(1).(stats_str).first_RT_date.date_ranges(10:end-1),Ys(10:end-1,:),'stacked');
% colors = distinguishable_colors(10); colors(4,:)=[];colors(4,:)=[];
% for i=1:4;  b(i).FaceColor=colors(i,:);  end
% xlabel('Years'); ylabel('Number of patients per clinic'); box off;
% legend(clients_ids,'Location','northwest')
% 
% Ys_ = Pal_Ys+Cur_Ys;
% ptyr = Pal_pt_per_yr+Cur_pt_per_yr;
% 
% figure('Color','w'); b = bar(model.client(1).(stats_str).first_RT_date.date_ranges(17:end-1),100*[Pal_pt_per_yr Cur_pt_per_yr]./ptyr,'stacked');
% colors = distinguishable_colors(10); colors(4,:)=[];colors(4,:)=[];
% b(1).FaceColor=colors(2,:); b(2).FaceColor=colors(1,:); 
% xlabel('Years'); ylabel('Number of patients per clinic'); box off;
% legend({'Palliative','Curative'},'Location','northwest')


%%%%%%%%%%%%%%%%
%%%Stats report
central_stats_report(model,stats_str)

 
stats_str1={'stats_cur','stats_pal'}; %,'stats_pal'

% Lost to follow up
disp('lost followup - incomplete follow up for 2 years')
for k=1:length(stats_str1)
    lost_follow = 0;
    disp(stats_str1{k})
    for i = 1:length(model.client)
        lost_follow = lost_follow+model.client(i).(stats_str1{k}).lost_followup;
    end
    disp(['Total lost follow up: ' num2str(lost_follow)])
end


%%%
category_var = {'histology','overall_stage','gender', 'smoking', 'n_stage', 't_stage','ecog','weightloss','laterality',...
    'tumour_loc','tumour_grade','marital_status','alcoholism','cardiac_comorbidity','technique','rt_intent','cause_of_death'};
%figure('Color','w');
pval_categorical = zeros(length(category_var),1);
for j=1:length(category_var)
    
    overall_num = zeros(size(model.client(1).(stats_str1{1}).(category_var{j}).num,1),length(stats_str1));
    Np = zeros(length(stats_str1),1);
    for k=1:length(stats_str1)
        disp(stats_str1{k})
        if isfield(model.client(i).(stats_str1{k}),category_var{j})
            disp(category_var{j})
            %subplot(4,5,j);
            stat=[]; 
            for i = 1:length(model.client)
                overall_num(:,k) = overall_num(:,k) + model.client(i).(stats_str1{k}).(category_var{j}).num;
                Np(k) = Np(k) + model.client(i).(stats_str1{k}).N;
            end
            disp(['Total missing: ' num2str(Np(k)-sum(overall_num(:,k))) ' (' num2str(100*((Np(k)-sum(overall_num(:,k)))/Np(k))) '%)'])
            disp(model.client(i).(stats_str1{1}).(category_var{j}).category')
            disp(overall_num(:,k)');
            disp([100*overall_num(:,k)'/Np(k)]);
            %             disp([model.client(i).client_name '-> Missing:' num2str(model.client(i).(stats_str1{k}).N-sum(model.client(i).(stats_str1{k}).(category_var{j}).num))...
            %                 ', N=' num2str(model.client(i).(stats_str1{k}).N) ' (' num2str(100*(sum(model.client(i).(stats_str1{k}).(category_var{j}).num))/model.client(i).(stats_str1{k}).N) '%)'])
            %             disp(model.client(i).(stats_str1{k}).(category_var{j}).category');
            %             disp(model.client(i).(stats_str1{k}).(category_var{j}).num');
            %             disp(100*model.client(i).(stats_str1{k}).(category_var{j}).num'/model.client(i).(stats_str1{k}).N);
            %             stat = [stat, [ model.client(i).(stats_str1{k}).(category_var{j}).num;model.client(i).(stats_str1{k}).N-sum(model.client(i).(stats_str1{k}).(category_var{j}).num)]];
            %
        end
        %bar(categorical([model.client(i).(stats_str).(category_var{j}).category;'missing']),stat)

    end
    
    if ismember(category_var{j},{'histology'}) && (length(Np)>1)
        E = Np(2)*(overall_num(1:end-1,1)/Np(1));
        O = overall_num(1:end-1,2);
    elseif ismember(category_var{j},{'ecog'})&& (length(Np)>1)
        
        overall_num(3,:) = sum(overall_num(3:end,:));
        overall_num(4:end,:)=[];
        E = Np(2)*(overall_num(:,1)/Np(1));
        O = overall_num(:,2);
    
% %     elseif ismember(category_var{j},{'cause_of_death'})
% %         E = (Np(2)-overall_num(end-1:end,2))*(overall_num(1:end-2,1)/(Np(1)-overall_num(end-1:end,1)));
% %         O = overall_num(1:end-2,2);
    else
        overall_num(overall_num(:,1)==0,:)=[];
        E = Np(2)*(overall_num(:,1)/Np(1));
        O = overall_num(:,2);
    end
    
    chisq=sum(((O-E).^2)./E);
    df=length(O)-1;
    pval_categorical(j) = chi2pdf(chisq,df);
    
end
pval_categorical
disp(['Categorical features not significantly different: ' category_var(pval_categorical>=0.001)])

dose_var = 'eqd2';
for k=1:length(stats_str1)
    stat = [];
    for i = 1:length(model.client)
        disp(stats_str1{k})
        disp(model.client(i).(stats_str1{k}).(dose_var).category');
        disp(model.client(i).(stats_str1{k}).(dose_var).num_category')
        stat = [stat, [model.client(i).(stats_str1{k}).(dose_var).num_category;model.client(i).(stats_str1{k}).N-sum(model.client(i).(stats_str1{k}).(dose_var).num_category)]];
    end
    figure('Color','w');bar(categorical([model.client(1).(stats_str).(dose_var).category;'missing']),stat)
    legend(clients_ids)
    sum(stat(:))
    disp(sum(stat,2))
    disp(sum(stat,2)/sum(stat(:)))
end

% continuous_var = {'fev', 'dlco', 'smoking_pack_yr', 'tumour_volume', 'age_start_rt',...
%     'presc_dose_ttl', 'survival', 'time_to_treat','time_to_cardiac'};
continuous_var = {'fev', 'dlco', 'smoking_pack_yr', 'tumour_volume', 'age_start_rt',...
    'presc_dose_ttl', 'eqd2', 'survival', 'time_to_treat'};
pval_continuous = zeros(length(continuous_var),1);
for j=1:length(continuous_var)
    disp(continuous_var{j})
    N=zeros(length(stats_str1),1);Np=zeros(size(N)); meansummation=zeros(size(N)); stdsum=zeros(size(N));
    for k=1:length(stats_str1)
        disp(stats_str1{k})
    
        if isfield(model.client(1).(stats_str1{k}),continuous_var{j})
            
             histlvl=0; cumhist=0;
            for i = 1:length(model.client)
                Np(k) = Np(k) + model.client(i).(stats_str1{k}).N;
                N(k) = N(k) + model.client(i).(stats_str1{k}).(continuous_var{j}).N;
                meansummation(k) = meansummation(k) + model.client(i).(stats_str1{k}).(continuous_var{j}).sum;
                stdsum(k) = stdsum(k) + model.client(i).(stats_str1{k}).(continuous_var{j}).stdsum;
                histlvl_c = model.client(i).(stats_str1{k}).(continuous_var{j}).CumHistLevels;
                cumhist_c = [0; diff(model.client(i).(stats_str1{k}).(continuous_var{j}).CumHist*model.client(i).(stats_str1{k}).N)];
                for m=1:length(histlvl_c)
                    if cumhist_c(m)>0
                        histlvl = sort([histlvl histlvl_c(m)]); ind=find(histlvl==histlvl_c(m),1,'first');
                        cumhist = [cumhist(1:ind-1) cumhist_c(m) cumhist(ind:end)];
                    end
                end
            end
            ind = find(double((cumsum(cumhist)/sum(cumhist))<=0.5),1,'last');
            disp(['Total values: ' num2str(N(k)) '. Missing->' num2str(Np(k)-N(k))]);
            disp(['Total mean: ' num2str(meansummation(k)/N(k))]);
            disp(['Total std: ' num2str(sqrt(stdsum(k)/(N(k)-1)))]);
            disp(['Total median: ' num2str(histlvl(ind))]);
            
        end
    end
    mu = meansummation./N; stdout = sqrt(stdsum./(N-1));
    DF = (((stdout(1)^2)/N(1) + (stdout(2)^2)/N(2))^2) / (((stdout(1)^2)/N(1))^2/(N(1)-1) + ((stdout(2)^2)/N(2))^2/(N(2)-1));
    v=round(DF);
    stdout_total = sqrt(((N(1)-1)*stdout(1)^2 + (N(2)-1)*stdout(2)^2)/(N(1)+N(2)-2));
    tstat = abs((mu(1) - mu(2)))/(stdout_total*sqrt((1/N(1) + 1/N(2))));
    %pval = betainc(v/(v+tstat^2),v/2,0.5)
    pval_continuous(j) = 1-tcdf(tstat,v);
end

disp(['Continuous features not significantly different: ' continuous_var(pval_continuous>=0.025)])

% for k=1:length(stats_str1)
%     disp(stats_str1{k})
%     
%     if isfield(model.client(i).(stats_str1{k}),'time_to_cardiac')
%         disp('Cardiac events:')
%         for i = 1:length(model.client)
%             disp([model.client(i).client_name '-> Missing:' num2str(model.client(i).(stats_str1{k}).N-sum(model.client(i).(stats_str1{k}).time_to_cardiac.N)) ', N=' num2str(model.client(i).(stats_str1{k}).N)...
%                 ' (' num2str(100*(model.client(i).(stats_str1{k}).N-model.client(i).(stats_str1{k}).time_to_cardiac.N)/model.client(i).(stats_str1{k}).N) '%)'])
%         end
%     end
% end

disp('Survival2yr:')
for k=1:length(stats_str1)
    disp(stats_str1{k})
    N=0;Np=0; histlvl=0; cumhist=0; meansummation = 0; stdsum=0;num_cens=0; num=0;
    for i = 1:length(model.client)
        disp(model.client(i).(stats_str1{k}).survival.surv2yr_cens)
        Np = Np + model.client(i).(stats_str1{k}).N;
        N = N + model.client(i).(stats_str1{k}).survival.N;
        num = num + model.client(i).(stats_str1{k}).survival.surv2yr;
        num_cens = num_cens + model.client(i).(stats_str1{k}).survival.surv2yr_cens;
        meansummation = meansummation + model.client(i).(stats_str1{k}).survival.sum;
        stdsum = stdsum + model.client(i).(stats_str1{k}).survival.stdsum;
        histlvl_c = model.client(i).(stats_str1{k}).survival.CumHistLevels;
        cumhist_c = [0; diff(model.client(i).(stats_str1{k}).survival.CumHist*model.client(i).(stats_str1{k}).N)];
        for j=1:length(histlvl_c)
            if cumhist_c(j)>0
            histlvl = sort([histlvl histlvl_c(j)]); ind=find(histlvl==histlvl_c(j),1,'first');
            cumhist = [cumhist(1:ind-1) cumhist_c(j) cumhist(ind:end)];
            end
        end
    end
    ind = find(double((cumsum(cumhist)/sum(cumhist))<=0.5),1,'last');
    disp(['Total proportion: Yes->' num2str(num) '. No->' num2str(N-num) '. Censored->' num2str(num_cens)]);
    disp(['Total mean: ' num2str(meansummation/N)]);
    disp(['Total std: ' num2str(sqrt(stdsum/(N-1)))]);
    disp(['Total median: ' num2str(histlvl(ind))]);
end


category_var = {'histology','overall_stage','gender', 'smoking', 'n_stage', 't_stage','ecog','weightloss','laterality',...
    'tumour_loc','tumour_grade','marital_status','alcoholism','cardiac_comorbidity','technique','rt_intent','cause_of_death'};
vari = {'overall_stage','gender','laterality','tumour_loc','ecog'};
vari = {'histology'};

stats_str = 'stats_cur';
for k = 1:length(vari)
    figure('Color','w','Position',[10 10 350 250]);
    stat=[];
    for j = 1:length(model.client)
        stat = [stat, [ model.client(j).(stats_str).(vari{k}).num;...
            model.client(j).(stats_str).N-sum(model.client(j).(stats_str).(vari{k}).num)]];
    end
    bar(categorical([model.client(j).(stats_str).(vari{k}).category;'missing']),stat)
    %legend(clients_ids,'Location','northwest');
end
 legend(clients_ids,'Location','northwest');

 
continuous_var = {'fev', 'dlco', 'smoking_pack_yr', 'tumour_volume', 'age_start_rt',...
    'presc_dose_ttl', 'survival', 'time_to_treat'};

vari = {'age_start_rt','survival','time_to_treat','tumour_volume'};
vari_str = {'Age at start of RT (years)','Survival (months)','Time from diagnosis to RT (days)','Tumour volume (cc)'};
stats_str = 'stats_pal';

for k = 1:length(vari)
    figure('Color','w','Position',[10 10 300 200]); hold on;
    for j = 1:length(model.client)
        x=model.client(j).(stats_str).(vari{k}).CumHistLevels;
        Y=model.client(j).(stats_str).(vari{k}).CumHist;
        plot(x, Y,'LineWidth',2.5) %yi is cdf 
    %     xi = 0:(max(x)/500):max(x);
    %     yi = interp1(x,Y,xi,'linear');
    %     PDF_on=diff([0 yi]);
    end
    %legend(clients_ids,'Location','southeast');
    %title(['Cumulative distribution'])
    xlabel(vari_str{k})
end
legend(clients_ids,'Location','southeast');