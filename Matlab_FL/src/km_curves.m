function [logrank_stat,p_value]=km_curves(KM, labels, model)


num_curve = length(KM);
for i=1:num_curve
    for j=1:size(KM(i).D,2)
        S=[]; se=[]; 
        time_inc = (0:(1/model.settings.km_periods):model.settings.km_limit);
        time_inc = time_inc(1:end-1)';
        if any(KM(i).D(:,j))       
            time_inc = time_inc(KM(i).D(:,j)>0);
            km_N_i = KM(i).N(KM(i).D(:,j)>0,j);
            km_D_i = KM(i).D(KM(i).D(:,j)>0,j);
            
            S = cumprod(1 - km_D_i./km_N_i);
            [mi mi_ind]=min(abs(time_inc-1))
            disp(S(mi_ind))
            se = S .* sqrt(cumsum(km_D_i ./ (km_N_i .* (km_N_i-km_D_i))));
            se(isnan(se))=0;
        end
        
        c1=S.*(numel(KM(i).N(KM(i).D(:,j)>0,j))+numel(KM(i).N(KM(i).D(:,j)>0,j)));
        if isempty(c1)
            lambda = 0;
        else
            c2=-(diff(log(c1(1:end-1)))./diff(time_inc(1:end-1)));
            lambda=mean(c2(c2~=0));
        end
        surv_struct(i,j).S = S; surv_struct(i,j).se = se; surv_struct(i,j).time_inc = time_inc;
        surv_struct(i,j).N = KM(i).N(:,j);
        surv_struct(i,j).D = KM(i).D(:,j);
        surv_struct(i,j).km_pts = KM(i).num_pts(j);
        surv_struct(i,j).ti = (0:(1/model.settings.km_periods):model.settings.km_limit);
        surv_struct(i,j).ti = surv_struct(i,j).ti(1:end-1)';
        surv_struct(i,j).lambda = lambda;
        
        
    end
end
    



logrank_stat=zeros(num_curve);
p_value=zeros(num_curve);
for j=1:num_curve
    for k=1:num_curve
        [logrank_stat(j,k), p_value(j,k)] = KM_logrank(surv_struct(j).N,surv_struct(j).D,...
                                                        surv_struct(k).N,surv_struct(k).D,0.05,...
                                                        [num2str(j) ' ' num2str(k)], [num2str(j) ' ' num2str(k)]);
    end
end



th = model.settings.survival_threshold/12;

f=figure('Color','w'); ax = axes; hold on;
col = 'rbgkmc'; alpha = 0.05; 
maxkm=[]; minkm=[];
% col = 'gbrkmc'; alpha = 0.05; 
for i=1:length(surv_struct)
    if ~isempty(surv_struct(i).S)
        zalpha = -norminv(alpha/2);
        halfwidth = zalpha*surv_struct(i).se;
        Flo = max(0, surv_struct(i).S-halfwidth);
        Flo(isnan(halfwidth)) = NaN; % max drops NaNs, put them back
        Fup = min(1, surv_struct(i).S+halfwidth);
        Fup(isnan(halfwidth)) = NaN; % max drops NaNs

        if ~(any(any(isnan([surv_struct(i).S Flo Fup]))) || length(surv_struct(i).S)==1)

            if model.settings.km_conf
                h=stairs(surv_struct(i).time_inc,[surv_struct(i).S Flo Fup]);
            else
                h=stairs(surv_struct(i).time_inc,[surv_struct(i).S]);
            end
            [~,j]=min(abs(surv_struct(i).time_inc-th)); 
            %text(th,surv_struct(i).S(j),[num2str(100*surv_struct(i).S(j),'%0.0f') '%'],'Color',col(i),'Fontsize',12)
            %plot(th,surv_struct(i).S(j),'r.','MarkerSize',20,'HandleVisibility','off')
            set(h(1),'Color',col(i),'Linewidth',3);
            if model.settings.km_conf
                set(h(2),'Color',0.6*ones(1,3),'HandleVisibility','off'); set(h(3),'Color',0.6*ones(1,3),'HandleVisibility','off');
            end
            set(ax,'XLim', [0 model.settings.km_limit], 'YLim',[0 1]);
        end
        maxkm = max([maxkm Flo' Fup']); minkm = min([minkm Flo' Fup']);
    end
end

if exist('h')

if all(contains(class(h),'Stair'))
    for j=1:length(surv_struct)
        labels(j) = {[labels{j} ', $N_{\mu}$=' num2str(surv_struct(j).km_pts)]};
%         labels(j) = {[labels{j} ', N=' num2str(surv_struct(j).km_pts)]};
    end
    legend(labels, 'Interpreter', 'latex');
%     legend(labels);
    xlabel('Time from start of radiotherapy (years)'); 
    ylabel('Cumulative probability of event'); %ylabel('Proportion of survival');
    ylim([max([minkm-.02 0]) min([maxkm+.02 1])])
end

end