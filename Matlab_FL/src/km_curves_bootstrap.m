function [logrank_stat,p_value, pv]=km_curves_bootstrap(KM, labels, model)

logrank_stat=0;
p_value=0;
ti = (0:(1/model.settings.km_periods):model.settings.km_limit);
ti = ti(1:end-1)';

surv = zeros(length(ti),size(KM(1).D,2)); 

num_curve = length(KM);
for i=1:num_curve
    for j=1:size(KM(i).D,2)
        S=1; se=[]; 
        time_inc = (0:(1/model.settings.km_periods):model.settings.km_limit);
        time_inc = time_inc(1:end-1)';
        if any(KM(i).D(:,j))       
            time_inc = time_inc(KM(i).D(:,j)>0);
            km_N_i = KM(i).N(KM(i).D(:,j)>0,j);
            km_D_i = KM(i).D(KM(i).D(:,j)>0,j);
            S = cumprod(1 - km_D_i./km_N_i);
            se = S .* sqrt(cumsum(km_D_i ./ (km_N_i .* (km_N_i-km_D_i))));
            se(isnan(se))=0;
        else
            time_inc=0;
        end
        
        c1=S.*(numel(KM(i).N(KM(i).D(:,j)>0,j))+numel(KM(i).N(KM(i).D(:,j)>0,j)));
        if isempty(c1)
            lambda = 0;
        else
            c2=-(diff(log(c1(1:end-1)))./diff(time_inc(1:end-1)));
            lambda=mean(c2(c2~=0));
        end
        if length(time_inc)<2
            Si = ones(size(ti));
%             Si = [interp1([0;time_inc],[1;S],ti,'next','extrap') interp1([0;time_inc],[1;S],ti,'previous','extrap')];
        else
            Si = [interp1(time_inc,S,ti,'next','extrap') interp1(time_inc,S,ti,'previous','extrap')];
        end
        Si2 = zeros(length(Si),1); 
        for t=1:length(Si)
            Si2(t) = Si(t,find(~isnan(Si(t,:)),1,'first'));
        end
        
        surv(:,j)=Si2;
        
        %stairs(time_inc,S)
        %stairs(ti,Si)
        %plot(time_inc,S)
%         surv(i).S
%         
%         surv(i).S = prctile(S,[2.5,50,97.5]); 
%         surv(i).se = prctile(se,[2.5,50,97.5]); 
%         surv(i).time_inc = time_inc;
%         surv(i).ti = (0:(1/model.settings.km_periods):model.settings.km_limit); surv_struct(i,j).ti = surv_struct(i,j).ti(1:end-1)';
%         surv(i).lambda = lambda;
        
        
    end
    
    ci_prctiles = [2.5 50 97.5];%[2.5 50 97.5]
    surv_struct(i).D = prctile(KM(i).D',ci_prctiles)';
    surv_struct(i).N = prctile(KM(i).N',ci_prctiles)';
    surv_struct(i).S = prctile(surv',ci_prctiles)';
    surv_struct(i).num = prctile(KM(i).num_pts',ci_prctiles)';
end

lim = length(ti);
logrank_stat=zeros(num_curve);
p_value=zeros(num_curve); pv = zeros(num_curve); 
perc5yr=zeros(num_curve);perc2yr=zeros(num_curve);
for j=1:num_curve
    for k=1:num_curve
        [logrank_stat(j,k), p_value(j,k)] = KM_logrank(surv_struct(j).N(1:lim,2),surv_struct(j).D(1:lim,2),...
                                                        surv_struct(k).N(1:lim,2),surv_struct(k).D(1:lim,2),0.05,...
                                                        [num2str(j) ' ' num2str(k)], [num2str(j) ' ' num2str(k)]);
        
        perc2yr(k) = surv_struct(k).S(abs(ti-2)==min(abs(ti-2)),2);
        perc2yr(j) = surv_struct(j).S(abs(ti-2)==min(abs(ti-2)),2);  
                                          
        p = calculatePvalueProportion(surv_struct(k).num(2), perc2yr(k),...
                                      surv_struct(j).num(2), perc2yr(j));
        pv(j,k) = p;
        disp(['For ' num2str(j) ',' num2str(k) ' group. 2-year survival proportions: p=' num2str(p)])
        OddsRatio = ((1-perc2yr(k))/perc2yr(k))/((1-perc2yr(j))/perc2yr(j));
        disp(['Odds Ratio 2-years -> ' num2str((1-perc2yr(k))/perc2yr(k)) ':' num2str((1-perc2yr(j))/perc2yr(j)) ' = ' num2str(OddsRatio)])


        perc5yr(k) = surv_struct(k).S(abs(ti-5)==min(abs(ti-5)),2);
        perc5yr(j) = surv_struct(j).S(abs(ti-5)==min(abs(ti-5)),2);

        p5 = calculatePvalueProportion(surv_struct(k).num(2), perc5yr(k),...
                                       surv_struct(j).num(2), perc5yr(j));
        disp(['For ' num2str(j) ',' num2str(k) ' group. 5-year survival proportions: p=' num2str(p5)])
        OddsRatio = ((1-perc5yr(k))/perc5yr(k))/((1-perc5yr(j))/perc5yr(j));
        disp(['Odds Ratio 5-years -> ' num2str((1-perc5yr(k))/perc5yr(k)) ':' num2str((1-perc5yr(j))/perc5yr(j)) ' = ' num2str(OddsRatio)])

        
    end
    
end


%p = calculatePvalueProportion(423, .512, 431, .592);

th = model.settings.survival_threshold/12;

f=figure('Color','w'); ax = axes; hold on;
col = [0.6350 0.0780 0.1840;
    0 0.4470 0.7410;
    0.4660 0.6740 0.1880];
% 
% col = [0 1 0;
%     1 0 0;
%     0.4660 0.6740 0.1880;
%     0.6350 0.0780 0.1840];

col = [0 1 0;
    0 0 1;
    1 0 0;
    0.4660 0.6740 0.1880;
    0 0.4470 0.7410;
    0.6350 0.0780 0.1840];


% if length(labels)>4
%     col = colormap(jet(length(labels)));
% end

alpha = 0.05;
maxkm=[]; minkm=[];
% col = 'gbrkmc'; alpha = 0.05;
for i=1:length(surv_struct)
    if ~isempty(surv_struct(i).S)
        
        if model.settings.km_conf
%             h=stairs(ti,[surv_struct(i).S(:,2) surv_struct(i).S(:,1) surv_struct(i).S(:,3)]);
            h=stairs(ti,surv_struct(i).S(:,2));
%             plot(ti,surv_struct(i).S(:,2),'.k','MarkerSize',20,'HandleVisibility','off')
            plot(ti(abs(ti-2)==min(abs(ti-2))),surv_struct(i).S(abs(ti-2)==min(abs(ti-2)),2),'.k','MarkerSize',15,'HandleVisibility','off')
            patch([ti; flipud(ti)], [surv_struct(i).S(:,1); flipud(surv_struct(i).S(:,3))], col(i,:), 'FaceAlpha',0.2, 'EdgeColor','none','HandleVisibility','off')

        else
            h=stairs(ti,surv_struct(i).S(:,2));
        end
        [~,j]=min(abs(ti-th));
        grid on;
        
        %text(th,surv_struct(i).S(j),[num2str(100*surv_struct(i).S(j),'%0.0f') '%'],'Color',col(i),'Fontsize',12)
        %plot(th,surv_struct(i).S(j),'r.','MarkerSize',20,'HandleVisibility','off')
        set(h(1),'Color',col(i,:),'Linewidth',3);
        
        disp(['2-year survival percentage for ' num2str(i) ' curve: ' num2str(surv_struct(i).S(abs(ti-2)==min(abs(ti-2)),2))])
        disp(['5-year survival percentage for ' num2str(i) ' curve: ' num2str(surv_struct(i).S(abs(ti-5)==min(abs(ti-5)),2))])
%         if model.settings.km_conf
%             set(h(2),'Color',0.6*ones(1,3),'HandleVisibility','off'); set(h(3),'Color',0.6*ones(1,3),'HandleVisibility','off');
%         end
        set(ax,'XLim', [0 model.settings.km_limit], 'YLim',[0 1]);
    end
    maxkm = max([maxkm surv_struct(i).S(:,1)' surv_struct(i).S(:,3)']); minkm = min([minkm surv_struct(i).S(:,1)' surv_struct(i).S(:,3)']);
end

if exist('h')
    
    if all(contains(class(h),'Stair'))
        for j=1:length(surv_struct)
            labels(j) = {[labels{j} ', $N_{\mu}$=' num2str(round(surv_struct(j).num(2))) ...
                ' [' num2str(round(surv_struct(j).num(1))) '-' num2str(round(surv_struct(j).num(3))) ']']};%surv_struct(j).km_pts
            %         labels(j) = {[labels{j} ', N=' num2str(surv_struct(j).km_pts)]};
        end
        legend(labels, 'Interpreter', 'latex');
        %     legend(labels);
        xlabel('Time from start of radiotherapy (years)');
        %ylabel('Cumulative probability of event'); %ylabel('Proportion of survival');
        ylabel('Proportion survival');
%         ylim([max([minkm-.02 0]) min([maxkm+.02 1])])
        ylim([0 1])
    end
    
end

