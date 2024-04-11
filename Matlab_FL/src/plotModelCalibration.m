function [f, b] = plotModelCalibration2(counts,calibBins)%posCount,allCount

sizeshrink = 0.005;


obs_prob_ci = squeeze(prctile(counts(:,1,:)./counts(:,2,:),[10,50,90],3));
cb = calibBins; 
cb(isnan(obs_prob_ci(:,2)))=[];
obs_prob_ci(any(isnan(obs_prob_ci), 2), :) = [];
medi_obs = obs_prob_ci(obs_prob_ci(:,2)>0,2);
medi_obs_cb = cb(obs_prob_ci(:,2)>0);

[b,bint,r,rint,stats] = regress(obs_prob_ci(:,2),[ones(length(cb),1) cb(:)]);
brier_score = sum((cb(:)-(obs_prob_ci(:,2))).^2)/length(cb);

% [b,bint,r,rint,stats] = regress(medi_obs,[ones(length(medi_obs_cb),1) medi_obs_cb(:)]);



f = figure('Color','w');
sb1 = subplot(2,1,2,'Position',[0.1 0.01 0.8 0.18]);
bar(calibBins+0.00,-sizeshrink*(median(counts(:,2,:),3,'omitnan')-median(counts(:,1,:),3,'omitnan')), 'barwidth', 0.01,'FaceColor','k'); %
hold on;
bar(calibBins,sizeshrink*(median(counts(:,1,:),3,'omitnan')), 'barwidth', 0.01,'FaceColor','k'); %
axis off;
set(gca,'xlim',[0 1])

sb2 = subplot(2,1,1,'Position',[0.1 0.275 0.8 0.7]);
plot(cb,obs_prob_ci(:,2),'.b','MarkerSize',20)
hold on;
plot([0 1],[0 1],'-k')

patch([cb'; flipud(cb')], [obs_prob_ci(:,1); flipud(obs_prob_ci(:,3))], 'b', 'FaceAlpha',0.2, 'EdgeColor','none','HandleVisibility','on')

box off;
xlabel('Predicted probability')
ylabel('Observed frequency')
set(gca,'xlim',[0 1],'ylim',[-0.05 1])
% annotation('textbox',[.12 .67 .3 .3],'String',{['Slope: ' num2str(b(2),'%.3f')],['Intercept: ' num2str(b(1),'%.3f')],['R^2: ' num2str(stats(1),'%.3f')], ['Brier: ' num2str(brier_score,'%.3f')]},'FitBoxToText','on');





obs_prob = mean(counts(:,1,:)./counts(:,2,:),3,'omitnan'); 
obs_prob_std = std(counts(:,1,:)./counts(:,2,:),[],3,'omitnan'); 

obs_prob_std(isnan(obs_prob)) = [];
cb = calibBins; cb(isnan(obs_prob))=[];
obs_prob(isnan(obs_prob)) = [];

[b,bint,r,rint,stats] = regress(obs_prob,[ones(length(cb),1) cb(:)]);
brier_score = sum((cb(:)-(obs_prob)).^2)/length(cb);

f = figure('Color','w');
sb1 = subplot(2,1,2,'Position',[0.1 0.01 0.8 0.18]); %[0.1 0.01 0.8 0.18]
bar(calibBins+0.00,-sizeshrink*(mean(counts(:,2,:),3,'omitnan')-mean(counts(:,1,:),3,'omitnan')), 'barwidth', 0.01,'FaceColor','k'); %
hold on;
bar(calibBins,sizeshrink*(mean(counts(:,1,:),3,'omitnan')), 'barwidth', 0.01,'FaceColor','k'); %
axis off;
set(gca,'xlim',[0 1])

sb2 = subplot(2,1,1,'Position',[0.1 0.275 0.8 0.7]);
plot(cb,obs_prob,'.b','MarkerSize',20, 'handlevisibility','off')
hold on;
%plot([0 1],[0 1],'-k')
% plot([0:0.01:1],b(1)+b(2)*[0:0.01:1],'--k')

patch([cb'; flipud(cb')], [obs_prob+obs_prob_std; flipud(obs_prob-obs_prob_std)], 'b', 'FaceAlpha',0.2, 'EdgeColor','none','HandleVisibility','off')


box off;
xlabel('Predicted probability')
ylabel('Observed frequency')
set(gca,'xlim',[0 1])
%annotation('textbox',[.12 .67 .3 .3],'String',{['Slope: ' num2str(b(2),'%.3f')],['Intercept: ' num2str(b(1),'%.3f')],['R^2: ' num2str(stats(1),'%.3f')], ['Brier: ' num2str(brier_score,'%.3f')]},'FitBoxToText','on');


end
