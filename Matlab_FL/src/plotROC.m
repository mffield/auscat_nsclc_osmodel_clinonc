function plotROC(AUC,TPR,FPR,clients_names)

figure('Color','w');
cols = 'rgbmk'; % use 'colorcube' colormap is more required
colors = distinguishable_colors(length(clients_names));
colors = distinguishable_colors(10);
colors(4,:)=[];colors(4,:)=[];

fpr_mu = squeeze(mean(FPR,2)); fpr_std = squeeze(std(FPR,[],2));
tpr_mu = squeeze(mean(TPR,2)); tpr_std = squeeze(std(TPR,[],2));
auc_mu = mean(AUC,2); auc_std = std(AUC,[],2);

lw = [2*ones(size(AUC,1)-1,1); 3];

for jj=1:size(AUC,1)-1
    if all(fpr_mu(:,jj)==0)
%         stairs(fpr_mu(:,jj),tpr_mu(:,jj),'Color',cols(jj),'linewidth',lw(jj)); hold on;
    else
        stairs(fpr_mu(:,jj),tpr_mu(:,jj),'Color',colors(jj,:),'linewidth',lw(jj)); hold on;
    end
    clients_names{jj} = [clients_names{jj} ' (AUC = ' num2str(auc_mu(jj),'%.3f') ')'];
end
stairs(fpr_mu(:,end),tpr_mu(:,end),'Color',[0 0 0],'linewidth',lw(end),'linestyle',':'); hold on;
clients_names{end} = [clients_names{end} ' (AUC = ' num2str(auc_mu(end),'%.3f') ')'];

title('ROC curve')
xlabel('False Positive Rate'); ylabel('True Positive Rate')
legend(clients_names,'Location','SouthEast')
end