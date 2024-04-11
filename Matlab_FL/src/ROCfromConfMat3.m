function [AUC, TPR, FPR]=ROCfromConfMat3(model)

AUC = zeros(length(model.client)+1,length(model.client(1).validation_result));
for k = 1:length(model.client(1).validation_result)
    
    CMTotal = zeros(size(model.client(1).validation_result(k).confusionMatrix));
    
    for i = 1:length(model.client)
        
        CM = model.client(i).validation_result(k).confusionMatrix;
        
        CMTotal = CMTotal + CM;        
        for t=1:size(CM,3)
            TPR(t,k,i) = CM(1,1,t)/sum(CM(1,:,t));
            FPR(t,k,i) = CM(2,1,t)/sum(CM(2,:,t));
        end
        AUC(i,k) = sum(abs(diff(FPR(2:end,k,i))).*TPR(3:end,k,i));
        
    end
    for t=1:size(CM,3)
        TPR(t,k,length(model.client)+1) = CMTotal(1,1,t)/sum(CMTotal(1,:,t));
        FPR(t,k,length(model.client)+1) = CMTotal(2,1,t)/sum(CMTotal(2,:,t));
    end
    AUC(length(model.client)+1,k) = sum(abs(diff(FPR(2:end,k,length(model.client)+1))).*TPR(3:end,k,length(model.client)+1));
    
end

end
