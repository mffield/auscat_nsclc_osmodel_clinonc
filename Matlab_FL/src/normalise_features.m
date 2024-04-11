function stats=normalise_features(model,Data,stats)


x = table2array(Data(:,model.features));
for i = 1:size(x,2)
    stats.localN(i) = size(x,1) - sum(isnan(x(:,i)));
    stats.localSum(i) = sum(x(~isnan(x(:,i)),i));
end

if isfield(model,'Mu')
    
    for i = 1:size(x,2)
        if any(isnan(x(:,i)))
            x(isnan(x(:,i)),i) = model.Mu(i);
        end
    end
    centredData = x - repmat(model.Mu(:)',size(x,1),1);
    stats.localCov = centredData'*centredData;
    stats.localStdSum = sum(centredData.^2);
end
