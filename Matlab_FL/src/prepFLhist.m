function FLHist = prepFLhist(model, cohort)

    FLHist=zeros(length(model.client(1).validation_result(1).(cohort).PHist),length(model.client),length(model.client(1).validation_result));
    for m = 1:length(model.client(1).validation_result)
        for c = 1:length(model.client)
            FLHist(:,c,m) = model.client(c).validation_result(m).(cohort).PHist';
        end
    end

end
