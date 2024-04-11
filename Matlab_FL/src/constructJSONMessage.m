function message = constructJSONMessage(model,parameters)

    if isempty(parameters)
        message = '{}';
    else
        for i=1:length(parameters)
            packet.(parameters{i}) = model.(parameters{i});
        end
        message = jsonencode(packet);
    end

end
