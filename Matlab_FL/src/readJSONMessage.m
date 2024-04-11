
function [model, received_parameters] = readJSONMessage(mode,model,message,client_index)

    if strcmp(mode,'central')
        % message from client to central server
        % remove the initial text from string
        %[~,jsonmsg]=strtok(message,'{');
        temp = jsondecode(message);
        if ~isempty(temp)
            msg_fields = fieldnames(temp);
            for i=1:length(msg_fields)
                model.client(client_index).(msg_fields{i}) = temp.(msg_fields{i});
            end
            received_parameters = fieldnames(model.client(client_index));
        else
            received_parameters=[];
        end
        %determine client index and insert array into appropriate field of
        %array?
        
%         agg_param_names = model.agg_param_names;
%         for j=find(ismember(agg_param_names,received_parameters))% 1:length(agg_param_names)
%             
%             model.(agg_param_names{j})(:,client_index) = model.client(client_index).(agg_param_names{j})(:);
% 
%         end
        
    elseif strcmp(mode,'local')
    
        % message from central server
        %[~,jsonmsg]=strtok(message,'{');
        model.central = jsondecode(message);
        received_parameters = fieldnames(model.central);
        for j=1:length(received_parameters)
            model.(received_parameters{j}) = model.central.(received_parameters{j});
        end
    end

end