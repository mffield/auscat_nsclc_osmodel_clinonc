% LocalLungModel.m

% This is intended to run at each data client and assumes the
% CentralLungDSS.m script is running and listening for clients.

% restoredefaultpath % this is useful for testing/debugging code locally

% Function can be compiled for distribution across the network.
% mcc -m LocalLungModel.m -d ./output  -a ./src -a '.\@MPIService'

function LocalLungModel(MPI_config,MPI_name)

config = readConfigFile(MPI_config);
mpi = initialise_client_comm(config,MPI_name);

if mpi.startedMPI
    prev_itr=-1;
    model.state='';
    while ~strcmp(model.state,'terminate_function')
        
        model = init_model('local',[],[],[],[],[]);
        model.client_name = config.clientName;
        
        [itr, model, ~, central_message] = Client_Receive_from_MPI(mpi,prev_itr,model);
        [model, ~] = readJSONMessage('local',model,central_message,model.alg_client_index);
        
        % Load data and preprocess
        [data, ~]= loadLungData(config, model);
       
        % Partition data for validation procedures
        [model, data] = initialize_validations(model, data, model.settings.validation, double(data.table.presc_dose_ttl>model.settings.curative_dose_threshold));
        
        while ~strcmp(model.state,'terminate_algorithm')
            if contains(central_message,'algorithm":"logreg')
                [var,model,mpi,itr]=logreg_local_update(mpi,itr,central_message,model,data);
            end
            % Sending back local message
            [itr,model,central_message]=sendMessage(mpi,model,itr,var);
            
        end
        disp('Terminating algorithm session and re-initialising.')
        [itr,model,central_message]=sendMessage(svc,model,itr,[]);
    end
    disp('Terminating function.')
    [itr,model,central_message]=sendMessage(svc,model,itr,[]);
    %         sendReplyToMaster(svc, itr, 'null'); prev_itr=itr;
    
end

end

    