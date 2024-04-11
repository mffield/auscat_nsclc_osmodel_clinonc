function model = init_model(mode,model,features,client_names,training_id,testing_id)

if strcmpi(mode,'central')
    model.state={''};

    model.d=length(features);% feature dimensions in simulated data test
    model.features = strjoin(features,',');
    model.max_itr=2;
    model.clients_names=client_names;
    model.training_centers_index = training_id;
    model.testing_centers_index = testing_id;
    model.training_centers = strjoin(model.clients_names(model.training_centers_index),',');
    model.testing_centers = strjoin(model.clients_names(model.testing_centers_index),',');
    number_of_training_centers = length(strsplit(model.training_centers,','));

    model.parameter_names = {''};
    model.setup_parameters = {'localSum','localCov','localN','Mu'}; % names of setup parameters
    
    model.validation={''};
    model.randseed=1;
    
    
    model.time_prop=zeros(2,1);
    model.rrss_prop=zeros(2,1);
    model.val_time_strat=0;
    model.localSum = zeros(model.d,number_of_training_centers); % mean imputation parameters
    model.localCov = zeros(model.d*model.d,number_of_training_centers); % mean imputation parameters
    model.localN = zeros(1,number_of_training_centers); % number of data points per feature
        
    model.Mu = zeros(model.d,1);
    model.Cov = zeros(model.d,model.d);
    
elseif strcmpi(mode,'local')
    
    if isempty(model)
        model.alg_client_index=1;
        model.training_centers = {' '};
        model.testing_centers = {' '};
        model.features={' '};
        model.setup_parameters = {'localSum','localCov','localN','Mu'};
        model.parameter_names = {''};
        model.state={''};
        model.target={''};
        model.target_threshold=0;
        model.validation={''};
        model.time_prop=zeros(2,1);
        model.rrss_prop=zeros(2,1);
        model.val_time_strat=0;
        model.randseed=1;


    else
        
        model.d = length(strsplit(features{1},','));
        model.parameter_names = {''};
        model.setup_parameters = {'localSum','localCov','localN','Mu'};

        model.localSum=zeros(model.d,1);
        model.localCov=zeros(model.d*model.d,1);
        model.localN=0;
        
        model.Mu=0;
    end
    
end
model.algorithm='gaussian';

end


