function model=initialise_network_data(mpi, model)

disp(['Starting algorithm: ' model.algorithm])

model.itr=1;
model.state='initialise';

central_message = constructJSONMessage(model,[{'algorithm'},{'state'},{'features'},...
    {'target'},{'settings'}]);
model_ = model;
[model, model.itr] = Broadcast_Receive_from_MPI(mpi.svc, central_message, mpi.num_clients, model.itr, model_);

% model.Mu = sum(model.localSum,2)./(sum(model.localN,2));
% model.Max = max(model.localMax,[],2);
% model.Min = min(model.localMin,[],2);

model.settings.client_names = {model.client.client_name};

model.central_stats = central_stats_collection(model,'stats');
model.central_stats_original=central_stats_collection(model,'stats_original');

central_message = constructJSONMessage(model,[{'algorithm'},{'state'},{'central_stats'},{'central_stats_original'}, {'settings'}]);
model_ = model;
[model, model.itr] = Broadcast_Receive_from_MPI(mpi.svc, central_message, mpi.num_clients, model.itr, model_);


model.central_stats = central_stats_collection(model,'stats');
model.central_stats_cur = central_stats_collection(model,'stats_cur');
model.central_stats_original=central_stats_collection(model,'stats_original');
model.Mu = model.central_stats_cur.Mu;


central_message = constructJSONMessage(model,[{'algorithm'},{'state'},{'central_stats_cur'},{'Mu'}]);
model_ = model;
[model, model.itr] = Broadcast_Receive_from_MPI(mpi.svc, central_message, mpi.num_clients, model.itr, model_);

model.central_stats = central_stats_collection(model,'stats');
model.central_stats_cur = central_stats_collection(model,'stats_cur');
model.Std = model.central_stats_cur.Std;


end