function [itr, model, received_parameters, central_message, mpi] = Client_Receive_from_MPI(mpi,prev_itr,model)

timeout=500;
time1=datevec(now());time2 = datevec(now());
while (str2double(readMasterItr(mpi.svc)) == prev_itr) && (etime(time2,time1)<timeout)
    pause(0.25)
    time2 = datevec(now());
end

if etime(time2,time1)>=timeout
    disp(['Message system timeout - no incoming central message: exceeded ' num2str(timeout) ' seconds'])
    mpi.timeouts = mpi.timeouts+1;
    itr=int64(str2double(readMasterItr(mpi.svc)));
    central_message=readMasterMessage(mpi.svc);
    disp(strcat('Last received Message from central algorithm; ',central_message))
else
    itr=int64(str2double(readMasterItr(mpi.svc)));
    central_message=readMasterMessage(mpi.svc);
    disp(strcat('Received Message from central algorithm; ',central_message))
end
[model, received_parameters] = readJSONMessage('local',model,central_message,model.alg_client_index);

end