function [itr,model,central_message,mpi]=sendMessage(mpi,model,itr,var)


local_message = constructJSONMessage(model,var);


% fid = fopen(['message' num2str(mpi.ClientIndex) '.dat'],'w');
% fprintf(fid,'%s\n',local_message);
% fclose(fid);


sendReplyToMaster(mpi.svc, itr, local_message); prev_itr=itr;
disp(strcat('Client ',num2str(mpi.ClientIndex),'  sent reply itr ',num2str(itr),' :  ',local_message));


retry = 0; success=false;
while retry < 5 && ~success
    try
        [itr, model, ~, central_message, mpi] = Client_Receive_from_MPI(mpi,prev_itr,model);
        success=true;
    catch err
        retry=retry+1;
        disp(['Encountered an exception relating to function call "Client_Receive_from_MPI()" for client: ' num2str(mpi.ClientIndex) '. Attempt: ' num2str(retry)])
        disp(err.message);
        for st = 1:size(err.stack,1)
            disp(err.stack(st).file);disp(err.stack(st).name);disp(err.stack(st).line);
        end
    end
end
if ~success
    error(['Error after ' num2str(retry) ' attempts recovering from exceptions thrown from "Client_Receive_from_MPI()" for client: ' num2str(mpi.ClientIndex)])
end

% Wait until central algorithm has an update



end