function MPI = ServerConnectionSetup(MPI_name,MPI_port,participating_clients,initialise)


%%% Process input arguments
number_of_participating_clients=length(strsplit(participating_clients,','));
number_of_clients=0;
disp(['Number of participating clients: ' num2str(number_of_participating_clients)])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Connect to distribution server
svc=MPIService(['http://localhost:' MPI_port '/' MPI_name '/MPIClient'],['http://localhost:' MPI_port '/' MPI_name '/MPIClient?wsdl']);
startedMPI = initialiseCentralMPI(svc);
if strcmp(startedMPI,'false')
    disp('MPI failed to start')
else
    disp('Initialisation of MPI successful')
    if initialise
        initClientInfo(svc);
    end
    
    if str2double(getNumberofCurrentClients(svc))~=number_of_participating_clients
        
        setAllowedClientsList(svc,participating_clients);
        setMasterOn(svc);
        disp('Wait for clients to connect to the server.')
        %%% Wait for each client to also connect to the server
        currentnumber=0;
        while currentnumber<number_of_participating_clients   %[str2double(getNumberofCurrentClients(svc)) number_of_participating_clients]
            pause(0.5);
            if str2double(getNumberofCurrentClients(svc)) > currentnumber
                currentnumber = str2double(getNumberofCurrentClients(svc));
                disp(currentnumber)
            end
        end
    end
    number_of_clients= str2double(getNumberofCurrentClients(svc));
    disp(currentnumber)
    disp('All clients have connected to the server')
    
end
MPI.svc=svc;
MPI.num_clients = number_of_clients;

end