function mpi = initialise_client_comm(config,MPI_name)

    mpi=struct;
    if ~isfield(config, 'MPIClientPort')
        config.MPIClientPort='80';
    end
    if ~isfield(config, 'MPIClientHost')
        config.MPIClientHost='localhost';
    end
%     mpi.svc=MPIService(['http://127.0.0.1:' config.MPIClientPort '/' MPI_name '/MPIClient'],...
%                    ['http://127.0.0.1:' config.MPIClientPort '/' MPI_name '/MPIClient?wsdl']);
    ['http://' config.MPIClientHost ':' config.MPIClientPort '/' MPI_name '/MPIClient']

    mpi.svc=MPIService(['http://' config.MPIClientHost ':' config.MPIClientPort '/' MPI_name '/MPIClient'],...
                   ['http://' config.MPIClientHost ':' config.MPIClientPort '/' MPI_name '/MPIClient?wsdl']);
               
    disp(['Client starting from ' config.clientName])

    mpi.startedMPI = initialiseClientMPI(mpi.svc);
    mpi.timeouts=0;
    if strcmp(mpi.startedMPI,'false')
        mpi.startedMPI=false;
        mpi.comm_message = 'MPI failed to start'; disp(mpi.comm_message)
    else
        disp('Initialisation of MPI successful. Waiting for central algorithm to start');
        while strcmp(isMasterOn(mpi.svc),'false');  pause(0.25);  end
        disp('Central algorithm is ON');
        
        mpi.ClientIndex=str2double(addClient(mpi.svc));
        if mpi.ClientIndex==-1
            mpi.comm_message = ['Cannot add Client: ' config.clientName]; disp(mpi.comm_message)
            mpi.startedMPI = false;
        else
            mpi.comm_message = strcat(['Client ' config.clientName ' added with index: '],num2str(mpi.ClientIndex)); disp(mpi.comm_message)
            mpi.startedMPI = true;
        end
        
    end

end
