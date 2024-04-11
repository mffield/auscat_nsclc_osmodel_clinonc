function xReturn = initialiseClientMPI(obj)
%initialiseClientMPI(obj)
%
%     Input:
%   
%     Output:
%       return = (boolean)

% Build up the argument lists.
values = { ...
   };
names = { ...
   };
types = { ...
   };

% Create the message, make the call, and convert the response into a variable.
soapMessage = createSoapMessage( ...
    'http://client.mpi.fl/', ...
    'initialiseClientMPI', ...
    values,names,types,'rpc');
response = callSoapService( ...
    obj.endpoint, ...
    '', ...
    soapMessage);
xReturn = parseSoapResponse(response);