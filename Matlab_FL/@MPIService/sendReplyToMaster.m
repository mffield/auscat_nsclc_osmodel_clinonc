function xReturn = sendReplyToMaster(obj,arg0,arg1)
%sendReplyToMaster(obj,arg0,arg1)
%
%     Input:
%       arg0 = (int)
%       arg1 = (string)
%   
%     Output:
%       return = (int)

% Build up the argument lists.
values = { ...
   arg0, ...
   arg1, ...
   };
names = { ...
   'arg0', ...
   'arg1', ...
   };
types = { ...
   '{http://www.w3.org/2001/XMLSchema}int', ...
   '{http://www.w3.org/2001/XMLSchema}string', ...
   };

% Create the message, make the call, and convert the response into a variable.
soapMessage = createSoapMessage( ...
    'http://client.mpi.fl/', ...
    'sendReplyToMaster', ...
    values,names,types,'rpc');
response = callSoapService( ...
    obj.endpoint, ...
    '', ...
    soapMessage);
xReturn = parseSoapResponse(response);
