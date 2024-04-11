function xReturn = setAllowedClientsList(obj,arg0)
%setAllowedClientsList(obj,arg0)
%
%     Input:
%       arg0 = (string)
%   
%     Output:
%       return = (int)

% Build up the argument lists.
values = { ...
   arg0, ...
   };
names = { ...
   'arg0', ...
   };
types = { ...
   '{http://www.w3.org/2001/XMLSchema}string', ...
   };

% Create the message, make the call, and convert the response into a variable.
soapMessage = createSoapMessage( ...
    'http://client.mpi.fl/', ...
    'setAllowedClientsList', ...
    values,names,types,'rpc');
response = callSoapService( ...
    obj.endpoint, ...
    '', ...
    soapMessage);
xReturn = parseSoapResponse(response);