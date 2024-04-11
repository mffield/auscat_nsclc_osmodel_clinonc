function obj = MPIService(endpoint,wsdl)

obj.endpoint = endpoint;
obj.wsdl = wsdl;

obj = class(obj,'MPIService');

