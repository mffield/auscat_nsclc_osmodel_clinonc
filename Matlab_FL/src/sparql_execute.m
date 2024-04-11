%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SPARQL query executor function
% 
% Description: This function will try to execute a SPARQL query (given by
% the 2nd input parameter) on a SPARQL endpoint (given by the 1st input
% parameter). When successful (HTTP response code is 200), the XML file
% will be parsed into a cell-array.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [header, Data, extra] = sparql_execute(endpoint, query, parseResult)

%create header statement only accepting xml as response
if parseResult==1
    header = http_createHeader('Accept','application/sparql-results+xml;charset=UTF-8');
%     header = http_createHeader('Accept','application/rdf+json;charset=UTF-8');
else
    header = http_createHeader('Accept','text/csv;charset=UTF-8');
end




%open the URL, including the query as a GET-parameter
[binaryData, extra] = urlread2(endpoint,'POST',['query=' urlencode(query)],header);


%if the HTML response is OK (200), parse the xml
if isempty(strfind(extra.firstHeaders.Response,'200'))==0 && parseResult
    %temporarily store XML result as file
    fid = fopen('myFile.xml', 'w');
    fprintf(fid,'%s',binaryData);
    fclose(fid);
    
    xmlstruct = xml2struct('myFile.xml');
    header=[];
    v = xmlstruct.sparql.head.variable;
    for i=1:length(v)
        header = [header, {v{i}.Attributes.name}];
    end
    
    result = xmlstruct.sparql.results.result;
    Data = repmat({nan(1)},size(result,2),length(header));
    for i=1:length(result)
        binding = result{i}(1).binding;
        for j=1:length(binding)
            if isfield(binding{j},'literal')
                value = binding{j}.literal.Text;
                variableName = binding{j}.Attributes.name;
                variableIndex = ismember(header,variableName);
                
                if isfield(binding{j}.literal,'Attributes')
                    if isfield(binding{j}.literal.Attributes,'datatype')
                        datatype = binding{j}.literal.Attributes.datatype;
                    else
                        datatype = '';
                    end
                    variableName = binding{j}.Attributes.name;
                    variableIndex = ismember(header,variableName);
                    if contains(datatype,'double')
                        value = {str2double(strrep(value,',','.'))};
                    elseif contains(datatype,'int')
                        value = {str2double(value)};
                    elseif contains(datatype,'date')
                        value = strrep(value,'T',' ');
                        try
                            value = strsplit(value,' ');
                            value = {datetime(value{1},'InputFormat','yyyy-MM-dd')};
                        catch err
                            disp(err.message)
                            disp(value)
                        end
                    elseif strcmp(datatype,'') 
                        if isnan(value)
                            value = {''};
                        else
                            value = {value};
                        end
                    end
                else
                    %disp(variableName)
                    if isnan(value)
                        value = {''};
                    else
                        value = {value};
                    end
                end
                Data(i,variableIndex) = value;
                
            elseif isfield(binding{j},'uri')
                
                value = binding{j}.uri.Text;
                variableName = binding{j}.Attributes.name;
                variableIndex = ismember(header,variableName);
                Data(i,variableIndex) = {value};
            end
        end
    end
    Data = cell2table(Data);
    Data.Properties.VariableNames = header;
    
else
   Data = 0; 
   header = 0;
   disp('No data returned by SPARQL query')
end



end