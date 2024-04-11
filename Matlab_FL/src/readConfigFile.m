%
function config=readConfigFile(filename)

fid = fopen(filename);

while (~feof(fid))
   
    line = fgetl(fid);
    parse_location = find(line=='=');
    if ~isempty(line)
    if line(1)~='#' % using # as comment indicator        
        if ~isempty(parse_location) 
            param = strtrim(line(1:parse_location-1));
            value = strtrim(line(parse_location+1:end));
            
            config.(param) = value;
        end
    end
    end
end

fclose(fid);
end