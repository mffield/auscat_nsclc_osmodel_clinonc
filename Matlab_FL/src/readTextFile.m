function output = readTextFile(filename)

output = [];
fid = fopen(filename);
disp(['Processing: ' filename])
if fid~=0
    while ~feof(fid)
        line = fgetl(fid);
        output = [output '   ' line '   '];
    end
else
    error(['Cannot open file ' filename '.'])
end
fclose(fid);

end


