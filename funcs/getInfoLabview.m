function info = getInfoLabview(path)

%%
f = fopen(path,'r');

s = (fread(f,'*char'))';
N = length(s);
s = double(s);
s(s == 32 | s == 9 | s== 10 | s == 13) =[]; % remove space, tab and newline char
s = char(s);
info = [];
mapv = strfind(s,',');
mapb = strfind(s,':');
mapc = strfind(s,'}');

i = 2; % first char is always '{'
mapv(end+1) = mapc(end);
mapv(mapv < i) = []; mapb(mapb < i) = []; mapc(mapc < i) = [];
while i < N
    field = getNextField(s(i:mapb(1)));
    if strcmp(s(i+numel(field)+3),'{') == 0 % if the next : is not followed by a {
        %grab the field
        [field,value] = getNextField(s(i:mapv(1)));
        info.(field) = value;
        % set the cursor after the 
        i = mapv(1)+1;
    else % there is a subfield 
        % routine to extract the subfields names and values
        sub = s(i:mapc(1)); % subfield name
        g = strfind(sub,'"');
        v = strfind(sub,',');
        subname = sub(g(1)+1:g(2)-1);  subname = strrep(subname,' ','_');
        g(1:2) = []; v(end+1) = length(sub);
        Nsubfields = length(v); % compute # subfields using # of ,
        ii = g(1);
        for k = 1:Nsubfields
            [field,value] = getNextField(sub(ii:v(k)));
            info.(subname).(field) = value;
            ii = v(k)+1;
        end
        i = i+length(sub)+1;
    end
    
    mapv(mapv < i) = []; mapb(mapb < i) = []; mapc(mapc < i) = [];
    if isempty(mapb)
        i = N;
    end
end
fclose(f);

end

function [field,value] = getNextField(s)
    g = strfind(s,'"');
    field = s(g(1)+1:g(2)-1);
    value = s(g(2)+2:end-1);
    field = strrep(field,' ','_');
    field = strrep(field,'?','_');
    field = strrep(field,'/','_');
    if ~isempty(str2num(value))
        value = str2num(value);
    end
    if ~isnan(str2double(field(1))) % if the first char is not a letter
        field = ['z_',field];
    end
end