function [wf,sim,wfDec,simDec] = batchSIM(data,doMask,doApod,dataDec,ki,k0,figID)

if ~exist('ki','var'); ki = 0.4; end
if ~exist('k0','var'); k0 = 0.8; end
if ~exist('figID','var'); figID = 0; end

 if mod(size(data,3),9+2) == 0
     Nph = [3 3 3];
     Ndiff = 1;
 elseif mod(size(data,3),15+2) == 0
     Nph = [5 5 5];
     Ndiff = 2;
     ki = ki/2;
 elseif mod(size(data,3),12+2) == 0
     Nph = [4 4 3 3]; 
     Ndiff = 1;
elseif mod(size(data,3),14+2) == 0
     Nph = [4 4 3 3]; 
     Ndiff = 1;
elseif mod(size(data,3),25+2) == 0
     Nph = [5 5 5 5 5]; 
     Ndiff = 2;
     ki = ki/2;
else
     Nph = [3 3 3];
     Ndiff = 1;
 end

nFrames = floor(size(data,3)/(sum(Nph)+2));
sy = size(data,1) - mod(size(data,1),2);
sx = size(data,2) - mod(size(data,2),2);
sim = zeros(2*sy,2*sx,size(data,4),nFrames); 
wf = zeros(sy,sx,size(data,4),nFrames);
if exist('dataDec','var')
    simDec = sim; wfDec = wf;
end

for f = 1:nFrames
    disp(['Processing frame # ',num2str(f)])
    fprintf(['[',repmat('-',[1,size(data,4)]),']\n']);
    fprintf('[');
    for k = 1:size(data,4)
        fprintf('-');
        im = data(:,:,(sum(Nph)+2)*(f-1)+(2:sum(Nph)+2),k);
        im = double(im(1:end-(mod(size(im,1),2)),1:end-(mod(size(im,2),2)),:));
        [sim(:,:,k,f),wf(:,:,k,f)]= getSIMauto(im,Nph,k0,ki,Ndiff,doMask,doApod,figID);
        if exist('dataDec','var')
            im = dataDec(:,:,(sum(Nph)+2)*(f-1)+(2:sum(Nph)+2),k);
            im = double(im(1:end-(mod(size(im,1),2)),1:end-(mod(size(im,2),2)),:));
            [simDec(:,:,k,f),wfDec(:,:,k,f)] = getSIMauto(im,Nph,k0,ki,Ndiff,doMask,doApod,figID);
        end
    end
    disp(']');
end