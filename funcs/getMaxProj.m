function maxProj = getMaxProj(data,cmap,T)

if nargin < 3; T = 0;end

data = linmap(data,0,1);
cmap = imresize(cmap,[size(data,3),3],'bilinear');

data(data < T) = T;
[dataMax,map] = max(data,[],3);
% map(dataMax == T) = 0;

for k = 1:max(map(:))
    
    map0 = cmap(k,:);
    map = [linspace(0,map0(1),255)',linspace(0,map0(2),255)',linspace(0,map0(3),255)'];
    slice = data(:,:,k);
%     slice = round(linmap(slice,0,1));
    
    sliceR = slice.*map0(1);
    sliceG = slice.*map0(2);
    sliceB = slice.*map0(3);
    
    RGB(:,:,:,k) = cat(3,sliceR,sliceG,sliceB);
end

maxProj = uint8(linmap(max(RGB,[],4),0,255));