% input : prefiltered absolute value image 
% output : most likely peak
function [x,y,ph,A] = getSIMPeak(Ik)

%%
N = 10;

% look for the first Nth maximum in the image
temp = abs(Ik);
for k = 1:N
    
    [~,ind] = max(temp(:));
    [cy,cx] = ind2sub(size(temp),ind);
    
    p(k).x = cx;
    p(k).y = cy;
    temp(cy,cx) = 0;
    
end

% cluster the maximums based on their peaks 
cim = zeros(size(Ik));
for k = 1:length(p)
    cim(p(k).y,p(k).x) = 1;
end

cim = imclose(cim,strel('diamond',2));
c = clusterIm(cim,1);


% find the largest cluster
cmax = c{1};
sz(1) = size(c{1},1); A(1) = abs(Ik(c{1}(2),c{1}(1)));
for k = 2:length(c)
   sz(k) = size(c{k},1); A(k) = abs(Ik(c{k}(2),c{k}(1)));
   if size(cmax,1) < size(c{k},1)
       cmax = c{k};
   end
end
if sum(sz) == numel(sz) % if all cluster have a size of 1
    % pick the brightest
    [~,ind] = max(A);
    cmax = c{ind};
end
% make mask from cluster
m = zeros(size(Ik));
xmin = 100000;
xmax = 0;
ymin = 100000;
ymax = 0;
for k = 1:size(cmax,1)
     m(cmax(k,2),cmax(k,1)) = 1;
    if ymin > cmax(k,2);ymin = cmax(k,2); end
    if ymax < cmax(k,2);ymax = cmax(k,2); end
    if xmin > cmax(k,1);xmin = cmax(k,1); end
    if xmax < cmax(k,1);xmax = cmax(k,1); end
    
end

imcrop = Ik(ymin:ymax,xmin:xmax).*m(ymin:ymax,xmin:xmax);
scale = 1;
temp = imresize(imcrop,scale,'bilinear');

[~,ind] = max(abs(temp(:)));
[y,x] = ind2sub(size(temp),ind);

ph = angle(temp(y,x));
A = abs(temp(y,x));
y = ymin-1 + y/scale;
x = xmin-1 + x/scale;