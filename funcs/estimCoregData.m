% input a stack
function [dx,dy,T,cx,x,cy,y] = estimCoregData(wf,loadCoreg,v1x,v2x,v1y,v2y)

if nargin < 2
    loadCoreg = 0;
end

if loadCoreg
    load('coreg.mat','dx','dy');
else
% estimate plane to plane shift
    dx = zeros(size(wf,3)-1,1); dy = dx;
    for k = 1:size(wf,3)-1
        [dx(k),dy(k)] = ccrShiftEstimation(wf(v1y:v2y,v1x:v2x,k),wf(v1y:v2y,v1x:v2x,k+1),5);
    end
    dx = cumsum(dx); dy = cumsum(dy);
end

% apply shift
for k = 1:size(wf,3)-1
    wf(:,:,k+1) = imtranslate(wf(:,:,k+1),[dx(k) dy(k)]);
end

% remove plane to plane non-overlaping pixels
[indy,indx] = find(sum(wf>0,3)  == 8);
cx = min(indx)+10; x = max(indx)-cx-12;
cy = min(indy)+10; y = max(indy)-cy-12;
wf2 = []; 
for k = 1:size(wf,3)
    wf2(:,:,k) = wf(cy:cy+y,cx:cx+x,k);
end
% normalize to the mean
T = (mean(mean(wf2,1),2))/mean(wf2(:));
%plotStack(wf2./T)