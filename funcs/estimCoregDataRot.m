% input a stack
function [dx,dy,T,cx,x,cy,y] = estimCoregDataRot(wf,loadCoreg,v1x,v2x,v1y,v2y)

if nargin < 2
    loadCoreg = 0;
end

if loadCoreg
    load('coreg.mat','dx','dy');
else
% estimate plane to plane shift
    dx = zeros(7,1); dy = dx;
    for k = 1:7
        [dx(k),dy(k)] = ccrShiftEstimation(wf(v1y:v2y,v1x:v2x,k),wf(v1y:v2y,v1x:v2x,k+1),5);
    end
    dx = cumsum(dx); dy = cumsum(dy);
end

% apply shift
for k = 1:7
    wf(:,:,k+1) = imtranslate(wf(:,:,k+1),[dx(k) dy(k)]);
end


% remove plane to plane non-overlaping pixels
[indy,indx] = find(sum(wf>0,3)  == 8);
cx = min(indx)+10; x = max(indx)-cx-12;
cy = min(indy)+10; y = max(indy)-cy-12;
wf2 = []; 
for k = 1:8
    wf2(:,:,k) = wf(cy:cy+y,cx:cx+x,k);
end

% estimate plane to plane rotation
for k = 1:7
   [dth(k),dx2(k),dy2(k)] = ccrRotEstimation(wf2(:,:,k),wf2(:,:,k+1),11,1);
end
dth = cumsum(dth);

% apply rotation
for k = 1:7
    wf2(:,:,k+1) = imrotate(wf2(:,:,k+1),dth(k),'bilinear','crop');
end
% normalize to the mean
T = (mean(mean(wf2,1),2))/mean(wf2(:));
plotStack(wf2./T)