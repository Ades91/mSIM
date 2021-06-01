% reorderPlanes raw camera images into x,y,z,t stack
function data = reorderPlanes(im,flip)

if nargin < 2
    flip = 1;
end

% flip camera 2
if flip
for k = 2:2:size(im,3)
    im(:,:,k) = fliplr(im(:,:,k));
end
end
%%

% off = 512;off2 = 0;
off = 450;off2 = 20;
data = zeros(600,off,size(im,3)/2,8,'uint16');
px = [45 560 1084 1575];
% px = (0:4).*512+1;
data(:,:,:,1) = im(:,px(1):px(1)+off-1,1:2:end-1);
data(:,:,:,3) = im(:,px(2):px(2)+off-1,1:2:end-1);
data(:,:,:,5) = im(:,px(3):px(3)+off-1,1:2:end-1);
data(:,:,:,7) = im(:,px(4):px(4)+off-1,1:2:end-1);
data(:,:,:,2) = im(:,off2+px(1):off2+px(1)+off-1,2:2:end);
data(:,:,:,4) = im(:,off2+px(2):off2+px(2)+off-1,2:2:end);
data(:,:,:,6) = im(:,off2+px(3):off2+px(3)+off-1,2:2:end);
data(:,:,:,8) = im(:,off2+px(4):off2+px(4)+off-1,2:2:end);
