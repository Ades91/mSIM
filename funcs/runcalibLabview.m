function cal = runcalibLabview

% select calibration file
cdir = cd;
try
    cd('D:\Adrien\')
end
[fname,fpath] = uigetfile({'*.*';'*.All'}); cd(cdir);
im = loadData([fpath,fname]);
% reorder planes
data = reorderPlanes(im);

% stack order : yxtc
%% 
tf=cell(7,1);

for ch=(1:7)

    img_seq_ch_fix = data(:,:,:,ch);
    img_seq_ch_moving = data(:,:,:,ch+1);

    [~,indcom]=max(squeeze(mean(max(img_seq_ch_fix))).*squeeze(mean(max(img_seq_ch_moving))));

    im_fix=medfilt2(img_seq_ch_fix(:,:,indcom),[3 3]);
    im_mov=medfilt2(img_seq_ch_moving(:,:,indcom),[3 3]);

    % extracting the center of gravities of the beads
    out=struct;
    out.ru=5;
    sys.bg=1;
    out.nh=0;
	cal.bgth=-20; % background threshold
	cal.logsize=2.1; % size of the Laplacian of Gaussian filter
	cal.alol=15;  % lower limit of the segments area
	cal.aupl=100; % upper limit of the segments area
    im_fix = im_fix - 100; im_fix(im_fix < 0) = 0;
    im_mov = im_mov - 100; im_mov(im_mov < 0) = 0;
    out=hriSegmentation(double(im_fix),cal.bgth,cal.logsize,out);
    out=hriFilterSegments(double(im_fix),cal.aupl,cal.alol,sys,out);
    cog_fix=[out.xoAbsolute+0.1*out.xoEstimate, out.yoAbsolute+0.1*out.yoEstimate];
    out=hriSegmentation(double(im_mov),cal.bgth,cal.logsize,out);
    out=hriFilterSegments(double(im_mov),cal.aupl,cal.alol,sys,out);
    cog_mov=[out.xoAbsolute+0.1*out.xoEstimate, out.yoAbsolute+0.1*out.yoEstimate];

    im_fix = im_fix(50:end-50,50:end-50);
    im_mov = im_mov(50:end-50,50:end-50);
    [my, mx] = ccrShiftEstimation(im_fix,im_mov,5); % ! Stefan switched mx and my

    cog_common_fix=zeros(1,2);
    cog_common_mov=zeros(1,2);

    %identify corresponding center of gravities
    pixel_tolerance=2;
    for k=1:size(cog_fix,1)
        d=(cog_mov-repmat(cog_fix(k,:)-[mx my],[size(cog_mov,1) 1]));
        d=sqrt(d(:,1).^2+d(:,2).^2);
        [p,ind]=min(d);
        if(p<pixel_tolerance)
            cog_common_fix(end+1,:)=cog_fix(k,:);
            cog_common_mov(end+1,:)=cog_mov(ind,:);
        end
    end

    cog_common_fix(1,:)=[];  %first element is irrelevant (0,0)
    cog_common_mov(1,:)=[];

    d=(cog_common_mov+repmat([mx my],[size(cog_common_mov,1) 1]))-cog_common_fix;
    disp(['avg error of COG coordinates using pure displacement at pixel level: ' num2str(mean(sqrt(d(:,1).^2+d(:,2).^2)))]);

    % affine transformation
    ccm=fliplr(cog_common_mov);
    ccf=fliplr(cog_common_fix);
    tf{ch}=cp2tform(ccm,ccf,'similarity');

    [x,y] = tformfwd(tf{ch},cog_common_mov(:,2),cog_common_mov(:,1));
end
     

%% compute crop ROI to keep only overlapping pixels
dataCoreg = zeros(size(data,1),size(data,2),size(data,4));
for ch = 1:8
    [~,z0]=max(squeeze(mean(max(data(:,:,:,ch)))));
    temp = data(:,:,z0,ch);
    if ch == 1
        dataCoreg(:,:,1) = temp;
    else
        for n = ch:-1:2
            T = tf{n-1}.tdata.T;
            temp = imtransform(temp,maketform('affine',T), 'XData', [1 size(temp,2)],...
                'YData', [1 size(temp,1)],'Size', size(temp));
        end
        dataCoreg(:,:,ch) = temp;
    end
end
% plotStack(dataCoreg,13)

mask = dataCoreg > 0; mask = sum(mask,3) == 8;

list = find(mask(:,round(end/2))==1);
cy0 = list(1); cy1 = list(end);
list = find(mask(round(end/2),:)==1);
cx0 = list(1); cx1 = list(end);
cal.tf = tf;
off = 5; % 5 pixels offset
cal.cy0 = cy0+off; cal.cy1 = cy1-off; cal.cx0 = cx0+off; cal.cx1 = cx1-off;

%%  check coregistration
T = tf{ch-1}.tdata.T;
im_mov_tr=imtransform(im_mov,maketform('affine',T), 'XData', [1 size(im_fix,2)],...
            'YData', [1 size(im_fix,1)],'Size', size(im_fix));
figure(1234);imshowpair(im_fix,im_mov_tr);set(gcf,'visible','on')
        
