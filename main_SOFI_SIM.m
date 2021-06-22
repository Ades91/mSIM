% 2D SOFI-SIM reconstruction
% load and format data, compute SOFI and reconstruct nth order SIM
% (see https://doi.org/10.1021/acsphotonics.1c00668)
addpath(genpath(['..',filesep,'ImDecorr',filesep,'funcs']));    % https://github.com/Ades91/ImDecorr
addpath(genpath(['..',filesep,'sofipackage']));                  % https://github.com/kgrussmayer/sofipackage
addpath(genpath(['..',filesep,'jsonlab-master']));              % https://github.com/fangq/jsonlab

addpath(genpath(cd));
%% custom data loading 
cdir = cd;
path = 'F:\20200612';
try
    cd(path)
end
[fname,pname] = uigetfile({'*.*;*.All'});
ind = strfind(fname,'.');
fnameC = fname(1:ind(end)-1); % cleaned filename
cd(cdir);
clear data; clear stack; clear tmp;
info = getInfoLabview(strrep([pname,fname],'.bin','.info'))
pause(0.01);
data = loadData([pname,fname]); 
% format data
disp('Reorder data')
Nph = ones(1,info.nA).*info.nPh
seqSize = sum(Nph);

% number of full sequences available
nSeq = floor(size(data,3)/(seqSize));
disp(['Number of SIM sequences: ',num2str(nSeq)])
stack = cell(seqSize,1);
for n = 1:seqSize; stack{n} = zeros(size(data,1),size(data,2),nSeq,'uint16'); end

count = 1; Nangles = 1; flipEvenPhases = 1;
if info.consecPattern == 0
    for n = 1:nSeq
        for an = 1:length(Nph)
            for k = 1:Nph(an)
                stackID = (an-1)*Nph(an) + k;
                if flipEvenPhases == 1
                    if mod(Nangles,2) == 1 % odd 
                        stack{stackID}(:,:,n) = data(:,:,count);
                        p = count;
                        stackID;
                    else
                        stack{stackID}(:,:,n) = data(:,:,Nangles*Nph(an)+1-k);
                        p = Nangles*Nph(an)+1-k;
                        stackID;
                    end
                else
                    stack{stackID}(:,:,n) = data(:,:,count);
                end
                count = count +1;
            end
            Nangles = Nangles +1;
        end
    end
else
    for n = 1:seqSize
        stack{n} = data(:,:,nSeq*(n-1)+1:nSeq*n);
    end
end
% clear data
%%
getDcorr(apodImRect(mean(stack{7},3),20))
%% optional pre-processing step: filter high-frequency noise
for n = 1:numel(stack)
    disp(['Stack #',num2str(n),'/',num2str(numel(stack))])
    stack{n} = filterImCirc(stack{n},0.5,0.6);
end
%% SIM only
SIMonlySTD = []; SIMonlyM = []; 
for k = 1:numel(stack)
    SIMonlyM(:,:,k) = mean(single(stack{k}),3);
    SIMonlySTD(:,:,k) = std(single(stack{k}),[],3);
end
% plotStack(SIMonlyM)
% pseudo WF image
pseudoWFM = mean(SIMonlyM,3);
pseudoWFSTD = mean(SIMonlySTD,3);
figure(34543);imagesc(pseudoWFSTD)

%% setup SOFI computation
settings = [];
% Cumulant calculation settings
settings.sys.orders = 2:3;
settings.sys.wsize = 500; % substack for analysis to avoid correlations from bleaching etc
settings.sys.sub = [];%[550];%10050%[20000]; % evaluate only first n frames (for quick preview)
settings.sys.start = [50]; %51; if empty, it takes all frames
settings.sys.end = [320]; % empty : use all frames
settings.sys.block = 1;
settings.sys.pxy = 97; % projected pixel size (in xy) [nm]  sofisetup = 96.0384, Hendriksetup = 104.8, AD-gut setup = 108
settings.sys.frames = size(stack{1},3);
settings.sys.jk = 0;
settings.sys.verbose = 0;

disp('Compute SOFIonly')
% rand perm allows to check if the subwindows are similar enough for
% SOFISIM
[sofionly_c, settings,stats] = sofiCumulantsAll(double(data(:,:,randperm(size(data,3)))), settings);
sofionly_c = padSofi(sofionly_c);
% clear data; pause(0.1);% release memory 

% compute SOFI
disp('Compute SOFI')
for k = 1:numel(stack)
    [sofi_c{k}, settings,stats] = sofiCumulantsAll(double(stack{k}), settings);
    sofi_c{k} = padSofi(sofi_c{k});
end

% Linearization settings
settings.dec.fwhm = 2.5;
settings.dec.iter = 8;
settings.dec.reconvolve = 0;
settings.dec.psfmodel = 'airy';
settings.dec.denoise = 0;

disp('Linearize')
for k = 1:numel(stack)
   sofi{k}=sofiLinearize(sofi_c{k},settings.dec.fwhm,[],settings.sys.orders,settings.dec.iter,settings);
end
disp('Linearize SOFIonly')
sofionly=sofiLinearize(sofionly_c,settings.dec.fwhm,[],settings.sys.orders,settings.dec.iter,settings);

%% intermediate result: display sofi images and estimate their resolution
o = 3;
figure(800*o);imagesc(mean(sofionly_c{o},3).^(1/o));title(['raw sofi',num2str(o)])
figure(801*o); imagesc(mean(sofionly{o},3));title(['sofi',num2str(o)])
imagescf(apodImRect(mean(sofionly_c{o},3),20),802*o);title(['rawsofi',num2str(o)])
imagescf(apodImRect(mean(sofionly{o},3),20),803*o);title(['sofi',num2str(o)])
getDcorr(gpuArray(apodImRect(mean(sofionly_c{o},3),20)))
getDcorr(gpuArray(apodImRect(mean(sofionly{o},3),20)))
%%
for n = settings.sys.orders
    sofionly{n} = mean(sofionly{n},3); 
    sofionly_c{n} = mean(sofionly_c{n},3);
end

% build SIM stack
SOFI_SIM = []; SOFI_SIMc = [];
for n = settings.sys.orders
    for k = 1:numel(stack)
        SOFI_SIM{n}(:,:,k) = mean(sofi{k}{n},3);
        SOFI_SIMc{n}(:,:,k) = mean(sofi_c{k}{n},3); 
    end
end

% plotStack(SOFI_SIM{2},342)
% plotStack(SOFI_SIM{3},343)
% figure(34534);imagesc(log(abs(fftshift(fftn(fftshift(apodImRect(SOFI_SIM3(:,:,8),20)))))+1));
%%
dk = 0.1; 
Ndiff = 1; k0 = 0.6;
pps = 80; ki = 0.46;pps/240/3;
doMask = 0; doApod = 0.05; % 'SIMonlyM' 'SIMonlySTD' 'SOFI_SIM2' 'SOFI_SIM3'
[K,A] = checkSIMpeak(double(SIMonlySTD),ki/Ndiff,dk,Nph,Ndiff) ;
ki = 0.52;median(K(:,1))*Ndiff
settings.sim.ki = ki; settings.sim.K = K; settings.sim.A = A;
%% compute SIM
% make data even before computation
Ndiff = 1; threshold = 0.1;
SIMonlySTD = double(SIMonlySTD(1:end-(mod(size(SIMonlySTD,1),2)),1:end-(mod(size(SIMonlySTD,2),2)),:));
SIMonlySTDSIM = getSIMauto(SIMonlySTD,Nph,k0,ki,Ndiff,doMask,doApod,500);
Ik = log(abs(fftshift(fftn(SIMonlySTDSIM)))+1);
SIMonlySTDSIM = resampleImage(SIMonlySTDSIM,ceil(sum(Ik(end/2,:) < threshold)/2)); pause(0.1)
Ndiffs = [1 2 3 4 5]; doApod = 0.01;
settings.sim.Ndiff = Ndiffs; settings.sim.doApod = doApod; settings.sys.doMask = 0;

SOFI_SIM_SIM = []; SOFI_SIMc_SIM = [];
for n = settings.sys.orders
    Ndiff = Ndiffs(n); 
    tmp = getSIMauto(SOFI_SIM{n},Nph,k0,ki/(size(SOFI_SIM{n},1)/size(SIMonlySTD,1)),Ndiff,doMask,doApod,(n==3)*1000);
    Ik = log(abs(fftshift(fftn(tmp)))+1);
    SOFI_SIM_SIM{n} = resampleImage(tmp,ceil(sum(Ik(end/2,:) < threshold)/2));
    % raw cumulants
    tmp = getSIMauto(SOFI_SIMc{n},Nph,k0,ki/(size(SOFI_SIMc{n},1)/size(SIMonlySTD,1)),Ndiff,doMask,doApod);
    Ik = log(abs(fftshift(fftn(tmp)))+1);
    SOFI_SIMc_SIM{n} = resampleImage(tmp,ceil(sum(Ik(end/2,:) < threshold)/2)); pause(0.1)
end

%% compute resolution and save data

resTh = [];
disp('Compute resolution')
resWF = 2*pps/getDcorr(gpuArray((pseudoWFSTD)),linspace(0,0.66,50)); kc = (2*pps)/resWF;
ppsSIM = pps/(size(SIMonlySTDSIM,1)/size(SIMonlySTD,1));
    resSIM = 2*ppsSIM/getDcorr(gpuArray(SIMonlySTDSIM));

% base theoretical resolution on sofi2 raw
ppsSOFIc = (pps/(size(sofionly_c{2},1)/size(SIMonlySTD,1)));
    resSOFI2raw= 2*ppsSOFIc/getDcorr(gpuArray(sofionly_c{2}));
kc = sqrt(2)*pps/resSOFI2raw;
resTh = 2*pps./[kc,(kc+ki)];
    
resSOFI = []; labelSOFI = [];
for n = settings.sys.orders
    ppsSOFIc(n) = (pps/(size(sofionly_c{n},1)/size(SIMonlySTD,1)));
        resSOFI(end+1) = 2*ppsSOFIc(n)/getDcorr(gpuArray(sofionly_c{n}));
    labelSOFI{end+1} = ['SOFIc',num2str(n)];
    resTh(end+1) = 2*pps/(sqrt(n)*kc);
    
    ppsSOFI(n) = (pps/(size(sofionly{n},1)/size(SIMonlySTD,1)));
        resSOFI(end+1) = 2*ppsSOFI(n)/getDcorr(gpuArray(sofionly{n}));
    labelSOFI{end+1} = ['SOFI',num2str(n)];
    resTh(end+1) = 2*pps/(n*kc);
end
disp('Compute resolution SOFI SIM')
resSOFISIM = []; labelSOFISIM = [];
for n = settings.sys.orders
    ppsSOFISIMc(n) = (pps/(size(SOFI_SIMc_SIM{n},1)/size(SIMonlySTD,1)));
        resSOFISIM(end+1) = 2*ppsSOFISIMc(n)/getDcorr(gpuArray(SOFI_SIMc_SIM{n}));
    labelSOFISIM{end+1} = ['SOFISIMc',num2str(n)];
    resTh(end+1) = 2*pps/(sqrt(n)*kc + n*ki); pause(0.1)
    
    ppsSOFISIM(n) = (pps/(size(SOFI_SIM_SIM{n},1)/size(SIMonlySTD,1)));
        resSOFISIM(end+1) = 2*ppsSOFISIM(n)/getDcorr(gpuArray(SOFI_SIM_SIM{n}));
    labelSOFISIM{end+1} = ['SOFISIM',num2str(n)];
    resTh(end+1) = 2*pps/((n)*kc + n*ki); pause(0.1)
end

%% plot resolution vs method
resSOFI(2:2:end)./resSOFI(1:2:end)
if sum(settings.sys.orders == 3)
    resArray = [resWF,resSIM, resSOFI,resSOFISIM  ];
    labelList = [{'WF','SIM'}  ,labelSOFI, labelSOFISIM];
end
settings.sim.resArray = resArray; settings.sim.labelList = labelList; settings.sim.kc = (2*pps)/resWF;

figure(34534);plot(resArray,'-ok','linewidth',1.5); hold on
plot(resTh,'--ob','LineWidth',1.5); hold off
ylabel('Resolution (nm)')
set(gca,'Xtick',1:numel(resArray))
set(gca,'XTickLabel',labelList)
set(gca,'xticklabelRotation',45)
ylim([25 300])
legend('Measured resolution','Theoretical resolution')
%% resolution plot
saveas(34534,[pname,fnameC,'_resplot'],'png');
saveas(34534,[pname,fnameC,'_resplot'],'epsc');
%% saving data
disp('Saving data')
writeTIFF(pseudoWFM,[pname,fnameC,'_pWFM'],[1/pps 1/pps])
writeTIFF(pseudoWFSTD,[pname,fnameC,'_pWFSTD'],[1/pps 1/pps])
writeTIFF(SIMonlySTDSIM,[pname,fnameC,'_SIMonly'],[1/ppsSIM 1/ppsSIM])
for n = settings.sys.orders
    writeTIFF(sofionly{n},[pname,fnameC,'_SOFI',num2str(n),'_only'],[1/ppsSOFI(n) 1/ppsSOFI(n)])
    writeTIFF(sofionly_c{n},[pname,fnameC,'_SOFI',num2str(n),'c_only'],[1/ppsSOFIc(n) 1/ppsSOFIc(n)])

    writeTIFF(SOFI_SIM_SIM{n},[pname,fnameC,'_SOFISIM',num2str(n)],[1/ppsSOFISIM(n) 1/ppsSOFISIM(n)])
    writeTIFF(SOFI_SIMc_SIM{n},[pname,fnameC,'_SOFISIM',num2str(n),'c'],[1/ppsSOFISIMc(n) 1/ppsSOFISIMc(n)])
end

savejson('',settings,[pname,fnameC,'_results.json']);
disp('Saving done')
%%
cmap = morgenstemning;% cmap = hot;

sz = size(SOFI_SIM_SIM{max(settings.sys.orders)}); resizeMethod = 'nearest';
in = [];
in{1} = imresize(pseudoWFSTD,sz,resizeMethod); in{2} = imresize(SIMonlySTDSIM,sz,resizeMethod);
for n = settings.sys.orders
   in{end+1} = imresize(sofionly{n},sz,resizeMethod); in{end}(in{end} < 0) = 0;
   in{end+1} = imresize(SOFI_SIM_SIM{n},sz,resizeMethod);  in{end}(in{end} < 0) = 0;
end
x0 = 0.55; x0 = [x0, x0+0.1]; montage = [];
x0 = round(x0.*sz(1)); xx = x0(1):x0(end); yy = xx;
montage = zeros(2*numel(xx),(numel(settings.sys.orders)+1)*numel(xx));
montage(1:end/2,1:numel(xx)) = linmap(in{1}(yy,xx),0,1);
montage(end/2+1:end,1:numel(xx)) = linmap(in{2}(yy,xx),0,1);
ystart = 1;
for n = 1:2:numel(in)-2
    ystart = ystart + numel(xx);
    montage(1:end/2,ystart:ystart+numel(xx)-1) = linmap(in{n+2}(yy,xx),0,1);
    montage(end/2+1:end,ystart:ystart+numel(xx)-1) = linmap(in{n+3}(yy,xx),0,1);
end

figure(34535);imagesc(montage);colormap(cmap);
%%
writeTIFF(montage,[pname,fnameC,'_montage'])

%% estimate drift from std
wLen = 1000;wf = [];
N = floor(size(data,3)/wLen)
for k = 1:N
    ind = (k-1)*wLen +1;
    wf(:,:,k) = apodImRect(medfilt2(std(single(data(:,:,ind:ind+wLen-1)),[],3),[2 2]),20);
end
plotStack(wf,5535353)
dx = []; dy = [];
for k = 1:size(wf,3)-1
    [dx(k),dy(k)] = ccrShiftEstimation(wf(:,:,k),wf(:,:,k+1),10);
end
dx = cumsum(dx); dy = cumsum(dy); dx = [0, dx]; dy = [0, dy];
t = linspace(1,size(data,3),numel(dx));
figure(345);plot(t,97.*dx); hold on;plot(t,97.*dy); hold off;legend('dx','dy')
ylabel('[nm]'); xlabel('Frame #'); title(['wLen: ',num2str(wLen),' frames'])
%% drift estimation plot
saveas(345,[pname,fnameC,'_driftplot'],'png');

%% phase-drift estimation
k0 = 0.55; r0 = 0.1;    

s = 2*info.nPh+(1:info.nPh); titleStr = 'All phases of angle #1'; legends = []; sname = '_phase_drift_all_phases';
% s = [1 8 15];  titleStr = ''; legends = {'Angle # 1','Angle #2','Angle #3'}; sname = '_phase_drift';
% 2','Angle # 3'};
ph = zeros(size(stack{1},3),numel(s)); phAvg = zeros(numel(s),1);
for ss = 1:numel(s)
        id = s(ss);
        tmp = std(single(stack{id}(:,:,1:10)),[],3);
%         tmp = single(stack{id}(:,:,1));
        cx = size(tmp,1)/2; cy = size(tmp,2)/2;
        x = linspace(-1,1,size(tmp,1)); [X,Y] = meshgrid(x); R = sqrt(X.^2 + Y.^2);
        OTFinv = (R./k0).*(R < k0);
        Ik = log(abs(fftshift(fftn(apodImRect(tmp,20))))+1);
        [p2x,p2y] = getSIMPeak(OTFinv.*Ik.*(Y > 0));    
        mask = sqrt((X - (p2x-cx)/cx).^2 + (Y- (p2y-cy)/cy).^2) < r0 | ...
       sqrt((X + (p2x-cx)/cx).^2 + (Y + (p2y-cy)/cy).^2) < r0;
    figure(101201);him = imagesc(OTFinv.*Ik); hold on;plot(p2x,p2y,'kx');hold off; 
    xlim([p2x-20,p2x+20]);ylim([p2y-20 p2y+20]);pause(0.01)
    for n = 1:size(stack{id},3)
        Ik = fftshift(fftn(apodImRect(single(stack{id}(:,:,n)),20)));
        ph(n,ss) = getSIMphase(Ik,p2x,p2y);
        set(him,'cdata',log(abs(Ik)+1)); pause(0.01);
    end
    phAvg(ss) = getSIMphase(fftshift(fftn(apodImRect(single(mean(single(stack{id}),3)),20))),p2x,p2y);
end

figure(234);plot(unwrap(ph(1:end,:))); hold on;
plot([1 info.N],[phAvg phAvg],'linewidth',1.5); hold off
xlabel('Frame number'); ylabel('Phase (rad)');
title(titleStr)
if ~isempty(legends); legend(legends,'location','southeast'); end
%%
saveas(234,[pname,fnameC,sname],'png');

%%
order = 2;
fluctuationMap = (std(SOFI_SIMc{order},[],3)./mean(SOFI_SIMc{order},3));
mask = mean(SOFI_SIMc{order},3) > mean(SOFI_SIMc{order}(:));

meanFluctuations = mean(fluctuationMap(mask));

figure(345325);imagesc(mask.*(fluctuationMap - min(fluctuationMap(mask)))+min(fluctuationMap(mask))); title(['Average signal fluctuations ',num2str(meanFluctuations*100,3),' %'])

%%
writeTIFF(fluctuationMap,[pname,fnameC,sname])

%% line profiles

order = 3;
figure(23454);imagesc(SOFI_SIM_SIM{order}); colormap(morgenstemning)
roi = imline(gca); roi = roi.getPosition; pts.point1 = roi(1,:); pts.point2 = roi(2,:);
refImgSz = size(SOFI_SIM_SIM{order},1); width = 10;
lpSOFISIM = getProfile(SOFI_SIM_SIM{order},pts,width); lpSOFISIM(lpSOFISIM <0) = 0;lpSOFISIM = linmap(lpSOFISIM,0,1);
lpWF = getProfile(imresize(pseudoWFSTD,[refImgSz, refImgSz]),pts,width); lpWF(lpWF <0) = 0;lpWF = linmap(lpWF,0,1);
lpSIM = getProfile(imresize(SIMonlySTDSIM,[refImgSz, refImgSz]),pts,width);lpSIM(lpSIM <0) = 0; lpSIM = linmap(lpSIM,0,1);
lpSOFIonly = getProfile(imresize(sofionly{order},[refImgSz, refImgSz]),pts,width);lpSOFIonly(lpSOFIonly <0) = 0; lpSOFIonly = linmap(lpSOFIonly,0,1);
figure(98765); plot(lpSOFISIM,'r','linewidth',1.5); hold on
plot(lpWF,'k','linewidth',1.5)
plot(lpSIM,'g','linewidth',1.5)
plot(lpSOFIonly,'b','linewidth',1.5); hold off
legend('SOFISIM','WF','SIM','SOFI')
fhwm = [getFWHM(lpWF) getFWHM(lpSIM) getFWHM(lpSOFIonly) getFWHM(lpSOFISIM)].*ppsSOFISIM(order)
%% correct phase drift
dxint = imresize(dx,[1 700],'bilinear'); dyint = imresize(dy,[1 700],'bilinear');
id = 8; r0 = 0.05;
tmp0 = apodImRect(mean(single(stack{id}(:,:,5)),3),20); cx = size(tmp0,1)/2; cy = size(tmp0,2)/2;
x = linspace(-1,1,size(tmp0,1)); [X,Y] = meshgrid(x); R = sqrt(X.^2 + Y.^2);
OTFinv = (R./k0).*(R < k0);
Ik = log(abs(fftshift(fftn(apodImRect(tmp0,20))))+1);
[p2x,p2y] = getSIMPeak(OTFinv.*Ik.*(Y > 0));
ki = sqrt(((p2x-cx)/cx).^2 + ((p2y-cy)/cy).^2);
th = atan2((p2y-cy)/cy,(p2x-cx)/cx);
T = 2/ki; % period in pixels
mask = sqrt((X - (p2x-cx)/cx).^2 + (Y- (p2y-cy)/cy).^2) < r0 | ...
       sqrt((X + (p2x-cx)/cx).^2 + (Y + (p2y-cy)/cy).^2) < r0;
% figure(ss);imagesc(OTFinv.*Ik); hold on;plot(p2x,p2y,'kx');hold off; pause(0.01)
ph = []; ph2 = [];
ph(1) = getSIMphase(fftshift(fftn(apodImRect(tmp0,20))),p2x,p2y);
ph2(1) = ph(1);
stack2 = zeros(size(tmp,1),size(tmp,2),size(stack{id},3));
stack2(:,:,1) = real(ifftshift(ifftn(ifftshift(mask.*fftshift(fftn(tmp0))))));
for n = 2:size(stack{id},3)
    tmp = apodImRect(single(stack{id}(:,:,n)),20);
    ph(n) = getSIMphase(fftshift(fftn(apodImRect(tmp,20))),p2x,p2y);
    Tn = (ph(n)-ph(1))*T/(2*pi);
     tmp = imtranslate(tmp,[Tn*cos(th) Tn*sin(th)]);
%     tmp = imtranslate(tmp,[dxint(n) dyint(n)]);
%     Ik = fftshift(fftn(fftshift(apodImRect(tmp,20)))); 
    stack2(:,:,n) = tmp;%real(ifftshift(ifftn(ifftshift(mask.*Ik))));
     ph2(n) = getSIMphase(fftshift(fftn(apodImRect(tmp,20))),p2x,p2y);
    
end
phaseDrift = (dx.*cos(th) + dy.*sin(th))*2*pi/T; % phase change originating from sample drift
N= size(stack{id},3);
figure(100);plot(linspace(1,N,N),unwrap(ph)); hold on; plot(linspace(0,N,N),unwrap(ph2)); hold off
xlabel('Frame number'); ylabel('Phase (rad)'); legend(['Angle #',num2str(ceil(id/7))],'Phase drift from sample drift')

plotStack(stack2)
%%
saveas(100,[pname,fnameC,'_phase_drift_vs_mechanical_drift'],'png');