% 2D SIM reconstruction
% Simple and robust reconstruction. Can handle any combination of number of
% phases and angles. Tools are provided to estimate the illumination
% spatial frequency. Reconstruction parameters are itteratively estimated
% from the data.

addpath(['.',filesep,'funcs'])
cdir = cd;
try
%     cd('D:\Measurements')
    cd('D:\Adrien\')
end
[fname,fpath] = uigetfile({'*.*';'*.All'});
cd(cdir);
% can load tiff stacks or binary (json file needed)
im = loadData([fpath,fname]);

plotStack(im,353453)

%% check if SIM peaks are properly found
im = im(1:min(size(im,1),size(im,2)),1:min(size(im,1),size(im,2)),:);
% check and make data even before computation
im = im(1:end-(mod(size(im,1),2)),1:end-(mod(size(im,2),2)),:);
% user defined parameters used in the reconstruction
pps = 97;           % pixel size in nanometer
ki = 2*pps/400;     % illumination spatial frequency
dk = 0.1;           % uncertainty of ki (search in ki +/- dk)
% number of phases for each angle
Nph = [5 5 5 5 5]; [7 7 7 7 7];[3 3 3 3 3];[3 3 3];[5 5 5]; [3 3 3 3 3]; [5 5 5];    [4 4 4 4]; [4 4 3 3];   [3 3 4 4];
Ndiff = 1;          % # of diffraction order to consider for reconstruction
% check if peaks are inside the searching area
K = checkSIMpeak(double(im),ki,dk,Nph,Ndiff);
ki = mean(K)
ki = ki(1);
T = 2*pps/ki        % average pattern period in nanometer

%% SIM reconstruction
k0 = 0.6;                   % input image normalized cutoff frequency
doMask = 0; doApod = 0.05;  % doMask: apply a circular mask of radius k0
                            % doApod: apodize the peaks intensity 
figID = 20;                 % figure number of the displayed results
% make data even before computation
im = double(im(1:end-(mod(size(im,1),2)),1:end-(mod(size(im,2),2)),:));
% getSIMauto estimate the reconstruction parameters directly from the image
% ki is just an input guess
[sim,wf] = getSIMauto(im,Nph,k0,ki,Ndiff,doMask,doApod,figID);


%% Real space resampling via Fourier space cropping
pps = 97;
% make sim square
s = padarray(sim,([max(size(sim)) max(size(sim))]-size(sim))./2,sim(1,1));
R2 = fftshift(fftn(fftshift(s)));
% R2(abs(R2)< 1) = 1.*10^4;

% cropX = floor(floor(max(sum(sum(R2,1)>1)/1.5,sum(sum(R2,2)>1)/1.5))/2);
crop = floor(sum(abs(R2(end/2,:))< 1)./2);
% crop = 190;
R3 = R2(crop+1:end-crop,crop+1:end-crop);
simResamp = real(ifftshift(ifftn(ifftshift(R3))));

ppsSIM = pps/(size(simResamp,1)/size(im,1));
% figure;imagesc(log(abs(R3)+1))
figure(4);imagesc(log(abs(R3)+1))
% Ik = fftshift(fftn(fftshift(wf)));
% Ik2 = padarray(Ik,(size(simResamp)-size(Ik))./2);

%% resolution estimate
addpath(['.',filesep,'ImDecorr']) % https://github.com/Ades91/ImDecorr

kcSIM = getDcorr(gpuArray(simResamp),linspace(0,1,50),10,1010);
kc = getDcorr(gpuArray(apodImRect(wf,20)));
resWF = 2*pps/kc
resSIM = 2*ppsSIM/kcSIM
resSIMtheoretical = 2*pps/(kc+ki)
%% saving
ind = strfind(fname,'.');
fnameC = fname(1:ind(end)-1);
writeTIFF(wf./max(wf(:)),[fpath,fnameC,'_wfMean'],[1/pps 1/pps]);
writeTIFF(simResamp./max(simResamp(:)),[fpath,fnameC,'_SIM'],[1/ppsSIM 1/ppsSIM]);
