function dataDec = deconvolveData(data,nIter)

%% load psf and resample it
psf = loadData('MP_psf.tif');
psf(:,:,end) = [];
addpath('C:\gitRepositories\wdt\calcOTF')

psfR = linmap(rescalePSF(linmap(psf,0,1),[size(psf,1)*0.108 size(psf,2)*0.108 size(psf,3)*0.2],...
            [size(data,1)*0.108 size(data,2)*0.108 8*0.35],[size(data,1) size(data,2) 8]),0,1);
psfR = gpuArray(linmap(psfR,0,1));
psf2 = gpuArray(psfR(:,:,5));

% deconvolve raw SIM frames
dataDec = zeros(size(data));
for f = 1:size(data,3)
    disp(['Processing frame # ',num2str(f)])
    ph = apodImRect(gpuArray(squeeze(data(:,:,f,:))),20);
    dataDec(:,:,f,:) = gather(deconvlucy(ph, psfR,nIter));
    T = squeeze(mean(mean(dataDec(:,:,f,:),1),2));
    % deconvolve the top and bottom plane in 2D
    dataDec(:,:,f,1) = gather(deconvlucy(ph(:,:,1),psf2,nIter));
    dataDec(:,:,f,8) = gather(deconvlucy(ph(:,:,8),psf2,nIter));
    dataDec(:,:,f,1) = T(1)*dataDec(:,:,f,1)/mean(mean(dataDec(:,:,f,1)));
    dataDec(:,:,f,8) = T(8)*dataDec(:,:,f,8)/mean(mean(dataDec(:,:,f,8)));
end
