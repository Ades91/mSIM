function saveRGBProj(wf,sim,fpath,fnameC,r,wfDec,simDec)

%% compute RBG max proj
projWF = []; projSIM = [];
if exist('wfDec','var')
    projWFDec = []; projSIMDec = [];
end
pps = 108; ppsSIM = pps/(size(sim,1)/size(wf,1));

cmap = flipud(isolum);
zrange = 2:8;
% r = [1 1 1 1]; [1.5 2.2 3 2];
for f = 1:size(sim,4)
    projWF(:,:,:,f)  = linmap(double(getMaxProj(wf(:,:,zrange,f) ,cmap)),0,r(1));
    projSIM(:,:,:,f) = linmap(double(getMaxProj(sim(:,:,zrange,f),cmap)),0,r(2));
    if exist('wfDec','var')
        projWFDec(:,:,:,f) = linmap(double(getMaxProj(wfDec(:,:,zrange,f) ,cmap)),0,r(3));
        projSIMDec(:,:,:,f) = linmap(double(getMaxProj(simDec(:,:,zrange,f),cmap)),0,r(4));
    end
end

vWF = mean(mean(mean(projWF(:,:,:,1))));
vSIM = mean(mean(mean(projSIM(:,:,:,1))));
if exist('wfDec','var')
    vWFDec = mean(mean(mean(projWFDec(:,:,:,1))));
    vSIMDec = mean(mean(mean(projSIMDec(:,:,:,1))));
end
for f = 2:size(sim,4)
    projWF(:,:,:,f)  = (vWF*double(projWF(:,:,:,f))/mean(mean(mean(double(projWF(:,:,:,f))))));
    projSIM(:,:,:,f)  = (vSIM*double(projSIM(:,:,:,f))/mean(mean(mean(double(projSIM(:,:,:,f))))));
    if exist('wfDec','var')
        projWFDec(:,:,:,f)  = (vWFDec*double(projWFDec(:,:,:,f))/mean(mean(mean(double(projWF(:,:,:,f))))));
        projSIMDec(:,:,:,f)  = (vSIMDec*double(projSIMDec(:,:,:,f))/mean(mean(mean(double(projSIMDec(:,:,:,f))))));
    end
end
% save RGB stack
ppsSIM = 0.108/2;
writeRGBTIFF(255.*projWF,[fpath,fnameC,'_maxproj_wf'],[1/pps 1/pps])
writeRGBTIFF(255.*projSIM,[fpath,fnameC,'_maxproj_SIM'],[1/ppsSIM 1/ppsSIM])
if exist('wfDec','var')
    writeRGBTIFF(255.*projWFDec,[fpath,fnameC,'_maxproj_wfDec'],[1/pps 1/pps])
    writeRGBTIFF(255.*projSIMDec,[fpath,fnameC,'_maxproj_SIMDec'],[1/ppsSIM 1/ppsSIM])
end
