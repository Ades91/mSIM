% input has the shape y,x,z,t

function [sim,simDec] = resampleSIM(sim,s,simDec)
%%
R2 = fftshift(fftn(fftshift(s)));
cropX = floor(sum(abs(R2(end/2,:))< 1)./2);
cropY = floor(sum(abs(R2(:,end/2))< 1)./2);
if cropX > cropY
    N = size(sim,2)-2*cropX;
    N = (size(sim,1)/size(sim,2))*N;
    cropY = round((size(sim,1)-N)/2);
else
    N = size(sim,1)-2*cropY;
    N = (size(sim,2)/size(sim,1))*N;
    cropX = round((size(sim,2)-N)/2);
end

%%
simResamp = zeros(size(sim,1)-2*cropY,size(sim,2)-2*cropX,size(sim,3),size(sim,4));
if exist('simDec','var')
    simDecResamp = zeros(size(sim,1)-2*cropY,size(sim,2)-2*cropX,size(sim,3),size(sim,4));
end

for f = 1:size(sim,4)
    for k = 1:size(sim,3)
    % make sim square
        s = sim(:,:,k,f);
        R2 = fftshift(fftn(fftshift(s)));
        R3 = R2(cropY+1:end-cropY,cropX+1:end-cropX);
        simResamp(:,:,k,f) = real(ifftshift(ifftn(ifftshift(R3))));
        if exist('simDec','var')
        	s = simDec(:,:,k,f);
            R2 = fftshift(fftn(fftshift(s)));
            R3 = R2(cropY+1:end-cropY,cropX+1:end-cropX);
            simDecResamp(:,:,k,f) = real(ifftshift(ifftn(ifftshift(R3))));
        end
            
    end
end

sim = simResamp;
if exist('simDec','var')
    simDec = simDecResamp;
end