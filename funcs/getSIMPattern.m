% Estimate sim pattern from data
% Author : Adrien Descloux, 2019 07 01
% inputs 
%   im      : input sim sequence 
%   Nph     : number of phases per angles (e.g. [5 5 5] means 3 angles with
%               5 phases each
%   k0      : normalized cut-off frequency (frequency support of the OTF)
%   ki      : normalized illumination frequency
%   sig     : size of the Gaussian mask used to filter the peaks
%   figID   : if supplied and >0 display results in figure(figID)

% outputs
%   pattern : estimated SIM pattern per angles and phases

function pattern = getSIMPattern(im,Nph,ki,sig,figID)

%im=im-mean(im,3);
ki=norm(ki./size(im(:,:,1)));

%% compute Fourier transform
Nangle = length(Nph);
I = zeros(size(im,1),size(im,2),Nangle,max(Nph));
id = 1;
for n = 1:Nangle
    for k = 1:Nph(n)
        temp = apodImRect(im(:,:,id),20); id = id+1;
        I(:,:,n,k) = fftshift(fftn(fftshift(temp)));
    end
end

%% compute ideal OTF
x = linspace(-0.5,0.5,size(im,2)); y = linspace(-0.5,0.5,size(im,1));
[X,Y] = meshgrid(x,y);
R = sqrt(X.^2 + Y.^2);
OTFinv = 1;% (1-(1-R./k0)).*(R<k0);
OTF = 1;%(1-R./k0).*(R<k0);

%% compute Fourier domain cross-correlations and extract peaks
p1x = zeros(length(Nph),size(I,4)); p1y = p1x;
p2x = p1x; p2y = p1x;
pattern = zeros(size(im)*2);
x = 1:size(im,2); y = 1:size(im,1);
[X,Y] = meshgrid(x,y);
count = 1;
for n = 1:length(Nph)
    for k = 2:size(I,4)
        cc = fftshift(fftn(fftshift(ifftshift(ifftn(ifftshift(OTFinv.*I(:,:,n,1)))).*...
            conj(ifftshift(ifftn(ifftshift(OTFinv.*I(:,:,n,k))))))));
        if k < 4
            mask = double(R > 0.7*ki & R < 1.3*ki);
        else
            mask = double(R > 0.7*2*ki & R < 1.3*2*ki);
        end
        
        % get SIM peak from the cross-correlation
        [cx,cy,~,~] = getSIMPeak(mask.*abs(cc));
        [dx,dy] = subPixelGauss(abs(cc(cy-1:cy+1,cx-1:cx+1)));
        p1x(n,k) = cx+dx; p1y(n,k) = cy+dy;
        r = sqrt((X-p1x(n,k)).^2 + (Y-p1y(n,k)).^2)  > 10;
        
        % mask previously found peak and look for the other one
        [cx,cy,~,~] = getSIMPeak(mask.*abs(cc).*r);
        [dx,dy] = subPixelGauss(abs(cc(cy-1:cy+1,cx-1:cx+1)));
        p2x(n,k) = cx+dx; p2y(n,k) = cy+dy;
        
        count = count +1;
    end
end
% pattern phase estimation from Wicker2013 => phases is the correct matrix M 
% for n = 1:Nangle
%     for k = 1:Nph(n)
%         phases(n,k) = getSIMphase(I(:,:,n,k).*OTF,round(p2x(n,3)),round(p2y(n,3)));
%     end
% %     phases(n,:) = sort(phases(n,:),'ascend');
% %     phases(n,:) = phases(n,:);
% end

%% generate pattern by filtering the images in Fourier space with two Gaussian
count = 1;
padPre=[ceil(size(im,1)/2),ceil(size(im,2)/2)];
padPost=[floor(size(im,1)/2),floor(size(im,2)/2)];
for n = 1:length(Nph)
    for k = 1:Nph(n)
        g = exp(-(sqrt((X-round(p1x(n,2))).^2 + (Y-round(p1y(n,2))).^2)/(sig)).^2);
        g = g + exp(-(sqrt((X-round(p2x(n,2))).^2 + (Y-round(p2y(n,2))).^2)/(sig)).^2);
        %g = g + exp(-(sqrt((X-floor(size(im,2)/2)-1).^2 + (Y-floor(size(im,1)/2)-1).^2)/(sig)).^2);
        tmp=padarray(padarray(g.*I(:,:,n,k),padPre,'pre'),padPost,'post'); % Sinc interpolation (we need to reconstruct on a twice finer grid)
        pattern(:,:,count)=linmap(real(ifftshift(ifftn(ifftshift(tmp)))),0,1);
        %pattern(:,:,count)=real(ifftshift(ifftn(ifftshift(tmp))));
        count = count +1;
    end
end
pattern=pattern(:,:,1:size(im,3));
pattern=pattern/(mean(pattern(:))*sum(Nph));

% Ndiff=1;
% %% generate pattern by direct rendering it in real space on 2x upsampled grid
% [X,Y] = meshgrid(-(size(im,2)-0.5)/2:0.5:(size(im,2)-0.5)/2,-(size(im,1)-0.5)/2:0.5:(size(im,1)-0.5)/2);
% [X,Y] = meshgrid(0:0.5:size(im,2)-1,0:0.5:size(im,1)-1);
% count = 1; pattern = zeros(size(X,1),size(X,2),size(im,3));
% Ny = size(im,1); Nx = size(im,2);
% cx = Nx/2 +1; cy = Ny/2 +1; % Fourier space center for even size images
% for n = 1:length(Nph)
%     for k = 1:Nph(n)
%         kx = pi*(p1x(n,2,1)-cx)/(Nx/2);
%         ky = pi*(p1y(n,2,1)-cy)/(Ny/2);
%         if ky < 0 % always pick ky < 0 for the pattern wavevector
%             kx = pi*(p1x(n,3,1)-cx)/(Nx/2);
%             ky = pi*(p1y(n,3,1)-cy)/(Ny/2);
%         end
%         if Ndiff == 1
%             temp = (cos(kx.*X + ky.*Y-phases(n,k))+1)/2;
%         else
%             temp = ((cos(kx.*X + ky.*Y-phases(n,k))+1)/2).^2;
%         end
%         tt = linmap(temp,0,1);
%         pattern(:,:,count)=tt/(mean(tt(:))*sum(Nph)); % Normalization such that the mean of each pattern is 1/#Patterns
%         count = count+1;
%     end
% end

%% plotting if figID > 0
if exist('figID','var')
    if figID > 0
        nn = ceil(sqrt(size(im,3)));
        figure(figID); count = 1;
        for n = 1:length(Nph)
            for k = 1:Nph(n)
                subplot(nn,nn,count)
                imagesc(log(abs(mask.*OTFinv.*I(:,:,n,k))+1)); hold on
                plot(p1x(n,2), p1y(n,2),'kx');
                plot(p2x(n,2), p2y(n,2),'kx'); hold off
                count = count +1;
            end
        end
        
    end
end
