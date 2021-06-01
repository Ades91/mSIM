function [sim,wf] = getSIM(im,Nph,k0,ki,dk,Ndiff,doApod,doMask,figID)

if ~exist('figID','var') || isempty(figID); figID = 0; end
if ~exist('doMask','var') || isempty(doMask); doMask = 0; end
if ~exist('doApod','var') || isempty(doApod); doApod = 0; end
if ~exist('Ndiff','var') || isempty(Ndiff); Ndiff = 1; end

%%
% init variables
Nangle = length(Nph);
im = double(im);

x = linspace(-1,1,size(im,2));
y = linspace(-1,1,size(im,1));
Ny = size(im,1); Nx = size(im,2);
[X,Y] = meshgrid(x,y);
R = sqrt(X.^2 + Y.^2);
th = rad2deg(atan2(Y,X));
CTF = R < ki/2;
psf = abs(ifftshift(ifftn(ifftshift(CTF)))).^2;
OTF = R < k0;

mask = R < k0;

%% peak finding
ph = [];     dir = [];
for m = 1:Ndiff
    OTF = (m.*ki+dk/2-R)./(m.*ki+dk/2); OTF(OTF<0) =0;
    OTFinv = (OTF > 0).*(1-OTF);
    maskPeak = R > (m.*ki-dk/2) & R < (m.*ki+dk/2);
    maskUp =  th > -170 & th < 10;

    id = 1;
    for n = 1:Nangle
        for k = 1:Nph(n)
%         id = (k-1)*Nangle+n
        temp = apodImRect(im(:,:,id),20);id = id+1;
        T = (fftshift(fftn(fftshift(temp))));
        
        Ik = maskPeak.*maskUp.*OTFinv.*T;
    
        [cx,cy,an,A] = getSIMPeak(Ik);
        [dx,dy] = subPixelGauss(abs(T(cy-1:cy+1,cx-1:cx+1)));
        
        ph(k,n,m) = an;
        p1x(k,n,m) = cx+dx;
        p1y(k,n,m) = cy+dy;
%         A1(k,n,m) = (A).^2./max(abs(T(:))).^2;
         
%         xx = (cx-Nx/2)/(Nx/2); yy = (cy-Ny/2)/(Ny/2);
        Ik = (maskPeak.*not(maskUp).*OTFinv.*T);
        [dx,dy] = subPixelGauss(abs(T(cy-1:cy+1,cx-1:cx+1)));
        [cx,cy,an,A] = getSIMPeak(Ik);
  
        p2x(k,n,m) = cx+dx;
        p2y(k,n,m) = cy+dy;
        
%         A2(k,n,m) = abs(A).^2./max(abs(T(:))).^2;
        
        end
        if figID
            if m == 1
                figure(figID+1);subplot(ceil(Nangle/ceil(sqrt(Nangle))),ceil(sqrt(Nangle)),n)
                imagesc(log(abs(T)+1)); hold on
            else
                hold on
                subplot(ceil(Nangle/ceil(sqrt(Nangle))),ceil(sqrt(Nangle)),n)
            end
            plot(p1x(k,n,m),p1y(k,n,m),'kx');
            plot(p2x(k,n,m),p2y(k,n,m),'kx');
        end
    end
end
%% compute fourier transform
I = zeros(size(im,1),size(im,2),Nangle,max(Nph));
id = 1;
for n = 1:Nangle
    for k = 1:Nph(n)
        temp = apodImRect(im(:,:,id),20); id = id+1;
        I(:,:,n,k) = fftshift(fftn(fftshift(mean(im(:)).*temp./(mean(temp(:))))));
%         I(:,:,n,k) = fftshift(fftn(fftshift(temp)));
    end
end

%% recovering I0 Ik+ Ik-
if sum(Nph==5) > 0 && Ndiff == 2
    Ns = 5;
else
    Ns = 3;
end

S = zeros(size(im,1),size(im,2),Nangle,Ns);
t = tic;
for n = 1:Nangle
    if Nph(n) == 3 || Nph(n) == 5
        A = 1- sqrt((p2x(1,n,1)-size(im,2)/2-1).^2 + (p2y(1,n,1)-size(im,1)/2-1).^2)/(size(im,1)*k0);
        M = [ones(Nph(n),1), A.*ones(Nph(n),1).*exp(1i.*ph(1:Nph(n),n,1)), A.*ones(Nph(n),1).*exp(-1i.*ph(1:Nph(n),n,1))];
        if Ndiff == 2 && Nph(n) == 5
            A2 = 1 - sqrt((p2x(1,n,2)-size(im,2)/2-1).^2 + (p2y(1,n,2)-size(im,1)/2-1).^2)/(size(im,1)*k0);
            M = [M, A2.*ones(Nph(n),1).*exp(1i.*ph(1:Nph(n),n,2)), A2.*ones(Nph(n),1).*exp(-1i.*ph(1:Nph(n),n,2))];
        else
%             M = M(1:3,1:3);
        end

            %%
            temp = I(:,:,n,1:size(M,1));
            temp = permute(temp,[4 1 2 3]); temp = reshape(temp,[size(M,1) size(I,1)*size(I,2)]);
            t2 = mldivide(M,temp);
            t2 = reshape(t2,[size(t2,1) size(I,1) size(I,2)]); 
            S(:,:,n,:) = permute(t2,[2 3 1]);
            %%
%             for xx = 1:size(I,2)
%                 for yy = 1:size(I,1)
%                     S(yy,xx,n,:) = M\squeeze(I(yy,xx,n,1:size(M,1)));
%                 end
%             end
    else
        A = 1- sqrt((p2x(1,n,1)-size(im,2)/2-1).^2 + (p2y(1,n,1)-size(im,1)/2-1).^2)/(size(im,1)*k0);
        M = [ones(Nph(n),1), A.*ones(Nph(n),1).*exp(1i.*ph(1:Nph(n),n,1)), A.*ones(Nph(n),1).*exp(-1i.*ph(1:Nph(n),n,1))];
        
        for xx = 1:size(I,2)
            for yy = 1:size(I,1)
                S(yy,xx,n,:) = mldivide(M,squeeze(I(yy,xx,n,1:Nph(n))));
            end
        end
    end
end
toc(t)

if figID
    Stemp = [];
    for k = 1:Nangle
        Stemp = cat(3,Stemp,squeeze(S(:,:,k,:)));
    end
    plotStack(log(abs(Stemp)+1),figID+2)
end


%% reconstruction
% compute the mean of all unmixed I0
im0 = zeros(size(S,1),size(S,2));
% for n = 1:Nangle
%       im0 = im0 + abs(ifftshift(ifftn(ifftshift(S(:,:,n,1)))));
% end
im0 = mean(im,3);
R2 = fftshift(fftn(fftshift(apodImRect(im0,20))))./(Nangle);

% extend Fourier space according to k0 and measured shift
% K = sqrt((p1x(:)-N/2).^2 + (p1y(:)-N/2).^2); K = max(K(p1x~= 0))/(N/2);
% extend Fourier space by a factor 2
if doMask
    R2 = R2.*mask;
end

R2 = padarray(R2,[Ny/2 Nx/2]);
map = double(R2 ~= 0);
map2 = map;

for m = 1:Ndiff
    for n =  1:Nangle
        c = 1;
        pady = (size(R2,1) - size(S,1))/2;
        padx = (size(R2,2) - size(S,2))/2;
    
        t1 = S(:,:,n,2+2.*(m-1));
        t2 = S(:,:,n,3+2.*(m-1));
    
        if doMask
            t1 = mask.*t1;
            t2 = mask.*t2;
        end
    
        if doApod > 0 
            rApod = doApod;
            alpha = 1;
            apodMask = (cos(pi.*R./rApod)+1)./2;
            apodMask(R > rApod) = 0;
            apodMask = 1 - alpha.*apodMask;
            t1 = apodMask.*t1;
            t2 = apodMask.*t2;
        end
    
        t1 = padarray(t1,[pady padx]);
        t2 = padarray(t2,[pady padx]);

        t1x = p1x(1,n,m)+padx; t2x = p2x(1,n,m)+padx;
        t1y = p1y(1,n,m)+pady; t2y = p2y(1,n,m)+pady;

        ph1 = angle(t1(round(t1y),round(t1x)));
        ph2 = angle(t2(round(t2y),round(t2x)));
    
        t1 = imtranslate(t1,-c.*([t1x-1 t1y-1]-[size(t1,2) size(t1,1)]/2));
        t2 = imtranslate(t2,-c.*([t2x-1 t2y-1]-[size(t2,2) size(t2,1)]/2));
        
        temp = exp(1i.*-ph1).*t1 + exp(1i.*-ph2).*t2;
        
        map = double(t1~=0) + double(t2~=0);
        temp = temp./sqrt(map);
        temp(isnan(temp)) = 0;
        
        map2 =  map2 + double(temp~=0);
        R2 = R2 + temp;
        
    end
end

R2 = R2./sqrt(map2);
R2(isnan(R2)) = 0;

sim = real(ifftshift(ifftn(ifftshift(R2))));
wf = mean(im,3);

if figID
    figure(figID);imagesc(wf);colormap(gray)
    figure(figID+30+ceil(100.*doApod));imagesc(sim);colormap(gray)
end
