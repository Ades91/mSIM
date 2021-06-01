function [sim,wf,param,pattern] = getSIMauto(im,Nph,k0,ki,Ndiff,doMask,doApod,figID)

if ~exist('figID','var') || isempty(figID); figID = 0; end
if ~exist('doMask','var') || isempty(doMask); doMask = 0; end
if ~exist('doApod','var') || isempty(doApod); doApod = 0; end
if ~exist('Ndiff','var') || isempty(Ndiff); Ndiff = 1; end

%% compute fourier transform
Nangle = length(Nph);
I = zeros(size(im,1),size(im,2),Nangle,max(Nph));
id = 1;
for n = 1:Nangle
    for k = 1:Nph(n)
%         temp = apodImRect(im(:,:,id),20); id = id+1;
        temp = im(:,:,id); id = id+1;
        I(:,:,n,k) = fftshift(fftn(fftshift(mean(im(:)).*temp./(mean(temp(:))))));
%         I(:,:,n,k) = fftshift(fftn(fftshift(temp)));
    end
end

%% recovering I0 Ik+ Ik- with unknown phase
if sum(Nph>=11) > 0 && Ndiff >= 5
    Ns = 11;
elseif sum(Nph>=9) > 0 && Ndiff >= 4
    Ns = 9; Ndiff = 4;
elseif sum(Nph>=7) > 0 && Ndiff >= 3
    Ns = 7; Ndiff = 3;
elseif sum(Nph>=5) > 0 && Ndiff >= 2
    Ns = 5; Ndiff = 2;
else
    Ns = 3; Ndiff = 1;
end

S = zeros(size(im,1),size(im,2),Nangle,Ns,class(I));
% t = tic;
for n = 1:Nangle
    if Nph(n) > 2 && mod(Nph(n),2) == 1
        ph = 0.01+linspace(0,2*pi,Nph(n)+1)'; ph(end) = [];
        
        A = 1;(1-ki/k0); % ideal OTF assumption
        M = [ones(Nph(n),1), A.*ones(Nph(n),1).*exp(1i.*ph(1:Nph(n))), A.*ones(Nph(n),1).*exp(-1i.*ph(1:Nph(n)))];
        % recursive construction of the unmixing matrix
        for k = 2:Ndiff
            if Nph(n) >= 2*k+1 % there are enough phases to perform the unmixing
                M = [M, ones(Nph(n),1).*exp(1i.*k.*ph(1:Nph(n))), ones(Nph(n),1).*exp(-1i.*k.*ph(1:Nph(n)))]; 
            end
        end

        %%
        temp = I(:,:,n,1:size(M,1));
        temp = permute(temp,[4 1 2 3]); temp = reshape(temp,[size(M,1) size(I,1)*size(I,2)]);
        t2 = mldivide(M,temp);
        t2 = reshape(t2,[size(t2,1) size(I,1) size(I,2)]); 
        S(:,:,n,:) = permute(t2,[2 3 1]);
        
    else
        ph = linspace(0,2*pi,Nph(n)+1)'; ph(end) = [];
        M = [ones(Nph(n),1), ones(Nph(n),1).*exp(1i.*ph(1:Nph(n))), ones(Nph(n),1).*exp(-1i.*ph(1:Nph(n)))];
        
        temp = I(:,:,n,1:Nph(n));
        temp = permute(temp,[4 1 2 3]); temp = reshape(temp,[size(M,1) size(I,1)*size(I,2)]);
        t2 = mldivide(M,temp);
        t2 = reshape(t2,[3 size(I,1) size(I,2)]); 
        S(:,:,n,:) = permute(t2,[2 3 1]);

    end
end
% toc(t)

% renormalize unmixed components to their noise level
x = (-size(im,2)/2:size(im,2)/2-1)./(size(im,2)/2);
y = (-size(im,1)/2:size(im,1)/2-1)./(size(im,1)/2);
% x = linspace(-1,1,size(im,2)); y = linspace(-1,1,size(im,1));
[X,Y] = meshgrid(x,y);
R = sqrt(X.^2 + Y.^2);
OTFinv = (1-(1-R./k0)).*(R<k0); OTF = (1-R./k0).*(R<k0);
mask = R < k0;
for n = 1:Nangle
   temp = squeeze(S(:,:,n,1));
   base =  mean(mean(abs(temp(not(mask)))));
   for k = 2:min(Nph(n),size(M,2))
       temp = S(:,:,n,k);
       w(n,k) = base/(mean(mean(abs(temp(not(mask))))));
       S(:,:,n,k) = S(:,:,n,k).*w(n,k);
   end
end
if figID
    Stemp = [];
    for k = 1:Nangle
        Stemp = cat(3,Stemp,squeeze(S(:,:,k,:)));
    end
%     plotStack(log(abs(Stemp)+1),figID+20)
end

%% compute componant cross-correlation
p2x = []; p2y = []; ph_offset = []; shifts = [];
cc = [];

for n = 1:length(Nph)
    for k = 2:3%size(S,4)
        cc = fftshift(fftn(fftshift(ifftshift(ifftn(ifftshift(OTFinv.*S(:,:,n,1)))).*...
            conj(ifftshift(ifftn(ifftshift(OTFinv.*S(:,:,n,k))))))));

        if k < 4
            mask = double(R >= 0.9*ki & R <= 1.1*ki);
        elseif k < 6
            mask = double(R > 0.8*2*ki & R < 1.2*2*ki);
        else
            mask = double(R > 0.8*3*ki & R < 1.2*3*ki);
        end
        [cx,cy,~,~] = getSIMPeak(mask.*abs(cc));
        [dx,dy] = subPixelGauss(abs(cc(cy-1:cy+1,cx-1:cx+1)));
        p2x(n,k) = cx+dx;
        p2y(n,k) = cy+dy;

%         ph_offset(n,k) = angle(S(cy,cx,n,k));
%         ph2(n,k)= angle(cc(cy,cx,n));
%         A(n) = abs(S(cy,cx,n,5))./abs(S(size(S,1)/2,size(S,2)/2,n,1));
    end
    [cx,cy,~,~] = size(S);cx = cx/2 +1; cy = cy/2 +1;
    for k = 2:Ndiff
        p2x(n,2*k) = k*(p2x(n,2)-cx) + cx;
        p2x(n,2*k+1) = k*(p2x(n,3)-cx) + cx;
        p2y(n,2*k) = k*(p2y(n,2)-cy) + cy;
        p2y(n,2*k+1) = k*(p2y(n,3)-cy) + cy;
    end
end

param.p2x = p2x; param.p2y = p2y;
param.angle = atan2(p2y(:,2)-cy,p2x(:,2)-cx);

%% pattern phase estimation from Wicker2013 
phases = [];
for n = 1:Nangle
    for k = 1:Nph(n)
        m = 1;
%         phases(n,k) = getSIMphase(I(:,:,n,k),round(p2x(n,3+(m-1)*2)),round(p2y(n,3+(m-1)*2)));
%     phases(n,:) = sort(phases(n,:),'ascend');
        r0 = 0.25;
        mask = sqrt((X - (p2x(n,3+(m-1)*2)-cx)/size(im,2)).^2 + (Y- (p2y(n,3+(m-1)*2)-cx)/size(im,1)).^2) < r0 | ...
            sqrt((X + (p2x(n,3+(m-1)*2)-cx)/size(im,2)).^2 + (Y + (p2y(n,3+(m-1)*2)-cx)/size(im,1)).^2) < r0;

        phases(k,n) = getSIMphase(I(:,:,n,k).*mask,round(p2x(n,3)),round(p2y(n,3)));

    end
%     phases(n,:) = sort(phases(n,:),'ascend');
        phases(:,n) = unwrap(phases(:,n));
      phases(:,n) = phases(:,n) - min(phases(:,n));
end

param.phases = phases;
phases
%%
S = zeros(size(im,1),size(im,2),Nangle,Ns,class(I));
% t = tic;
for n = 1:Nangle
    if Nph(n) > 2 && mod(Nph(n),2) == 1
        ph = phases(:,n);
        
        A = 1;(1-ki/k0); % ideal OTF assumption
        M = [ones(Nph(n),1), A.*ones(Nph(n),1).*exp(1i.*ph(1:Nph(n))), A.*ones(Nph(n),1).*exp(-1i.*ph(1:Nph(n)))];
        % recursive construction of the unmixing matrix
        for k = 2:Ndiff
            if Nph(n) >= 2*k+1 % there are enough phases to perform the unmixing
                M = [M, ones(Nph(n),1).*exp(1i.*k.*ph(1:Nph(n))), ones(Nph(n),1).*exp(-1i.*k.*ph(1:Nph(n)))]; 
            end
        end

        %%
        temp = I(:,:,n,1:size(M,1));
        temp = permute(temp,[4 1 2 3]); temp = reshape(temp,[size(M,1) size(I,1)*size(I,2)]);
        t2 = mldivide(M,temp);
        t2 = reshape(t2,[size(t2,1) size(I,1) size(I,2)]); 
        S(:,:,n,:) = permute(t2,[2 3 1]);
        

    else
        ph = phases(n,:);
        M = [ones(Nph(n),1), ones(Nph(n),1).*exp(1i.*ph(1:Nph(n))), ones(Nph(n),1).*exp(-1i.*ph(1:Nph(n)))];
        
        temp = I(:,:,n,1:Nph(n));
        temp = permute(temp,[4 1 2 3]); temp = reshape(temp,[size(M,1) size(I,1)*size(I,2)]);
        t2 = mldivide(M,temp);
        t2 = reshape(t2,[3 size(I,1) size(I,2)]); 
        S(:,:,n,:) = permute(t2,[2 3 1]);

    end
end
% toc(t)
% plotStack(log(abs(S)+1),figID+20)
% renormalize unmixed components to their noise level
x = linspace(-1,1,size(im,2)); y = linspace(-1,1,size(im,1));
[X,Y] = meshgrid(x,y);
R = sqrt(X.^2 + Y.^2);
TH = atan2(Y,X);
mask = R < k0;
for n = 1:Nangle
   temp = squeeze(S(:,:,n,1));
   base =  mean(mean(abs(temp(not(mask)))));
   for k = 2:min(Nph(n),size(M,2))
       temp = S(:,:,n,k);
       w(n,k) = base/(mean(mean(abs(temp(not(mask))))));
       S(:,:,n,k) = S(:,:,n,k).*w(n,k);
   end
end

if figID
    Stemp = [];
    for k = 1:Nangle
        Stemp = cat(3,Stemp,squeeze(S(:,:,k,:)));
    end
    plotStack(log(abs(Stemp)+1),figID+21)
end


%% shifts = sqrt((p2y-size(im,1)/2-1).^2 + (p2x-size(im,2)/2-1).^2);
if nargout == 4
    [Ny,Nx,~] = size(I);
    pattern = zeros(2*size(I,1),2*size(I,2),sum(Nph));
    for n = 1:Nangle
        for k = 1:Nph(n)
            temp = zeros(size(I,1),size(I,2));
            temp(round(p2y(n,2)),round(p2x(n,2))) = I(round(p2y(n,2)),round(p2x(n,2)),n,k);
%             temp(round(p2y(n,3)),round(p2x(n,3))) = I(round(p2y(n,3)),round(p2x(n,3)),n,k);
            temp = padarray(temp,[Ny/2 Nx/2]);
            pattern(:,:,k+Nph(n)*(n-1)) = linmap(real(ifftshift(ifftn(ifftshift(temp)))),0,1);
        end
    end
end

%% reconstruction
% compute the mean of all unmixed I0
im0 = zeros(size(S,1),size(S,2));
for n = 1:Nangle
      im0 = im0 + abs(ifftshift(ifftn(ifftshift(S(:,:,n,1)))));
end
% im0 = mean(im,3);
R2 = fftshift(fftn(fftshift(apodImRect(im0,20))))./(Nangle);


mask = R < k0;
if doMask
    R2 = R2.*mask;
end
Ny = size(im,1); Nx = size(im,2);

% extend Fourier space according to k0 and measured shift
K = sqrt((p2x(:)-Nx/2).^2 + (p2y(:)-Ny/2).^2); K = max(K(p2x(:)~= 0));
R2 = padarray(R2,round(1.5.*[K K]));

% map = double(R2 ~= 0);
map0 = double(R2~=0);

[~,ind] = max(abs(R2(:)));
[y0,x0] = ind2sub(size(R2),ind);

for m = 1:Ndiff
    for n =  1:Nangle
        
        pady = (size(R2,1) - size(S,1))/2;
        padx = (size(R2,2) - size(S,2))/2;

        t1 = S(:,:,n,2+2.*(m-1));
        t2 = S(:,:,n,3+2.*(m-1));
    
        if doMask
            %maskTh = TH > mod(param.angle(n)-0.1+pi,2*pi)-pi & TH < mod(param.angle(n)+0.1+pi,2*pi)-pi;
            % t1 = (conv2(single(mask),single(maskTh),'same')>0).*t1;
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

        cx1 = p2x(n,2+2.*(m-1)); cy1 = p2y(n,2+2.*(m-1));
        cx2 = p2x(n,3+2.*(m-1)); cy2 = p2y(n,3+2.*(m-1));
        
        cx1 = (cx1-size(im,2)/2-1);
        cy1 = (cy1-size(im,1)/2-1);
        cx2 = (cx2-size(im,2)/2-1);
        cy2 = (cy2-size(im,1)/2-1);
        
%         ph1 = m.*(ph_offset(n));
        
        t1 = imtranslate(t1,[cx1 cy1]);
        t2 = imtranslate(t2,[cx2 cy2]);
        ph1 = angle(t1(y0,x0));
        ph2 = angle(t2(y0,x0));
        

        temp = exp(1i.*-ph1).*t1 + exp(1i.*-ph2).*t2;
%         figure;imagesc(abs(temp))
        map = double(t1~=0) + double(t2~=0);
        temp = temp./sqrt(map);
        temp(isnan(temp)) = 0;
        map0 = map0 + double(temp~=0) ;
        R2 = R2 + temp;
        
        k(m,n,1) = 2*sqrt(cx1.^2 + cy1.^2)./size(im,1); k(m,n,2) = 2*sqrt(cx2.^2 + cy2.^2)./size(im,1);
        kph(m,n,1) = ph1; kph(m,n,2) = ph2;
        an(m,n,1) = atan2(cy1,cx1); an(m,n,2) = atan2(cy2,cx2);
    end
end

% store reconstruction parameters
param.k = k; param.kph = kph; param.an = an;

R2 = R2./sqrt(map0);
R2(isnan(R2)) = 0;
sim = real(ifftshift(ifftn(ifftshift(R2))));

wf = mean(im,3);

if figID
    figure(figID+30+ceil(100.*doApod))
    p = get(figID+30+ceil(100.*doApod),'position');
    set(figID+30+ceil(100.*doApod),'position',[p(1) p(2) 1000 p(4)]);
    subplot(121);imagesc(sim);colormap(gray);title('SIM')
    subplot(122);imagesc(im0);colormap(gray);title('Pseudo WF')
end
