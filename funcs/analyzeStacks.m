
% load im

cdir = cd;
try
cd('D:\Measurements')
end
[fname,fpath] = uigetfile({'*.*';'*.All'});
cd(cdir);

im = loadData([fpath,fname]);

im(:,:,1) = []; 

plotStack(im)

%% check and make data is even before computation
im = im(1:end-(mod(size(im,1),2)),1:end-(mod(size(im,2),2)),:);

ki = 0.4;
dk = 0.2;
Nph = [3 3 3];
% init variables
im = double(im);

x = linspace(-1,1,size(im,2));
y = linspace(-1,1,size(im,1));
Ny = size(im,1); Nx = size(im,2);
[X,Y] = meshgrid(x,y);

R = sqrt(X.^2 + Y.^2);
% OTF = fftshift(fftn(fftshift(psf))); 
OTF = (ki+dk/2-R)./(ki+dk/2); OTF(OTF<0) =0;
OTFinv = (OTF > 0).*(1-OTF);

Nangle = length(Nph);
K = [];

%% peak finding
    maskPeak = R > (ki-dk/2) & R < (ki+dk/2);
    ph = []; dir = [];

for k = 1:size(im,3)
%         id = (k-1)*Nangle+n
        temp = apodImRect(im(:,:,k),20);
        m(k) = mean(temp(:))./mean(im(:));
        s(k) = std(temp(:));
        temp = temp.*s(k);
        T = fftshift(fftn(fftshift(temp)));
        
        Ik = maskPeak.*OTFinv.*T;
        
        [cx,cy,an,A] = getSIMPeak(Ik);
        
        ph(k) = an;
        dir(k) = atan2(cy-Ny/2,cx-Nx/2);
        A1(k) = (A);%.^2./max(abs(T(:))).^2;
        
        p1x(k) = cx;
        p1y(k) = cy;
        
        xx = (cx-Nx/2)/(Nx/2); yy = (cy-Ny/2)/(Ny/2);
        mask2 = (X-xx).^2 + (Y-yy).^2 < 0.2^2;
        Ik = maskPeak.*not(mask2).*OTFinv.*T;
        [cx,cy,an,A] = getSIMPeak(Ik);
        
        p2x(k) = cx;
        p2y(k) = cy;
        A2(k) = abs(A).^2./max(abs(T(:))).^2;
        

%         figure(345234534+1);subplot(ceil(Nangle/ceil(sqrt(Nangle))),ceil(sqrt(Nangle)),n)
%         imagesc(log(abs(T)+1)); hold on
% %         imagesc(angle(T)); hold on
%         plot(p1x(k,n),p1y(k,n),'kx');
%         plot(p2x(k,n),p2y(k,n),'kx');hold off
        
        K(k) = sqrt((p1x(k)-Nx/2).^2 + (p1y(k)-Ny/2).^2)/((Nx+Ny)/4);
end


%%

figure(3543)
subplot(221);plot(m)
subplot(222);plot(A1)
subplot(223);plot(s)
subplot(224);plot(ph)