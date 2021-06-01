function [K,A1] = checkSIMpeak(im,ki,dk,Nph,Ndiff)

% init variables
Nangle = length(Nph);
im = double(im);

x = linspace(-1,1,size(im,2)+1);x(end) = [];
y = linspace(-1,1,size(im,1)+1);y(end) = [];
Ny = size(im,1); Nx = size(im,2);
[X,Y] = meshgrid(x,y);

R = sqrt(X.^2 + Y.^2);

Nangle = length(Nph);
K = [];

%% peak finding
figure(345234534+1);
plot(0,0);delete(gca)
hold on
for m = 1:Ndiff
    % OTF = fftshift(fftn(fftshift(psf))); 
    OTF = (m.*ki+dk/2-R)./(m.*ki+dk/2); OTF(OTF<0) =0;
    OTFinv = (OTF > 0).*(1-OTF);
        maskPeak = R > (m.*ki-dk/2) & R < (m.*ki+dk/2);
    ph = []; dir = [];
    id = 1;
for n = 1:Nangle
    for k = 1%:Nph(n)
        id = (n-1)*Nph(n)+k
        temp = apodImRect(im(:,:,id),20);%id = id+1;
        T = fftshift(fftn(fftshift(temp)));
        
        Ik = maskPeak.*OTFinv.*T;
    
        [cx,cy,an,A] = getSIMPeak(OTFinv.*Ik);
        
        ph(k,n,m) = getSIMphase(Ik,cx,cy);
        dir(k,n) = atan2(cy-Ny/2,cx-Nx/2);
        
        p1x(k,n,m) = cx;
        p1y(k,n,m) = cy;
        A1(k,n,m) = (A)./max(abs(T(:)));
        
        xx = (cx-Nx/2)/(Nx/2); yy = (cy-Ny/2)/(Ny/2);
        mask2 = (X-xx).^2 + (Y-yy).^2 < 0.2^2;
        Ik = maskPeak.*not(mask2).*OTFinv.*T;
        [cx,cy,an,A] = getSIMPeak(Ik);
        
        p2x(k,n,m) = cx;
        p2y(k,n,m) = cy;
        A2(k,n,m) = abs(A)./max(abs(T(:)));
        

    end
    
        th = linspace(0,2*pi,100);
        rmax = m.*ki+dk/2;
        rmin = m.*ki-dk/2;
        if m == 1
            
            subplot(ceil(Nangle/ceil(sqrt(Nangle))),ceil(sqrt(Nangle)),n);
            imagesc(y,x,log(abs(T)+1)); hold on
            plot(rmin.*cos(th),rmin.*sin(th),'w');
            plot(rmax.*cos(th),rmax.*sin(th),'w');
            title(['A : ',num2str(mean(A1(:,n,m)))])
        else
            subplot(ceil(Nangle/ceil(sqrt(Nangle))),ceil(sqrt(Nangle)),n)
            hold on
            plot(rmin.*cos(th),rmin.*sin(th),'--w');
            plot(rmax.*cos(th),rmax.*sin(th),'--w');
        end

        
        plot(linmap(p1x(k,n,m),1,Nx,-1,1),linmap(p1y(k,n,m),1,Ny,-1,1),'kx');
        plot(linmap(p2x(k,n,m),1,Nx,-1,1),linmap(p2y(k,n,m),1,Ny,-1,1),'kx');

        kx= linmap(p1x(k,n,m),1,Nx,-1,1);
        ky= linmap(p1y(k,n,m),1,Nx,-1,1);
        K(n,m) = sqrt(kx.^2 + ky.^2);
end
end

        hold off
        
%         figure(234532);plot(ph(:,:,1))
%         
%         p1 = ph(:,1,1);
