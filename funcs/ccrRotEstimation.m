function [dth,x,y] = ccrRotEstimation(im1,im2,Na,amax)
off = 10;
im1 = im1(off:end-off,off:end-off);
im1 = (im1-mean(im1(:)))/std(im1(:)); s1 = sum(abs(im1(:)).^2);

a = linspace(-amax,amax,Na);
dx = -2:2; dy = -2:2; % center of rotation shift

for k = 1:length(a)
    for xx = 1:length(dx)
        for yy = 1:length(dy)
            temp = im2;
            temp = imtranslate(temp,[dx(xx) dy(yy)],'bilinear'); % translate first
            temp = imrotate(temp,a(k),'bilinear','crop'); % then rotate
            temp = temp(off:end-off,off:end-off);
            temp = (temp-mean(temp(:)))/std(temp(:));
            s2 = sum(abs(temp(:)).^2);
%             cc(xx,yy,k) = sum(sum(im1.*temp))./sqrt(s1*s2);
            t = corrcoef(im1,temp); cc(xx,yy,k) = t(2,1);
        end
    end
end

[~,ind] = max(cc(:));
[cy,cx,cz] = ind2sub(size(cc),ind);
dth = a(cz);
x = dx(cx);
y = dy(cy);