% this function implements the gaussian approximation described in
% Fisher2001_A comparison of algorithms for subpixel peak detection

function [dx,dy] = subPixelGauss(data)

% data is a 3x3 matrix centered on the coarse peak finding
x0 = 2; y0 = 2;

dx = 0.5*(log(data(y0,x0-1)) - log(data(y0,x0+1)))./...
    (log(data(y0,x0-1)) - 2*log(data(y0,x0)) + log(data(y0,x0+1)));
if isnan(dx)
    dx = 0;
end
dy = 0.5*(log(data(y0-1,x0)) - log(data(y0+1,x0)))./...
    (log(data(y0-1,x0)) - 2*log(data(y0,x0)) + log(data(y0+1,x0)));
if isnan(dy)
    dy = 0;
end
