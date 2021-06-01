% this function implements the method described in
% Wicker2013_Non-iterative determination of pattern phase in structured illumination 
% microscopy using auto-correlations in Fourier space

% input, Fourier reweigthed transform of raw image, cx,cy, position of the peak in pixels
function ph = getSIMphase(I,cx,cy)

D = ifftshift(ifftn(ifftshift(fftshift(fftn(fftshift(I))).*conj(fftshift(fftn(fftshift(I)))))));

ph= angle(D(cy,cx));



