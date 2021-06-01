function im = resampleImage(im,crop)

Ik = fftshift(fftn(fftshift(apodImRect(im,20))));

im = real(ifftshift(ifftn(ifftshift(Ik(crop+1:end-crop,crop+1:end-crop)))));