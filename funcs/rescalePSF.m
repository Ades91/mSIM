function psf = rescalePSF(psf,FOVpsf,FOVdat,Ndat)
% return the 3D psf rescaled and resampled to the data Field of view

if length(FOVpsf) == 2 % expand the inputs to 3D
    FOVpsf(3) = FOVpsf(2); FOVpsf(2) = FOVpsf(1);
    FOVdat(3) = FOVdat(2); FOVdat(2) = FOVdat(1);
    Ndat(3) = Ndat(2); Ndat(2) = Ndat(1);
end
% number of pixels
Npsf(1) = size(psf,1);
Npsf(2) = size(psf,2);
Npsf(3) = size(psf,3);

%%%%%%%%%% FOV = [FOVxy FOVz] %% Experimental FOV
%%%%%%%%%% NE  = [NxE   NzE ] %% Experimental number of pixels

% rescale/resample xy and z in a separate way

% rescale dim 1
if FOVpsf(1) <= FOVdat(1)
    disp('Scale the psf down (XY)')
    % scales psf down to match sample size
    scaling_x = FOVpsf(1)*Ndat(1)/FOVdat(1);
    
    % resample the psf
	[y x z] = ndgrid(linspace(1,size(psf,1),scaling_x),...
                 linspace(1,size(psf,2),size(psf,2)),...
                 linspace(1,size(psf,3),size(psf,3)));
    psfI = interp3(abs(psf),x,y,z);
%     psfP = interp3(angle(psf),x,y,z);

    % then pad everything with 0 to match meas size
    padVal = 0;round(abs(mean(psfI(:))));
    psfI = padarray(psfI,[round((Ndat(1)-scaling_x)/2) 0 0],padVal);
%     psfP = padarray(psfP,[round((Ndat(1)-scaling_x)/2) 0 0],padVal);

else % need to crop the psf
    disp('Scale the psf up (XY)')
    offset = (Npsf(1)*FOVdat(1)/FOVpsf(1));
    offset = round((Npsf(1)-offset)/2);
    gridX = offset:Npsf(1)-offset+1;
    
    temp = psf(gridX,:,:);
    psf = temp;

end

psf = psfI;%.*exp(1i.*psfP);

% rescale dim 2
if FOVpsf(2) <= FOVdat(2)
    disp('Scale the psf down (XY)')
    % scales psf down to match sample size
    scaling_y = FOVpsf(2)*Ndat(2)/FOVdat(2);
    
    % resample the psf
    [y x z] = ndgrid(linspace(1,size(psf,1),size(psf,1)),...
                 linspace(1,size(psf,2),scaling_y),...
                 linspace(1,size(psf,3),size(psf,3)));
    psfI = interp3(abs(psf),x,y,z);
%     psfP = interp3(angle(psf),x,y,z);

    % then pad everything with 0 to match meas size
    padVal = 0;round(abs(mean(psfI(:))));
    psfI = padarray(psfI,[0 round((Ndat(2)-scaling_y)/2) 0],padVal);
%     psfP = padarray(psfP,[0 round((Ndat(2)-scaling_y)/2) 0],padVal);

else % need to crop the psf
    disp('Scale the psf up (XY)')
    offset = (Npsf(2)*FOVdat(2)/FOVpsf(2));
    offset = round((Npsf(2)-offset)/2);
    gridX = offset:Npsf(2)-offset+1;
    
    temp = psf(:,gridX,:);
    psf = temp;

end

psf = psfI;%.*exp(1i.*psfP);

% repeat the same procedure for Z
if FOVpsf(3) <= FOVdat(3)
    
    disp('Scale the psf down (Z)')
    % scales psf down to match sample size
    scaling_z = FOVpsf(3)*Ndat(3)/FOVdat(3);
    
    % resample the psf
    [y x z] = ndgrid(linspace(1,size(psf,1),size(psf,1)),...
                 linspace(1,size(psf,2),size(psf,2)),...
                 linspace(1,size(psf,3),scaling_z));
    psfI = interp3(abs(psf),x,y,z);
%     psfP = interp3(angle(psf),x,y,z);

    % then pad everything with 0 to match meas size
    psfI = padarray(psfI,[0 0 round((Ndat(3)-scaling_z)/2)],round(abs(mean(psfI(:)))));
%     psfP = padarray(psfP,[0 0 round((Ndat(3)-scaling_z)/2)],round(abs(mean(psfP(:)))));
       
else % need to crop the psf
    
    disp('Scale the psf up (Z)')
    offset = (Npsf(3)*FOVdat(3)/FOVpsf(3));
    offset = round((Npsf(3)-offset)/2);
    gridZ = offset:Npsf(3)-offset+1;
    
    temp = psf(:,:,gridZ);
    psf = temp;
 
end

psf = psfI;%.*exp(1i.*psfP);

% and correct rounding error
[y x z] = ndgrid(linspace(1,size(psf,1),Ndat(1)),...
                 linspace(1,size(psf,2),Ndat(2)),...
                 linspace(1,size(psf,3),Ndat(3)));
psfI = interp3(abs(psf),x,y,z);
% psfP = interp3(angle(psf),x,y,z);

psf = (psfI);%.*exp(1i.*psfP);
