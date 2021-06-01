
function writeTIFF(data, filename,res)

if nargin == 2
    res = [1 1]; % save generic resolution if none provided
end

%% save to tiff while keeping floating value

% writes data as a multi-channel TIFF with single prec. float pixels
filename = [filename '.tif'];

    for k = 1:size(data,3)
        t = Tiff(filename, 'a');
        tagstruct.ImageLength = size(data, 1);
        tagstruct.ImageWidth = size(data, 2);
        tagstruct.Compression = Tiff.Compression.None;
        tagstruct.XResolution = res(1);
        tagstruct.YResolution = res(2); % YResolution encode z-res
        tagstruct.ResolutionUnit = 1;
        tagstruct.SampleFormat = Tiff.SampleFormat.IEEEFP;
        tagstruct.Photometric = 1;
        tagstruct.BitsPerSample =  32;                        % float data
        tagstruct.SamplesPerPixel = 1;
        tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
        tagstruct.ImageDescription = 'unit=micron';
        t.setTag(tagstruct);
        t.write(single(data(:,:,k)));
        t.close();
    end

end