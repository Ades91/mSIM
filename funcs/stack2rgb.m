function proj = stack2rgb(stack,cmap)



cmap = imresize(cmap,[size(stack,3) 3],'bilinear');
[~,map] = max(stack,[],3);

projR = zeros(size(stack,1),size(stack,2));
projG = zeros(size(stack,1),size(stack,2));
projB = zeros(size(stack,1),size(stack,2));

for k = 1:max(map(:))
    mask = map == k;
    slice = stack(:,:,k);
   projR(mask) = projR(mask) + slice(mask).*cmap(k,1);
   projG(mask) = projG(mask) + slice(mask).*cmap(k,2);
   projB(mask) = projB(mask) + slice(mask).*cmap(k,3);
    
end

proj = cat(3,projR,projG,projB);