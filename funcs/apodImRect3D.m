% apodize the input image in with a cosine edge mask over a length of N
% pixels
function out = apodImRect3D(in,N)

Nx = min(size(in,1),size(in,2));
Nz = size(in,3);

x = abs(linspace(-Nx/2,Nx/2,Nx));
map = x > Nx/2 - N(1);
z = abs(linspace(-Nz/2,Nz/2,Nz));
mapZ = z > Nz/2 - N(end);

val = mean(mean(in(:)));

d = (-abs(x)- mean(-abs(x(map)))).*map;
d = linmap(d,-pi/2,pi/2);
d(not(map)) = pi/2;
mask = (sin(d)+1)/2;

dz = (-abs(z)- mean(-abs(z(mapZ)))).*mapZ;
dz = linmap(dz,-pi/2,pi/2);
dz(not(mapZ)) = pi/2;
temp = (sin(dz)+1)/2;

maskZ = zeros(1,1,length(temp));
maskZ(1,1,:) = temp;

% make it 2D
if size(in,1) > size(in,2)
    mask = mask.*imresize(mask',[size(in,1) 1],'bilinear');
elseif size(in,1) < size(in,2)
    mask = imresize(mask,[1 size(in,2)],'bilinear').*mask';
else
    mask = mask.*mask';
end
	
mask = repmat(mask,[1 1 length(temp)]).*repmat(maskZ,[Nx,Nx,1]);

out = (in-val).*mask + val;