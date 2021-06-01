% coregister and crop the data using dx and dy 
function data2 = coregData(data,dx,dy,T,cx,x,cy,y)

data2 = zeros(y+1,x+1,size(data,3),8);
for k = 1:8
    disp(['Plane # ',num2str(k)])
    for n = 1:size(data,3)
        if k ~= 8
            data(:,:,n,k+1)= imtranslate(data(:,:,n,k+1),[dx(k),dy(k)]);
        end
        data2(:,:,n,k) = data(cy:cy+y,cx:cx+x,n,k)./T(k);
    end 
end
