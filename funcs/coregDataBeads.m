% coregister and crop the data using dx and dy 
function data2 = coregDataBeads(data,cal)
y = cal.cy1-cal.cy0+1;
x = cal.cx1-cal.cx0+1;
data2 = zeros(y,x,size(data,3),8);
for ch = 1:8
    for t = 1:size(data,3)
        temp = data(:,:,t,ch);
        for n = ch:-1:2
            T = cal.tf{n-1}.tdata.T;
            temp = imtransform(temp,maketform('affine',T), 'XData', [1 size(temp,2)],...
                'YData', [1 size(temp,1)],'Size', size(temp));
        end
        data2(:,:,t,ch) = temp(cal.cy0:cal.cy1,cal.cx0:cal.cx1);
    end
end
