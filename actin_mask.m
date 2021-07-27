function [mask Ccell_cent Carea Car Ccirc Cang coords]=actin_mask(I)

sigma=6;hsize=[12,12];
h=fspecial('gaussian',hsize,sigma);
I2=imfilter(I,h,'replicate');
%figure(1); imshow(I2)
BW=edge(I2,'sobel',0.001);
%figure(2);imshow(BW)

se=strel('disk',6,0);
mask=imdilate(BW,se);
mask=imfill(mask,'holes');
mask=imerode(mask,se);
%figure(3);imshow(mask)

cc=bwconncomp(mask);
r=regionprops(cc,'Centroid','Area','Perimeter','MajorAxisLength','MinorAxisLength','PixelList','Orientation');

for i=1:cc.NumObjects
    area(i)=r(i).Area;
end

[~,ind]=max(area);

for i=1:cc.NumObjects
    if i~=ind
        mask(cc.PixelIdxList{i})=0;
    end
end

Carea=area(ind);
Ccell_cent=r(ind).Centroid;
Car=r(ind).MajorAxisLength/r(ind).MinorAxisLength;
Ccirc=(r(ind).Area*4*pi)/(r(ind).Perimeter^2);
Cang=r(ind).Orientation;
x=r(ind).PixelList(1,1);
y=r(ind).PixelList(1,2);
coords=bwtraceboundary(mask,[y x],'N');

% imshow(image)
% hold on
% plot(coords(:,2),coords(:,1),'m')