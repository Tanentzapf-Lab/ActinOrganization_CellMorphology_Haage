function []=actin_extravaganza_max(file)

folder2=dir([file '/*.tif']);
for i=2:length(folder2)
    temp=imread(fullfile([file '/' folder2(i).name]));
    images(:,:,i-1)=temp(:,:,1); %% change 1 here to pick R=1, G=2
end
image=max(images,[],3); clear images
image=im2double(image);
image=image(10:end-10,10:end-10);
[ys xs]=size(image);
se=strel('disk',12);
im2=imtophat(image,se);
[mask Ccell_cent Carea Car Ccirc Cang coords]=actin_mask(image);

w=32;
x=w+1:w/2:xs-w;
y=w+1:w/2:ys-w;

sigma=w/2;
hsize=[w*2,w*2];
h=fspecial('gaussian',hsize,sigma);

w1=cos(linspace(-pi/2,pi/2,w*2));
w2=cos(linspace(-pi/2,pi/2,w*2));
window=w1'*w2;

l=5;
figure(1);imshow(image);hold on
plot(coords(:,2),coords(:,1),'w')
stdev=zeros(length(y),length(x));
theta=zeros(length(y),length(x));
ar=zeros(length(y),length(x));
ecc=zeros(length(y),length(x));
for i=1:length(x)
    for j=1:length(y)
        I=inpolygon(x(i),y(j),coords(:,2),coords(:,1));
        if I==1
            idx(j,i)=1;
        temp=im2(y(j)-w:y(j)+w-1,x(i)-w:x(i)+w-1);
        temp2=temp.*h;
        temp3=fft2(temp2);
        temp4=abs(fftshift(temp3));
        temp4=mat2gray(log(temp4+1));
        temp5=(temp4>0.07);
        temp5=bwareaopen(temp5,2);
        temp5=imfill(temp5,'holes');
        cc=bwconncomp(temp5,4);
        r=regionprops(cc,'Orientation','Eccentricity','MinorAxisLength','MajorAxisLength','Area');
        for k=1:cc.NumObjects
            areatemp(k)=r(k).Area;
        end
        [~,ind]=max(areatemp);clear areatemp
        theta(j,i)=90+r(ind).Orientation;
        ecc(j,i)=r(ind).Eccentricity;
        ar(j,i)=r(ind).MajorAxisLength./r(ind).MinorAxisLength;
        stdev(j,i)=std(temp4(:));
        line([x(i)-l*cos(theta(j,i)*pi/180),x(i)+l*cos(theta(j,i)*pi/180)],...
          [y(j)+l*sin(theta(j,i)*pi/180),y(j)-l*sin(theta(j,i)*pi/180)],'Color','g','LineWidth',1)
        end
    end
end

theta2=theta(:);
stdev2=stdev(:);
ar2=ar(:);
ecc2=ecc(:);
theta2=theta2(theta2~=0);
stdev2=stdev2(stdev2~=0);
ecc2=ecc2(ecc2~=0);
ar2=ar2(ar2~=0);

saveas(gcf,sprintf([folder2(2).name(1:end-4) '_directors.fig']));
close all
save(sprintf([folder2(2).name(1:end-4) 'actin_data.mat']),'coords','Cang','Ccell_cent','mask',...
    'Carea','Car','Ccirc','theta','stdev','ar','ecc','theta2','stdev2','ar2',...
    'ecc2','x','y');
