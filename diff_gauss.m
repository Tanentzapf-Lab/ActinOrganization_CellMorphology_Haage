function dogImg=diff_gauss(img)

%img=imadjust(img);
% if not high res
sigma1 =  1;
sigma2 = 10;

hsize1 = [2,2];
hsize2 = [20,20];

% if high res
% sigma1 =  4;
% sigma2 = 15;
% 
% hsize1 = [14, 14];
% hsize2 = [40,40];

h1 = fspecial('gaussian', hsize1, sigma1);
h2 = fspecial('gaussian', hsize2, sigma2);

gauss1 = imfilter(img,h1,'replicate');
gauss2 = imfilter(img,h2,'replicate');

dogImg = gauss1 - gauss2;