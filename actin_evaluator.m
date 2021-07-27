clear all;close all
scale=0.155; %% change this first!
ARcutoff=4;
th=10;

foldername=uigetdir;
folder=dir([foldername '/*.mat']);
nf=length(folder);
cell_area=zeros(nf,1);
cell_circ=zeros(nf,1);
cell_angle=zeros(nf,1);
cell_AR=zeros(nf,1);
cell_irreg=zeros(nf,1);

actin_angle_adj=cell(nf,1);
actin_mean_AR=zeros(nf,1);
actin_AR_inner=zeros(nf,1);
actin_AR_outer=zeros(nf,1);
actin_mean_angle_adj=zeros(nf,1);
idx=cell(nf,1);
fib=zeros(nf,1);
fib_inner=zeros(nf,1);
fib_outer=zeros(nf,1);
filenames=cell(nf,1);
S=zeros(nf,1);
for i=1:nf
    load(fullfile([foldername '/' folder(i).name]));
    
    se=strel('disk',th);
    mask2=imerode(mask,se);
    cc=bwconncomp(mask2);
    r=regionprops(cc,'PixelList');
    coords2=bwtraceboundary(mask2,[r(1).PixelList(1,2),r(1).PixelList(1,1)],'N');
    clear se mask2 cc r
    idx{i}=zeros(length(y),length(x));
    for k=1:length(x)
        for j=1:length(y)
            I=inpolygon(x(k),y(j),coords(:,2),coords(:,1));
            if I==1
                idx{i}(j,k)=1;
            end
            I=inpolygon(x(k),y(j),coords2(:,2),coords2(:,1));
            if I==1
                idx{i}(j,k)=2;
            end
        end
    end
    
    filenames{i}=folder(i).name;
    %CELL DATA
    if Cang<0
        Cang=180+Cang;
    end
    cell_angle(i)=Cang;
    cell_area(i)=Carea*scale^2;
    cell_circ(i)=Ccirc;
    cell_AR(i)=Car;
    cc1=bwconncomp(mask);
    r1=regionprops(cc1,'Perimeter','ConvexImage');
    ci=r1.ConvexImage;
    cc2=bwconncomp(ci);
    r2=regionprops(cc2,'Perimeter');
    cell_irreg(i)=r1.Perimeter/r2.Perimeter;
    
    %%ACTIN DATA
    actin_mean_angle_adj(i)=mean(theta2-cell_angle(i));
    fib(i)=length(ar2(ar2>ARcutoff))/length(ar2);
    actin_mean_AR(i)=mean(ar2);
    AR=ar(:);IDX=idx{i}(:);
    AR_inner=ar(IDX==2);
    AR_outer=ar(IDX==1);
    actin_AR_inner(i)=mean(AR_inner);
    actin_AR_outer(i)=mean(AR_outer);
    fib_inner(i)=length(AR_inner(AR_inner>ARcutoff))/length(AR_inner);
    fib_outer(i)=length(AR_outer(AR_outer>ARcutoff))/length(AR_outer);
    actin_angle_adj{i}=theta2-cell_angle(i);
    actin_mean_angle_adj(i)=mean(actin_angle_adj{i}(ar2>ARcutoff));
    Sk=zeros(length(theta2),1);
    for j=1:length(theta2)
        Sk(j)=cos(pi/180*(theta2(j)-cell_angle(i)));
    end
    S(i)=2*(mean(Sk(ar2>ARcutoff))-1/2);
    clear Sk
end


results=[filenames num2cell(cell_area) num2cell(cell_AR) num2cell(cell_circ) num2cell(cell_irreg) num2cell(fib.*100) num2cell(fib_inner.*100) num2cell(fib_outer.*100) num2cell(actin_mean_AR) num2cell(actin_AR_inner) num2cell(actin_AR_outer) num2cell(actin_mean_angle_adj) num2cell(S)];
rownames={'Name';'Cell Area';'Cell aspect ratio';'Cell circularity';'Cell irregularity';'Percent fibrousness';'Percent fibrousness inner';'Percent fibrousness outer';'Mean actin aspect ratio';'Mean actin aspect ratio inner';'Mean actin aspect ratio outer';'Mean adjusted actin angle';'Order Paramter'}; 
data=[rownames'; results];

prompt = {'Input file name'};dlg_title = 'Data Output File Name';num_lines = 1;def = {'filename'};
ofilename=char(inputdlg(prompt,dlg_title,num_lines,def));
xlswrite(ofilename,data)


save(sprintf([ofilename '.mat']),'actin_mean_AR',...
    'actin_angle_adj',...
    'cell_irreg','cell_area','cell_circ','cell_AR')
