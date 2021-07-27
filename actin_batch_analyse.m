clear all;close all
foldername=uigetdir;
folder=dir([foldername '/*frames*']);
nf=length(folder);

for i=1:nf
    tic
    file=fullfile([foldername '/' folder(i).name]);
    actin_extravaganza(file);
    toc
end