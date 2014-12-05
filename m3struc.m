function [] = m3struc(nreso,nbig,nsmall,airbubble,simg,sft) %,sbf)

% author: Bo Xin (bxin@lsst.org)
%         Large Synoptic Survey Telescope, Tucson, AZ 85719

% Terms of use is at the end of this code

% Description: synthesize M3 surface map with crow's feet, and evaluate
% structure function
%
% input:
% nreso: surface map resolution
% nbig:  number of "big" crow's feet included
% nsmall:  number of "small" crow's feet included
% airbubble:  0: no air bubble included
%             1: with air bubbles included
% simg: switch, 0: don't synthesize the image, 
%                   load it from previously saved mat file
%               1: synthezie the image, and save it into mat file
% sft:  switch, 0: don't run the FFT-based structure function calculation,
%                   load it from previously saved mat file 
%               1: run the FFT-based structure function calculation, and
%                  save it into mat file
% sbf is the same thing, except bf stands for Brute Force 
%           (this not used for now)
%
% output file names are like
%       sfm3_1006_25_300_0_ft.mat
%       sfm3_1006_25_300_0_bf.mat
%       sfm3_1006_0_0_0_ml.mat
% the variables inside are like
%       sfft, spft, sfbf, spbf, sfsp, spsp
%       (ft for Fourier Transform, bf for Brute Force,sp for SPecification,
%       sf for Structure Function, sp for SeParation)
% examples:
%       m3struc(1006,0,0,0,1,1,1);
%       m3struc(1006,10,100,600,1,1,1);
%       

OR1=4202.5; 
IR1=2533; 
OR3=2533; 
IR3=533.5;

%use these number for clear aperture to define mask
OR1C=4180;
IR1C=2558;
OR3C=2508;
IR3C=550;

%% calculate the structure function

filename=sprintf('m1m3/sfm3_%d_%d_%d_%d_img.mat',nreso,nbig,nsmall,airbubble);

if simg

    aftAOM3File='internalData/SOML141019/M3 141019 -22modes M1FA -FC.h5';
    [mlaftAOM3,centerRow,centerCol,pixelSize] = readH5Map(aftAOM3File);
    
    % we need to define the clear aperture
    [rows, cols] = size(mlaftAOM3);
    xVec = 1:cols;
    xVec = (xVec - centerCol) * pixelSize; % m
    yVec = 1:rows;
    yVec = (yVec - centerRow) * pixelSize; % m, increasing upward
    [x,y] = meshgrid(xVec,yVec); % rows x cols
    xnorm=x*1000/OR3;
    ynorm=y*1000/OR3;
    %     nMir=max(imgm,imgn);
    mirr=sqrt(xnorm.^2+ynorm.^2);
    
    mlaftAOM3=flipud(mlaftAOM3);
    idx=mirr<OR3C/OR3&mirr>IR3C/OR3;
    mlaftAOM3(~idx)=nan; %this masks off everything outside of the clear aperture
    mlaftAOM3=mlaftAOM3*1e3; %from um to nm
    
    newS=mlaftAOM3;
    if (cols==nreso && nbig==0 && nsmall == 0  && airbubble==0)
        % do nothing, use mlaftAOM3 as the map
    else
        % add crow's feet
        if nbig~=0
            [newS,idx,xnorm,ynorm]=generate_crows_feet(newS,nreso,xnorm,ynorm,OR3,OR3C,IR3C,nbig,50,100,2500,airbubble);
        end
        if nsmall~=0
            [newS,idx,xnorm,ynorm]=generate_crows_feet(newS,nreso,xnorm,ynorm,OR3,OR3C,IR3C,nsmall,5,50,1000,airbubble);
        end
    end
    save(filename,'newS','idx','xnorm','ynorm');
else
    load(filename);
end
cols=size(newS,2);

figure(1);
clf;
imagesc(newS);axis square xy;colorbar;
title(sprintf('RMS inside CA: %5.2fnm',rms(newS(~isnan(newS)))));

% write the surface map into a fits file 
% OR=2533;
% npixel=size(xnorm,1);
% pixelSize=(xnorm(1,2)-xnorm(1,1))*OR;%in mm
% header.pixSize=pixelSize;
% header.centrCol=npixel/2+0.5;
% header.centrRow=npixel/2+0.5;
% if airbubble
%     fits_write(sprintf('m1m3/m3_%d_wBubble.fits',npixel),newS,header);
% else
%     fits_write(sprintf('m1m3/m3_%d.fits',npixel),newS,header);
% end

%% spml and sfml can only be loaded from, e.g., Excel files
mldata=xlsread('internalData/SOML141019/M3 141019 -22modes M1FA -FC_masked ICA=1.1 OCA=5.016.xls');
spml=mldata(:,1);
sfml=mldata(:,2);

load('m1m3/sfm3_1006_0_0_0_ft.mat');
spmlft=spft;
sfmlft=sfft;

%% ft: fourier transform

filename=sprintf('m1m3/sfm3_%d_%d_%d_%d_ft.mat',nreso,nbig,nsmall,airbubble);
if sft
    d = 2;        
    sep = 1:d:round(cols*.99); %cut off the bad tails on the low frequency end
    spft = sep.*(OR3*2/1000)./cols; %separation in meter

    [sfft, ~]=strucFuncFFT(idx,newS,sep);
    sfft=sfft*2; % wavefront phase instead of surface
    save(filename,'sfft','spft');
else
    load(filename);
end

%% bf: brute force
% no need to calculate this most of the time, since we've verified that
% with larger and larger grid, this converges to FFT SF.

% filename=sprintf('m1m3/sfm3_%d_%d_%d_%d_bf.mat',nreso,nbig,nsmall,airbubble);
% if sbf
%     oversample=2;
%     
%     sfbf=strucFuncSurf(newS,oversample);
%     sfbf=sfbf(1:round(size(sfbf,1)*.99)); %cut off the bad tails on the low frequency end
%     sfbf=sfbf*2; % wavefront phase instead of surface
% 
%     spbf=1:size(sfbf,1);
%     pixelSize=(xnorm(1,2)-xnorm(1,1))*OR3/1000;
%     spbf=spbf*pixelSize*oversample;
%     
%     save(filename,'sfbf','spbf');
% else
%     load(filename);
% end

%% LSST specification
spsp = 0.01:0.1:5.1;
r0 = 1.51; %---- Fried parameter in m
theta = 11;   %-- Scattering at small scale nm------------
w = 500;     
D = 5.016;  

sfsp= 2.*theta.^2+(w./2./pi).^2*6.88*(spsp./r0).^(5./3).*(1-0.975.*(spsp./D).^(1/3));
sfsp = sfsp.^(1/2);

%% plot them
figure(2);clf;
loglog(spsp,sfsp,'-r');grid on;
hold on;
loglog(spml,sfml,'--g');
loglog(spmlft,sfmlft,'--b');
loglog(spft,sfft,'-k');
% loglog(spbf,sfbf,'-m');

% legend({'LSST specification','from SOML','FFT Calculation','Brute Force Calculation'},'Location','SouthEast','FontSize',12);
legend({'LSST specification','SOML Reduced M3 map (SOML calculation)',...
    'SOML Reduced M3 map (our calculation)','M3 Synthetic Map (our calculation)'},'Location','SouthEast','FontSize',12);
xlabel('Separation (m)');ylabel('Structure Function (sqrt(D)) (nm)');
% text(0.01,30,'LSST M3','FontSize',20);

end

%% Terms of use
% This subroutine is copyrighted in the name of its author(s) and the
% affiliated institutions. 
% It is free and you are welcome to use it for your research.  It is 
% distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY.  
% We would appreciate acknowledgement when the subroutine is used.  
% You can cite the subroutine as
% Author, computer code title of the file,
% (https://github.com/bxin/m1m3crowsfeet),2014