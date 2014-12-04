function [newS,idxmask,xnorm,ynorm] = generate_crows_feet(S,nreso,xnorm,ynorm,...
    OR,ORC,IRC,nfeet,sigmal,sigmah,depth,airbubble)

% author: Bo Xin (bxin@lsst.org)
%         Large Synoptic Survey Telescope, Tucson, AZ 85719

% Terms of use is at the end of this code

% Description: generate given number of crow's feet and add them to a
% existing surface S.

% input:
% S:      the existing surface, in our case, it is the M3 reduced map from
%           SOML
% nreso: surface map resolution
% xnorm and ynorm: normalzied coordinate arrays on which S is defined.
% OR:    outer radius, in mm
% ORC:   clear outer radius, in mm
% IRC:   clear inner radius, in mm
% nfeet: number of crow's feet to be added
%[sigmal,sigmah] is the range for standard deviation of the Gaussian along
%               the semi-major axis (along 
%               crow's feet wing span) Random numbers are used.
% depth: depth of the center hole
% airbubble:  0: no air bubble included
%             1: with air bubbles included

%% remake the surface map on a higher resolution grid if needed.
npixel0=size(xnorm,1);
if nreso~=size(S,1)
    [xnorm1,ynorm1]=meshgrid(-1:2/(nreso-1):1);
    newS=interp2(xnorm,ynorm,S,xnorm1,ynorm1);
    S=newS;
    xnorm=xnorm1;
    ynorm=ynorm1;
end

mirr=sqrt(xnorm.^2+ynorm.^2);
idxmask=mirr<ORC/OR&mirr>IRC/OR;

npixel=size(xnorm,1);
pixelSize=(xnorm(1,2)-xnorm(1,1))*OR;%in mm

%% determine the center coordinates
s = RandStream('mcg16807','Seed',nfeet);
RandStream.setGlobalStream(s);

%make the probability along r linear with r. and range from 0.5 to 1 from
%IRC to ORC.
% we need to calculate volumes of cones to throw out events
% r is the ratio of random numbers that would be kept.
h=0.5*IRC/ORC/(1-IRC/ORC);
r=1-(1./3*(ORC^2*(h+0.5)-IRC^2*h)-IRC^2*0.5)/((ORC^2-IRC^2));
scalef=4*ORC*ORC/(pi*(ORC*ORC-IRC*IRC))/r; 
ngen=round(scalef*nfeet);
x0=rand(ngen,1)*2-1;
y0=rand(ngen,1)*2-1;
r0=sqrt(x0.*x0+y0.*y0);
idx=r0<ORC/OR&r0>IRC/OR&rand(ngen,1)<0.5+(r0-IRC/ORC)*0.5/(1-IRC/ORC);
x0=x0(idx);
y0=y0(idx);
nreal=size(x0,1);

%% create stamps

%1.first 2D elliptical Gaussian that represents the wing
A1=300*ones(nreal,1);% in nm
sig1major=sigmal+(sigmah-sigmal)*rand(nreal,1); %in mm
sig1minor=0.75*ones(nreal,1); %in mm
%2. second elliptical Gaussian that represents the center hole (not the
%air bubble), this is the center part as seen on the Fizeau test fringes
A2=depth*ones(nreal,1);
sig2major=3.3*ones(nreal,1); %in mm
sig2minor=1.2*ones(nreal,1); %in mm
%3. air bubble (light trap)
rbubble=(5+4*rand(nreal,1))*1e-1; %0.1-0.9mm
%4. angle from the tangent
theta0=-atan(y0./x0);%+pi/180*20*(2*rand(nreal,1)-1);
theta1=theta0+pi/180*(25+10*(2*rand(nreal,1)-1));
theta2=theta0-pi/180*(25+10*(2*rand(nreal,1)-1));

% use a txt to override the x0 and y0
overlayFile='m1m3/m3overlay.txt';
m3ol=load(overlayFile);
olsize=min(nreal,size(m3ol,1));
y0(1:olsize)=(m3ol(1:olsize,1)-npixel0/2)/(npixel0/2);
x0(1:olsize)=(m3ol(1:olsize,2)-npixel0/2)/(npixel0/2);
theta1(1:olsize)=pi/2-atan((m3ol(1:olsize,3)-m3ol(1:olsize,1))./(m3ol(1:olsize,4)-m3ol(1:olsize,2)));
theta2(1:olsize)=pi/2-atan((m3ol(1:olsize,5)-m3ol(1:olsize,1))./(m3ol(1:olsize,6)-m3ol(1:olsize,2)));

% this is for making good-resulution subS plots
% pixelSize=0.5; %unit is mm, 
% A1(1)=300;
% sig1major(1)=200; %one representative size of a relatively big one in the big category
% sig1minor(1)=0.75;
% A2(1)=2500;
% sig2major(1)=3.3;
% sig2minor(1)=1.2;
% theta1(1)=-pi/7;
% theta2(1)=pi/7;

subn=ceil(sigmah/pixelSize*10); %5 sigma on each side
if mod(subn,2)==1 %to simplify debugging, we don't like subn to be odd number
    subn=subn+1;
end
subx0=subn/2+0.5;
suby0=subx0;
[subxg,subyg]=meshgrid((1:subn)-subx0,(1:subn)-suby0);

x0=round((x0+1)/2*npixel);
y0=round((y0+1)/2*npixel);

newS=S;
for i=1:nreal
    subx=(subxg+rand())*pixelSize; %crow's feet do not always center at pixel centers
    suby=(subyg+rand())*pixelSize;
    %make the pair of deeper wings
    Iwing1=A1(i)*exp(-0.5*( ((subx*cos(theta1(i))+suby*sin(theta1(i)))/sig1major(i)).^2+...
        ((-subx*sin(theta1(i))+suby*cos(theta1(i)))/sig1minor(i)).^2 ))+ ...
        A2(i)*exp(-0.5*( ((subx*cos(theta1(i))+suby*sin(theta1(i)))/sig2major(i)).^2+...
        ((-subx*sin(theta1(i))+suby*cos(theta1(i)))/sig2minor(i)).^2 ));
    Iwing2=A1(i)*exp(-0.5*( ((subx*cos(theta2(i))+suby*sin(theta2(i)))/sig1major(i)).^2+...
        ((-subx*sin(theta2(i))+suby*cos(theta2(i)))/sig1minor(i)).^2 ))+ ...
        A2(i)*exp(-0.5*( ((subx*cos(theta2(i))+suby*sin(theta2(i)))/sig2major(i)).^2+...
        ((-subx*sin(theta2(i))+suby*cos(theta2(i)))/sig2minor(i)).^2 ));
    subS=max(Iwing1,Iwing2);
    if airbubble
        Ibubble=max(0,rbubble(i)^2-subx.^2-suby.^2);
        subS=max(subS, sqrt(Ibubble)*1e6);
    end

    xstart=x0(i)-subn/2;xs=1;
    xend=x0(i)+subn/2-1;xe=subn;
    ystart=y0(i)-subn/2;ys=1;
    yend=y0(i)+subn/2-1;ye=subn;
    if xstart<1
        xs=xs+1-xstart;
        xstart=1;
    elseif xend>npixel
        xe=xe-(xend-npixel);
        xend=npixel;
    end
    if ystart<1
        ys=ys+1-ystart;
        ystart=1;
    elseif yend>npixel
        ye=ye-(yend-npixel);
        yend=npixel;
    end
    newS(xstart:xend,ystart:yend)=newS(xstart:xend,ystart:yend)-subS(xs:xe,ys:ye);
    %     a=subS;m=300;a(a>1000)=0;imagesc((1:m)*.1,(1:m)*.1,extractArray(a,m));axis square xy;xlabel('mm');ylabel('mm');
    % imagesc((1:2000)*.05,(1:2000)*.05,subS);axis square xy;colorbar;xlabel('cm');ylabel('cm');
    %     a=newS;a(a<-1000)=0;imagesc(a);colorbar;axis square;
    
    % for making contour plot for the example image
    %     sideD=10; % in cm
    %     np=sideD*10/pixelSize; %number of pixels on plot
    %     np=2*np;
    %     b=extractArray(subS,np);
    %     hdiff=20*300; %from center to edge, if without center hole
    %     [subx,suby]=meshgrid(-hdiff+hdiff/np:2*hdiff/np:hdiff-hdiff/np);
    %     subr=(subx.*subx+suby.*suby)/(hdiff);
    %     xx=pixelSize/10:pixelSize/10:sideD;
    %     b=b';
    %     b1=circshift(b,[30 -15]);
    %     b2=circshift(b,[-30 15]);
    %     b=max(b1,b2);
    %     c=-b-subr;
    %     zmin=min(min(c));zmax=max(max(c));zind=zmin:300:zmax;
    %     contourf(xx,xx,extractArray(c,np/2),zind);axis xy square;

end

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

