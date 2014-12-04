function sf1d=sf1dcomp(sf2d,sep)

% original author: Brent Ellerbroek (brente@tmt.org)
% modified and adopted for LSST use by:
%        Bo Xin (bxin@lsst.org) and Ming Liang (mliang@noao.edu)

% Terms of use is at the end of this code

% Description: evaluate the 1d structure function at given separation
%
% input:
% sep:   vector of separations at which to compute the 1-d structure
%        function
% sf2d:  2d structure function of size 2n by 2n, with the origin located at
%        the point (n+1,n+1)
% output:
% sf1d:  1-d structure function evaluated at the separations sep

n=size(sf2d,1)/2;
x=-1/2 : 1/(2*n) : (1/2-1/(2*n));
[x,y]=meshgrid(x);
kap=sqrt(x.^2+y.^2);

sfdom = (sf2d ~= 0);
sfdom(n+1,n+1)=1;
sfhat=fftshift(fft2(fftshift(sf2d)));
sfdomhat=fftshift(fft2(fftshift(sfdom)));

sf1d=zeros(size(sep));
for j=1:length(sep)
    r=sep(j);
    wf=besselj(0,2*pi*r*kap);
    sf1d(j)=real(sum(sum(wf.*sfhat))/sum(sum(wf.*sfdomhat)));
end

end

%% Terms of use
% This subroutine is copyrighted in the name of its author(s) and the 
% Large Synoptic Survey Telescope, under management of 
% the Association of Universities for Research in Astronomy (AURA).
% If it is modified from code written by others from a different institution, 
% we have obtained appropriate permission to use and distribute the code.
% It is free and you are welcome to use it for your research.  It is 
% distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY.  
% We would appreciate acknowledgement when the subroutine is used.  
% You can cite the subroutine as
% Author, computer code Title of the file (https://github.com/bxin/m1m3crowsfeet), 
% Large Synoptic Survey Telescope, Tucson, Arizona, 2014

