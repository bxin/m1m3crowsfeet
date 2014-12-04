function sf2d=sf2dcomp(mask,S)

% original author: Brent Ellerbroek (brente@tmt.org)
%        Thirty Meter Telescope Observatory, Pasadena, CA 91125
% modified and adopted for LSST use by:
%        Bo Xin (bxin@lsst.org) and Ming Liang (mliang@noao.edu)
%         Large Synoptic Survey Telescope, Tucson, AZ 85719

% Terms of use is at the end of this code

% Description: calculate the 2d structure function of an optical surface.

% input
% mask:   n by n, [0,1]-valued function defining the aperture, with origin at
%        the point (n/2+1,n/2+1)
% S:     n by n surface profile (must be zero outside of the aperture)

% output:
% sf2d:  2d structure function of size 2n by 2n, with the origin located at
%        the point (n+1,n+1)

n=size(mask,1);
mask2=zeros(2*n);
mask2((n/2+1):(2*n-n/2),(n/2+1):(2*n-n/2))=mask;
S2=zeros(2*n);
S2((n/2+1):(2*n-n/2),(n/2+1):(2*n-n/2))=S;
Ssq2=S2.*S2.*mask2;

maskhat=fftshift(fft2(fftshift(mask2)));
Shat=fftshift(fft2(fftshift(S2)));
Ssqhat=fftshift(fft2(fftshift(Ssq2)));

maskmaskstar=maskhat.*conj(maskhat);
SSstar=Shat.*conj(Shat);
Ssqmaskstar=Ssqhat.*conj(maskhat);

maskacf=fftshift(ifft2(fftshift(maskmaskstar)));
Sacf=fftshift(ifft2(fftshift(SSstar)));
Smaskccf=fftshift(ifft2(fftshift(Ssqmaskstar)));

ind=[1 2*n:-1:2];
numer=Smaskccf+Smaskccf(ind,ind)-2*Sacf;
denom=maskacf;

numer=real(numer);
denom=real(denom);
ind=find(denom > 1d-6*sum(sum(mask)));
sf2d=zeros(2*n);
sf2d(ind)=numer(ind)./denom(ind);

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
