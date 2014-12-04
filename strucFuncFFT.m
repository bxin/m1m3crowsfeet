function [sf1d,sf2d]=strucFuncFFT(mask,S,sep)

% author: Brent Ellerbroek (brente@tmt.org)
% modified and adopted for LSST use by:
%        Bo Xin (bxin@lsst.org) and Ming Liang (liang@noao.edu)

% Terms of use is at the end of this code

% input:
% mask:  n by n, [0,1]-valued function defining the aperture, with origin at
%        the point (n/2+1,n/2+1)
% S:     n by n surface map (outside of the aperture will be zeroed out by mask)
% sep:   vector of separations at which to compute the 1-d structure
%        function

% output:
% sf1d:  1-d structure function evaluated at the separations sep  
% sf2d:  2d structure function of size 2n by 2n, with the origin located at
%        the point (n+1,n+1)

fprintf('strucFuncFFT started at %s\n',datestr(datetime));

S(isnan(S))=0;
S=S.*mask;

sf2d=sf2dcomp(mask,S);
sf1d=sf1dcomp(sf2d,sep);

%our structure function is defined as sqrt(D)
sf1d = abs((sf1d).^0.5);

fprintf('strucFuncFFT ended at %s\n',datestr(datetime));

end

%% Terms of use
% This subroutine is copyrighted in the name of its author(s) and the 
% Large Synoptic Survey Telescope, under management of 
% the Association of Universities for Research in Astronomy (AURA).
% In the case that the original author of the code is from a different institution, 
% we have obtained appropriate permission to use and distribute the code.
% It is free and you are welcome to use it for your research.  It is 
% distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY.  
% We would appreciate acknowledgement when the subroutine is used.  
% You can cite subroutines as
%
% Author, computer code Title of the file (https://github.com/bxin/m1m3crowsfeet), 
% Large Synoptic Survey Telescope, Tucson, Arizona, 2014
