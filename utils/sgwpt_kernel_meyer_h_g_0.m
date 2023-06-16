%--------------------------------------------------------------------------
% File: sgwpt_kernel_meyer_h_g_0.m
%
% Goal: evaluates meyer h0 and g0 functions
%
% Use: r = sgwpt_kernel_meyer_h_g_0(x,kerneltype,lamx)
%
% Inputs: x - variable in [0, lmax]
%         kerneltype - type of kernel
%         lmax - upper bound on spectrum
%
% Outputs: r - the Meyer filters h_0 and g_0 whitout dilation
%
% Recalls: find.m, abs.m
%
% Authors: IM Bulai, S Saliani
% Date last modified: September, 2022
%
% This file is part of the SGWPT toolbox
% (Spectral Graph Wavelet Packet Transform Toolbox)
% Copyright (C) 2022, Iulia Martina Bulai and Sandra Saliani.
%
% The SGWPT toolbox is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% The SGWPT toolbox is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with the SGWPT toolbox.  If not, see <http://www.gnu.org/licenses/>.
%--------------------------------------------------------------------------
function r=sgwpt_kernel_meyer_h_g_0(x,kerneltype,lmax)
a=lmax/3;
l1=a;
l2=2*a;
l3=2*a+a;
v=@(x) x.^4.*(35-84*x+70*x.^2-20*x.^3) ;
r1ind=find(x>=0 & x<l1);
r2ind=find(x>=l1 & x<l2);
r3ind=find(x>=l2 & x<=l3);
% as we initialize r with zero, computed function will implicitly be zero
% for all x not in one of the three regions defined above
r=zeros(size(x));
switch kerneltype
    case 'h0'
        r(r1ind)=1;
        r(r2ind)=cos((pi/2)*v(abs(x(r2ind))/l1-1));
        r(r3ind)=0;
    case 'g0'
        r(r1ind)=0;
        r(r2ind)=sin((pi/2)*v(abs(x(r2ind))/l1-1));
        r(r3ind)= 1;
    otherwise
        error(sprintf('unknown kernel type %s',kerneltype));
end