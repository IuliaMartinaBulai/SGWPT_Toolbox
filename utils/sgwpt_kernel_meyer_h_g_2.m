%--------------------------------------------------------------------------
% File: sgwpt_kernel_meyer_h_g_2.m
%
% Goal: evaluates meyer h_2 and g_2 functions
%
% Use: r = sgwpt_kernel_meyer_h_g_2(x,kerneltype,lamx)
%
% Inputs: x - variable in [0, lmax]
%         kerneltype - type of kernel
%         lmax - upper bound on spectrum
%
% Outputs: r - the Meyer filters h_2 and g_2 with dilation 4
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
function r = sgwpt_kernel_meyer_h_g_2(x,kerneltype,lmax)
a = lmax/3;
D = 3*a/2;
l1 = a/4;
l2 = a/2;
l3 = a;
l4 = 5/4*a;
l5 = 7/4*a;
l6 = 2*a;
l7 = 5/2*a;
l8 = 11/4*a;
l9 = lmax;
v = @(x) x.^4.*(35-84*x+70*x.^2-20*x.^3) ;
r1ind = find(x>=0 & x<l1);
r2ind = find(x>=l1 & x<l2);
r3ind = find(x>=l2 & x<l3);
r4ind = find(x>=l3 & x<l4);
r5ind = find(x>=l4 & x<l5);
r6ind = find(x>=l5 & x<l6);
r7ind = find(x>=l6 & x<l7);
r8ind = find(x>=l7 & x<l8);
r9ind = find(x>=l8 & x<=l9);
% as we initialize r with zero, computed function will implicitly be zero for
% all x not in one of the three regions defined above
r = zeros(size(x));
switch kerneltype
    case 'h2'
        r(r1ind) = 1;
        r(r2ind) = cos((pi/2)*v(4*abs(x(r2ind))/a-1));
        r(r3ind) = 0;
        r(r4ind) = sin((pi/2)*v(4*(abs(x(r4ind))-D/2)/a-1));
        r(r5ind) = 1;
        r(r6ind) = cos((pi/2)*v(4*(abs(x(r6ind))-D)/a-1));
        r(r7ind) = 0;
        r(r8ind) = sin((pi/2)*v(4*(abs(x(r8ind))-3*D/2)/a-1));
        r(r9ind) = 1;
    case 'g2'
        r(r1ind) = 0;
        r(r2ind) = sin((pi/2)*v(4*abs(x(r2ind))/a-1));
        r(r3ind) = 1;
        r(r4ind) = cos((pi/2)*v(4*(abs(x(r4ind))-D/2)/a-1));
        r(r5ind) = 0;
        r(r6ind) = sin((pi/2)*v(4*(abs(x(r6ind))-D)/a-1));
        r(r7ind) = 1;
        r(r8ind) = cos((pi/2)*v(4*(abs(x(r8ind))-3*D/2)/a-1));
        r(r9ind) = 0;
    otherwise
        error(sprintf('unknown kernel type %s',kerneltype));
end