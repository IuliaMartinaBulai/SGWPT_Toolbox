%--------------------------------------------------------------------------
% File: sgwpt_kernel_abspline3_g.m
%
% Goal: evaluates spline functions
%
% Use: sgwpt_kernel_abspline3_g(x,lmax,sc)
%
% Inputs: x - variable in [0, lmax]
%         lmax - upper bound on spectrum
%         sc - the dilation coefficient
%
% Outputs: r - the spline filters
%
% Recalls: find.m
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
function r = sgwpt_kernel_abspline3_g(x,lmax,sc)
s = 1/lmax;
l1 = 1/(sc*s);
l2 = 2/(sc*s);
r1ind = find(x>=0 & x<=l1);
r2ind = find(x>l1 & x<=l2);
r3ind = find(x>l2);
v = @(x) (-5+11*x-6*x.^2+x.^3) ;
% as we initialize r with zero, computed function will implicitly be zero for
% all x not in one of the three regions defined above
r = zeros(size(x));
r(r1ind) = (sc*s*x(r1ind)).^2;
r(r2ind) = v(sc*s*x(r2ind));
r(r3ind) = 4./(sc*s*x(r3ind)).^2;
end