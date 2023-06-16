%--------------------------------------------------------------------------
% File: sgwpt_cheby_product.m
%
% Goal: Chebyshev coefficients for product of two polynomials
%
% Use: d = sgwpt_cheby_product(c1,c2)
%
% Inputs:  c1 - Chebyshev coefficients for pi(x) = sum c1(1+k) T_k(x) ; 0<=K<=M
%          c2 - Chebyshev coefficients for pj(x) = sum c2(1+k) T_k(x) ; 0<=K<=M
%
% Outputs:  d - Chebyshev coefficients for pi(x)*pj(x) = sum d(1+k) T_k(x) ; 0<=k<=2*M
%
% Recalls: lenght.m, sqrt.m, sgwpt_filter_design.m, sgwt_cheby_coeff.m,
%          sgwt_cheby_product.m
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
function d = sgwpt_cheby_product(c1,c2)
M = numel(c1)-1;
cp1 = c1;
cp1(1) = .5*c1(1);
cp2 = c2;
cp2(1) = .5*c2(1);
dp = zeros(1,2*M+1);
for m=0:(2*M)
    if (m==0)
        dp(1+m)=dp(1+m)+.5*cp1(1)*cp2(1);
        for i=0:M
            dp(1+m)=dp(1+m)+.5*cp1(1+i)*cp2(1+i);
        end
    elseif (m<=M)
        for i=0:m
            dp(1+m)=dp(1+m)+.5*cp1(1+i)*cp2(1+m-i);
        end
        for i=0:(M-m)
            dp(1+m)=dp(1+m)+.5*cp1(1+i)*cp2(1+i+m);
        end
        for i=m:M
            dp(1+m)=dp(1+m)+.5*cp1(1+i)*cp2(1+i-m);
        end
    else % M<m<=2*M
        for i=(m-M):M
            dp(1+m)=dp(1+m)+.5*cp1(1+i)*cp2(1+m-i);
        end
    end
end
d=dp;
d(1)=2*dp(1);
