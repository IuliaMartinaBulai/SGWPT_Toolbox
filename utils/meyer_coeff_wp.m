%--------------------------------------------------------------------------
% File: meyer_coeff_wp.m
%
% Goal: compute the Chebyshev polinomial approximation coefficients for Meyer filters
%
% Use: [c_0, dp1, dp2] = meyer_coeff_wp(m, lmax)
%
% Inputs:  m - order of polynomial approximation
%         lmax - upper bound on spectrum
%
% Outputs: c_0 - Chebyshev coefficients for h0 and g0
%          dp1 - Chebyshev coefficients for  h1h0, g1h0, h1g0, g1g0
%          dp2 - Chebyshev coefficients for g2h1h0, h2g1h0, g2g1h0, h2h1g0,
%          g2h1g0, h2g1g0, g2g1g0
%
% Recalls: lenght.m, sqrt.m, sgwpt_filter_design.m, sgwt_cheby_coeff.m,
%          sgwpt_cheby_product.m
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
function [c_0, dp1, dp2] = meyer_coeff_wp(m, lmax)
% Design filters for transform
fprintf('Designing transform in spectral domain\n');
% Meyer filters
designtype_0 = 'meyer_h_g_0';
designtype_1 = 'meyer_h_g_1';
designtype_2 = 'meyer_h_g_2';
[g_0,t] = sgwpt_filter_design(lmax,'designtype',designtype_0);
[g_1,t] = sgwpt_filter_design(lmax,'designtype',designtype_1);
[g_2,t] = sgwpt_filter_design(lmax,'designtype',designtype_2);
arange = [0 lmax];
% Chebyshev polynomial approximation
fprintf('Computing Chebyshev polynomials of order %g for fast transform \n',m);
% compute Chebyshev coefficients for h0 and g0
for k = 1:numel(g_0)
    c_0{k} = sgwt_cheby_coeff(g_0{k},m,m+1,arange);
end
% compute Chebyshev coefficients for h_1 and g_1
for k = 1:numel(g_1)
    c_1{k} = sgwt_cheby_coeff(g_1{k},m,m+1,arange);
end
% compute Chebyshev coefficients for h2 and g2
for k = 1:numel(g_2)
    c_2{k} = sgwt_cheby_coeff(g_2{k},m,m+1,arange);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Filtri al secondo livello  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% h1h0 = h1.*h0;
% g1h0 = g1.*h0;
% h1g0 = h1.*g0;
% g1g0 = g1.*g0;
fprintf('Computing Chebyshev polynomials of order %g for h1h0, g1h0, h1g0, g1g0  \n',2*m);
for k = 1:numel(g_0)
    for j = 1:numel(g_1)
        dp1{k,j} = sgwpt_cheby_product(c_0{k},c_1{j});
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Filtri al terzo livello  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% h2h1h0 = h2.*h1.*h0;
% g2h1h0 = g2.*h1.*h0;
% h2g1h0 = h2.*g1.*h0;
% g2g1h0 = g2.*g1.*h0;
% h2h1g0 = h2.*h1.*g0;
% g2h1g0 = g2.*h1.*g0;
% h2g1g0 = h2.*g1.*g0;
% g2g1g0 = g2.*g1.*g0;
fprintf('Computing Chebyshev polynomials of order %g for g2h1h0, h2g1h0, g2g1h0, h2h1g0, g2h1g0, h2g1g0, g2g1g0  \n',4*m);
for k = 1: numel(c_2)
    for j = 1: numel(dp1)
        dp2{k,j} = sgwpt_cheby_product(c_2{k},dp1{j});
    end
end
end