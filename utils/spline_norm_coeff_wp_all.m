%--------------------------------------------------------------------------
% File: spline_norm_coeff_wp_all.m
%
% Goal: Exact wavelets and Chebyshev polinomial approximation coefficients
% for spline normalized filters
%
% Use: [c_0_all, dp1_all, dp2_all]= spline_norm_coeff_wp_all(m, lmax)
%
% Inputs:  m - order of polynomial approximation
%         lmax - upper bound on spectrum
%
% Outputs: c_0_all - Chebyshev coefficients for h0 and g0
%          dp1_all - Chebyshev coefficients for  h1h0, g1h0, h1g0, g1g0
%          dp2_all - Chebyshev coefficients for g2h1h0, h2g1h0, g2g1h0, h2h1g0,
%          g2h1g0, h2g1g0, g2g1g0
%
% Recalls: lenght.m, sqrt.m, sgwpt_filter_design_wp.m, sgwt_cheby_coeff.m,
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
function [c_0_all, dp1_all, dp2_all] = spline_norm_coeff_wp_all(m, lmax)
% Design filters for transform
fprintf('Designing transform in spectral domain\n');
% spline filters
designtype_0 = 'abspline3_h_g_0';
designtype_1 = 'abspline3_h_g_1';
designtype_2 = 'abspline3_h_g_2';
[g_0,t] = sgwpt_filter_design(lmax,'designtype',designtype_0);
[g_1,t] = sgwpt_filter_design(lmax,'designtype',designtype_1);
[g_2,t] = sgwpt_filter_design(lmax,'designtype',designtype_2);
% Chebyshev polynomial approximation
fprintf('Computing Chebyshev polynomials of order %g for fast transform \n',m);
arange = [0 lmax];
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
% h1h=h1.*h;
% g1h=g1.*h;
% h1g=h1.*g;
% g1g=g1.*g;
fprintf('Computing Chebyshev polynomials of order %g for h1h0, g1h0, h1g0, g1g0  \n',2*m);
for k = 1:numel(g_0)
    for j = 1:numel(g_1)
        dp1{k,j} = sgwpt_cheby_product(c_0{k},c_1{j});
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Filtri al terzo livello  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% h2h1h=h2.*h1.*h;
% g2h1h=g2.*h1.*h;
% h2g1h=h2.*g1.*h;
% g2g1h=g2.*g1.*h;
% h2h1g=h2.*h1.*g;
% g2h1g=g2.*h1.*g;
% h2g1g=h2.*g1.*g;
% g2g1g=g2.*g1.*g;
fprintf('Computing Chebyshev polynomials of order %g for g2h1h0, h2g1h0, g2g1h0, h2h1g0, g2h1g0, h2g1g0, g2g1g0  \n',4*m);
for k = 1: numel(c_2)
    for j = 1: numel(dp1)
        dp2{k,j} = sgwpt_cheby_product(c_2{k},dp1{j});
    end
end
% compute Chebyshev coefficients for h0 and g0
for k=1:numel(g_0)
    c_0_all{k}=sgwt_cheby_coeff(g_0{k},m,m+1,arange)./sqrt(3); % da cancelare
end
% compute Chebyshev coefficients for h_1 and g_1
for k=1:numel(g_1)
    c_1_all{k}=sgwt_cheby_coeff(g_1{k},m,m+1,arange)./sqrt(3);
end
% compute Chebyshev coefficients for h2 and g2
for k=1:numel(g_2)
    c_2_all{k}=sgwt_cheby_coeff(g_2{k},m,m+1,arange)./sqrt(3);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Filtri al secondo livello  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% h1h=h1.*h;
% g1h=g1.*h;
% h1g=h1.*g;
% g1g=g1.*g;
fprintf('Computing Chebyshev polynomials of order %g for h1h, g1h, h1g, g1g  \n',2*m);
for k=1:numel(c_1_all)
    for j=1:numel(g_1)
        dp1_all{k,j}=sgwpt_cheby_product(c_0{k},c_1{j})./sqrt(3);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Filtri al terzo livello  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% h2h1h=h2.*h1.*h;
% g2h1h=g2.*h1.*h;
% h2g1h=h2.*g1.*h;
% g2g1h=g2.*g1.*h;
% h2h1g=h2.*h1.*g;
% g2h1g=g2.*h1.*g;
% h2g1g=h2.*g1.*g;
% g2g1g=g2.*g1.*g;
fprintf('Computing Chebyshev polynomials of order %g for g2h1h, h2g1h, g2g1h, h2h1g, g2h1g, h2g1g, g2g1g  \n',4*m);
for k=1: numel(c_2_all)
    for j=1: numel(dp1_all)
        dp2_all{k,j}=sgwpt_cheby_product(c_2{k},dp1{j})./sqrt(3);
    end
end
end