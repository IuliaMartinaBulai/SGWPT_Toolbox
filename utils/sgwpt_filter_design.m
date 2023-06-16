%--------------------------------------------------------------------------
% File: sgwpt_filter_design.m
%
% Goal: Return list of scaled wavelet kernels
%
% Use: [g,t] = sgwpt_filter_design(lmax,varargin)
%
% Inputs:  lmax - upper bound on spectrum
%
% Outputs: g - h and g functions handles (g{1} is h and g{2} is g)
%          t - values of scale parameters used for wavelet kernels
% Selectable parameters :  designtype
%
% Recalls: argselectAssign.m, argselectCheck.m, sgwpt_kernel_meyer_h_g_i.m
%          with i = 0,1,2, sgwpt_kernel_abspline3_g.m
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
function [g,t] = sgwpt_filter_design(lmax,varargin)
control_params={'designtype','abspline3','lpfactor',20,...
    'a',2,...
    'b',2,...
    't1',1,...
    't2',2,...
    };
argselectAssign(control_params);
argselectCheck(control_params,varargin);
argselectAssign(varargin);
t = 1;
K = lmax;
lmin = lmax/K;
switch designtype
    case 'meyer_h_g_0'
        g{1} = @(x) sgwpt_kernel_meyer_h_g_0(t*x,'h0',lmax);
        g{2} = @(x) sgwpt_kernel_meyer_h_g_0(t*x,'g0',lmax);
    case 'meyer_h_g_1'
        g{1} = @(x) sgwpt_kernel_meyer_h_g_1(t*x,'h1',lmax);
        g{2} = @(x) sgwpt_kernel_meyer_h_g_1(t*x,'g1',lmax);
    case 'meyer_h_g_2'
        g{1} = @(x) sgwpt_kernel_meyer_h_g_2(t*x,'h2',lmax);
        g{2} = @(x) sgwpt_kernel_meyer_h_g_2(t*x,'g2',lmax);
    case 'abspline3_h_g_0'
        sc = 1;
        lfac = (lmax*0.6*lmin);
        h = @(x) exp(-x.^4);
        g{1} = @(x) h(sc*x./lfac) ;
        g{2} = @(x) sgwpt_kernel_abspline3_g(x,lmax,1)
    case 'abspline3_h_g_1'
        sc = 2;
        lfac = (lmax*0.6*lmin);
        h = @(x) exp(-x.^4);
        g{1} = @(x) h((sc*x)./lfac);
        g{2} = @(x) sgwpt_kernel_abspline3_g(x,lmax,2)
    case 'abspline3_h_g_2'
        sc = 4;
        lfac = (lmax*0.6*lmin);
        h = @(x) exp(-x.^4);
        g{1} = @(x) h((sc*x)./lfac);
        g{2} = @(x) sgwpt_kernel_abspline3_g(x,lmax,4)
    otherwise
        error('Unknown design type');
end