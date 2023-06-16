%--------------------------------------------------------------------------
% File: rmse_run_rfMRI_100_norm_meyer.m
%
% Goal: Compute the RMSE for rfMRI_100 for Shannon entropy best basis for
% Meyer filters
%
% It is recomended a computing cluster or at least the use of paralel
% computing (with parfor 3 hours on a personal computer)
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
%%
clc
clear all
tic
load('sgwpt_best_basis_run_rfMRI_100_norm_meyer')
tic
parpool(4)
parfor i = 1:1:100
    tau = 0;
    rmse = 0;
    [tau, best_basis_thresholded] = sgwpt_thresholding(i, best_basis_T1);
    
    % substract tau from the coefficients
    best_basis_thresholded_new = cellfun(@(x) max(0,x - tau)-max(0,-x-tau),...
        best_basis_thresholded, 'un', 0);
    
    % inverse using Chebyshev polynomial of Laplacian applied to best basis
    % thresholded
    r4_thresholded = sgwt_inverse(best_basis_thresholded_new,L,coef_bb_T1,arange);
    tau_T1(i) = tau;
    rmse = sqrt(immse(rfMRI_norm, r4_thresholded));
    rmse_T1(i) = rmse;
end
toc 
tic
[Min, index_min] = min(rmse_T1);
[tau, best_basis_thresholded] = sgwpt_thresholding(index_min, best_basis_T1);

% substract tau from the coefficients
best_basis_thresholded_new = cellfun(@(x) max(0,x - tau)-max(0,-x-tau),...
    best_basis_thresholded, 'un', 0);

% inverse using Chebyshev polynomial of Laplacian applied to best basis
% thresholded
r4_thresholded = sgwt_inverse(best_basis_thresholded_new,L,coef_bb_T1,arange);
toc
clear Ad;
clear D;
cd ./data
filename = ('rmse_run_rfMRI_100_norm_meyer');
save(filename);
toc
