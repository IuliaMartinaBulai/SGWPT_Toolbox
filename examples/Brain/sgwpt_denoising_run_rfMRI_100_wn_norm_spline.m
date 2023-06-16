%--------------------------------------------------------------------------
% File:  sgwpt_denoising_run_rfMRI_100_wn_norm_spline.m
%
% Goal: Main file for the denoising of the normalized noisy signal 
% for the Shannon best basis using spline filters (for the remaining 
% entropy measures decomment the code)
%
% Authors: IM Bulai, S Saliani
% Date last modified: September, 2022
%
% This file is part of the SGWPT toolbox
% (Spectral Graph Wavelet Packet Transform Toolboxm toolbox)
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
load('sgwpt_best_basis_run_rfMRI_100_wn_norm_spline')

prc_range_1 = 54;
prc_range_2 = 80;
prc_range_3 = 50;
fprintf('Thresholding for the best basis with Shannon energy\n');
% thresholding from the best basis with Shannon
[tau_T1_54, best_basis_T1_thresholded_54] = sgwpt_thresholding(prc_range_1, best_basis_T1);
[tau_T1_80, best_basis_T1_thresholded_80] = sgwpt_thresholding(prc_range_2, best_basis_T1);
[tau_T1_50, best_basis_T1_thresholded_50] = sgwpt_thresholding(prc_range_3, best_basis_T1);
% substract tau from the coefficients
best_basis_T1_thresholded_54_new = cellfun(@(x) max(0,x - tau_T1_54)-max(0,-x-tau_T1_54),...
    best_basis_T1_thresholded_54, 'un', 0);
best_basis_T1_thresholded_80_new = cellfun(@(x) max(0,x - tau_T1_80)-max(0,-x-tau_T1_80),...
    best_basis_T1_thresholded_80, 'un', 0);
best_basis_T1_thresholded_50_new = cellfun(@(x) max(0,x - tau_T1_50)-max(0,-x-tau_T1_50),...
    best_basis_T1_thresholded_50, 'un', 0);
fprintf('Thresholding for the best basis with threshold energy\n');
% thresholding from the best basis with threshold
[tau_T2_54, best_basis_T2_thresholded_54] = sgwpt_thresholding(prc_range_1, best_basis_T2);
[tau_T2_80, best_basis_T2_thresholded_80] = sgwpt_thresholding(prc_range_2, best_basis_T2);
[tau_T2_50, best_basis_T2_thresholded_50] = sgwpt_thresholding(prc_range_3, best_basis_T2);
% substract tau from the coefficients
best_basis_T2_thresholded_54_new = cellfun(@(x) max(0,x - tau_T2_54)-max(0,-x-tau_T2_54),...
    best_basis_T2_thresholded_54, 'un', 0);
best_basis_T2_thresholded_80_new = cellfun(@(x) max(0,x - tau_T2_80)-max(0,-x-tau_T2_80),...
    best_basis_T2_thresholded_80, 'un', 0);
best_basis_T2_thresholded_50_new = cellfun(@(x) max(0,x - tau_T2_50)-max(0,-x-tau_T2_50),...
    best_basis_T2_thresholded_50, 'un', 0);
fprintf('Thresholding for the best basis with norm energy\n');
% thresholding from the best basis with norm
[tau_T3_54, best_basis_T3_thresholded_54] = sgwpt_thresholding(prc_range_1, best_basis_T3);
[tau_T3_80, best_basis_T3_thresholded_80] = sgwpt_thresholding(prc_range_2, best_basis_T3);
[tau_T3_50, best_basis_T3_thresholded_50] = sgwpt_thresholding(prc_range_3, best_basis_T3);
% substract tau from the coefficients
best_basis_T3_thresholded_54_new = cellfun(@(x) max(0,x - tau_T3_54)-max(0,-x-tau_T3_54),...
    best_basis_T3_thresholded_54, 'un', 0);
best_basis_T3_thresholded_80_new = cellfun(@(x) max(0,x - tau_T3_80)-max(0,-x-tau_T3_80),...
    best_basis_T3_thresholded_80, 'un', 0);
best_basis_T3_thresholded_50_new = cellfun(@(x) max(0,x - tau_T3_50)-max(0,-x-tau_T3_50),...
    best_basis_T3_thresholded_50, 'un', 0);
fprintf('Thresholding for the best basis with log energy\n');
% thresholding from the best basis with log energy
[tau_T4_54, best_basis_T4_thresholded_54] = sgwpt_thresholding(prc_range_1, best_basis_T4);
[tau_T4_80, best_basis_T4_thresholded_80] = sgwpt_thresholding(prc_range_2, best_basis_T4);
[tau_T4_50, best_basis_T4_thresholded_50] = sgwpt_thresholding(prc_range_3, best_basis_T4);
% substract tau from the coefficients
best_basis_T4_thresholded_54_new = cellfun(@(x) max(0,x - tau_T4_54)-max(0,-x-tau_T4_54),...
    best_basis_T4_thresholded_54, 'un', 0);
best_basis_T4_thresholded_80_new = cellfun(@(x) max(0,x - tau_T4_80)-max(0,-x-tau_T4_80),...
    best_basis_T4_thresholded_80, 'un', 0);
best_basis_T4_thresholded_50_new = cellfun(@(x) max(0,x - tau_T4_50)-max(0,-x-tau_T4_50),...
    best_basis_T4_thresholded_50, 'un', 0);
fprintf('Thresholding for the best basis with sure energy\n');
% thresholding from the best basis with sure
[tau_T5_54, best_basis_T5_thresholded_54] = sgwpt_thresholding(prc_range_1, best_basis_T5);
[tau_T5_80, best_basis_T5_thresholded_80] = sgwpt_thresholding(prc_range_2, best_basis_T5);
[tau_T5_50, best_basis_T5_thresholded_50] = sgwpt_thresholding(prc_range_3, best_basis_T5);
% substract tau from the coefficients
best_basis_T5_thresholded_54_new = cellfun(@(x) max(0,x - tau_T5_54)-max(0,-x-tau_T5_54),...
    best_basis_T5_thresholded_54, 'un', 0);
best_basis_T5_thresholded_80_new = cellfun(@(x) max(0,x - tau_T5_80)-max(0,-x-tau_T5_80),...
    best_basis_T5_thresholded_80, 'un', 0);
best_basis_T5_thresholded_50_new = cellfun(@(x) max(0,x - tau_T5_50)-max(0,-x-tau_T5_50),...
    best_basis_T5_thresholded_50, 'un', 0);

% thresholding from the best basis with pnorm
% [tau_T6_54, best_basis_T6_thresholded_54] = sgwpt_thresholding(prc_range_1, best_basis_T6)
% [tau_T6_80, best_basis_T6_thresholded_80] = sgwpt_thresholding(prc_range_2, best_basis_T6)
% [tau_T6_50, best_basis_T6_thresholded_50] = sgwpt_thresholding(prc_range_3, best_basis_T6)
% substract tau from the coefficients
% best_basis_T6_thresholded_54_new = cellfun(@(x) max(0,x - tau_T6_54)-max(0,-x-tau_T6_54),...
%     best_basis_T6_thresholded_54, 'un', 0);
% best_basis_T6_thresholded_80_new = cellfun(@(x) max(0,x - tau_T6_80)-max(0,-x-tau_T6_80),...
%     best_basis_T6_thresholded_80, 'un', 0);
% best_basis_T6_thresholded_50_new = cellfun(@(x) max(0,x - tau_T6_50)-max(0,-x-tau_T6_50),...
%     best_basis_T6_thresholded_50, 'un', 0);
toc

tic
fprintf('Recover the denoised signal from the best basis with Shannon energy (54)\n');
% inverse using Chebyshev polynomial of Laplacian applied to best basis
% thresholded
r4_T1_thresholded_54 = sgwt_inverse(best_basis_T1_thresholded_54_new,L,coef_bb_T1,arange);
% r4_T2_thresholded_54 = sgwt_inverse(best_basis_T2_thresholded_54_new,L,coef_bb_T2,arange);
% r4_T3_thresholded_54 = sgwt_inverse(best_basis_T3_thresholded_54_new,L,coef_bb_T3,arange);
% r4_T4_thresholded_54 = sgwt_inverse(best_basis_T4_thresholded_54_new,L,coef_bb_T4,arange);
% r4_T5_thresholded_54 = sgwt_inverse(best_basis_T5_thresholded_54_new,L,coef_bb_T5,arange);
% %r4_T6_thresholded_54 = sgwt_inverse(best_basis_T6_thresholded_54_new,L,coef_bb_T6,arange);

toc
tic
fprintf('Recover the denoised signal from the best basis with Shannon energy (80)\n');
r4_T1_thresholded_80 = sgwt_inverse(best_basis_T1_thresholded_80_new,L,coef_bb_T1,arange);
% r4_T2_thresholded_80 = sgwt_inverse(best_basis_T2_thresholded_80_new,L,coef_bb_T2,arange);
% r4_T3_thresholded_80 = sgwt_inverse(best_basis_T3_thresholded_80_new,L,coef_bb_T3,arange);
% r4_T4_thresholded_80 = sgwt_inverse(best_basis_T4_thresholded_80_new,L,coef_bb_T4,arange);
% r4_T5_thresholded_80 = sgwt_inverse(best_basis_T5_thresholded_80_new,L,coef_bb_T5,arange);
% %r4_T6_thresholded_80 = sgwt_inverse(best_basis_T6_thresholded_80_new,L,coef_bb_T6,arange);
toc
tic
fprintf('Recover the denoised signal from the best basis with Shannon energy (50)\n');
r4_T1_thresholded_50 = sgwt_inverse(best_basis_T1_thresholded_50_new,L,coef_bb_T1,arange);
% r4_T2_thresholded_50 = sgwt_inverse(best_basis_T2_thresholded_50_new,L,coef_bb_T2,arange);
% r4_T3_thresholded_50 = sgwt_inverse(best_basis_T3_thresholded_50_new,L,coef_bb_T3,arange);
% r4_T4_thresholded_50 = sgwt_inverse(best_basis_T4_thresholded_50_new,L,coef_bb_T4,arange);
% r4_T5_thresholded_50 = sgwt_inverse(best_basis_T5_thresholded_50_new,L,coef_bb_T5,arange);
% %r4_T6_thresholded_50 = sgwt_inverse(best_basis_T6_thresholded_50_new,L,coef_bb_T6,arange);
toc

clear Ad;
clear D;
cd ../data
filename = ('sgwpt_denoising_run_rfMRI_100_wn_norm_spline');
save(filename);

toc
