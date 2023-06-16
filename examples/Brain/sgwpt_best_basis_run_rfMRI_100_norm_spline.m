%--------------------------------------------------------------------------
% File:  sgwpt_best_basis_run_rfMRI_100_norm_spline.m
%
% Goal: Main file for the reconstruction of the normalized signal rfMRI_100 
% with lev0, lev1, lev2, and the five best basis using spline filters
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
m = 50;
Ad = importdata('A_wb.mat');
% number of vertices 
N = length(Ad);
D = sum(Ad);
% Normalized symmetric matrix 
Dn = 1./sqrt(D);
Dn = spdiags(Dn.',0,spalloc(N,N,N));
L = speye(N)-Dn*Ad*Dn;
lmax=sgwt_rough_lmax(L);
fprintf('%g\n',lmax);
arange = [0 lmax];
tic
fprintf('Compute Chebyshev polynomial coefficients for spline filters \n');
[c_0, dp1, dp2] = spline_norm_coeff_wp(m, lmax);
toc

tic
fprintf('Load and normalize signal \n');
% load the eigenvectors
load('WB_Utr_1000_short.mat')
% the rfMRI_100 was taken from CorrectedGSR_Lambda9_102513_session_rfMRI_206_smoothing6;
load('rfMRI_100.mat')

Utr = eigenvector;
rfMRI = rfMRI_100-Utr'*rfMRI_100*Utr;
rfMRI_norm = (rfMRI)./norm(rfMRI,2); % the signal normalized as anjali

toc 

tic
fprintf('Chebyshev polynomial of Laplacian applied to vector for pi(x)\n');
cplv = sgwt_cheby_op(rfMRI_norm,L,c_0,arange);

fprintf('Chebyshev polynomial of Laplacian applied to vector for pi(x)*pj(x)\n');
dplv1 = sgwt_cheby_op(rfMRI_norm,L,dp1,arange);

fprintf('Chebyshev polynomial of Laplacian applied to vector for pi(x)*pj(x))*pi(x) \n');
dplv2 = sgwt_cheby_op(rfMRI_norm,L,dp2,arange);
toc
T1 = 'shannon';
val1 = 0;
T2 = 'threshold';
val2 = sqrt(2*log(N))/sqrt(N);
T3 = 'norm';
val3 = 2;
T4 = 'log energy';
val4 = 0;
T5 = 'sure';
val5 = sqrt(2*log(N))/sqrt(N);
% pnorm was deleted from Matlab 
% T6 = 'pnorm';
% val6 = 0.5;
tic
fprintf('Apply the best basis algorithm for Shannon energy \n');
% best basis with 'Shannon'
[best_basis_T1, coef_bb_T1, e1_T1, e2_T1, e3_T1, e4_T1] = ...
    best_basis_frame(N,c_0,dp1,dp2,cplv,dplv1,dplv2,T1,val1);
fprintf('Apply the best basis algorithm for threshold energy \n');
% best basis with 'threshold'
[best_basis_T2, coef_bb_T2, e1_T2, e2_T2, e3_T2, e4_T2] = ...
    best_basis_frame(N,c_0,dp1,dp2,cplv,dplv1,dplv2,T2,val2);
fprintf('Apply the best basis algorithm for norm energy \n');
% best basis with 'norm'
[best_basis_T3, coef_bb_T3, e1_T3, e2_T3, e3_T3, e4_T3] = ...
    best_basis_frame(N,c_0,dp1,dp2,cplv,dplv1,dplv2,T3,val3);
fprintf('Apply the best basis algorithm for log energy \n');
% best basis with 'log energy'
[best_basis_T4, coef_bb_T4, e1_T4, e2_T4, e3_T4, e4_T4] = ...
    best_basis_frame(N,c_0,dp1,dp2,cplv,dplv1,dplv2,T4,val4);
fprintf('Apply the best basis algorithm for sure energy \n');
% best basis with 'sure'
[best_basis_T5, coef_bb_T5, e1_T5, e2_T5, e3_T5, e4_T5] =...
    best_basis_frame(N,c_0,dp1,dp2,cplv,dplv1,dplv2,T5,val5);
% fprintf('Apply the best basis algorithm for pnorm energy \n');
% % best basis with 'pnorm' 
% [best_basis_T6, coef_bb_T6, e1_T6, e2_T6, e3_T6, e4_T6] = ...
%     best_basis_frame(N,c_0,dp1,dp2,cplv,dplv1,dplv2,T6,val6);
toc
tic
fprintf('Recover the signal for pi(x) \n');
% inverse using 'Chebyshev polynomial of Laplacian applied to vector for pi(x)
r1 = sgwt_inverse(cplv,L,c_0,arange);
fprintf('Recover the signal for pi(x)*pj(x) \n');
% inverse using Chebyshev polynomial of Laplacian applied to vector for pi(x)*pj(x)
r2 = sgwt_inverse(dplv1,L,dp1,arange);
fprintf('Recover the signal for pi(x)*pj(x))*pi(x) \n');
% inverse using Chebyshev polynomial of Laplacian applied to vector for pi(x)*pj(x))*pi(x)
r3 = sgwt_inverse(dplv2,L,dp2,arange);
tic
fprintf('Recover the signal for the Shannon best basis \n');
% inverse using Chebyshev polynomial of Laplacian applied to best basis
r4_T1 = sgwt_inverse(best_basis_T1,L,coef_bb_T1,arange);
fprintf('Recover the signal for the threshold best basis \n');
r4_T2 = sgwt_inverse(best_basis_T2,L,coef_bb_T2,arange);
fprintf('Recover the signal for the norm best basis \n');
r4_T3 = sgwt_inverse(best_basis_T3,L,coef_bb_T3,arange);
fprintf('Recover the signal for the log energy best basis \n');
r4_T4 = sgwt_inverse(best_basis_T4,L,coef_bb_T4,arange);
fprintf('Recover the signal for the sure best basis \n');
r4_T5 = sgwt_inverse(best_basis_T5,L,coef_bb_T5,arange);
%r4_T6 = sgwt_inverse(best_basis_T6,L,coef_bb_T6,arange);
toc
clear Ad;
clear D;
cd ../data
filename = ('sgwpt_best_basis_run_rfMRI_100_norm_spline');
save(filename);
toc