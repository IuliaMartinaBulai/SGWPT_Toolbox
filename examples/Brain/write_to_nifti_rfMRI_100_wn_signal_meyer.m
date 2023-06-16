%--------------------------------------------------------------------------
% File: write_to_nifti_rfMRI_100_wn_signal_meyer.m
%
% Goal: File to write to nifti a signal for the case rfMRI_100_wn_signal_meyer
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
clc
clear all
load('sgwpt_denoising_run_rfMRI_100_wn_norm_meyer') 
load('indices_wb.mat')
filename_in = '102513_T1w_acpc_dc_restore_1.25.nii';

% original signal
xsignal = rfMRI_noise_norm; 
filename = 'rfMRI_100_wn_meyer.nii';
WriteToNifti(filename_in, xsignal,indices_wb, filename)

% signal recontructed using lev0 coefficients
xsignal = r1; 
filename = 'r1_rfMRI_100_wn_meyer.nii';
WriteToNifti(filename_in, xsignal,indices_wb, filename)

% signal recontructed using lev1 coefficients
xsignal = r2; 
filename = 'r2_rfMRI_100_wn_meyer.nii';
WriteToNifti(filename_in, xsignal,indices_wb, filename)

% signal recontructed using lev2 coefficients
xsignal = r3; 
filename = 'r3_rfMRI_100_wn_meyer.nii';
WriteToNifti(filename_in, xsignal,indices_wb, filename)

% signal recontructed using best basis Shannon coefficients
xsignal = r4_T1; 
filename = 'r4_rfMRI_100_wn_meyer_shannon.nii';
WriteToNifti(filename_in, xsignal,indices_wb, filename)

% signal recontructed using 80% cut best basis Shannon coefficients 
xsignal = r4_T1_thresholded_80; 
filename = 'r4_rfMRI_100_wn_meyer_shannon_thresholded_80.nii';
WriteToNifti(filename_in, xsignal,indices_wb, filename)

% signal recontructed using 50% cut best basis Shannon coefficients 
xsignal = r4_T1_thresholded_50; 
filename = 'r4_rfMRI_100_wn_meyer_shannon_thresholded_50.nii';
WriteToNifti(filename_in, xsignal,indices_wb, filename)

% signal recontructed using 33% cut best basis Shannon coefficients 
xsignal = r4_T1_thresholded_33; 
filename = 'r4_rfMRI_100_wn_meyer_shannon_thresholded_33.nii';
WriteToNifti(filename_in, xsignal,indices_wb, filename)
