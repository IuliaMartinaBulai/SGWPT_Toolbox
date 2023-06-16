%--------------------------------------------------------------------------
% File: SNR_coef.m
%
% Goal: 
%
% Use: [tau, SNR_coef] = SNR_coef(arange, L, best_basis, coef_bb, signal, noise)
%
% Inputs:  arange - 
%          L -
%          best_basis - input vector
%          coef_bb - 
%          signal - 
%          noise -
%
% Outputs: tau - threshold value
%          SNR_coef - 
%
% Recalls: cell2mat.m, reshape.m, sort.m, max.m, thresholding_coef.m,
% sgwt_inverse.m, SNR.m
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
function [tau, SNR_coef] = SNR_coef(arange,L,best_basis,coef_bb,signal,noise)

a = cell2mat(best_basis);
best_basis_col = reshape(a ,[],1);
% compute tau (the threshold value)
best_basis_ord = sort(best_basis_col);
tau = [0:1000:max(best_basis_ord)]';
for i=1:size(tau)
    best_basis_t_coef = thresholding_coef(tau(i),best_basis);
    r_thresholded_coef = sgwt_inverse(best_basis_t_coef,L,coef_bb,arange);
    SNR_t_coef(i) = SNR(r_thresholded_coef, signal-noise);
end

SNR_coef = SNR_t_coef';