%--------------------------------------------------------------------------
% File: sgwpt_thresholding:
%
% Goal: Given a percentage of coefficient to be tresholded and a vector,
%       find the thresholding coefficient and the tresholded vector
%
% Use: [tau, best_basis_thresholded] = sgwpt_thresholding(prc_range, best_basis)
%
% Inputs: best_basis - input vector
%         prc_range - percentage of thresholding
%
% Outputs: tau - threshold value
%          best_basis_thresholded - thresholded vector
%
% Recalls: cell2mat.m, reshape.m, prctile.m, abs.m
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
function [tau, best_basis_thresholded] = sgwpt_thresholding(prc_range, best_basis)
% convert the cell to mat and then to a collumn vector
a = cell2mat(best_basis);
best_basis_col = reshape(a ,[],1);
% compute tau (the threshold value)
tau = prctile(abs(best_basis_col),prc_range);
% do the thresholding
for i=1:numel(best_basis)
    best_basis_thresholded{i} = best_basis{i};
    best_basis_thresholded{i}( abs(best_basis_thresholded{i}) <= tau ) = 0;
end
end