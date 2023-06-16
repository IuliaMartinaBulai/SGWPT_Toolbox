%--------------------------------------------------------------------------
% File:  best_basis_frame.m
%
% Goal: obtaine the best basis frame for a given basis
%
% Use: [best_basis, coef_bb, e1, e2, e3, e4] = best_basis_frame(N, c, dp1,...
%                                              dp2, cplv, dplv1, dplv2, T, val)
%
% Inputs:  c_0 - Chebyshev coefficients for h0 and g0
%          dp1 - Chebyshev coefficients for  h1h0, g1h0, h1g0, g1g0
%          dp2 - Chebyshev coefficients for g2h1h0, h2g1h0, g2g1h0, h2h1g0,
%          g2h1g0, h2g1g0, g2g1g0
%          cplv - Chebyshev coefficients of Laplacian applied to vector for
%                 h0 and g0
%          dplv1 - Chebyshev coefficients of Laplacian applied to vector
%                  for h1h0, g1h0, h1g0, g1g0
%          dplv2 - Chebyshev coefficients of Laplacian applied to vector for
%                  g2h1h0, h2g1h0, g2g1h0, h2h1g0, g2h1g0, h2g1g0, g2g1g0
%          T - type of entropy T1 = 'shannon'; T2 = 'threshold'; T3 = 'norm';
%              T4 = 'log energy'; T5 = 'sure'; T6 = 'pnorm';
%          val - value associated to the chosen entropy
%
% Outputs: best_basis - the best basis obtained in terms of the chosen
%                       entropy measure
%         coef_bb - the coefficients ordered used to recontruct the signal
%                   from best_basis
%         e1 - returns the entropy specified by T1 of the vector in input
%         e2 - returns the entropy specified by T2 of the vector in input
%         e3 - returns the entropy specified by T3 of the vector in input
%         e4 - returns the entropy specified by T4 of the vector in input
%
% Recalls: cat.m, wentropy.m, numel.m, cell2mat.m, iscell.m
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
function [best_basis, coef_bb, e1, e2, e3, e4] = ...
    best_basis_frame(N,c,dp1,dp2,cplv,dplv1,dplv2,T,val)
f1 = {};
f2 = {};
f3 = {};
dp2_temp = [];
coef_ric_1 = [];
coef_ric_2 = [];
for i = 1: numel(dplv1)
    % entropy on the concatenation of vectors
    f1 = {cat(1, dplv2{2*i-1}, dplv2{2*i})};
    dp2_temp = {cat(1, dp2{2*i-1}, dp2{2*i})};
    e1(i) = wentropy(f1{:},T,val);
    e2(i) = wentropy(dplv1{1,i},T,val);
    if e1(i)<e2(i)
        f2{i} = f1;
        coef_ric_1{i} = dp2_temp{1,1}; % coef_ric_1{i}= dp2_temp;
    else
        f2{i} = dplv1{1,i};
        coef_ric_1{i} = dp1{i};
    end
end
for i = 1:numel(f2)
    if iscell(f2{i})
        f2{i} = cell2mat(f2{i});
    end
end
for i = 1:numel(cplv)
    % entropy on the matrix
    f3 = {cat(1, f2{2*i-1}, f2{2*i})};
    coef_ric_2_temp = {cat(1, coef_ric_1{2*i-1}, coef_ric_1{2*i})};
    e3(i) = wentropy(f3{:},T, val);
    e4(i) = wentropy(cplv{1,i},T,val);
    if e3(i)<e4(i)
        bb{i} = f3;
        coef_ric_2{i} = coef_ric_2_temp{1,1}; % coef_ric_2{i}=coef_ric_2_temp;
    else
        bb{i} = cplv{1,i};
        coef_ric_2{i} = c{i};
    end
end
for i = 1:numel(bb)
    if iscell(bb{i})
        bb{i} = cell2mat(bb{i});
    end
end
best_basis = {};
coef_bb = {};
for i = 1:numel(bb)
    best_basis_1 = {};
    coef_1 = {};
    index = numel(bb{i})/N;
    if index~=1
        for j = 1:index
            best_basis_1{j} = bb{i}(N*j-(N-1):N*j);
            coef_1{j} = coef_ric_2{1,i}(j,:);
        end
    else
        best_basis_1 = bb{i};
        coef_1 = coef_ric_2{1,i};
    end
    best_basis = [best_basis,best_basis_1];
    coef_bb = [coef_bb,coef_1];
end
[best_basis, coef_bb, e1, e2, e3, e4];
end
