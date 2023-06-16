%--------------------------------------------------------------------------
% File: example_toronto_graph.m
%
% Goal: compare the SGWPT with the existing methods using as example
% denoising signal on Toronto graph
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
close all
% Example Toronto road graph with traffic volume data as in Irion
% upload the graph and signal with noise
load 'Toronto_SNR7.mat' 
% upload the graph and the signal
G1 = load('Toronto.mat')

A = G.W;
xy = G.xy;
N=size(A,1);
x=xy(:,1);
y=xy(:,2);
Gr = graph(A);
p = plot(Gr,'XData',xy(:,1),'YData',xy(:,2),'MarkerSize',3);
title('Toronto road network')
L=sgwt_laplacian(A,'opt','normalized');
fprintf('Measuring largest eigenvalue, lmax = ');
lmax=sgwt_rough_lmax(L);
fprintf('%g\n',lmax);
arange=[0 lmax];
% the signal with noise
f = G.f;
% the noise
noise = G.f-G1.G.f;
% plot the graph, signal without noise
figure
Gr.Nodes.Value = f-noise;
p_2 = plot(Gr,'XData',xy(:,1),'YData',xy(:,2),'MarkerSize',3)
p_2.NodeCData = Gr.Nodes.Value;
colormap jet
caxis([min(f-noise) max(f-noise)]);
title('signal');
hcb=colorbar('location','north');
cxt=get(hcb,'Xtick');
cxt=[cxt(1),0,cxt(end)];
cpos=get(hcb,'Position');
cpos(4)=.02; % make colorbar thinner
set(hcb,'Position',cpos);
set(hcb,'Position',[.25 .91 .6 .02]);
 saveas(gcf, 'toronto_graph', 'epsc');

% plot the graph, signal+noise
figure
g = f;
Gr.Nodes.Value = g;
p_2 = plot(Gr,'XData',xy(:,1),'YData',xy(:,2),'MarkerSize',3)
p_2.NodeCData = Gr.Nodes.Value;
colormap jet
caxis([min(f-noise) max(f-noise)]);
title('signal+noise');
hcb=colorbar('location','north');
cxt=get(hcb,'Xtick');
cxt=[cxt(1),0,cxt(end)];
cpos=get(hcb,'Position');
cpos(4)=.02; % make colorbar thinner
set(hcb,'Position',cpos);
set(hcb,'Position',[.25 .91 .6 .02]);
 saveas(gcf, 'toronto_graph_noise', 'epsc');

arange=[0 lmax];

%% Chebyshev polynomial approximationf
m=100; % order of polynomial approximation
fprintf('Compute Chebyshev polynomial coefficients for spline filters \n');
[c_0, dp1, dp2] = spline_norm_coeff_wp(m, lmax);

fprintf('Compute Chebyshev polynomial coefficients for spline filters for all levels \n');
[c_0_all, dp1_all, dp2_all] = spline_norm_coeff_wp_all(m, lmax);

%% 
tic
signal = g;

fprintf('Chebyshev polynomial of Laplacian applied to vector for pi(x)\n');
cplv=sgwt_cheby_op(signal,L,c_0,arange);

fprintf('Chebyshev polynomial of Laplacian applied to vector for pi(x)*pj(x)\n');
dplv1=sgwt_cheby_op(signal,L,dp1,arange);

fprintf('Chebyshev polynomial of Laplacian applied to vector for (pi(x)*pj(x))*pi(x) \n');
dplv2=sgwt_cheby_op(signal,L,dp2,arange);

T1 = 'shannon';
val1 = 1;
[best_basis_T1, coef_bb_T1, e1_T1, e2_T1, e3_T1, e4_T1] = ...
    best_basis_frame(N,c_0,dp1,dp2,cplv,dplv1,dplv2,T1,val1);

coef_all_lev = [c_0_all,dp1_all(1,:),dp1_all(2,:),dp2_all(1,:),dp2_all(2,:)];

all_lev = sgwt_cheby_op(signal,L,coef_all_lev,arange);

tau =  16000;   

[cplv_t] = thresholding_coef(tau, cplv);
[dplv1_t] = thresholding_coef(tau, dplv1);
[dplv2_t] = thresholding_coef(tau, dplv2);
[best_basis_T1_t] = thresholding_coef(tau, best_basis_T1);

% inverse using 'Chebyshev polynomial of Laplacian applied to vector for pi(x)
r1 = sgwt_inverse(cplv,L,c_0,arange);
r1_t = sgwt_inverse(cplv_t,L,c_0,arange);
% inverse using Chebyshev polynomial of Laplacian applied to vector for pi(x)*pj(x)
r2 = sgwt_inverse(dplv1,L,dp1,arange);
r2_t = sgwt_inverse(dplv1_t,L,dp1,arange);
% inverse using Chebyshev polynomial of Laplacian applied to vector for (pi(x)*pj(x))*pi(x)
r3 = sgwt_inverse(dplv2,L,dp2,arange);
r3_t = sgwt_inverse(dplv2_t,L,dp2,arange);
% inverse using Chebyshev polynomial of Laplacian applied to best basis
r4_T1 = sgwt_inverse(best_basis_T1,L,coef_bb_T1,arange);
r4_T1_t = sgwt_inverse(best_basis_T1_t,L,coef_bb_T1,arange);
% inverse using Chebyshev polynomial of Laplacian applied to all levels
r_all_lev = sgwt_inverse(all_lev,L,coef_all_lev,arange);

figure; plot(signal); hold on; plot(f);

% plot the graph, reconstructed signal using BB1 with Shannon
figure
Gr.Nodes.Value = r4_T1;
p_2 = plot(Gr,'XData',xy(:,1),'YData',xy(:,2),'MarkerSize',3)
p_2.NodeCData = Gr.Nodes.Value;
colormap jet
caxis([min(f-noise) max(f-noise)]);
title('BB with Shannon');
hcb=colorbar('location','north');
cxt=get(hcb,'Xtick');
cxt=[cxt(1),0,cxt(end)];
cpos=get(hcb,'Position');
cpos(4)=.02; % make colorbar thinner
set(hcb,'Position',cpos);
set(hcb,'Position',[.25 .91 .6 .02]);
%saveas(gcf, 'toronto_graph_BB', 'epsc');

% plot the graph, reconstructed signal using lev_0
figure
Gr.Nodes.Value = r1;
p_2 = plot(Gr,'XData',xy(:,1),'YData',xy(:,2),'MarkerSize',3)
p_2.NodeCData = Gr.Nodes.Value;
colormap jet
caxis([min(f-noise) max(f-noise)]);
title('lev_0');
hcb=colorbar('location','north');
cxt=get(hcb,'Xtick');
cxt=[cxt(1),0,cxt(end)];
cpos=get(hcb,'Position');
cpos(4)=.02; % make colorbar thinner
set(hcb,'Position',cpos);
set(hcb,'Position',[.25 .91 .6 .02]);
%saveas(gcf, 'toronto_graph_lev0', 'epsc');

% plot the graph, reconstructed signal using lev_1
figure
Gr.Nodes.Value = r2;
p_2 = plot(Gr,'XData',xy(:,1),'YData',xy(:,2),'MarkerSize',3)
p_2.NodeCData = Gr.Nodes.Value;
colormap jet
caxis([min(f-noise) max(f-noise)]);
title('lev_1');
hcb=colorbar('location','north');
cxt=get(hcb,'Xtick');
cxt=[cxt(1),0,cxt(end)];
cpos=get(hcb,'Position');
cpos(4)=.02; % make colorbar thinner
set(hcb,'Position',cpos);
set(hcb,'Position',[.25 .91 .6 .02]);
%saveas(gcf, 'toronto_graph_lev1', 'epsc');

% plot the graph, reconstructed signal using lev_2
figure
Gr.Nodes.Value = r3;
p_2 = plot(Gr,'XData',xy(:,1),'YData',xy(:,2),'MarkerSize',3)
p_2.NodeCData = Gr.Nodes.Value;
colormap jet
caxis([min(f-noise) max(f-noise)]);
title('lev_2');
hcb=colorbar('location','north');
cxt=get(hcb,'Xtick');
cxt=[cxt(1),0,cxt(end)];
cpos=get(hcb,'Position');
cpos(4)=.02; % make colorbar thinner
set(hcb,'Position',cpos);
set(hcb,'Position',[.25 .91 .6 .02]);
%aveas(gcf, 'toronto_graph_lev2', 'epsc');

% plot the graph, reconstructed signal using all lev
figure
Gr.Nodes.Value = r_all_lev;
p_2 = plot(Gr,'XData',xy(:,1),'YData',xy(:,2),'MarkerSize',3)
p_2.NodeCData = Gr.Nodes.Value;
colormap jet
caxis([min(f-noise) max(f-noise)]);
title('all lev');
hcb=colorbar('location','north');
cxt=get(hcb,'Xtick');
cxt=[cxt(1),0,cxt(end)];
cpos=get(hcb,'Position');
cpos(4)=.02; % make colorbar thinner
set(hcb,'Position',cpos);
set(hcb,'Position',[.25 .91 .6 .02]);
% saveas(gcf, 'toronto_graph_all_lev', 'epsc');

%%
% plot the graph, reconstructed signal using BB1 with Shannon
figure
Gr.Nodes.Value = r4_T1_t;
p_2 = plot(Gr,'XData',xy(:,1),'YData',xy(:,2),'MarkerSize',3)
p_2.NodeCData = Gr.Nodes.Value;
colormap jet
caxis([min(f-noise) max(f-noise)]);
title('BB with Shannon and thresholding');
hcb=colorbar('location','north');
cxt=get(hcb,'Xtick');
cxt=[cxt(1),0,cxt(end)];
cpos=get(hcb,'Position');
cpos(4)=.02; % make colorbar thinner
set(hcb,'Position',cpos);
set(hcb,'Position',[.25 .91 .6 .02]);
saveas(gcf, 'toronto_graph_BB_thresh', 'epsc');

% plot the graph, reconstructed signal using lev_0
figure
Gr.Nodes.Value = r1_t;
p_2 = plot(Gr,'XData',xy(:,1),'YData',xy(:,2),'MarkerSize',3)
p_2.NodeCData = Gr.Nodes.Value;
colormap jet
caxis([min(f-noise) max(f-noise)]);
title('lev_0 and thresholding');
hcb=colorbar('location','north');
cxt=get(hcb,'Xtick');
cxt=[cxt(1),0,cxt(end)];
cpos=get(hcb,'Position');
cpos(4)=.02; % make colorbar thinner
set(hcb,'Position',cpos);
set(hcb,'Position',[.25 .91 .6 .02]);
saveas(gcf, 'toronto_graph_lev0_thresh', 'epsc');

% plot the graph, reconstructed signal using lev_1
figure
Gr.Nodes.Value = r2_t;
p_2 = plot(Gr,'XData',xy(:,1),'YData',xy(:,2),'MarkerSize',3)
p_2.NodeCData = Gr.Nodes.Value;
colormap jet
caxis([min(f-noise) max(f-noise)]);
title('lev_1 and thresholding');
hcb=colorbar('location','north');
cxt=get(hcb,'Xtick');
cxt=[cxt(1),0,cxt(end)];
cpos=get(hcb,'Position');
cpos(4)=.02; % make colorbar thinner
set(hcb,'Position',cpos);
set(hcb,'Position',[.25 .91 .6 .02]);
saveas(gcf, 'toronto_graph_lev1_thresh', 'epsc');

% plot the graph, reconstructed signal using lev_2
figure
Gr.Nodes.Value = r3_t;
p_2 = plot(Gr,'XData',xy(:,1),'YData',xy(:,2),'MarkerSize',3)
p_2.NodeCData = Gr.Nodes.Value;
colormap jet
caxis([min(f-noise) max(f-noise)]);
title('lev_2 and thresholding');
hcb=colorbar('location','north');
cxt=get(hcb,'Xtick');
cxt=[cxt(1),0,cxt(end)];
cpos=get(hcb,'Position');
cpos(4)=.02; % make colorbar thinner
set(hcb,'Position',cpos);
set(hcb,'Position',[.25 .91 .6 .02]);
saveas(gcf, 'toronto_graph_lev2_thresh', 'epsc');
toc

%% SGWT
tic
Nscales=5;
fprintf('%g\n',lmax);
fprintf('Designing transform in spectral domain\n');

designtype='abspline3';

[g_sgwt,t]=sgwt_filter_design(lmax,Nscales,'designtype',designtype);

% compute Chebyshev coefficients for p(x) = sum c(1+k) T_k(x); 0<=K<=M
for k=1:numel(g_sgwt)
    c_sgwt{k}=sgwt_cheby_coeff(g_sgwt{k},m,m+1,arange);
end

cplv_sgwt=sgwt_cheby_op(signal,L,c_sgwt,arange);
r_sgwt = sgwt_inverse(cplv_sgwt,L,c_sgwt,arange);

[cplv_sgwt_t] = thresholding_coef(tau, cplv_sgwt);
r_sgwt_t = sgwt_inverse(cplv_sgwt_t,L,c_sgwt,arange);
% plot the graph, reconstructed signal using BB1 with Shannon
figure
Gr.Nodes.Value = r_sgwt;
p_2 = plot(Gr,'XData',xy(:,1),'YData',xy(:,2),'MarkerSize',3)
p_2.NodeCData = Gr.Nodes.Value;
colormap jet
caxis([min(f-noise) max(f-noise)]);
title('SGWT');
hcb=colorbar('location','north');
cxt=get(hcb,'Xtick');
cxt=[cxt(1),0,cxt(end)];
cpos=get(hcb,'Position');
cpos(4)=.02; % make colorbar thinner
set(hcb,'Position',cpos);
set(hcb,'Position',[.25 .91 .6 .02]);

figure
Gr.Nodes.Value = r_sgwt_t;
p_2 = plot(Gr,'XData',xy(:,1),'YData',xy(:,2),'MarkerSize',3)
p_2.NodeCData = Gr.Nodes.Value;

colormap jet
caxis([min(f-noise) max(f-noise)]);
title('SGWT thresholded');
hcb=colorbar('location','north');
cxt=get(hcb,'Xtick');
cxt=[cxt(1),0,cxt(end)];
cpos=get(hcb,'Position');
cpos(4)=.02; % make colorbar thinner
set(hcb,'Position',cpos);
set(hcb,'Position',[.25 .91 .6 .02]);
saveas(gcf, 'toronto_graph_sgwt_thresh', 'epsc');

toc
%%
tic;[r1_ord, SNR_r1_coef] = SNR_coef(arange,L,cplv,c_0,signal,noise); toc;
tic;[r2_ord, SNR_r2_coef] = SNR_coef(arange,L,dplv1,dp1,signal,noise);toc;
tic;[r3_ord, SNR_r3_coef] = SNR_coef(arange,L,dplv2,dp2,signal,noise);toc;
tic;[r4_ord, SNR_r4_coef] = SNR_coef(arange,L,best_basis_T1,coef_bb_T1,signal,noise);toc;
tic;[r5_ord, SNR_r5_coef] = SNR_coef(arange,L,cplv_sgwt,c_sgwt,signal,noise);toc;
tic;[r6_ord, SNR_r6_coef] = SNR_coef(arange,L,all_lev,coef_all_lev,signal,noise);toc;
%
figure;
plot(r1_ord, SNR_r1_coef,'LineWidth',2);  hold on; plot(r2_ord, SNR_r2_coef,'LineWidth',2);
plot(r3_ord, SNR_r3_coef,'LineWidth',2); plot(r4_ord, SNR_r4_coef,'LineWidth',2);
plot(r5_ord, SNR_r5_coef,'LineWidth',2); %plot(r6_ord, SNR_r6_coef,'LineWidth',2);
yline( 8.96,'LineWidth',2);
legend('lev_0 (SGWPT)', 'lev_1 (SGWPT)', 'lev_2 (SGWPT)', 'BB (SGWPT)',...
    'SGWT','HGLET Best Basis', 'Location', 'Best')
ylabel('SNR');
%title('SNR');
xlabel('Threshold value');

 saveas(gcf, 'snr', 'epsc');
% saveas(gcf, 'snr', 'png');


[max_r1,pos_r1]=max(SNR_r1_coef);
[max_r2,pos_r2]=max(SNR_r2_coef);
[max_r3,pos_r3]=max(SNR_r3_coef);
[max_r4,pos_r4]=max(SNR_r4_coef);
[max_r5,pos_r5]=max(SNR_r5_coef);
[max_r6,pos_r6]=max(SNR_r6_coef);

max_1 = max(max( max( max(max_r1,max_r2), max(max_r3,max_r4) ), max_r5),max_r6);