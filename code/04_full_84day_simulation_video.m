% =============================================================================
% 04_full_84day_simulation_video.m
% Project : Hippocampal CA1 Scaffold Remodeling (BSE 639 - Tissue Engineering, Group 3, IIT Indore)
% Full 84-day model: 22^3 grid, 1000 Gaussian depots, astrocytes (AST), endothelial cells (ECL), Wnt/BMP schedule, daily slice images, MP4 video.
% Requires: MATLAB R2020b or newer, Parallel Computing Toolbox and a CUDA GPU (gpuArray). Windows path: uses %USERPROFILE%.
% Output  : Day_001.png ... Day_084.png and FullScaffold_84days_WntBMP_GPU_Final.mp4 in <Downloads>/tissue_images
% =============================================================================

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OPTIMIZED 84-DAY CA1 SIMULATION
% - AST/ECL corrected according to depot ratio
% - N_global grows correctly
% - Scaffold persists
% - Daily images + final video
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc; clear; close all;
%% ----------------- CUDA FORWARD COMPATIBILITY --------------------------
parallel.gpu.enableCUDAForwardCompatibility(true);
%% ----------------- SETTINGS --------------------------------------------
Ng = 22;                
L  = 1.0;
[x,y,z] = meshgrid(linspace(0,L,Ng));
dx = x(1,2,1)-x(1,1,1);
t_final = 84;          
substeps_per_day = 10;  
dt_day = 1.0/substeps_per_day; 
nDays = t_final;
%% ----------------- BIOPHYSICAL PARAMETERS ------------------------------
D_T_min = 0.004; D_T_max = 0.02;
D_m_min = 0.003; D_m_max = 0.015;
D_s_min = 0.001; D_s_max = 0.006;
D_E_min = 0.002; D_E_max = 0.01;
D_P_min = 0.002; D_P_max = 0.008;
D_ast_min = 0.0008; D_ast_max = 0.006;
D_ecl_min = 0.001; D_ecl_max = 0.007;
R0=1.0; mu=0.25; lambdaT=0.10;
alpha_P=0.3; alpha_A=0.25; beta_M=0.20; gamma_T=0.40;
k_d=0.05; % reduced scaffold degradation
k_E=0.10; k_mE=0.05;
r_P=0.50; P_max=0.8; k_diff=1.0; % increased neuron growth
delta_N=0.01;
r_ast = 0.20; r_ecl = 0.25; k_ast = 0.10; k_ecl = 0.12;
k_w_p = 0.08;
BMP_base = 0.02; BMP_final = 0.18;
bmp_ramp_start_day = 60; bmp_ramp_end_day = 84;
%% ----------------- DEPOTS (1000 Gaussians) -----------------------------
numDepots = 1000; sigma = 0.02*L;
rng(10); depots = L*rand(numDepots,3);
amps_T  = 0.5 + 1.5 * rand(numDepots,1);
amps_M  = 0.7 + 1.6 * rand(numDepots,1);
amps_S  = 0.3 + 1.2 * rand(numDepots,1);
amps_E  = 0.4 + 1.0 * rand(numDepots,1);
amps_P  = 0.4 + 1.0 * rand(numDepots,1);
amps_AST = 0.8 + 1.2 * rand(numDepots,1);
amps_ECL = 0.8 + 1.2 * rand(numDepots,1);
AST_n = round(numDepots*(1/61));
ECL_n = round(numDepots*(10/61));
perm = randperm(numDepots);
AST_idx = perm(1:AST_n);
ECL_idx = perm(AST_n+1:AST_n+ECL_n);
% GPU arrays
DepotField_all = gpuArray.zeros(Ng,Ng,Ng);
DepotField_T   = gpuArray.zeros(Ng,Ng,Ng);
DepotField_M   = gpuArray.zeros(Ng,Ng,Ng);
DepotField_S   = gpuArray.zeros(Ng,Ng,Ng);
DepotField_E   = gpuArray.zeros(Ng,Ng,Ng);
DepotField_P   = gpuArray.zeros(Ng,Ng,Ng);
DepotField_AST = gpuArray.zeros(Ng,Ng,Ng);
DepotField_ECL = gpuArray.zeros(Ng,Ng,Ng);
fprintf('Building %d Gaussian depots on GPU...\n', numDepots);
for i = 1:numDepots
   dx_i = gpuArray(x - depots(i,1));
   dy_i = gpuArray(y - depots(i,2));
   dz_i = gpuArray(z - depots(i,3));
   g = exp(-(dx_i.^2 + dy_i.^2 + dz_i.^2)/(2*sigma^2));
   DepotField_all = DepotField_all + g;
   DepotField_T   = DepotField_T   + amps_T(i)  * g;
   DepotField_M   = DepotField_M   + amps_M(i)  * g;
   DepotField_S   = DepotField_S   + amps_S(i)  * g;
   DepotField_E   = DepotField_E   + amps_E(i)  * g;
   DepotField_P   = DepotField_P   + amps_P(i)  * g;
   if ismember(i, AST_idx), DepotField_AST = DepotField_AST + amps_AST(i)*g; end
   if ismember(i, ECL_idx), DepotField_ECL = DepotField_ECL + amps_ECL(i)*g; end
end
% Normalize depot fields relative to global max to preserve differences
maxDepot = max(DepotField_all(:));
DepotField_T   = DepotField_T / maxDepot;
DepotField_M   = DepotField_M / maxDepot;
DepotField_S   = DepotField_S / maxDepot;
DepotField_E   = DepotField_E / maxDepot;
DepotField_P   = DepotField_P / maxDepot;
DepotField_AST = DepotField_AST / maxDepot;
DepotField_ECL = DepotField_ECL / maxDepot;
DepotField     = DepotField_all / maxDepot;
%% ----------------- DIFFUSION MATRICES ---------------------------------
D_T_matrix   = D_T_min + (D_T_max - D_T_min).*DepotField_T;
D_m_matrix   = D_m_min + (D_m_max - D_m_min).*DepotField_M;
D_s_matrix   = D_s_min + (D_s_max - D_s_min).*DepotField_S;
D_E_matrix   = D_E_min + (D_E_max - D_E_min).*DepotField_E;
D_P_matrix   = D_P_min + (D_P_max - D_P_min).*DepotField_P;
D_ast_matrix = D_ast_min + (D_ast_max - D_ast_min).*DepotField_AST;
D_ecl_matrix = D_ecl_min + (D_ecl_max - D_ecl_min).*DepotField_ECL;
%% ----------------- Wnt field ------------------------------------------
W = exp(-((x-0.5).^2 + (y-0.5).^2 + (z-0.5).^2)/0.06);
W = W ./ max(W(:));
W = gpuArray(W);
%% ----------------- INITIAL CONDITIONS ---------------------------------
T = gpuArray.zeros(Ng,Ng,Ng);
M = gpuArray(0.4 + 0.05*randn(Ng,Ng,Ng)); M(M<0)=0;
S = gpuArray.ones(Ng,Ng,Ng);
E = gpuArray.zeros(Ng,Ng,Ng);
P = gpuArray(0.05 + 0.2*exp(-((x-0.5).^2 + (y-0.5).^2 + (z-0.5).^2)/0.08));
P(P>P_max) = P_max;
AST = max(P./50,1e-6);
ECL = max(P./5,1e-6);
N_global = 0;
%% ----------------- IMAGE SETTINGS -----------------------------
imageFolder = fullfile(getenv('USERPROFILE'),'Downloads','tissue_images');
if ~exist(imageFolder,'dir'), mkdir(imageFolder); end
slice_pos = L*[0.1 0.3 0.5 0.7 0.9];
%% ----------------- MAIN SIMULATION LOOP -----------------------------
fprintf('Starting 84-day simulation...\n');
tic;
for day = 1:nDays
   for sub = 1:substeps_per_day
       t = day-1 + (sub-0.5)/substeps_per_day; % current time in days
      
       % BMP schedule
       if t < bmp_ramp_start_day
           B_scale = BMP_base;
       elseif t <= bmp_ramp_end_day
           frac = (t - bmp_ramp_start_day)/(bmp_ramp_end_day - bmp_ramp_start_day);
           B_scale = BMP_base + frac*(BMP_final - BMP_base);
       else
           B_scale = BMP_final;
       end
      
       % Laplacians
       LapT  = del2(T, dx, dx, dx);
       LapM  = del2(M, dx, dx, dx);
       LapS  = del2(S, dx, dx, dx);
       LapE  = del2(E, dx, dx, dx);
       LapP  = del2(P, dx, dx, dx);
       LapAST = del2(AST, dx, dx, dx);
       LapECL = del2(ECL, dx, dx, dx);
      
       % --- TIMP
       [dTx,dTy,dTz] = gradient(T, dx);
       [dDTx,dDTy,dDTz] = gradient(D_T_matrix, dx);
       diffusion_T = D_T_matrix.*LapT + dDTx.*dTx + dDTy.*dTy + dDTz.*dTz;
       Source_T = R0 * exp(-0.25*t).*DepotField_all;
       T = T + dt_day*(diffusion_T + Source_T - lambdaT.*T); T(T<0)=0;
      
       % --- MMP
       [dMx,dMy,dMz] = gradient(M, dx);
       [dDm_x,dDm_y,dDm_z] = gradient(D_m_matrix, dx);
       diffusion_M = D_m_matrix.*LapM + dDm_x.*dMx + dDm_y.*dMy + dDm_z.*dMz;
       M = M + dt_day*(diffusion_M + alpha_P.*P + alpha_A - beta_M.*M - gamma_T.*(T.*M));
       M(M<0)=0;
      
       % --- Scaffold (with reduced degradation)
       [dSx,dSy,dSz] = gradient(S, dx);
       [dDs_x,dDs_y,dDs_z] = gradient(D_s_matrix, dx);
       S = S + dt_day*(D_s_matrix.*LapS + dDs_x.*dSx + dDs_y.*dSy + dDs_z.*dSz - k_d.*(min(M,1).*S));
       S = min(max(S,0),1);
      
       % --- ECM
       [dEx,dEy,dEz] = gradient(E, dx);
       [dDe_x,dDe_y,dDe_z] = gradient(D_E_matrix, dx);
       E = E + dt_day*(D_E_matrix.*LapE + dDe_x.*dEx + dDe_y.*dEy + dDe_z.*dEz + k_E.*P - k_mE.*(M.*E));
       E(E<0)=0;
      
       % --- NPC
       [dPx,dPy,dPz] = gradient(P, dx);
       [dDp_x,dDp_y,dDp_z] = gradient(D_P_matrix, dx);
       P = P + dt_day*(D_P_matrix.*LapP + dDp_x.*dPx + dDp_y.*dPy + dDp_z.*dPz + ...
           r_P.*P.*(1-P./P_max) + k_w_p.*W - 0.5*B_scale.*P - k_diff.*(S.*P));
       P = min(max(P,0),P_max);
      
       % --- AST (driven by depot amplitude)
       [dASTx,dASTy,dASTz] = gradient(AST, dx);
       [dDa_x,dDa_y,dDa_z] = gradient(D_ast_matrix, dx);
       AST = AST + dt_day*(D_ast_matrix.*LapAST + dDa_x.*dASTx + dDa_y.*dASTy + dDa_z.*dASTz + r_ast.*DepotField_AST - k_ast.*(M.*AST));
       AST(AST<0)=0;
      
       % --- ECL (driven by depot amplitude)
       [dECLx,dECLy,dECLz] = gradient(ECL, dx);
       [dDecl_x,dDecl_y,dDecl_z] = gradient(D_ecl_matrix, dx);
       ECL = ECL + dt_day*(D_ecl_matrix.*LapECL + dDecl_x.*dECLx + dDecl_y.*dECLy + dDecl_z.*dECLz + r_ecl.*DepotField_ECL - k_ecl.*(M.*ECL));
       ECL(ECL<0)=0;
      
       % --- Neurons (scaled for Ng^3)
       SP = S.*P; Integral_SP = sum(SP(:))*dx^3;
       N_global = N_global + dt_day*(k_diff*Integral_SP*(Ng^3) - delta_N*N_global);
       N_global = max(N_global,0);
   end
  
   % --- Save daily image
   fig = figure('Color','w','Visible','off');
   tlay = tiledlayout(3,3,'TileSpacing','compact','Padding','compact');
   vars = {T,M,S,E,P,S.*P,AST,ECL,DepotField_all};
   names = {'TIMP','MMP','Scaffold','ECM','NPC','Neurons','AST','ECL','Depots'};
   clims = {[0 1.5],[0 2],[0 1],[0 1.5],[0 P_max],[0 1],[0 1],[0 1],[0 1]};
   for i=1:9
       nexttile;
       h = slice(x,y,z,gather(vars{i}), slice_pos, slice_pos, slice_pos);
       if ~isempty(h), set(h,'EdgeColor','none','FaceAlpha',0.25); end
       colormap(turbo); caxis(clims{i}); colorbar('Location','eastoutside'); title(names{i});
       axis([0 L 0 L 0 L]); axis vis3d; view([30 30]);
   end
   sgtitle(sprintf('Day %d – CA1 Model', day),'FontWeight','bold');
   drawnow;
   saveas(fig, fullfile(imageFolder, sprintf('Day_%03d.png', day)));
   close(fig);
  
   % --- Progress
   elapsed = toc;
   fprintf('Day %d/%d saved. N_global=%.3f, Elapsed: %.1f s\n', day, nDays, N_global, elapsed);
end
%% ----------------- CREATE VIDEO FROM IMAGES -----------------------------
imageFiles = dir(fullfile(imageFolder,'Day_*.png'));
videoFile = fullfile(imageFolder,'FullScaffold_84days_WntBMP_GPU_Final.mp4');
v = VideoWriter(videoFile,'MPEG-4'); v.FrameRate = 6; open(v);
for i = 1:length(imageFiles)
   img = imread(fullfile(imageFolder,imageFiles(i).name));
   writeVideo(v,img);
end
close(v);
fprintf('Final video saved: %s\n', videoFile);
