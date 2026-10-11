% =============================================================================
% 02_spatial_surface_3d.m
% Project : Hippocampal CA1 Scaffold Remodeling (BSE 639 - Tissue Engineering, Group 3, IIT Indore)
% Spatial 3D PDE model (25^3 grid, 25 Gaussian depots, 12 time units = weeks). Plots mid-plane (z = L/2) surfaces at week 12.
% Requires: MATLAB R2020b or newer (turbo colormap, tiledlayout). Runs on CPU.
% Output  : One figure with 5 surface plots (TIMP, MMP, scaffold, ECM, NPC)
% =============================================================================

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FULL COUPLED CA1 MODEL – 5 SPATIAL PARAMETERS (S, M, T, E, P)
% (Global N(t) plot removed as it is non-spatial)
%
% With 25 Gaussian diffusion depots:
%   - Used explicitly for TIMP source
%   - Used to modulate MMP diffusion (D_m(x))
%
% Equations:
%
% 1) TIMP Release (T):
%    ∂T/∂t = D_T ∇²T + R0 e^{-μt} * DepotField(x) - λ_T T
%
% 2) MMP Dynamics (M):
%    ∂M/∂t = D_m(x) ∇²M + α_P P + α_A - β_M M - γ_T T M
%
% 3) Scaffold Degradation (S):
%    ∂S/∂t = D_s ∇²S - k_d M S
%
% 4) ECM Deposition (E):
%    ∂E/∂t = k_E P - k_mE M E
%
% 5) NPC Proliferation & Differentiation (P):
%    ∂P/∂t = r_P P (1 - P/P_max) - k_diff S P + η W(x)
%
% 6) Neuronal Maturation (N): (Only computed, not plotted)
%    dN/dt = k_diff ∫ S(x,t) P(x,t) dV - δ_N N
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc; clear; close all;
%% ----------------- GRID & TIME -----------------------------------------
Ng = 25;                         % number of grid points per dimension
L  = 1.0;                        % normalized length (mm)
[x,y,z] = meshgrid(linspace(0,L,Ng), ...
                  linspace(0,L,Ng), ...
                  linspace(0,L,Ng));
dx = x(1,2,1) - x(1,1,1);        % spatial step
t_final = 12;                    % weeks
dt      = 0.01;                  % time step (stable for D~0.01)
Nt      = round(t_final/dt);     % number of time steps
weeks      = 1:12;
week_times = weeks;
nWeeks     = numel(weeks);
% Snapshots for each variable at weekly times (kept for completeness, though not plotted here)
T_snap = cell(1,nWeeks);
M_snap = cell(1,nWeeks);
S_snap = cell(1,nWeeks);
E_snap = cell(1,nWeeks);
P_snap = cell(1,nWeeks);
N_hist = zeros(Nt,1);            % global neuron count vs time
t_vec  = dt*(1:Nt);              % time vector
%% ----------------- PARAMETERS (meaning from slides) --------------------
% TIMP (Inhibitor) Release (T)
D_T     = 0.01;    % D_T: effective diffusion of TIMP
R0      = 1.0;     % R0: initial TIMP release rate
mu      = 0.25;    % μ: release decay rate constant
lambdaT = 0.10;    % λ_T: TIMP degradation
% MMP Dynamics (M)
% baseline diffusion will be modulated spatially via 25 depots
D_m_background = 0.003;   % base diffusion away from depots
D_m_depot_max  = 0.015;   % higher diffusion near depots (more porous)
alpha_P = 0.3;            % α_P: MMP secretion by NPCs
alpha_A = 0.25;           % α_A: baseline MMP (AD inflammatory)
beta_M  = 0.20;           % β_M: enzyme decay
gamma_T = 0.40;           % γ_T: TIMP inhibition term
% Scaffold Degradation (S)
D_s  = 0.002;      % D_s: effective diffusion (fragment smoothing)
k_d  = 0.15;       % k_d: enzymatic degradation constant
% ECM Deposition (E)
k_E  = 0.10;       % k_E: ECM produced by NPCs
k_mE = 0.05;       % k_mE: ECM degraded by MMPs
% NPC Proliferation & Differentiation (P)
r_P    = 0.50;     % r_P: proliferation rate
P_max  = 0.8;      % carrying capacity for NPC density
k_diff = 0.05;     % k_diff: differentiation ∝ scaffold loss
eta    = 0.25;     % η: Wnt morphogen strength (neuronal bias)
% Neuronal Maturation (N)
delta_N = 0.01;    % δ_N: neuron death rate (AD conditions)
%% ----------------- 25 GAUSSIAN DEPOTS (CORE PART) ----------------------
numDepots = 25;
sigma     = 0.04 * L;     % spatial spread of each depot
rng(10);                  % reproducible positions
depots = L * rand(numDepots,3);   % [0,1]^3 positions
% DepotField: sum of 25 Gaussians – used for TIMP source
DepotField = zeros(Ng,Ng,Ng);
for i = 1:numDepots
   dx_i = x - depots(i,1);
   dy_i = y - depots(i,2);
   dz_i = z - depots(i,3);
   DepotField = DepotField + exp(-(dx_i.^2 + dy_i.^2 + dz_i.^2)/(2*sigma^2));
end
DepotField = DepotField ./ max(DepotField(:));   % normalize [0,1]
% D_m(x): spatially varying MMP diffusion coefficient
D_m_matrix = D_m_background + (D_m_depot_max - D_m_background) .* DepotField;
% Near depots -> D_m ≈ D_m_depot_max
% Away from depots -> D_m ≈ D_m_background
%% ----------------- Wnt Morphogen Field W(x) ----------------------------
% W(x,t) ~ static spatial field here, higher in CA1 niche centre
W = exp(-((x-0.5).^2 + (y-0.5).^2 + (z-0.5).^2)/0.06);
W = W ./ max(W(:));   % normalized [0,1]
%% ----------------- INITIAL CONDITIONS ----------------------------------
T = zeros(Ng,Ng,Ng);        % TIMP: zero initially, released over time
% MMP: slightly elevated baseline, some small noise
M = 0.8 + 0.05*randn(Ng,Ng,Ng);
M(M<0) = 0;
% Scaffold: intact (=1)
S = ones(Ng,Ng,Ng);
% ECM: none initially
E = zeros(Ng,Ng,Ng);
% NPCs: small central cluster
P = 0.05 + 0.20*exp(-((x-0.5).^2 + (y-0.5).^2 + (z-0.5).^2)/0.08);
P(P>P_max) = P_max;
% Global neuron population N(t)
N_global = 0.0;
%% ----------------- TIME LOOP (EXPLICIT EULER) --------------------------
wk_idx = 1;
for n = 1:Nt
   t = n*dt;
   % ---------- Laplacians (Neumann approx via del2) ----------
   LapT = del2(T,dx,dx,dx);
   LapM = del2(M,dx,dx,dx);
   LapS = del2(S,dx,dx,dx);
   % ---------- 1) TIMP PDE (with depot field) ----------
   % ∂T/∂t = D_T ∇²T + R0 e^{-μt} * DepotField - λ_T T
   Source_T = R0 * exp(-mu*t) .* DepotField;
   dT = D_T*LapT + Source_T - lambdaT*T;
   T  = T + dt*dT;
   T(T<0) = 0;
   % ---------- 2) MMP PDE (with variable D_m(x)) ----------
   % ∂M/∂t = D_m(x) ∇²M + α_P P + α_A - β_M M - γ_T T M
   dM = D_m_matrix .* LapM ...
      + alpha_P*P ...
      + alpha_A ...
      - beta_M*M ...
      - gamma_T*(T.*M);
   M  = M + dt*dM;
   M(M<0) = 0;
   % ---------- 3) Scaffold PDE ----------
   % ∂S/∂t = D_s ∇²S - k_d M S
   dS = D_s*LapS - k_d*(M.*S);
   S  = S + dt*dS;
   S(S<0) = 0; S(S>1) = 1;
   % ---------- 4) ECM PDE ----------
   % ∂E/∂t = k_E P - k_mE M E
   dE = k_E*P - k_mE*(M.*E);
   E  = E + dt*dE;
   E(E<0) = 0;
   % ---------- 5) NPC PDE ----------
   % ∂P/∂t = r_P P (1 - P/P_max) - k_diff S P + η W(x)
   dP = r_P*P.*(1 - P/P_max) - k_diff*(S.*P) + eta*W;
   P  = P + dt*dP;
   P(P<0) = 0;  P(P>P_max) = P_max;
   % ---------- 6) Neuronal Maturation ODE (N is still calculated) ----------
   % dN/dt = k_diff ∫ S(x,t)P(x,t) dV - δ_N N
   SP = S .* P;
   Integral_SP = sum(SP(:)) * dx^3;
   dN = k_diff*Integral_SP - delta_N*N_global;
   N_global = N_global + dt*dN;
   if N_global < 0, N_global = 0; end
   N_hist(n) = N_global;
   % ---------- Save weekly snapshots (retained for future use) ----------
   if wk_idx <= nWeeks && abs(t - week_times(wk_idx)) < dt/2
       T_snap{wk_idx} = T;
       M_snap{wk_idx} = M;
       S_snap{wk_idx} = S;
       E_snap{wk_idx} = E;
       P_snap{wk_idx} = P;
       wk_idx = wk_idx + 1;
   end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PLOTTING BLOCKS – CONSOLIDATED 3D SURFACE PLOTS (Z=L/2 SLICE)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define the slice index for the mid-plane (Z=L/2)
zslice = round(Ng/2);
% Get the coordinates for the surface plot (X-Y plane)
X_plot = squeeze(x(:,:,zslice));
Y_plot = squeeze(y(:,:,zslice));
% Create the figure with 5 subplots arranged in 3 rows, 2 columns
figure('Color','w','Position',[50 50 1000 900]);
tiledlayout(3,2,'TileSpacing','compact','Padding','tight');
colormap(turbo); % Set a visually strong colormap for the surfaces
% A) TIMP (T) - Final concentration (Week 12)
plot_variable_surface(X_plot, Y_plot, T, zslice, ...
   '1. TIMP Concentration (T)', 'T Value');
% B) MMP (M) - Final concentration (Week 12)
plot_variable_surface(X_plot, Y_plot, M, zslice, ...
   '2. MMP Concentration (M)', 'M Value');
% C) Scaffold (S) - Final fraction (Week 12)
plot_variable_surface(X_plot, Y_plot, S, zslice, ...
   '3. Scaffold Fraction (S)', 'S Fraction');
% D) ECM (E) - Final concentration (Week 12)
plot_variable_surface(X_plot, Y_plot, E, zslice, ...
   '4. ECM Density (E)', 'E Density');
% E) NPC (P) - Final density (Week 12)
plot_variable_surface(X_plot, Y_plot, P, zslice, ...
   '5. NPC Density (P)', 'P Density');
% The 6th tile remains empty.
sgtitle('CA1 Model Results: Spatial Distribution of 5 Variables (Z=L/2 Surface at Week 12)','FontWeight','heavy','FontSize',14);

%% ----------------- LOCAL FUNCTIONS (must be at the end of a script file) -----------------
% --- Helper function for surface plotting ---
function plot_variable_surface(X_p, Y_p, Data3D, slice_idx, title_str, z_label_str)
   % Extract the 2D slice data for the surface height
   Z_slice = squeeze(Data3D(:,:,slice_idx));
  
   nexttile;
   surf(X_p, Y_p, Z_slice,'EdgeColor','none');
  
   % Get min/max for color scaling based on the entire 3D field
   caxis([min(Data3D(:)) max(Data3D(:))]);
  
   % Customize plot view and labels
   view(45, 30); % Set a consistent 3D view
   title(title_str,'FontWeight','bold');
   xlabel('x (mm)'); ylabel('y (mm)'); zlabel(z_label_str);
   grid on; box on;
   % Add a color bar
   cb = colorbar;
   ylabel(cb, z_label_str, 'FontWeight','bold');
end
