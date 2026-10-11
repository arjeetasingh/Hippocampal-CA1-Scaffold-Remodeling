% =============================================================================
% 01_dashboard_dynamics.m
% Project : Hippocampal CA1 Scaffold Remodeling (BSE 639 - Tissue Engineering, Group 3, IIT Indore)
% Lumped (non-spatial) 12-week model: NPCs, neurons, scaffold, MMP/TIMP, Wnt/BMP, BDNF, laminin.
% Requires: MATLAB R2019b or newer (tiledlayout)
% Output  : One 3x2 dashboard figure; optional PNG export at the bottom of the script
% =============================================================================

%% HIPPOCAMPAL CA1 SCAFFOLD – DASHBOARD DYNAMICS (FINAL)
% Project: BSE 639 (Tissue Engineering), Group 3, IIT Indore
% Description:
% 12-week simulation showing NPC proliferation, neuronal maturation,
% scaffold degradation, enzyme regulation, and morphogen interactions.
clc; clear; close all;
weeks = 0:1:12;  % time in weeks
%% --- Biological Parameters ---
% NPC and neuron growth
rP = 0.45;          % NPC proliferation rate per week
deltaN = 0.18;      % Maturation rate (NPC → Neuron)
carrying_capacity = 1;
% Scaffold degradation
S0 = 1;             % Initial scaffold integrity (normalized)
kd = 0.08;          % Degradation constant
Ds = 0.01;          % Diffusion-degradation factor
% MMP/TIMP system
M0 = 0.2; T0 = 0.5;
alphaP = 0.25; alphaA = 0.15;
betaM = 0.12; gammaT = 0.5;
R0 = 0.3; lambdaT = 0.05;
% Morphogens
Wnt0 = 0.9;  % Wnt activator (CHIR99021)
BMP0 = 0.4;  % BMP (differentiation signal)
BDNF0 = 0.3; % Synaptogenesis / neuronal support
Laminin0 = 0.8; % Laminin-511 for progenitor adhesion
%% --- Preallocate Arrays ---
P = zeros(size(weeks));
N = zeros(size(weeks));
S = zeros(size(weeks));
M = zeros(size(weeks));
T = zeros(size(weeks));
Wnt = zeros(size(weeks));
BMP = zeros(size(weeks));
BDNF = zeros(size(weeks));
Laminin = zeros(size(weeks));
%% --- Initial Conditions ---
P(1) = 0.05;
N(1) = 0;
S(1) = S0;
M(1) = M0;
T(1) = T0;
Wnt(1) = Wnt0;
BMP(1) = BMP0;
BDNF(1) = BDNF0;
Laminin(1) = Laminin0;
%% --- Week-wise Simulation Loop ---
for t = 2:length(weeks)
   % NPC proliferation (logistic + Wnt/BMP modulation)
   dP = rP * P(t-1) * (1 - P(t-1)/carrying_capacity) * (Wnt(t-1) - BMP(t-1));
  
   % Neuronal maturation (BDNF-driven)
   dN = deltaN * P(t-1) * (BDNF(t-1) + 0.2);
  
   % Scaffold degradation (MMP-mediated)
   dS = -kd * S(t-1) * (M(t-1) / (T(t-1) + 0.1));
  
   % MMP and TIMP dynamics
   dM = alphaP * P(t-1) + alphaA * 0.02 - betaM * M(t-1) - gammaT * T(t-1) * M(t-1);
   dT = R0 * exp(-lambdaT * weeks(t-1)) - 0.1 * M(t-1) * T(t-1);
   % Morphogen evolution
   Wnt(t) = Wnt(t-1) - 0.02 * (t-1)/12; % slight decline
   BMP(t) = BMP(t-1) + 0.01 * (t-1)/12; % gradual rise
   BDNF(t) = BDNF(t-1) + 0.015 * (t-1)/12; % increases with time
   Laminin(t) = Laminin(t-1) - 0.005 * (t-1)/12; % slight decay
   % Update biological states
   P(t) = max(P(t-1) + dP, 0);
   N(t) = N(t-1) + dN;
   S(t) = max(S(t-1) + dS, 0);
   M(t) = max(M(t-1) + dM, 0);
   T(t) = max(T(t-1) + dT, 0);
end
%% --- DASHBOARD PLOTS ---
figure('Color','w','Position',[100 100 1200 700]);
tiledlayout(3,2,'Padding','compact','TileSpacing','compact');
% 1. NPCs vs Neurons
nexttile;
plot(weeks,P,'-g','LineWidth',2); hold on;
plot(weeks,N,'-b','LineWidth',2);
xlabel('Weeks'); ylabel('Normalized Cell Count');
title('NPC Proliferation & Neuronal Maturation');
legend('NPCs','Neurons','Location','best'); grid on;
% 2. Scaffold Degradation
nexttile;
plot(weeks,S,'-r','LineWidth',2);
xlabel('Weeks'); ylabel('Remaining Scaffold Fraction');
title('Scaffold Degradation (MMP-Controlled)');
grid on;
% 3. MMP vs TIMP
nexttile;
plot(weeks,M,'--r','LineWidth',2); hold on;
plot(weeks,T,'--b','LineWidth',2);
xlabel('Weeks'); ylabel('Enzyme/Regulator Level');
title('MMP–TIMP Enzymatic Dynamics');
legend('MMP','TIMP','Location','best'); grid on;
% 4. Wnt vs BMP
nexttile;
plot(weeks,Wnt,'-g','LineWidth',2); hold on;
plot(weeks,BMP,'-m','LineWidth',2);
xlabel('Weeks'); ylabel('Signal Level');
title('Wnt (CHIR99021) vs BMP Pathway');
legend('Wnt','BMP','Location','best'); grid on;
% 5. BDNF & Laminin
nexttile;
plot(weeks,BDNF,'-c','LineWidth',2); hold on;
plot(weeks,Laminin,'-k','LineWidth',2);
xlabel('Weeks'); ylabel('Relative Level');
title('BDNF & Laminin-511 (Neuronal Support)');
legend('BDNF','Laminin-511','Location','best'); grid on;
% 6. Combined System Overlay
nexttile;

plot(weeks,P,'-g','LineWidth',1.5); hold on;
plot(weeks,N,'-b','LineWidth',1.5);
plot(weeks,S,'-r','LineWidth',1.5);
plot(weeks,Wnt,'--g','LineWidth',1);
plot(weeks,BMP,'--m','LineWidth',1);
plot(weeks,BDNF,'--c','LineWidth',1);
xlabel('Weeks'); ylabel('Normalized Units');
title('System-Level Interaction Summary');
legend({'NPC','Neuron','Scaffold','Wnt','BMP','BDNF'},'Location','bestoutside');
grid on;
sgtitle('Hippocampal CA1 Scaffold Remodeling Dashboard (12 Weeks)','FontWeight','bold');
%% --- Export (optional) ---
% exportgraphics(gcf,'CA1_Scaffold_Dashboard.png','Resolution',300);
