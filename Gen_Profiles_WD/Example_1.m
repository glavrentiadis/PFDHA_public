%% Example 1
% Gerenrate rupture profiles for 40km rupture with 1.5m average slip
close all
clear 

%profile parameters
srl = 40;
dx = 0.2; 
thres2accept = [0.02,0.005];
nprof = 10;

%amplitude model parameters
Bk0 = 1.5*40;
%eq specific parameters
KC = 0.006907;
Np = 1;
%global relationships
% KC = 10^(-2.031 -1.009*(log10(srl) - 1.6));
% Np = 1.236;
st_fa = 0.265;

%phase derivative parameters
mu = -125.6637;
s = 28.9977;
%global relationships
% mu = -10^(2.097 + (log10(srl) - 1.6));
% s = 10^(1.493 + 0.996*(log10(srl) - 1.6))

%disp. array
dist_array = (0:dx:srl)';
%generate profiles
[dist_gen,disp_gen,iter2pass] = CreateDispProfLogistPhaseDeriv(dist_array,Bk0,KC,Np,st_fa,mu,s,nprof,thres2accept);

%plot profiles
figid = figure;
hl = plot(dist_gen,disp_gen);
xlabel('Along strike dist. (m)')
ylabel('Displacement (m)')
% ylimits = ylim(); 
% ylim([0,ylimits(2)]);

