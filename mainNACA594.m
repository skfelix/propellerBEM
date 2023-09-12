%% MAIN PROP - ALL BETA
clc; clear; fclose all; 

% Add folders to path  
folder = fileparts(which(mfilename)); 
addpath(genpath(folder));
%%%% Theodorsen - NACA Report 594 (1946) - imperial units
load('NACA594_data.mat')
propData = NACA594.propD;
P_max = propData.perf.beta15.Wdot;
n = 1000/60;
D = propData.geom.D;
B = propData.geom.B;
h = 0;
unit_system = propData.geom.unit_system;

inputs.P_max = P_max; inputs.n = n; inputs.D = D; inputs.B = B; inputs.h = h; inputs.unit_system = unit_system;

temp_offset = 0; % 10 for ISA+10, or -10 for ISA-10 for example
rho = ISA(h,0,unit_system);
polars  = 1;
x        = propData.geom.x;
betadeg    = propData.geom.betadeg;

Cp_max = P_max/(rho*n^3*D^5);

run = 0;
% Performance
if run 
    % BETA75 = 15 deg
    V = 60:5:160;
    beta_ref_deg = propData.perf.beta15.beta_ref_deg+2;
    x_beta_ref = propData.perf.beta15.x_beta_ref;
    an1   = criglerPerf(inputs,propData.geom,polars,V,beta_ref_deg,x_beta_ref);
    polars = 0;
    an2   = adkinsPerf(inputs,propData.geom,polars,V,beta_ref_deg,x_beta_ref);
    % BETA75 = 30 deg
    V = 110:10:270;
    beta_ref_deg = propData.perf.beta30.beta_ref_deg+2;
    x_beta_ref = propData.perf.beta30.x_beta_ref;
    polars = 1;
    an3   = criglerPerf(inputs,propData.geom,polars,V,beta_ref_deg,x_beta_ref);
    polars = 0;
    an4   = adkinsPerf(inputs,propData.geom,polars,V,beta_ref_deg,x_beta_ref);
    % BETA75 = 45 deg
    V = 270:10:440;
    beta_ref_deg = propData.perf.beta45.beta_ref_deg+2;
    x_beta_ref = propData.perf.beta45.x_beta_ref;
    polars = 1;
    an5  = criglerPerf(inputs,propData.geom,polars,V,beta_ref_deg,x_beta_ref);
    polars = 0;
    an6  = adkinsPerf(inputs,propData.geom,polars,V,beta_ref_deg,x_beta_ref);
    save('results/NACA594_results.mat','an1','an2','an3','an4','an5','an6');
else
    load('results/NACA594_results.mat','an1','an2','an3','an4','an5','an6');
end

%% Plots

nfig = 1;
lgd= {
    '$\beta_{75}=15^\circ$ - Crigler'
    '$\beta_{75}=30^\circ$ - Crigler'
    '$\beta_{75}=45^\circ$ - Crigler'
    '$\beta_{75}=15^\circ$ - Adkins'
    '$\beta_{75}=30^\circ$ - Adkins'
    '$\beta_{75}=45^\circ$ - Adkins'
    'Experimental'
    };
title = strcat('NACA TN 594');
set(groot,'defaultAxesTickLabelInterpreter','latex'); 


% Perf
% nfig = nfig+1;
plotVelocity = 0;
figure(nfig);
perf1 = {an1,an3,an5};
perf2 = {an2,an4,an6};
group = 1;
plotPerf(perf1,title,lgd,group,plotVelocity)
group = 2;
plotPerf(perf2,title,lgd,group,plotVelocity)
hold on
subplot(211)
plot(propData.perf.beta15.J,propData.perf.beta15.Ct,'ok','MarkerFaceColor','k')
plot(propData.perf.beta30.J,propData.perf.beta30.Ct,'ok','MarkerFaceColor','k')
plot(propData.perf.beta45.J,propData.perf.beta45.Ct,'ok','MarkerFaceColor','k')
subplot(212)
plot(propData.perf.beta15.J,propData.perf.beta15.Cp,'ok','MarkerFaceColor','k')
plot(propData.perf.beta30.J,propData.perf.beta30.Cp,'ok','MarkerFaceColor','k')
plot(propData.perf.beta45.J,propData.perf.beta45.Cp,'ok','MarkerFaceColor','k')