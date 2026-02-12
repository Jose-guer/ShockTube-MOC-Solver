% This code solves the shock tube problem using the method of
% characteristics. The test time at the end wall is computed.

clc;
clear;
close all;

%%% Constants
torr = 101325/760;
atm = 101325;
Psi = 101235/14.7;
kPa = 10^3;
in2m = 0.0254;
Ru = 8.314; %kJ/kg-K

%%%Initial Conditions
P4 = 220*Psi;
P1 = 80*torr; %P1 needs to be in pascals

T1 = 293; %Initial temperature - set for all gases in tube
T4 = 293;

%name of gas mixture files
% gas4  = 'helium'; % 'HeAr' & 'N2He' & 'helium'
gas4 = 'helium';
gas1  = 'air'; %'nitrogen' % 'air' % 'argon' & 'venus' & 'arox' & 'ArH2O2' & 'HeH2N2'
gases = {gas1 gas4};
plot_xt = 1;
eta = 1; %Diaphragm efficienciy, applied to Mach number

L4 = 2.215;
L1 = 11.25;

L  = [L1, L4];
P  = [P1, P4];
T  = [T1, T4];

%%% Volume of regions
D_tube = 5.5 * in2m;
Ac = pi/4*D_tube^2;

VH = Ac*L4; %Vol driver [m^3]
VL = Ac*L1; % Vol driven [m^3]
Vtot = VL + VH; % [m^3]

[M5, T5, P5, U5, MW1, gas1, T2, P2, rho2, M2, U2, Us1, M1, testtime] = ShockTube_MOC_Solve(P,T,L,gases,eta,plot_xt);

%% Display Information
% clc;

disp('*** Initial States ***')
fprintf('P4 = %.3f [atm] \n',P(2)/101325)
fprintf('P1 = %.3f [torr] \n',P(1)*760/101325)
disp(['gas4: ',gas4])
disp(['gas1: ',gas1])
fprintf('\n')

disp('*** Incident Shock ***')
fprintf('M1 = %.3f \n',M1)
fprintf('S = %0.3f [m/s] \n',Us1)
fprintf('P2/P1 = %0.3f \n',P2/P1)
fprintf('T2/T1 = %0.3f \n \n',T2/T1)

disp('*** Post-shock state ***')
fprintf('M2 = %.3f \n',M2)
fprintf('U2 = %.3f \n',U2)
fprintf('T2 = %.3f [K] \n',T2)
fprintf('P2 = %.3f [atm] \n',P2/101325)
fprintf('R2 = %.4f [kg/m^3] \n \n',rho2)

disp('*** Post-reflected-shock state ***')
fprintf('M5 = %.3f \n',M5)
fprintf('T5 = %.3f  [K]\n',T5)
fprintf('P5 = %.3f [atm] \n',P5/101325)
fprintf('test time = %0.4f (ms) \n\n',testtime * 1000)


disp('*** End State ***')
P_tube = (P4*VH + P1*VL)/Vtot;
fprintf('Tube pressure = %0.4f [atm] \n',P_tube/101325)


fsave = '/Users/joseguerrero/Desktop/';
fname = 'XTdiagram.png';
exportgraphics(gcf,[fsave,fname],'resolution',300)
