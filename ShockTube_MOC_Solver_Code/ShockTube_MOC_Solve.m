function [M5, T5, P5, Us5, MW1, gas1, T2, P2, rho2, M2, U2, Us1, Ms1,testtime] = ShockTube_MOC_Solve(P,T,L,gases,eta,plot_xt)
% This script solves the standard shock tube problem using the method of
% characteristics (MOC)
%
%     Driver              Driven
%|---------------|--------------------------------|
%|  P4,T4,L4     |            P1,T1,L1            |
%|---------------|--------------------------------|
%
% Regions:
% 1 : Driven section 
% 2 : Post shock state
% 3 : Post expansion wave state
% 4 : Driven section 
% 5 : Test gas (post refelected shock state)
%
% Inputs -   P : Initial pressures [P1,P4], Pa
%            T : Initial temperatures [T1, T4], K
%            L : Shock tube geometry [L1, L4], m
%
% Outputs -  M5 : Reflected shock Mach number
%            T5 : Temperature in region 5, K
%            P5 : Pressure in region 5, Pa

% The code is organized into 3 sections:

% *** Section 1 ***

% This section determines the incident and reflected shocks speed and 
% Mach number given the initial driver and driven pressures. 

% *** Section 2 ***

% This section implements the method of characteristics to calculate
% the trajectory of the head of the reflected expansion fan created at
% the driver/driven interface. 

% *** Section 3 ***

% An XT-diagram is plot using the EW trajectory and computed shock speeds.

% References

% [1] John D. Anderson. Modern Compressible Flow: With a Historical
% Perspective

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Universal gas constant
Ru = 8.314;

P1  = P(1);
P4  = P(2);

T1  = T(1);
T4  = T(2);

L1  = L(1);
L4  = L(2);

gas1  = char(gases(1));
gas4  = char(gases(2));

%function for calculating speed of sound (m/s)
a = @(g,R,T) sqrt(g*R*T);  %g = gamma, R = gas constant, J/k-kg, T = temp, K

%gammm - Cp/Cv
[g1, MW1]  = getgamma(T1,gas1);
[g4, MW4]  = getgamma(T4,gas4);

%speed of sound
a1  = a(g1,Ru/MW1,T1);     %m/s
a4  = a(g4,Ru/MW4,T4);     %m/s

%density, ideal gas
rho1 = P1/(Ru/MW1*T1);     %kg/m^3
rho4 = P4/(Ru/MW4*T4);     %kg/m^3

%% Section 1

% P4/P1 relation written in terms of Ms1. 'x' is Ms1 in this equation.
fun1 = @(x) (1 + 2*g1/(g1+1)*(x^2-1))/(1-(g4-1)/(g1+1)*a1/a4*(x^2-1)/x)^(2*g4/(g4-1)) - P4/P1;
options=optimset('Display','off');
Ms1 = fsolve(fun1,2,options);  
Ms1 = Ms1*eta;  %Incident Shock Mach Number * efficiency
Us1 = Ms1*a1;   %Incident Shock Velocity (lab frame)

%P2/P1 relation as a function of Ms1. 'x' is Ms1 in this equation.
P2P1 = @(g,x) 1 + 2*g/(g+1)*(x^2-1);
P2 = P2P1(g1,Ms1)*P1;  %calculate P2

%T2/T1 relation as a function of Ms1. 'x' is Ms1 in this equation.
T2T1 = @(g,x) (1 + 2*g/(g+1)*(x^2-1)) * (2 + (g-1)*x^2) / (g+1) / x^2;
T2 = T2T1(g1,Ms1)*T1;   %calculate T2

%rho2/rho1 relation as a function of Ms1. 'x' is Ms1 in this equation.
rho2rho1 = @(g,x) (g+1)*x^2 / (2 + (g-1)*x^2);
rho2 = rho2rho1(g1,Ms1)*rho1;

%Calculate other properties for region2 (shocked gas)
g2 = getgamma(T2,gas1);
a2 = a(g2,Ru/MW1,T2);
U2 = Us1 * (1 - rho1/rho2); %shock induced gas velocity (lab frame)
M2 = U2/a2;

% Calculate the properties for region 5 after the shock has been reflected
% from the end of the shocktube
fun2 = @(MR) MR/(MR^2-1)- Ms1/(Ms1^2-1)*(1 + 2*(g2-1)/(g2+1)^2*(Ms1^2-1)*(g2+1/Ms1^2))^0.5;
M5 = fsolve(fun2,Ms1,options); %reflected shock Mach number
T5 = T2T1(g2,M5)*T2;
P5 = P2P1(g2,M5)*P2;
g5 = getgamma(T5,gas1);
a5 = a(g5,Ru/MW1,T5);
U5_rel = M5*a2; 
Us5 = U5_rel - U2; %reflected shock velocity (lab frame)

%% Section 2

%Calculate travel times for shock front
Us1 = Ms1*a1;
ts1 = L1/Us1;   %time it takes incident shock to travel driven section length

%Caluclates travel time for the contact surface (CS).
tc1 = L1/U2;    %time it takes CS to travel driven section length

%Calculates travel time for expansion fan head reach the driver section end
%wall.
tx1 = L4/a4; %time it takes expansion fan head to travel driver section length

%The following analysis is performed assuming that the expansion fan
%originates at (x,t) = (0,0) and is eventually translated to (x,t) = (L4,0)

%Find a3 and u3 for calculating J+ along reflected expansion fan head
P3  = P2;
T3  = (P3/P4)^((g4-1)/g4)*T4; %isentropic expansion from (4) -> (3)
g3 = getgamma(T3,gas4);
a3  = a(g3,Ru/MW4,T3);
U3 = U2;

%Riemann invariants
a  = @(Jp,Jn,g) (g-1)/4*(Jp-Jn);
u  = @(Jp,Jn) 1/2*(Jp+Jn);
Jplus  = @(u,a,g) u+2*a/(g-1);
Jminus = @(u,a,g) u-2*a/(g-1);

%Number of points used for MOC
n = 100;

% In driver/driven region, velocity ranges linearly from 0 (end wall) to U3
% (contact surface).  The local speed of sound, a_fan, is determined from J+
% since it is contstant throughout the simple region

u_fan = linspace(0,U3,n);
a_fan = a4-(g4-1)/2.*u_fan; %J+ = constant, u4 = 0

%Step through velocities to determine J+ and J-, and the local u and a
%along the reflected expansion fan head
for i = 1:length(u_fan)
    Jp(i) = Jplus(U3,a3,g3);  %constant in this region
    Jn(i) = Jminus(u_fan(i),a_fan(i),g4);

    ahead(i) = a(Jp(i),Jn(i),g4);
    uhead(i) = u(Jp(i),Jn(i));
end

%calculate slopes for the c+ and c- characteristics.
c_plus  = uhead + ahead;  % dx/dt = u + a for a c+ chracteristic
c_minus = u_fan - a_fan;  % dx/dt = u - a for a c- chracteristic

% invert the slopes so they work in x-t space, and average the c+ slopes
% between the two points being used in the evaluation
c_plus = ( 1./c_plus(1:end-1) + 1./c_plus(2:end) ) / 2;
c_minus = 1./(c_minus(2:end));

% Find intersection of c- and c+ characteristics, using the average slope
% between the current and previous c+ points being evaluated.

% c+ char. : y2 - y1 = c+ (x2 - x1)
% c- char. : y2 = c- * x2

%combine line equations: x2 = ( y1 - (c+ * x1) ) / (c- - c+)

% xhead(1) = (c_plus(1)*L4+tx1)/(c_minus(1)-c_plus(1));
xhead(1) = -L4;
yhead(1) = c_minus(1)*xhead(1);

for i = 2:length(c_plus)
    xhead(i) = (yhead(i-1)-c_plus(i)*xhead(i-1))/(c_minus(i)-c_plus(i));
    yhead(i) = c_minus(i)*xhead(i);
end

%translate coordinate system to that the driver section 
% end wall is at x = 0, the diaphragm is at x = L4, and the driven section
% end wall is at x = L1 + L4

xhead = xhead+L4;

%% Section 3 - Plots

% Determine intersection points and times:
% Two lines (y = m*x + b) intersect at x = (b2 - b1)/(m1 - m2)

% Contact surface - Expansion wave
bhead = yhead(end) - c_plus(end)*xhead(end); %y-int EW head
bcs = -L4/U2;                                %y-int contact surface

x0 = (bhead - bcs)/(1/U2 - c_plus(end));
t0 = (x0-L4)/U2;

% Contact surface - Reflected Shock
brefl = ts1 + (L1+L4)/Us5;                  %y-int reflected shock
x1 = (brefl - bcs)/(1/U2 + 1/Us5);
t1 = (x1-L4)/U2;

% Head of expansion wave - Reflected Shock
x2 = (brefl - bhead)/(c_plus(end) + 1/Us5);
t2 = (x2-xhead(end))*c_plus(end)+yhead(end);

% Determine if the contact surface or EW reach the shock first
tvec = [t1, t2]; %[cs, exp-wave]
xvec = [x1,x2];

% Take the minimum time, i.e. the first collision
[t_int,ind] = min(tvec);
x_int = xvec(ind);

% Reflected EW back to end wall
trw = (L1 + L4 - x_int)/a5;
testtime = (L1 + L4 - x_int)/Us5 + (L1 + L4 - x_int)/a5;

%plot x-t diagram upto the point of intersection
if ind == 1 %contact surface
    x_cs = x_int;
    x_ew = (t_int-yhead(end))/c_plus(end)+xhead(end);
else
    x_cs = t_int*U2 + L4;
    x_ew = x_int;
end

%Set limits on plot

violet = [138,43,226]/255;

if plot_xt == 1

    yhead = yhead * 1000;
    ts1 = ts1 * 1000;
    t_int = t_int * 1000;
    trw = trw * 1000;
    tx1 = tx1 * 1000;

    figure
    hold on

    LW = 2;

    %Incident shock
    line([L4 L4+L1],[0 ts1],'linestyle','-','color',[1 0 0],'linewidth',LW)

    %Reflected shock
    line([L4+L1, x_int],[ts1,t_int],'linestyle','--','color','r','linewidth',LW)

    %Contact Surface
    line([L4 x_cs],[0 t_int],'linestyle','--','color',[0 0 0],...
        'linewidth',LW) %contact surface between driver/driven

    %Reflected EW from Contact Surface
    line([x_int,L1+L4],[t_int,t_int + trw],'linestyle','-.','color','b','linewidth',LW)

    %Expansion fans
    line([L4 0],[0 tx1],'linestyle','--','color',[0 0 1],...
        'linewidth',LW) %driver/driven expansion fan

    %Reflected head for driver/driven
    line(xhead,yhead,'linestyle','--','color',violet,...
        'linewidth',LW)

    %Extend of reflected EW head to shock tube end wall
    line([xhead(end) x_ew],[yhead(end) t_int],...
        'linestyle','-','color',violet,'linewidth',LW)

    %c- lines
    for i = 1:10:length(xhead)
        line([L4 xhead(i)],[0 yhead(i)],'linestyle','-',...
            'color',[0 0 0],'linewidth',0.25)
    end

    %last c- line
    line([L4 xhead(end)],[0 yhead(end)],'linestyle','-',...
        'color',[0 0 0],'linewidth',0.25)


    hold off
    xlabel('distance (m)')
    ylabel('time (ms)')

    str_Incident = ['Incident Shock: $M_\mathrm{S} =$ ',num2str(round(Ms1,2)), ', $U_\mathrm{S} =$ ', num2str(round(Us1,2)) ' m/s' ];
    str_Refl_S   = ['Reflected Shock: $M_\mathrm{R} =$ ',num2str(round(M5,2)), ', $U_\mathrm{R} =$ ', num2str(round(Us5,2)) ' m/s' ];
    str_Refl_EW  = ['Reflected Expansion Wave: $a_5 =$ ',num2str(round(a5,2)), ' m/s' ];
    str_CS       = ['Contact Surface: $U_2 =$ ',num2str(round(U2,2)), ' m/s' ];
    str_EW_head  = ['Expansion Wave Head: $a_4 =$ ',num2str(round(a4,2)), ' m/s' ];
    str_ref_EW_head  = ['Refl. Expansion Wave Head: $U_\mathrm{head} =$ ',num2str(round(1./c_plus(end),2)), ' m/s' ];
    legend(str_Incident, str_Refl_S, str_CS, str_Refl_EW, str_EW_head,'', str_ref_EW_head, 'C- lines',...
        'location','northwest','fontsize',20)
    set(gcf,'Position',[289 167 1076 669])

end


end
