function [gamma, MW, enthalpy, cp_mix, mu]=getgamma(T,mixture_file)
%%%This function takes in a temperature and a mixture file containing NASA9
%%%polynomial data used to determine cp, s, and h. It returns the specific
%%%heat ratio gamma, molecular weight MW [kg/mol], enthalpy [MJ/kg], cp_mix
%%%[J/mol-K] and viscosity [units]

% mixture file must have [x], [poly], [mw], and [visco] defined such that:
%x = [molefrac1 molefrac2...]
%poly = [polynomial(poly#,Trange,component)]
%mw = [component1 component2 ...] %kg/mol
%polynomial fit is A*T^-2+B*T^-1+C+D*T+E*T^2+F*T^3+G*T^4
%visco = [c To etao]

%%%NASA9 Polynomial Data Documentation
% The NASA thermo data file format was documented in:
%
% Sanford Gordon and Bonnie J. McBride, "Computer Program for Calculation of
% Complex Chemical Equilibrium Compositions and Applications: I. Analysis",
% NASA Reference Publication 1311, October 1994.
%
% Bonnie J. McBride and Sanford Gordon, "Computer Program for Calculation of
% Complex Chemical Equilibrium Compositions and Applications: II. Users Manual
% and Program Description", NASA Reference Publication 1311, June 1996.
%
% The equations below for nondimensional specific heat, enthalpy, and
% entropy, are given in Sanford and Bonnie (1994).  Eqs. 4.6-4.8 are the
% "old" NASA format, and Eqs. 4.9-4.11 are the "new" NASA format as discussed
% in this file.
%
% Eq. 4.6: Cp0/R = a1 + a2*T + a3*T^2 + a4*T^3 + a5*T^4
% Eq. 4.7: H0/RT = a1 + a2/2*T + a3/3*T^2 + a4/4*T^3 + a5/5*T^4 + a6/T
% Eq. 4.8: S0/R = a1*ln(T) + a2*T + a3/2*T^2 + a4/3*T^3 + a5/4*T^4 + a7
%
% Eq. 4.9: Cp0/R = a1*T^-2 + a2*T^-1 + a3 + a4*T + a5*T^2 + a6*T^3 + a7*T^4
% Eq. 4.10: H0/RT = -a1*T^-2 + a2*T^-1*ln(T) + a3 + a4*T/2 + a5*T^2/3 +
%                       a6*T^3/4 + a7*T^4/5 + b1/T
% Eq. 4.11: S0/R = -a1*T^-2/2 - a2*T^-1 + a3*ln(T) + a4*T + a5*T^2/2 +
%                     a6*T^3/6 + a7*T^4/4 + b2


eval(mixture_file);
R = 8.314; %J/mol-K
Tref = 298.15; %K

c    = visco(1);
To   = visco(2);
etao = visco(3);

h = @(coeff,T,R) -coeff(1,1)*R*T^-1 +...
    coeff(1,2)*R*log(T) +...
    coeff(1,3)*R*T +...
    coeff(1,4)*R*(T^2)/2 +...
    coeff(1,5)*R*(T^3)/3 +...
    coeff(1,6)*R*(T^4)/4 +...
    coeff(1,7)*R*(T^5)/5 +...
    coeff(1,8)*R;

for i = 1:length(x) % per species in mixture
    
    poly(:,:,i) = poly(:,:,i);
    
    if T<1000
        cp(i)=R*(poly(1,1,i)*T^-2+poly(1,2,i)*T^-1+poly(1,3,i)+poly(1,4,i)*T+poly(1,5,i)*T^2+poly(1,6,i)*T^3+poly(1,7,i)*T^4);
        enthalpy(i) = h(poly(1,:,i),T,R)-h(poly(1,:,i),Tref,R);
    elseif T<6000
        cp(i)=R*(poly(2,1,i)*T^-2+poly(2,2,i)*T^-1+poly(2,3,i)+poly(2,4,i)*T+poly(2,5,i)*T^2+poly(2,6,i)*T^3+poly(2,7,i)*T^4);
        enthalpy(i) = h(poly(2,:,i),T,R)-h(poly(1,:,i),Tref,R);
    else
        cp(i)=R*(poly(3,1,i)*T^-2+poly(3,2,i)*T^-1+poly(3,3,i)+poly(3,4,i)*T+poly(3,5,i)*T^2+poly(3,6,i)*T^3+poly(3,7,i)*T^4);
        enthalpy(i) = h(poly(3,:,i),T,R)-h(poly(1,:,i),Tref,R);
    end
    
    %{
    if T<1000
        enthalpy = h(poly(1,:,i),T,R)-h(poly(1,:,i),Tref,R);
    elseif T<6000
        h1000 = h(poly(1,:,i),1000,R)-h(poly(1,:,i),Tref,R);
        h6000 = h(poly(2,:,i),T,R)-h(poly(2,:,i),1000,R);
        enthalpy = h1000+h6000;
    else
        h1000 = h(poly(1,:,i),1000,R)-h(poly(1,:,i),Tref,R);
        h6000 = h(poly(2,:,i),6000,R)-h(poly(2,:,i),1000,R);
        hHIGH = h(poly(3,:,i),T,R)-h(poly(3,:,i),6000,R);
        enthalpy(i) = hHIGH+h1000+h6000;
    end
    %}
    
end

MW = sum(mw.*x);
enthalpy = enthalpy./1e6;
enthalpy = sum(enthalpy.*x)./MW; %MJ/kg

cp_mix = sum(cp.*x);  %J/mol-K
gamma=cp_mix/(cp_mix-R);

mu = etao*(To+c)/(T+c)*(T/To)^(3/2);

clear x poly mw visco

end