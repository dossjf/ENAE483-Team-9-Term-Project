%Function: vary_dv_isp.
%Purpose: For a given X value (fraction of total Delta V provided by the
%first stage), provides a first pass mass estimation of each stage.

%Inputs: X, isp1, isp2
    %X: Unitless - Range 0 to 1 (inclusive). Fraction of total mission DV
    %provided by first stage.
    %isp1: seconds - Range 269 to 366 seconds (inclusive). Specific Impulse
    %of first stage.
    %isp2: seconds - Range 269 to 366 seconds (inclusive). Specific Impulse
    %of second stage.

%Ouputs:
    %m1: kilograms - Range unbounded. Total first stage gross mass. 
    %m2: kilograms - Range unbounded. Total second stage gross mass. 
    %mpr1: kilograms - Range unbounded. Total first stage propellant mass. 
    %mpr2: kilograms - Range unbounded. Total first stage propellant mass. 
    %m_total: kilograms - Range unbounded. Total launch vehicle gross mass. 
    %deltaV1_fraction: Unitless - Range 0 to 1 inclusive. X.

%Authors: Christian Waidner, James Doss, Arnav Kalotra, Anthony Huynh,
%Fouad Ayoub

function [m1, m2, mpr1, mpr2, m_total, deltaV1_fraction] = vary_dv_isp(X, isp1, isp2)
    delt1 = 0.05;
    delt2 = 0.08;
    mpl = 63000; % kg
    vtot = 9.8e3; % m/s
    g = 9.81; % m/s^2
    
    syms min1 min2 mpr1 mpr2
    eq1 = delt1 == min1 / (min1 + mpr1 + min2 + mpr2 + mpl);
    eq2 = delt2 == min2 / (min2 + mpl + mpr2);
    eq3 = -isp1*g*log((min1 + min2 + mpl + mpr2)/(min1 + mpr1 + min2 + mpl + mpr2)) == X * vtot;
    eq4 = -isp2*g*log((min2 + mpl)/(min2 + mpl + mpr2)) == (1-X) * vtot;
    
    [min1, min2, mpr1, mpr2] = vpasolve([eq1 eq2 eq3 eq4], [min1 min2 mpr1 mpr2]);
    min1 = double(min1);
    min2 = double(min2);
    mpr1 = double(mpr1);
    mpr2 = double(mpr2);
    m2 = min2 + mpr2;
    m1 = min1 + mpr1;
    m_total = m1 + m2 + mpl;
    deltaV1_fraction = X;
end