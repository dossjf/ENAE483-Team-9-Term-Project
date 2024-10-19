%% Author: Anthony Huynh
%% ENAE483 Principles of Space System Design
%% Team Project 1 - Individual Submission 1

clc; close all; clear; format shortG;

%--------------------------------------------------------------------------
% 1.2a Individual Sample Mass Trends
% Author: Anthony Huynh

           % LOX/LCH4 LOX/LH2 LOX/RP1 Solid Storables
prop_matrix = [3.6      6.03    2.72    NaN     2.67; % Oxidizer:Fuel Ratio
              327      366     311     269     285; % Isp (s)
              2.26     1.86    1.92    4.5     1.75; % Thrust 1st Stage (MN)
              0.745    0.099   0.061   2.94    0.067; % Thrust 2nd Stage (MN)
              2.4      2.4     3.7     6.6     1.5; % Engine Exhaust Diameter 1st Stage (m)
              1.5      2.15    0.92    2.34    1.13; % Engine Exhaust Diameter 2nd Stage (m)
              35.16    20.64   25.8    10.5    15.7; % Chamber Pressure 1st Stage (MPa)
              10.1     4.2     6.77    5       14.7; % Chamber Pressure 2nd Stage (MPa)
              34.34    78      37      16      26.2; % Nozzle area ratio sea-level 1st Stage
              45       84      14.5    56      81.3]; % Nozzle area ratio sea-level 2nd Stage

% Initializes result array for each combination
% Columns: m1, m2, mpr1, mpr2, gross mass, deltaV fraction, stage1_cost, stage2_cost, combined_cost
results_LOXLH2_LOXLH2 = NaN(length(0:0.001:1),9); 
results_LOXLH2_LOXRP1 = NaN(length(0:0.001:1),9); 
results_LOXLH2_Solid = NaN(length(0:0.001:1),9); 
results_LOXLH2_Storables = NaN(length(0:0.001:1),9); 
results_LOXLH2_LOXLCH4 = NaN(length(0:0.001:1),9);
 
i = 1; % index to traverse the rows in the result matrices
for X = 0:0.001:1 % Iterates through values of X between 0 and 1 with step size of 0.001
    [results_LOXLH2_LOXLH2(i,1), results_LOXLH2_LOXLH2(i,2), results_LOXLH2_LOXLH2(i,3), results_LOXLH2_LOXLH2(i,4), results_LOXLH2_LOXLH2(i,5), results_LOXLH2_LOXLH2(i,6)] = vary_dv_isp(X, prop_matrix(2,2), prop_matrix(2,2));
    if (~isnan(results_LOXLH2_LOXLH2(i,1))) % checks if the entry is not NaN (meaning there is actual data)
        [results_LOXLH2_LOXLH2(i,7), results_LOXLH2_LOXLH2(i,8), results_LOXLH2_LOXLH2(i,9)] = cost_calculator(results_LOXLH2_LOXLH2(i,1)-results_LOXLH2_LOXLH2(i,3),results_LOXLH2_LOXLH2(i,2)-results_LOXLH2_LOXLH2(i,4));
    end
    [results_LOXLH2_LOXRP1(i,1), results_LOXLH2_LOXRP1(i,2), results_LOXLH2_LOXRP1(i,3), results_LOXLH2_LOXRP1(i,4), results_LOXLH2_LOXRP1(i,5), results_LOXLH2_LOXRP1(i,6)] = vary_dv_isp(X, prop_matrix(2,2), prop_matrix(2,3));
    if (~isnan(results_LOXLH2_LOXRP1(i,1))) % checks if the entry is not NaN (meaning there is actual data)
        [results_LOXLH2_LOXRP1(i,7), results_LOXLH2_LOXRP1(i,8), results_LOXLH2_LOXRP1(i,9)] = cost_calculator(results_LOXLH2_LOXRP1(i,1)-results_LOXLH2_LOXRP1(i,3),results_LOXLH2_LOXRP1(i,2)-results_LOXLH2_LOXRP1(i,4));
    end
    [results_LOXLH2_Solid(i,1), results_LOXLH2_Solid(i,2), results_LOXLH2_Solid(i,3), results_LOXLH2_Solid(i,4), results_LOXLH2_Solid(i,5), results_LOXLH2_Solid(i,6)] = vary_dv_isp(X, prop_matrix(2,2), prop_matrix(2,4));
    if (~isnan(results_LOXLH2_Solid(i,1))) % checks if the entry is not NaN (meaning there is actual data)
        [results_LOXLH2_Solid(i,7), results_LOXLH2_Solid(i,8), results_LOXLH2_Solid(i,9)] = cost_calculator(results_LOXLH2_Solid(i,1)-results_LOXLH2_Solid(i,3),results_LOXLH2_Solid(i,2)-results_LOXLH2_Solid(i,4));
    end
    [results_LOXLH2_Storables(i,1), results_LOXLH2_Storables(i,2), results_LOXLH2_Storables(i,3), results_LOXLH2_Storables(i,4), results_LOXLH2_Storables(i,5), results_LOXLH2_Storables(i,6)] = vary_dv_isp(X, prop_matrix(2,2), prop_matrix(2,5));
    if (~isnan(results_LOXLH2_Storables(i,1))) % checks if the entry is not NaN (meaning there is actual data)
        [results_LOXLH2_Storables(i,7), results_LOXLH2_Storables(i,8), results_LOXLH2_Storables(i,9)] = cost_calculator(results_LOXLH2_Storables(i,1)-results_LOXLH2_Storables(i,3),results_LOXLH2_Storables(i,2)-results_LOXLH2_Storables(i,4));
    end
    [results_LOXLH2_LOXLCH4(i,1), results_LOXLH2_LOXLCH4(i,2), results_LOXLH2_LOXLCH4(i,3), results_LOXLH2_LOXLCH4(i,4), results_LOXLH2_LOXLCH4(i,5), results_LOXLH2_LOXLCH4(i,6)] = vary_dv_isp(X, prop_matrix(2,2), prop_matrix(2,1));
    if (~isnan(results_LOXLH2_LOXLCH4(i,1))) % checks if the entry is not NaN (meaning there is actual data)
        [results_LOXLH2_LOXLCH4(i,7), results_LOXLH2_LOXLCH4(i,8), results_LOXLH2_LOXLCH4(i,9)] = cost_calculator(results_LOXLH2_LOXLCH4(i,1)-results_LOXLH2_LOXLCH4(i,3),results_LOXLH2_LOXLCH4(i,2)-results_LOXLH2_LOXLCH4(i,4));
    end
    i = i + 1;
end

% Calculates the minimum LV gross mass and corresponding row number in matrix for all five combinations
[minmass_LOXLH2_LOXLCH4, index_minmass1] = min(results_LOXLH2_LOXLCH4(:,5));
[minmass_LOXLH2_LOXLH2, index_minmass2] = min(results_LOXLH2_LOXLH2(:,5));
[minmass_LOXLH2_LOXRP1, index_minmass3] = min(results_LOXLH2_LOXRP1(:,5));
[minmass_LOXLH2_Solid, index_minmass4] = min(results_LOXLH2_Solid(:,5));
[minmass_LOXLH2_Storables, index_minmass5] = min(results_LOXLH2_Storables(:,5));

% Displays the Minimized Mass Solution in Metric Tons
disp("Stage 1 LOX/LH2 and Stage 2 LOX/LCH4 Minimized Mass Solution: " + round(minmass_LOXLH2_LOXLCH4/1000,2) + " tons at X = " + results_LOXLH2_LOXLCH4(index_minmass1,6));
disp("Stage 1 LOX/LH2 and Stage 2 LOX/LH2 Minimized Mass Solution: " + round(minmass_LOXLH2_LOXLH2/1000,2) + " tons at X = " + results_LOXLH2_LOXLCH4(index_minmass2,6));
disp("Stage 1 LOX/LH2 and Stage 2 LOX/RP1 Minimized Mass Solution: " + round(minmass_LOXLH2_LOXRP1/1000,2) + " tons at X = " + results_LOXLH2_LOXRP1(index_minmass3,6));
disp("Stage 1 LOX/LH2 and Stage 2 Solid Minimized Mass Solution: " + round(minmass_LOXLH2_Solid/1000,2) + " tons at X = " + results_LOXLH2_Solid(index_minmass4,6));
disp("Stage 1 LOX/LH2 and Stage 2 Storables Minimized Mass Solution: " + round(minmass_LOXLH2_Storables/1000,2) + " tons at X = " + results_LOXLH2_Storables(index_minmass5,6) + newline);

% Plots a graph of the First Stage Mass, Second Stage Mass, and Gross Mass in metric tons as trends of X 
% Selected Stage 1 LOX/LH2 and Stage 2 LOX/LH2 for Sample Mass Trends
figure()
hold on
% X bound selected to be 100 = 0:0.001:0.099 for better view of data
plot(results_LOXLH2_LOXLH2(100:1001,6),results_LOXLH2_LOXLH2(100:1001,1)/100, 'LineWidth', 2); % Corresponding First Stage Mass in Metric Tons
plot(results_LOXLH2_LOXLH2(100:1001,6),results_LOXLH2_LOXLH2(100:1001,2)/100, 'LineWidth', 2); % Corresponding Second Stage Mass in Metric Tons
plot(results_LOXLH2_LOXLH2(100:1001,6),results_LOXLH2_LOXLH2(100:1001,5)/100, 'LineWidth', 2); % Corresponding Gross Mass in Metric Tons
plot(results_LOXLH2_LOXLH2(index_minmass2,6),results_LOXLH2_LOXLH2(index_minmass2,5)/1000,"k.", "MarkerSize", 30); % Minimized Solution
ax = gca;
ax.YAxis.Exponent = 0;
set(gca,'FontSize',25, 'FontName', 'Courier')
title("Stage 1 LOX/LH2 and Stage 2 LOX/LH2 Mass Optimization")
xlabel("First Stage \DeltaV Fraction (X)")
xlim([0.1 1]);
ylabel("Mass (Metric Tons)")
legend("First Stage Mass", "Second Stage Mass", "Gross Mass", "Minimized Solution");
hold off

%--------------------------------------------------------------------------
% 1.3a) Individual Sample Cost Trends 
% Author: Anthony Huynh

% Calculates the minimum program cost solutions and corresponding row number in matrix for all five combinations
[mincost_LOXLH2_LOXLCH4, index_mincost1] = min(results_LOXLH2_LOXLCH4(:,9));
[mincost_LOXLH2_LOXLH2, index_mincost2] = min(results_LOXLH2_LOXLH2(:,9));
[mincost_LOXLH2_LOXRP1, index_mincost3] = min(results_LOXLH2_LOXRP1(:,9));
[mincost_LOXLH2_Solid, index_mincost4] = min(results_LOXLH2_Solid(:,9));
[mincost_LOXLH2_Storables, index_mincost5] = min(results_LOXLH2_Storables(:,9));

% Displays the Minimized Cost Solution in $M2024
disp("Stage 1 LOX/LH2 and Stage 2 LOX/LCH4 Minimized Cost Solution: " + round(mincost_LOXLH2_LOXLCH4,2) + " ($M2024) at X = " + results_LOXLH2_LOXLCH4(index_mincost1,6));
disp("Stage 1 LOX/LH2 and Stage 2 LOX/LH2 Minimized Cost Solution: " + round(mincost_LOXLH2_LOXLH2,2) + " ($M2024) at X = " + results_LOXLH2_LOXLH2(index_mincost2,6));
disp("Stage 1 LOX/LH2 and Stage 2 LOX/RP1 Minimized Cost Solution: " + round(mincost_LOXLH2_LOXRP1,2) + " ($M2024) at X = " + results_LOXLH2_LOXRP1(index_mincost3,6));
disp("Stage 1 LOX/LH2 and Stage 2 Solid Minimized Cost Solution: " + round(mincost_LOXLH2_Solid,2) + " ($M2024) at X = " + results_LOXLH2_Solid(index_mincost4,6));
disp("Stage 1 LOX/LH2 and Stage 2 Storables Minimized Cost Solution: " + round(mincost_LOXLH2_Storables,2) + " ($M2024) at X = " + results_LOXLH2_Storables(index_mincost5,6));

% Plots a graph of the First Stage Cost, Second Stage Cost, and Combined LV Cost in metric tons as trends of X 
% Selected Stage 1 LOX/LH2 and Stage 2 LOX/LH2 for Sample Cost Trends

figure()
hold on
% X bound selected to be 100 = 0:0.001:0.099 for better view of data
plot(results_LOXLH2_LOXLH2(100:1001,6),results_LOXLH2_LOXLH2(100:1001,7), 'LineWidth', 2); % Corresponding First Stage Cost in $M2024
plot(results_LOXLH2_LOXLH2(100:1001,6),results_LOXLH2_LOXLH2(100:1001,8), 'LineWidth', 2); % Corresponding Second Stage Cost in $M2024
plot(results_LOXLH2_LOXLH2(100:1001,6),results_LOXLH2_LOXLH2(100:1001,9), 'LineWidth', 2); % Corresponding Gross Vehicle Cost in $M2024
plot(results_LOXLH2_LOXLH2(index_mincost2,6),results_LOXLH2_LOXLH2(index_mincost2,9),"k.", "MarkerSize", 30); % Minimized Solution
ax = gca;
ax.YAxis.Exponent = 0;
set(gca,'FontSize',25, 'FontName', 'Courier')
title("Stage 1 LOX/LH2 and Stage 2 LOX/LH2 Cost Optimization")
xlabel("First Stage \DeltaV Fraction (X)")
xlim([0.1 1]);
ylabel("Total Cost ($M2024)")
legend("First Stage Cost", "Second Stage Cost", "Launch Vehicle Cost", "Minimized Solution");
hold off

%---------------------------Helper Functions-------------------------------
% Authors: Fouad Ayoub, James Doss, Anthony Huynh, Arnav Kalotra, Christian Waidner
% Function: vary_dv_isp_mass
% Description: Takes fraction of total Delta V provided by first stage (X)
% to calculate the corresponding mass estimations for a given combination of Isp
% values for stage 1 and 2 propellants
% Inputs: 
    % deltaV1 Fraction (X, unitless)
    % Stage 1 Propellant Isp (s)
    % Stage 2 Propellent Isp (s)
% Outputs: 
    % Stage 1 Mass (kg)
    % Stage 2 Mass (kg)
    % Stage 1 Propellant Mass (kg)
    % Stage 2 Propellant Mass (kg)
    % Gross Mass (kg)
    % deltaV1 Fraction (X, unitless) 
function [m1, m2, mpr1, mpr2, m_total, deltaV1_fraction] = vary_dv_isp(X, isp1, isp2)
% Givens:
delt1 = 0.05; % inert mass fraction Stage 1 (NOT stage inert mass fraction)
delt2 = 0.08; % inert mass fraction Stage 2 (NOT stage inert mass fraction)
mpl = 63000; % kg
vtot = 9.8e3; % m/s
g = 9.81; % m/s^2

%Solution:
syms min1 min2 mpr1 mpr2
eq1 = delt1 == min1 / (min1 + mpr1 + min2 + mpr2 + mpl);
eq2 = delt2 == min2 / (min2 + mpl + mpr2);
eq3 = -isp1*g*log((min1 + min2 + mpl + mpr2)/(min1 + mpr1 + min2 + mpl + mpr2)) == X * vtot;
eq4 = -isp2*g*log((min2 + mpl)/(min2 + mpl + mpr2)) == (1-X) * vtot;
[min1, min2, mpr1, mpr2] = vpasolve([eq1 eq2 eq3 eq4], [min1 min2 mpr1 mpr2]);
if (~isempty(min1)) % If statement to check for non-NaN results from vpasolve (this is due to asymptote that yields no results)
    min1 = double(min1);
    min2 = double(min2);
    mpr1 = double(mpr1);
    mpr2 = double(mpr2);
    m2 = min2 + mpr2; % Stage mass 2 defined as inert + propellant ONLY
    m1 = min1 + mpr1; % Stage mass 1 defined as inert + propellant ONLY
    m_total = m1 + m2 + mpl; % Also referred to as m0
    deltaV1_fraction = X;
else % Else returns values as NaN to ignore for outside calculations
    m1 = nan; m2 = nan; mpr1 = nan; mpr2 = nan; m_total = nan; deltaV1_fraction = nan;
end
end

% Author: Anthony Huynh
% Function: cost_calculator
% Description: Takes the inert masses of stage 1 and 2 to calculate the
% corresponding stage and combined program costs
% Inputs: 
    % Stage 1 Inert Mass (kg)
    % Stage 2 Inert Mass (kg)
% Outputs: 
    % Total Stage 1 Cost ($M2024)
    % Total Stage 2 Cost ($M2024)
    % Program Cost ($M2024)
function [stage1_cost, stage2_cost, combined_cost] = cost_calculator(min1, min2)
% Givens:
P = log(0.85)/log(2);
flights = 10;

% Solution:
stage1_NRcost = 12.73*(min1)^0.55; % NR = non-recurring cost
stage2_NRcost = 12.73*(min2)^0.55;
stage1_fupc = 0.3024*(min1)^0.662; % fupc = first unit production cost
stage2_fupc = 0.3024*(min2)^0.662;
stage1_Rcost = 0;
stage2_Rcost = 0;
for unit_number=1:flights
    stage1_Rcost = stage1_Rcost + stage1_fupc*unit_number^P;
    stage2_Rcost = stage2_Rcost + stage2_fupc*unit_number^P;
end
stage1_cost = stage1_NRcost + stage1_Rcost;
stage2_cost = stage2_NRcost + stage2_Rcost;
combined_cost = (stage1_cost + stage2_cost)*1.025; % 2023-2024 inflation rate assumed to be 2.5% from lecture notes
end