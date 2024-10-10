%Function: SampleMassTrend_RP1FirstStage
%Purpose: Plots mass trends of a LOX/RP1 first stage. Varies second stage
%propellant type as well as first stage delta-V fraction.

%Inputs:

%Ouputs: [min_LCH4,min_LH2,min_RP1,min_Solid,min_Storable]
    %min_Y_Second: unitless - Range 0 to 1 inclusive. Returns the X value
    %that minimizes total gross launch vehicle mass for a given RP1 first
    %stage and Y second stage.

%Authors: James Doss, using prop_matrix from Anthony.

import vary_dv_isp.*

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
N = 100; %Total number of steps in X.
X_vals = linspace(0,1,N); %First Stage DV Ratio to be plotted. Not all values of X will generate a solution capable of reaching orbit in which case the results array will hold NaN
RP1_First_CH4_Second_results = NaN(length(X_vals),6);
RP1_First_H2_Second_results = NaN(length(X_vals),6);
RP1_First_RP1_Second_results = NaN(length(X_vals),6);
RP1_First_Solid_Second_results = NaN(length(X_vals),6);
RP1_First_Storable_Second_results = NaN(length(X_vals),6);

for i = 1:length(X_vals) %Generates the results.
    if(~isempty(vary_dv_isp(X_vals(i),prop_matrix(2,3),prop_matrix(2,1))))
        [RP1_First_CH4_Second_results(i,1),RP1_First_CH4_Second_results(i,2),RP1_First_CH4_Second_results(i,3),RP1_First_CH4_Second_results(i,4),RP1_First_CH4_Second_results(i,5),RP1_First_CH4_Second_results(i,6)] = vary_dv_isp(X_vals(i),prop_matrix(2,3),prop_matrix(2,1));
    end
    if(~isempty(vary_dv_isp(X_vals(i),prop_matrix(2,3),prop_matrix(2,2))))
        [RP1_First_H2_Second_results(i,1),RP1_First_H2_Second_results(i,2),RP1_First_H2_Second_results(i,3),RP1_First_H2_Second_results(i,4),RP1_First_H2_Second_results(i,5),RP1_First_H2_Second_results(i,6)] = vary_dv_isp(X_vals(i),prop_matrix(2,3),prop_matrix(2,2));
    end
    if(~isempty(vary_dv_isp(X_vals(i),prop_matrix(2,3),prop_matrix(2,3))))
        [RP1_First_RP1_Second_results(i,1),RP1_First_RP1_Second_results(i,2),RP1_First_RP1_Second_results(i,3),RP1_First_RP1_Second_results(i,4),RP1_First_RP1_Second_results(i,5),RP1_First_RP1_Second_results(i,6)] = vary_dv_isp(X_vals(i),prop_matrix(2,3),prop_matrix(2,3));
    end
    if(~isempty(vary_dv_isp(X_vals(i),prop_matrix(2,3),prop_matrix(2,4))))
        [RP1_First_Solid_Second_results(i,1),RP1_First_Solid_Second_results(i,2),RP1_First_Solid_Second_results(i,3),RP1_First_Solid_Second_results(i,4),RP1_First_Solid_Second_results(i,5),RP1_First_Solid_Second_results(i,6)] = vary_dv_isp(X_vals(i),prop_matrix(2,3),prop_matrix(2,4));
    end
    if(~isempty(vary_dv_isp(X_vals(i),prop_matrix(2,3),prop_matrix(2,5))))
        [RP1_First_Storable_Second_results(i,1),RP1_First_Storable_Second_results(i,2),RP1_First_Storable_Second_results(i,3),RP1_First_Storable_Second_results(i,4),RP1_First_Storable_Second_results(i,5),RP1_First_Storable_Second_results(i,6)] = vary_dv_isp(X_vals(i),prop_matrix(2,3),prop_matrix(2,5));
    end
end

[min_LCH4,min_LCH4_Index] = min(RP1_First_CH4_Second_results(:,5));
[min_LH2,min_LH2_Index] = min(RP1_First_H2_Second_results(:,5));
[min_RP1,min_RP1_Index] = min(RP1_First_RP1_Second_results(:,5));
[min_Solid,min_Solid_Index] = min(RP1_First_Solid_Second_results(:,5));
[min_Storable,min_Storable_Index] = min(RP1_First_Storable_Second_results(:,5));

disp("Minimum LCH4 2nd Stage Gross Mass (kg): " + min_LCH4 + ". Occuring at X = " + X_vals(min_LCH4_Index));
disp("Minimum LH2 2nd Stage Gross Mass (kg): " + min_LH2 + ". Occuring at X = " + X_vals(min_LH2_Index));
disp("Minimum RP1 2nd Stage Gross Mass (kg): " + min_RP1 + ". Occuring at X = " + X_vals(min_RP1_Index));
disp("Minimum Solid 2nd Stage Gross Mass (kg): " + min_Solid + ". Occuring at X = " + X_vals(min_Solid_Index));
disp("Minimum Storable 2nd Stage Gross Mass (kg): " + min_Storable + ". Occuring at X = " + X_vals(min_Storable_Index));

figure
hold on 
grid on
title("Mass Trends: RP1 First Stage")
xlabel("First Stage DV Fraction")
ylabel("Gross Vehicle Mass (kg)")
plot(RP1_First_CH4_Second_results(:,6),RP1_First_CH4_Second_results(:,5))
plot(RP1_First_H2_Second_results(:,6),RP1_First_H2_Second_results(:,5))
plot(RP1_First_RP1_Second_results(:,6),RP1_First_RP1_Second_results(:,5))
plot(RP1_First_Solid_Second_results(:,6),RP1_First_Solid_Second_results(:,5))
plot(RP1_First_Storable_Second_results(:,6),RP1_First_Storable_Second_results(:,5))
legend("LCH4 Second Stage","LH2 Second Stage","RP1 Second Stage","Solid Second Stage","Storable Second Stage")






