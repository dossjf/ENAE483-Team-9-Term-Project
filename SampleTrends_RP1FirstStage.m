%Function: SampleTrends_RP1FirstStage
%Purpose: Plots mass and cost trends of a LOX/RP1 first stage. Varies second stage
%propellant type as well as first stage delta-V fraction.

%Inputs:

%Ouputs: [min_LCH4,min_LH2,min_RP1,min_Solid,min_Storable]
    %min_Y_Second: unitless - Range 0 to 1 inclusive. Returns the X value
    %that minimizes total gross launch vehicle mass for a given RP1 first
    %stage and Y second stage.

%Authors: James Doss
clc; close all; clear
%--------------------------------------------------------------------------
%1.2a) Sample Mass Trends - James Doss
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
    RP1_First_CH4_Second_results = NaN(length(X_vals),8);
    RP1_First_H2_Second_results = NaN(length(X_vals),8);
    RP1_First_RP1_Second_results = NaN(length(X_vals),8);
    RP1_First_Solid_Second_results = NaN(length(X_vals),8);
    RP1_First_Storable_Second_results = NaN(length(X_vals),8);

    for i = 1:length(X_vals) %Generates the results. 
        if(~isempty(vary_dv_isp(X_vals(i),prop_matrix(2,3),prop_matrix(2,1)))) %If the return of vary_dv_isp is empty (i.e. the design is not capable of reaching orbital velocity), the output is not saved and the NaNs are used instead to prevent errors in graphing.
            [RP1_First_CH4_Second_results(i,1),RP1_First_CH4_Second_results(i,2),RP1_First_CH4_Second_results(i,3),RP1_First_CH4_Second_results(i,4),RP1_First_CH4_Second_results(i,5),RP1_First_CH4_Second_results(i,6),RP1_First_CH4_Second_results(i,7),RP1_First_CH4_Second_results(i,8)] = vary_dv_isp(X_vals(i),prop_matrix(2,3),prop_matrix(2,1));
        end
        if(~isempty(vary_dv_isp(X_vals(i),prop_matrix(2,3),prop_matrix(2,2))))
            [RP1_First_H2_Second_results(i,1),RP1_First_H2_Second_results(i,2),RP1_First_H2_Second_results(i,3),RP1_First_H2_Second_results(i,4),RP1_First_H2_Second_results(i,5),RP1_First_H2_Second_results(i,6),RP1_First_H2_Second_results(i,7),RP1_First_H2_Second_results(i,8)] = vary_dv_isp(X_vals(i),prop_matrix(2,3),prop_matrix(2,2));
        end
        if(~isempty(vary_dv_isp(X_vals(i),prop_matrix(2,3),prop_matrix(2,3))))
            [RP1_First_RP1_Second_results(i,1),RP1_First_RP1_Second_results(i,2),RP1_First_RP1_Second_results(i,3),RP1_First_RP1_Second_results(i,4),RP1_First_RP1_Second_results(i,5),RP1_First_RP1_Second_results(i,6),RP1_First_RP1_Second_results(i,7),RP1_First_RP1_Second_results(i,8)] = vary_dv_isp(X_vals(i),prop_matrix(2,3),prop_matrix(2,3));
        end
        if(~isempty(vary_dv_isp(X_vals(i),prop_matrix(2,3),prop_matrix(2,4))))
            [RP1_First_Solid_Second_results(i,1),RP1_First_Solid_Second_results(i,2),RP1_First_Solid_Second_results(i,3),RP1_First_Solid_Second_results(i,4),RP1_First_Solid_Second_results(i,5),RP1_First_Solid_Second_results(i,6),RP1_First_Solid_Second_results(i,7),RP1_First_Solid_Second_results(i,8)] = vary_dv_isp(X_vals(i),prop_matrix(2,3),prop_matrix(2,4));
        end
        if(~isempty(vary_dv_isp(X_vals(i),prop_matrix(2,3),prop_matrix(2,5))))
            [RP1_First_Storable_Second_results(i,1),RP1_First_Storable_Second_results(i,2),RP1_First_Storable_Second_results(i,3),RP1_First_Storable_Second_results(i,4),RP1_First_Storable_Second_results(i,5),RP1_First_Storable_Second_results(i,6),RP1_First_Storable_Second_results(i,7),RP1_First_Storable_Second_results(i,8)] = vary_dv_isp(X_vals(i),prop_matrix(2,3),prop_matrix(2,5));
        end
    end
    
    [min_LCH4,min_LCH4_Index] = min(RP1_First_CH4_Second_results(:,5)); %Finding the minimums is done simply with the "min" function. LERPing could be used for low values of N to improve estimates but I run this at 100-1000 steps.
    [min_LH2,min_LH2_Index] = min(RP1_First_H2_Second_results(:,5));
    [min_RP1,min_RP1_Index] = min(RP1_First_RP1_Second_results(:,5));
    [min_Solid,min_Solid_Index] = min(RP1_First_Solid_Second_results(:,5));
    [min_Storable,min_Storable_Index] = min(RP1_First_Storable_Second_results(:,5));
    
    disp("Minimum LCH4 2nd Stage Gross Mass (kg): " + min_LCH4 + ". Occuring at X = " + X_vals(min_LCH4_Index));
    disp("Minimum LH2 2nd Stage Gross Mass (kg): " + min_LH2 + ". Occuring at X = " + X_vals(min_LH2_Index));
    disp("Minimum RP1 2nd Stage Gross Mass (kg): " + min_RP1 + ". Occuring at X = " + X_vals(min_RP1_Index));
    disp("Minimum Solid 2nd Stage Gross Mass (kg): " + min_Solid + ". Occuring at X = " + X_vals(min_Solid_Index));
    disp("Minimum Storable 2nd Stage Gross Mass (kg): " + min_Storable + ". Occuring at X = " + X_vals(min_Storable_Index));
    
    figure %Plots the mass trends using a fairly formulatic approach.
    hold on 
    grid on
    title("Mass Trends: RP1/LOX First Stage, LCH4/LOX Second Stage")
    xlabel("First Stage DV Fraction (X)")
    ylabel("Mass (Metric Tons)")
    plot(RP1_First_CH4_Second_results(:,6),RP1_First_CH4_Second_results(:,1)/1000)
    plot(RP1_First_CH4_Second_results(:,6),RP1_First_CH4_Second_results(:,2)/1000)
    plot(RP1_First_CH4_Second_results(:,6),RP1_First_CH4_Second_results(:,5)/1000)
    plot(X_vals(min_LCH4_Index),min_LCH4/1000,'O')
    legend("Gross First Stage Mass","Gross Second Stage Mass","Gross Vehicle Mass")
    
    figure
    hold on 
    grid on
    title("Mass Trends: RP1/LOX First Stage, LH2/LOX Second Stage")
    xlabel("First Stage DV Fraction (X)")
    ylabel("Mass (Metric Tons)")
    plot(RP1_First_H2_Second_results(:,6),RP1_First_H2_Second_results(:,1)/1000)
    plot(RP1_First_H2_Second_results(:,6),RP1_First_H2_Second_results(:,2)/1000)
    plot(RP1_First_H2_Second_results(:,6),RP1_First_H2_Second_results(:,5)/1000)
    plot(X_vals(min_LH2_Index),min_LH2/1000,'O')
    legend("Gross First Stage Mass","Gross Second Stage Mass","Gross Vehicle Mass")
    
    figure
    hold on 
    grid on
    title("Mass Trends: RP1/LOX First Stage, RP1/LOX Second Stage")
    xlabel("First Stage DV Fraction (X)")
    ylabel("Mass (Metric Tons)")
    plot(RP1_First_RP1_Second_results(:,6),RP1_First_RP1_Second_results(:,1)/1000)
    plot(RP1_First_RP1_Second_results(:,6),RP1_First_RP1_Second_results(:,2)/1000)
    plot(RP1_First_RP1_Second_results(:,6),RP1_First_RP1_Second_results(:,5)/1000)
    plot(X_vals(min_RP1_Index),min_RP1/1000,'O')
    legend("Gross First Stage Mass","Gross Second Stage Mass","Gross Vehicle Mass")
    
    figure
    hold on 
    grid on
    title("Mass Trends: RP1/LOX First Stage, Solid Second Stage")
    xlabel("First Stage DV Fraction (X)")
    ylabel("Mass (Metric Tons)")
    plot(RP1_First_Solid_Second_results(:,6),RP1_First_Solid_Second_results(:,1)/1000)
    plot(RP1_First_Solid_Second_results(:,6),RP1_First_Solid_Second_results(:,2)/1000)
    plot(RP1_First_Solid_Second_results(:,6),RP1_First_Solid_Second_results(:,5)/1000)
    plot(X_vals(min_Solid_Index),min_Solid/1000,'O')
    legend("Gross First Stage Mass","Gross Second Stage Mass","Gross Vehicle Mass")
    
    figure
    hold on 
    grid on
    title("Mass Trends: RP1/LOX First Stage, Storable Second Stage")
    xlabel("First Stage DV Fraction (X)")
    ylabel("Mass (Metric Tons)")
    plot(RP1_First_Storable_Second_results(:,6),RP1_First_Storable_Second_results(:,1)/1000)
    plot(RP1_First_Storable_Second_results(:,6),RP1_First_Storable_Second_results(:,2)/1000)
    plot(RP1_First_Storable_Second_results(:,6),RP1_First_Storable_Second_results(:,5)/1000)
    plot(X_vals(min_Storable_Index),min_Storable/1000,'O')
    legend("Gross First Stage Mass","Gross Second Stage Mass","Gross Vehicle Mass")
%--------------------------------------------------------------------------
%1.3a) Sample Cost Trends - James Doss
    RP1CH4Costs = zeros(length(X_vals),3); %Creating arrays to hold costing information returned by apply_SLVCM
    RP1LH2Costs = zeros(length(X_vals),3);
    RP1RP1Costs = zeros(length(X_vals),3);
    RP1SolidCosts = zeros(length(X_vals),3);
    RP1StorableCosts = zeros(length(X_vals),3);
    for i = 1:length(X_vals) %applying SLVCM for each data point we made from mass and storing in the array. SLVCM with return a NaN when passed in NaNs so it automatically works at those edge cases where a design cannot reach orbit.
        [RP1CH4Costs(i,1),RP1CH4Costs(i,2),RP1CH4Costs(i,3)] = apply_SLVCM(RP1_First_CH4_Second_results(i,7),RP1_First_CH4_Second_results(i,8));
        [RP1LH2Costs(i,1),RP1LH2Costs(i,2),RP1LH2Costs(i,3)] = apply_SLVCM(RP1_First_H2_Second_results(i,7),RP1_First_H2_Second_results(i,8));
        [RP1RP1Costs(i,1),RP1RP1Costs(i,2),RP1RP1Costs(i,3)] = apply_SLVCM(RP1_First_RP1_Second_results(i,7),RP1_First_RP1_Second_results(i,8));
        [RP1SolidCosts(i,1),RP1SolidCosts(i,2),RP1SolidCosts(i,3)] = apply_SLVCM(RP1_First_Solid_Second_results(i,7),RP1_First_Solid_Second_results(i,8));
        [RP1StorableCosts(i,1),RP1StorableCosts(i,2),RP1StorableCosts(i,3)] = apply_SLVCM(RP1_First_Storable_Second_results(i,7),RP1_First_Storable_Second_results(i,8));
    end
    
    [RP1CH4Min,RP1CH4Min_Index] = min(RP1CH4Costs(:,3));
    [RP1LH2Min,RP1LH2Min_Index] = min(RP1LH2Costs(:,3));
    [RP1RP1Min,RP1RP1Min_Index] = min(RP1RP1Costs(:,3));
    [RP1SolidMin,RP1SolidMin_Index] = min(RP1SolidCosts(:,3));
    [RP1StorableMin,RP1StorableMin_Index] = min(RP1StorableCosts(:,3));
    disp("Minimum Total Program Cost (LCH4 2nd Stage) ($M2024)" + RP1CH4Min + ". Occuring at X = " + X_vals(RP1CH4Min_Index));
    disp("Minimum Total Program Cost (LH2 2nd Stage) ($M2024)" + RP1LH2Min + ". Occuring at X = " + X_vals(RP1LH2Min_Index));
    disp("Minimum Total Program Cost (RP1 2nd Stage) ($M2024)" + RP1RP1Min + ". Occuring at X = " + X_vals(RP1RP1Min_Index));
    disp("Minimum Total Program Cost (Solid 2nd Stage) ($M2024)" + RP1SolidMin + ". Occuring at X = " + X_vals(RP1SolidMin_Index));
    disp("Minimum Total Program Cost (Storable 2nd Stage) ($M2024)" + RP1StorableMin + ". Occuring at X = " + X_vals(RP1StorableMin_Index));
    
    figure
    hold on 
    grid on
    title("Cost Trends: RP1/LOX First Stage, LCH4/LOX Second Stage")
    xlabel("First Stage DV Fraction (X)")
    ylabel("Total Cost ($M2024)")
    plot(RP1_First_CH4_Second_results(:,6),RP1CH4Costs(:,1))
    plot(RP1_First_CH4_Second_results(:,6),RP1CH4Costs(:,2))
    plot(RP1_First_CH4_Second_results(:,6),RP1CH4Costs(:,3))
    plot(X_vals(RP1CH4Min_Index),RP1CH4Min,'O')
    legend("First Stage Production and NRE","Second Stage Production and NRE","Total Program Cost")
    
    figure
    hold on 
    grid on
    title("Cost Trends: RP1/LOX First Stage, LH2/LOX Second Stage")
    xlabel("First Stage DV Fraction (X)")
    ylabel("Total Cost ($M2024)")
    plot(RP1_First_H2_Second_results(:,6),RP1LH2Costs(:,1))
    plot(RP1_First_H2_Second_results(:,6),RP1LH2Costs(:,2))
    plot(RP1_First_H2_Second_results(:,6),RP1LH2Costs(:,3))
    plot(X_vals(RP1LH2Min_Index),RP1LH2Min,'O')
    legend("First Stage Production and NRE","Second Stage Production and NRE","Total Program Cost")
    
    figure
    hold on 
    grid on
    title("Cost Trends: RP1/LOX First Stage, RP1/LOX Second Stage")
    xlabel("First Stage DV Fraction (X)")
    ylabel("Total Cost ($M2024)")
    plot(RP1_First_RP1_Second_results(:,6),RP1RP1Costs(:,1))
    plot(RP1_First_RP1_Second_results(:,6),RP1RP1Costs(:,2))
    plot(RP1_First_RP1_Second_results(:,6),RP1RP1Costs(:,3))
    plot(X_vals(RP1RP1Min_Index),RP1RP1Min,'O')
    legend("First Stage Production and NRE","Second Stage Production and NRE","Total Program Cost")
    
    figure
    hold on 
    grid on
    title("Cost Trends: RP1/LOX First Stage, Solid Second Stage")
    xlabel("First Stage DV Fraction (X)")
    ylabel("Total Cost ($M2024)")
    plot(RP1_First_Solid_Second_results(:,6),RP1SolidCosts(:,1))
    plot(RP1_First_Solid_Second_results(:,6),RP1SolidCosts(:,2))
    plot(RP1_First_Solid_Second_results(:,6),RP1SolidCosts(:,3))
    plot(X_vals(RP1SolidMin_Index),RP1SolidMin,'O')
    legend("First Stage Production and NRE","Second Stage Production and NRE","Total Program Cost")
    
    figure
    hold on 
    grid on
    title("Cost Trends: RP1/LOX First Stage, Storable Second Stage")
    xlabel("First Stage DV Fraction (X)")
    ylabel("Total Cost ($M2024)")
    plot(RP1_First_Storable_Second_results(:,6),RP1StorableCosts(:,1))
    plot(RP1_First_Storable_Second_results(:,6),RP1StorableCosts(:,2))
    plot(RP1_First_Storable_Second_results(:,6),RP1StorableCosts(:,3))
    plot(X_vals(RP1StorableMin_Index),RP1StorableMin,'O')
    legend("First Stage Production and NRE","Second Stage Production and NRE","Total Program Cost")

%--------------------------- HELPER FUNCTIONS -----------------------------
%Function: apply_SLVCM.
%Purpose: For a given set of inert masses, a learning curve and a total
%number of vehicles to produce, the function returns the 2024 adjusted
%total costs of the first stage (production+NRE), second stage, and both
%first and second stage together. Uses the NASA SLVCM.
%Inputs: m_in1, m_in2, learningCurve, numberFlights
    %m_in1: kilograms - Range unbounded. Stage 1 Inert Mass.
    %m_in2: kilograms - Range unbounded. Stage 2 Inert Mass.
%Ouputs: stageOneTotalCost,stageTwoTotalCost,combinedTotalCost
    %stageOneTotalCost: millions of dollars (2024) - Range unbounded. Total cost of first
        %stage production and NRE.
    %stageTwoTotalCost: millions of dollars (2024) - Range unbounded. Total
        %cost of second stage production and NRE.
    %combinedTotalCost: millions of dollars (2024) - Range unbounded. Total
        %cost of program.
%Authors: James Doss
function [stageOneTotalCost, stageTwoTotalCost, combinedTotalCost] = apply_SLVCM(m_in1, m_in2)
    learningCurve = 0.85; %Given in problem statement.
    numberFlights = 10; %Also given in problem statement.
    p = log(learningCurve)/log(2); %Find p value for the given learning curve.
    StageNR_A = 12.73; %Non recurring coefficients for costing a launch vehicle stage.
    StageNR_B = 0.55;
    NRE_1 = Dollars2023to2024(StageNR_A*m_in1^StageNR_B); %Calculate non-recurring costs. Notably this only works with inert stage masses in KILOGRAMS!
    NRE_2 = Dollars2023to2024(StageNR_A*m_in2^StageNR_B);
    Stage_A = 0.3024; %First unit costs of each stage.
    Stage_B = 0.662;
    Stage_1_Unit1 = Stage_A*m_in1^Stage_B; %C_1 of stage one and two.
    Stage_2_Unit1 = Stage_A*m_in2^Stage_B;
    stage1Costs = zeros(numberFlights,1); %Creating an array to store costs data.
    stage2Costs = zeros(numberFlights,1);
    stage1Costs(1) = Stage_1_Unit1; %First entry of each stage cost array is obviously the cost of the first unit.
    stage2Costs(1) = Stage_2_Unit1;
    for n = 2:numberFlights %For the remaining units, apply the costing formula and populate the array.
        Stage_1_Unit_n = Stage_1_Unit1*n^p;
        Stage_2_Unit_n = Stage_2_Unit1*n^p;
        stage1Costs(n) = Stage_1_Unit_n;
        stage2Costs(n) = Stage_2_Unit_n;
    end
    stageOneTotalCost = Dollars2023to2024(sum(stage1Costs)+NRE_1); %Sum all the stage costs and add the non-reccuring costs. Convert to 2024 dollars.
    stageTwoTotalCost = Dollars2023to2024(sum(stage2Costs)+NRE_2);
    combinedTotalCost = stageOneTotalCost+stageTwoTotalCost; %Since both stage costs are already in 2024 dollars we dont need to convert here.
end
%--------------------------------------------------------------------------
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
function [m1, m2, mpr1, mpr2, m_total, deltaV1_fraction, min1, min2] = vary_dv_isp(X, isp1, isp2)
    delt1 = 0.05; %Defining mass fractions. As given by problem statement.
    delt2 = 0.08;
    mpl = 63000; % kg - Required payload to orbit per problem statement.
    vtot = 9.8e3; % m/s - Required vehicle delta v per problem statement.
    g = 9.81; % m/s^2
    
    syms min1 min2 mpr1 mpr2 %Setting up 4 equations to solve as a system. The first two come from the definition of the mass fractions involved and the second two come from the required delta V.
    eq1 = delt1 == min1 / (min1 + mpr1 + min2 + mpr2 + mpl);
    eq2 = delt2 == min2 / (min2 + mpl + mpr2);
    eq3 = -isp1*g*log((min1 + min2 + mpl + mpr2)/(min1 + mpr1 + min2 + mpl + mpr2)) == X * vtot;
    eq4 = -isp2*g*log((min2 + mpl)/(min2 + mpl + mpr2)) == (1-X) * vtot;
    
    [min1, min2, mpr1, mpr2] = vpasolve([eq1 eq2 eq3 eq4], [min1 min2 mpr1 mpr2]); %Using VPAsolve to find the values and returning.
    min1 = double(min1);
    min2 = double(min2);
    mpr1 = double(mpr1);
    mpr2 = double(mpr2);
    m2 = min2 + mpr2;
    m1 = min1 + mpr1;
    m_total = m1 + m2 + mpl;
    deltaV1_fraction = X;
end
%--------------------------------------------------------------------------
%Function: Dollars2023to2024
%Purpose: For a given amount of money in 2023, calculate the present value
%in 2024 using a 2.5% inflation rate.
%Inputs: input
    % input: dollars - Range unbounded. 2023 dollars.
%Ouputs: output
    % output: dollars - Range unbounded. 2024 dollars.
%Authors: James Doss
function output = Dollars2023to2024(input)
    inflationRate = 0.025; %Nominal NASA (PPP) 2023-2024 inflation rate was 2.5%.
    output = (1+inflationRate)*input;
end