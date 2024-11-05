clear; clc; close all;

% LOX/LCH4, LOX/LH2, LOX/RP1, Solid, N2O4:UDMH (Storables)
% OX/FUEL RATIO
% SPECIFIC IMPULSE, SEA LVL (S)
% STAGE 1 THRUST (MN)
% STAGE 2 THRUST (MN)
% STAGE 1 Pe
% STAGE 2 Pe
% STAGE 1 P0
% STAGE 2 P0
% STAGE 1 NOZZLE AREA RATIO 
% STAGE 2 NOZZLE AREA RATIO 
% DENSITY Oxidizer (kg/m^3)
% DENSITY Fuel (kg/m^3)
% MER MASS COEFFICIENT Oxidizer M(tank)<kg> = c * M(prop)<kg>
% MER MASS COEFFICIENT Fuel M(tank)<kg> = c * M(prop)<kg>
% MER MASS COEFFICIENT Cryogenic Insulation Oxidizer M(insul)<kg> = c * A(tank)<m^2>
% MER MASS COEFFICIENT Cryogenic Insulation Fuel M(insul)<kg> = c * A(tank)<m^2>
prop_matrix = [3.6        6.03        2.72        NaN         2.67; 
              327         366         311         269         285;
              2.26        1.86        1.92        4.5         1.75;
              0.745       0.099       0.061       2.94        0.067;
              2.4         2.4         3.7         6.6         1.5;
              1.5         2.15        0.92        2.34        1.13;
              35.16       20.64       25.8        10.5        15.7;
              10.1        4.2         6.77        5           14.7;
              34.34       78          37          16          26.2;  
              45          84          14.5        56          81.3;
              1140        1140        1140        0           1442  
              423         71          820         1680        791;
              12.16/1140  12.16/1140  12.16/1140  0           12.16/1442; 
              12.16/423   9.09/71     12.16/820   12.16/1680  12.16/791; 
              1.123       1.123       1.123       0           0;
              1.123       2.88        0           0           0 
              ];

flights = 10;

%% Section 1
% the function vary_dv_isp finds mass values given isps of two stages and
% the dv fraction.
%% 1.2b
% the function gross_mass_min automates finding the minimum mass solution.
%% 1.3b
% the function min_total_cost automates the calculation of the minimum
% total cost solution.
%% 1.5
% now we will calculate the gross min mass for each combination of
% propellant, then put them into gross_mass_min_table
% shape (<gross min mass, dv_frac, m1, m2, mpr1, mpr2, cost>, second stage,
% first stage)
for ind1=1:5
    for ind2=1:5
        % find all relevant mass values and dv fraction for this prop mixture
        [gross_mass_min_table(1, ind2, ind1), ... % gross minimum mass
        gross_mass_min_table(2, ind2, ind1), ... % dv fraction
        gross_mass_min_table(3, ind2, ind1), ... % m1
        gross_mass_min_table(4, ind2, ind1), ... % m2
        gross_mass_min_table(5, ind2, ind1), ... % mpr1
        gross_mass_min_table(6, ind2, ind1)] ... % mpr2
        = gross_mass_min(ind1, ind2, prop_matrix);
        % find the total cost estimate (and discard the stage 1 and stage 2
        % individual costs)
        [~, ~, gross_mass_min_table(7, ind2, ind1)] = total_cost_estimate(gross_mass_min_table(3, ind2, ind1) - gross_mass_min_table(5, ind2, ind1), gross_mass_min_table(4, ind2, ind1) - gross_mass_min_table(6, ind2, ind1), flights);
    end
end
% convert gross mass to mt
gross_mass_min_table(1, :, :) = gross_mass_min_table(1, :, :)/1000;
% convert total cost estimate to ($B)
gross_mass_min_table(7, :, :) = gross_mass_min_table(7, :, :)/1000;
%%
% next, we will calculate the minimum program cost solution for each
% combination of propellant and organize it in min_program_cost_table
% shape (<min cost, dv_frac, gross_mass>, second stage, first stage)
min_program_cost_table = zeros(3, 5, 5);
% rows: second stage
% columns: first stage
for ind1=1:5
    for ind2=1:5
        isp1 = prop_matrix(2,ind1);
        isp2 = prop_matrix(2,ind2);
        [min_program_cost_table(1, ind2, ind1), ... % minimum cost
        min_program_cost_table(2, ind2, ind1), ... % dv_frac
        min_program_cost_table(3, ind2, ind1)] ... % gross mass
        = min_total_cost(isp1, isp2, flights);
    end
end
% convert mass to mt
min_program_cost_table(3, :, :) = min_program_cost_table(3, :, :)/1000;
% convert total cost estimate to ($B)
min_program_cost_table(1, :, :) = min_program_cost_table(1, :, :)/1000;
%%
% now we will make a table for each second stage/first stage combination
prop_names = ["LOX/LCH4", "LOX/LH2", "LOX/RP1", "Solid", "Storables"];
for ind2=1:5
    fprintf("------------------------------------------------------------------------\n")
    fprintf("Second stage prop. %s\n", prop_names(ind2));
    fprintf("------------------------------------------------------------------------\n")
    fprintf("\t\t\t\t\t\t\t%s\t%s\t\t%s\t\t%s\t\t%s\n", prop_names);
    fprintf("Minimum LV gross mass (t)\t%0.2f\t\t%0.2f\t\t%0.2f\t\t%0.2f\t\t%0.2f\n", gross_mass_min_table(1, ind2, :))
    fprintf("Min. LV mass sln. DV-frac \t%0.3f\t\t%0.3f\t\t%0.3f\t\t%0.3f\t\t%0.3f\n", gross_mass_min_table(2, ind2, :))
    fprintf("Min. LV mass sln. cost ($B)\t%0.3f\t\t%0.3f\t\t%0.3f\t\t%0.3f\t\t%0.3f\n", gross_mass_min_table(7, ind2, :))
    fprintf("------------------------------------------------------------------------\n")
    fprintf("Minimum cost sln. ($B2024) \t%0.3f\t\t%0.3f\t\t%0.3f\t\t%0.3f\t\t%0.3f\n", min_program_cost_table(1, ind2, :))
    fprintf("Min. cost sln. DV-frac \t\t%0.3f\t\t%0.3f\t\t%0.3f\t\t%0.3f\t\t%0.3f\n", min_program_cost_table(2, ind2, :))
    fprintf("Min. cost sln. gross mass(t)%0.2f\t\t%0.2f\t\t%0.2f\t\t%0.2f\t\t%0.2f\n", min_program_cost_table(3, ind2, :))
    fprintf("\n\n")
end
%%
% now creating the 5x5 matrix for minimum gross mass solution
fprintf("\t\t\t%s\t%s\t\t%s\t\t%s\t\t%s\n", prop_names);
for ind2=1:5
    if ind2 == 1 || ind2 == 5
        fprintf("%s\t%0.2f\t\t%0.2f\t\t%0.2f\t\t%0.2f\t\t%0.2f\n", prop_names(ind2), gross_mass_min_table(1, ind2, :))
    else
        fprintf("%s\t\t%0.2f\t\t%0.2f\t\t%0.2f\t\t%0.2f\t\t%0.2f\n", prop_names(ind2), gross_mass_min_table(1, ind2, :))
    end
end
%%
% now creating the 5x5 matrix for minimum cost solution
fprintf("\t\t\t%s\t%s\t\t%s\t\t%s\t\t%s\n", prop_names);
for ind2=1:5
    if ind2 == 1 || ind2 == 5
        fprintf("%s\t%0.3f\t\t%0.3f\t\t%0.3f\t\t%0.3f\t\t%0.3f\n", prop_names(ind2), min_program_cost_table(1, ind2, :))
    else
        fprintf("%s\t\t%0.3f\t\t%0.3f\t\t%0.3f\t\t%0.3f\t\t%0.3f\n", prop_names(ind2), min_program_cost_table(1, ind2, :))
    end
end

%(<gross min mass, dv_frac, m1, m2, mpr1, mpr2, cost>, second stage,
% first stage)

gross_mass_min_table(1, :, :) = gross_mass_min_table(1, :, :)*1000; % return to normal from mt
min_program_cost_table(3, :, :) = min_program_cost_table(3, :, :)*1000; % return to normal from mt
min_program_cost_table(1, :, :) = min_program_cost_table(1, :, :)*1000; % return to normal from $B2024

%% Section 2

% Pass 1 assuming cylindrical tanks with constant diameter

% now we will calculate the mass estimates for each subsystem of each stage
% and put them into a new table with shape:
% (<propellant, propellant tanks, propellant tank insulation, 
% engines, thrust structure, casing, gimbals, avionics, wiring, payload 
% fairing, intertank fairing, interstate fairing, aft fairing>, second stage, first stage)
minimized_mass_subsystem_table_stage1 = zeros(5, 5, 13);
minimized_mass_subsystem_table_stage2 = zeros(5, 5, 13);
% rows: second stage
% columns: first stage
% 3D dimension: subsystem mass

D_list = 5.2:0.1:10;

final_results = zeros(5,5);
for ind2 = 1:5
    for ind1 = 1:5
        % Initialization of useful parameters
        m0 = gross_mass_min_table(1,ind2,ind1);
        m1 = gross_mass_min_table(3,ind2,ind1);
        m2 = gross_mass_min_table(4,ind2,ind1);
        mpr1 = gross_mass_min_table(5,ind2,ind1);
        mpr2 = gross_mass_min_table(6,ind2,ind1);

        minimized_mass_subsystem_table_stage1(ind2,ind1,1) = mpr1; % propellant mass 1
        minimized_mass_subsystem_table_stage2(ind2,ind1,1) = mpr2; % propellant mass 2

        [mpr1_tank, mpr1_ox, mpr1_fuel] = calc_m_tank_pr(prop_matrix, mpr1, ind1);
        [mpr2_tank, mpr2_ox, mpr2_fuel] = calc_m_tank_pr(prop_matrix, mpr2, ind2);
        minimized_mass_subsystem_table_stage1(ind2,ind1,2) = mpr1_tank; % propellant tank mass 1 
        minimized_mass_subsystem_table_stage2(ind2,ind1,2) = mpr2_tank; % propellant tank mass 2
        
        optimized_eval_values_stage1 = NaN(length(D_list),13);
        optimized_eval_values_stage2 = NaN(length(D_list),13);

        i = 1;
        optimized_eval_values_diameter_set = zeros(length(D_list),30);
        for D=D_list % for loop to iterate through different L and D values
            % Radius
            r = D/2;

            [L] = calc_L(mpr1_ox, mpr2_ox, mpr1_fuel, mpr2_fuel, ind1, ind2, r, prop_matrix);
            % L = [L_ox1, L_ox2, L_fuel1, L_fuel2, L_entire_tanks]
            % L_entire_tanks = L_ox1 + L_ox2 + L_fuel1 + L_fuel2 + 4*D
           
            minimized_mass_subsystem_table_stage1(ind2,ind1,3) = calc_m_insul_pr(prop_matrix, L(1), L(3), ind1, r); % tank 1 insulation mass 
            minimized_mass_subsystem_table_stage2(ind2,ind1,3) = calc_m_insul_pr(prop_matrix, L(2), L(4), ind2, r); % tank 2 insulation mass
            
            [m_engines1, m_engines2, m_thrust_structures1, m_thrust_structures2, m_gimbals1, m_gimbals2] = calc_mass_required_engines_thrust_structure_gimbals(m0, m1, mpr1, mpr2, prop_matrix, ind1, ind2);
            minimized_mass_subsystem_table_stage1(ind2,ind1,4) = m_engines1; % total mass for stage 1 engines
            minimized_mass_subsystem_table_stage2(ind2,ind1,4) = m_engines2; % total mass for stage 2 engines
    
            minimized_mass_subsystem_table_stage1(ind2,ind1,5) = m_thrust_structures1; % total thrust structure mass 1
            minimized_mass_subsystem_table_stage2(ind2,ind1,5) = m_thrust_structures2; % total thrust structure mass 2
    
            minimized_mass_subsystem_table_stage1(ind2,ind1,6) = calc_m_casing(mpr1, ind1); % solids casing mass
            minimized_mass_subsystem_table_stage2(ind2,ind1,6) = calc_m_casing(mpr2, ind2); % solids casing mass
    
            minimized_mass_subsystem_table_stage1(ind2,ind1,7) = m_gimbals1; % total gimbal mass 1
            minimized_mass_subsystem_table_stage2(ind2,ind1,7) = m_gimbals2; % total gimbal mass 2
    
            m02_gross = gross_mass_min_table(4,ind2,ind1) + gross_mass_min_table(6,ind2,ind1);
            m01_gross = gross_mass_min_table(3,ind2,ind1) + gross_mass_min_table(5,ind2,ind1);
            minimized_mass_subsystem_table_stage1(ind2,ind1,8) = 10*m01_gross^0.361; % avionics mass
            minimized_mass_subsystem_table_stage2(ind2,ind1,8) = 10*m02_gross^0.361; % avionics mass
            minimized_mass_subsystem_table_stage1(ind2,ind1,9) = 1.058*sqrt(m01_gross)*(L(1) + L(3))^0.25; % wiring mass
            minimized_mass_subsystem_table_stage2(ind2,ind1,9) = 1.058*sqrt(m02_gross)*(L(2) + L(4))^0.25; % wiring mass

            [m_fairing_array, h_fairing_total] = calc_mass_height_fairing(D, r); %m_fairing_array = [M_fair_payload_total,  M_fair_intertank1, M_fair_interstage, M_fair_intertank2, M_fair_aft]
            minimized_mass_subsystem_table_stage1(ind2,ind1,10) = 0; % payload fairing considered as stage 2
            minimized_mass_subsystem_table_stage2(ind2,ind1,10) = m_fairing_array(1); % payload fairing
            minimized_mass_subsystem_table_stage1(ind2,ind1,11) = m_fairing_array(2); % intertank1 fairing mass
            minimized_mass_subsystem_table_stage2(ind2,ind1,11) = m_fairing_array(4); % intertank1 fairing mass
            minimized_mass_subsystem_table_stage1(ind2,ind1,12) = m_fairing_array(3); % interstage fairing mass
            minimized_mass_subsystem_table_stage2(ind2,ind1,12) = 0; % interstage fairing considered as stage 1
            minimized_mass_subsystem_table_stage1(ind2,ind1,13) = m_fairing_array(5); % aft fairing mass
            minimized_mass_subsystem_table_stage2(ind2,ind1,13) = 0; % aft fairing considered as stage 1
            
            optimized_eval_values = zeros(length(D_list),29);
            if ~any(isnan(L), 'all') % Only evaluates if all L values are mathematically possible solutions
                if (L(5) + h_fairing_total)/D <= 12 % Checks if L/D meets specs 
                    %disp(prop_names(ind1) + " + " + prop_names(ind2));
                    %disp("Total Height: " + L(5) + h_fairing_total + " | L/D: " + (L(5) + h_fairing_total)/D)         
                    stage1 = reshape(minimized_mass_subsystem_table_stage1(ind2,ind1,1:13), 1, 13);
                    stage2 = reshape(minimized_mass_subsystem_table_stage2(ind2,ind1,1:13), 1, 13);
                    optimized_eval_values(i, :) = [sum(minimized_mass_subsystem_table_stage1(ind2,ind1,1:13)) + sum(minimized_mass_subsystem_table_stage2(ind2,ind1,1:13)), (L(5) + h_fairing_total)/D, (L(5) + h_fairing_total), stage1, stage2];
                    test = 1;
                    % disp(optimized_eval_values(i, :))
                    % [new total inert mass of all subsystems, L/D ratio, total height, stage 1 subsystem mass results(13 long), stage2 subsystem mass results(13 long)]
                end
            end
            disp("test1")
            min1 = m1 - mpr1;
            min2 = m2 - mpr2;
            disp(optimized_eval_values(i, :))
            %disp(optimized_eval_values(1,1));
            [max_total_inert_mass, max_index] = max(optimized_eval_values(:,1));
            max_inert_mass_margin = ((max_total_inert_mass/(min1 + min2)) - 1)*100;
            optimized_eval_values_diameter_set(i, :) = [max_inert_mass_margin, optimized_eval_values(max_index,1:29)];
            i = i + 1;
        end
        min1 = m1 - mpr1;
        min2 = m2 - mpr2;
        [max_total_inert_mass, max_index] = max(optimized_eval_values_diameter_set(:,1));
        max_inert_mass_margin = ((max_total_inert_mass/(min1 + min2)) - 1)*100;
        final_results(ind2,ind1) = max_inert_mass_margin;
    end
end

% shape (<propellant, propellant tanks, propellant tank insulation, 
% engines, thrust structure, casing, gimbals, avionics, wiring, 
% payload fairing, intertank fairing, interstage fairing, aft fairing>,
% second stage, first stage)

% Helper Functions
function [m_tank_pr, mpr_ox, mpr_fuel] = calc_m_tank_pr(prop_matrix, mpr, ind)
        if ~isnan(prop_matrix(1,ind)) % tests if there is an ox/fuel ratio
            mpr_fuel = mpr/(prop_matrix(1,ind) + 1);
            mpr_ox = prop_matrix(1,ind)*mpr_fuel;
            m_tank_pr = mpr_ox*prop_matrix(13,ind) + mpr_fuel*prop_matrix(14,ind);
        else % for case of Solid where there is no oxidizer
            mpr_ox = 0;
            mpr_fuel = mpr;
            m_tank_pr = mpr_fuel*prop_matrix(14,ind);
        end
end

function mass = fair_mass(area)
    % Fairing mass: 4.95*[A_fairing]^1.15
    mass = 4.95*area^1.15;
end

function L = calc_L(mpr1_ox, mpr2_ox, mpr1_fuel, mpr2_fuel, ind1, ind2, r, prop_matrix)
    % Determine the length required of each tank
    D = 2*r;
    L_ox1 = 0; L_ox2 = 0; % default set to 0
    if prop_matrix(11,ind1) ~= 0 % tests if there is an oxidizer, if not then cryogenic tank insulation mass for oxidizer 1 is default 0
        Vpr1_ox = mpr1_ox/prop_matrix(11,ind1);
        if Vpr1_ox < (4/3)*pi*r^3
            L_ox1 = NaN;
        else
            L_ox1 = (Vpr1_ox - (4/3)*pi*r^3)/(pi*r^2);
        end
    end
    Vpr1_fuel = mpr1_fuel/prop_matrix(12,ind1);
    if Vpr1_fuel < (4/3)*pi*r^3
        L_fuel1 = NaN;
    else
        L_fuel1 = (Vpr1_fuel - (4/3)*pi*r^3)/(pi*r^2);
    end
    if prop_matrix(11,ind2) ~= 0 % tests if there is an oxidizer, if not then cryogenic tank insulation mass for oxidizer 2 is default 0
        Vpr2_ox = mpr2_ox/prop_matrix(11,ind2);
        if Vpr2_ox < (4/3)*pi*r^3
            L_ox2 = NaN;
        else
            L_ox2 = (Vpr2_ox - (4/3)*pi*r^3)/(pi*r^2);
        end
    end
    Vpr2_fuel = mpr2_fuel/prop_matrix(12,ind2);
    if Vpr2_fuel < (4/3)*pi*r^3
        L_fuel2 = NaN;
    else
        L_fuel2 = (Vpr2_fuel - (4/3)*pi*r^3)/(pi*r^2);
    end
    L_entire_tanks = L_ox1 + L_ox2 + L_fuel1 + L_fuel2 + 4*D;
    L = [L_ox1, L_ox2, L_fuel1, L_fuel2, L_entire_tanks];
end

function m_insul_pr = calc_m_insul_pr(prop_matrix, L_ox, L_fuel, ind, r)
        m_insul_pr_ox = 0; m_insul_pr_fuel = 0; % default
        if prop_matrix(11,ind) ~= 0 % tests if there is an oxidizer, if not then cryogenic tank insulation mass for oxidizer 1 is default 0
            SA_ox = 4*pi*r^2 + 2*pi*r*L_ox;
            m_insul_pr_ox = prop_matrix(15,ind)*SA_ox;
        end
        if prop_matrix(16,ind) ~= 0 % tests if cryogenic tank insulation is required for fuel 1 
            SA_fuel = 4*pi*r^2 + 2*pi*r*L_fuel;
            m_insul_pr_fuel = prop_matrix(16,ind)*SA_fuel;
        end
        m_insul_pr = m_insul_pr_ox + m_insul_pr_fuel;
end

function m_engine = engine_mass(mpr, thrust_N, ind, Ae_At)
    if ind ~= 4 % determine if engine is solid
        m_engine = (7.81 * 10^-4)*thrust_N + (3.37 * 10^-5)*thrust_N*sqrt(Ae_At) + 59;
    else
        m_engine = 0.135*mpr;
    end
end

function [m_engines1, m_engines2, m_thrust_structures1, m_thrust_structures2, m_gimbals1, m_gimbals2] = calc_mass_required_engines_thrust_structure_gimbals(m0, m1, mpr1, mpr2, prop_matrix, ind1, ind2)
    thrust_req1= 1.3*m0;
    thrust_per_engine1 = prop_matrix(3, ind1)*10^6;
    num_engines1 = ceil(thrust_req1 / thrust_per_engine1);
    Ae_At1 = prop_matrix(9,ind1); % Exit-Throat Area Ratio 1
    P0_1 = prop_matrix(7,ind1);
    m_engines1 = num_engines1*engine_mass(mpr1, thrust_per_engine1, ind1, Ae_At1);
    m_thrust_structures1 = num_engines1*(2.55*10^-4 * thrust_per_engine1);
    m_gimbals1 = num_engines1*237.8*(thrust_per_engine1/P0_1)^0.9375;

    thrust_req2 = 0.76*(m0-m1);
    thrust_per_engine2 = prop_matrix(4, ind2)*10^6;
    num_engines2 = ceil(thrust_req2 / thrust_per_engine2);
    Ae_At2 = prop_matrix(10,ind2); % Exit-Throat Area Ratio 2
    P0_2 = prop_matrix(8,ind2);
    m_engines2 = num_engines2*engine_mass(mpr2, thrust_per_engine2, ind2, Ae_At2);
    m_thrust_structures2 = num_engines2*(2.55*10^-4 * thrust_per_engine2);
    m_gimbals2 = num_engines2*237.8*(thrust_per_engine2/P0_2)^0.9375;
end

function m_casing = calc_m_casing(mpr, ind)
    m_casing = 0; % default
    if ind == 4
        m_casing = 0.135*mpr;
    end
end

function [m_fairing, h_fairing_total] = calc_mass_height_fairing(d, r)
    % Calculate all fairing areas
    h_fair_intertank = 0.4*d;
    h_fair_interstage = 1*d;
    h_fair_aft = 0.5*d;
    
    vary_N_results = zeros(length(1:1:10),6);
    % holds [M_fair_payload_total,  M_fair_intertank1, M_fair_interstage, M_fair_intertank2, M_fair_aft, h_fair_total]
    i = 1;
    for N = 1:0.1:10
        V_fair_payload_min = pi*(2.6^2)*13; 
        h_fair_cone = ceil((V_fair_payload_min/(pi*r^2))/(N + 1/3)*10)/10;
        h_fair_cylinder = 3*h_fair_cone;
    
        A_fair_cone = pi*r*sqrt(r^2 + h_fair_cone^2);
        A_fair_cylinder = pi*d*h_fair_cylinder;
        A_fair_intertank1 = pi*d*h_fair_intertank;
        A_fair_interstage = 1.5*pi*d*h_fair_interstage;
        A_fair_intertank2 = pi*d*h_fair_intertank;
        A_fair_aft = pi*d*h_fair_aft;
    
        % Calculate all fairing masses
        M_fair_payload_total = fair_mass(A_fair_cone) + fair_mass(A_fair_cylinder);
        M_fair_intertank1 = fair_mass(A_fair_intertank1);
        M_fair_interstage = fair_mass(A_fair_interstage);
        M_fair_intertank2 = fair_mass(A_fair_intertank2);
        M_fair_aft = fair_mass(A_fair_aft);
        m_fairing = [M_fair_payload_total,  M_fair_intertank1, M_fair_interstage, M_fair_intertank2, M_fair_aft];
        h_fairing_total = h_fair_cone + h_fair_cylinder + 2*h_fair_intertank + h_fair_interstage + h_fair_aft;
        vary_N_results(i, :) = [m_fairing(1:5), h_fairing_total]; 
        i = i + 1;
    end
    [min_h_fairing_total, min_index] = min(vary_N_results(:,6));
    m_fairing = vary_N_results(min_index, 1:5);
    h_fairing_total = min_h_fairing_total;
end

% min_total_cost: finds the minimum total cost given the two stage isps and
% the number of flights
% Parameters:
% isp1: isp of the first stage propellant
% isp2: isp of the second stage propellant
% flights: number of flights
% Returns:
% min_cost: the minimized gross cost, in $2024
% dv_frac: relevant dv fraction for the minimized gross cost
% gross_mass: relevant gross mass for the minimized gross cost
function [min_cost, dv_frac, gross_mass] = min_total_cost(isp1, isp2, flights)
mpl = 63000;
dvgrid = 0:0.001:1;
% store cost with respect to dv fraction, shape = (<cost, gross_mass>,
% length dvgrid)
cost_wrt_dv = zeros(2, length(dvgrid));
for ind=1:length(dvgrid)
    % calculate the relevant masses given prop combinations and dv_frac
    [m1, m2, mpr1, mpr2] = vary_dv_isp(dvgrid(ind), isp1, isp2);
    % store the cost estimate for the prop combination and dv_frac
    [~, ~, cost_wrt_dv(1, ind)] = total_cost_estimate(m1-mpr1, m2-mpr2, flights);
    % calculate gross mass of this solution
    cost_wrt_dv(2, ind) = m1 + m2 + mpl;
end
% find the minimum cost along cost_wrt_dv
[min_cost, min_cost_ind]= min(cost_wrt_dv(1, :));
% find the relevant dv fraction
dv_frac = dvgrid(min_cost_ind);
% find the relevant gross mass
gross_mass = cost_wrt_dv(2, min_cost_ind);
end

% total_cost_estimate: calculate the total cost estimate using the NASA
% approximation model. adjusts for inflation from 2023 to 2024
% Parameters:
% mi1: inert mass of stage 1
% mi2: inert mass of stage 2
% flights: number of flights
% Returns:
% stage1_costs: stage 1 costs including NR and production of all units, in $2024
% stage2_costs: stage 2 costs including NR and production of all units, in $2024
% cost: total cost estimation, in $2024
function [stage1_costs, stage2_costs, cost] = total_cost_estimate(mi1, mi2, flights)
% curve defined in problem statement and inflation factor from $2023 to
% $2024
curve = 0.85;
p = log(curve)/log(2);
infl = 1.025;
% total NR costs for both the first and second stage
NRcost1 = infl*(12.73*(mi1)^0.55);
NRcost2 = infl*(12.73*(mi2)^0.55);
% production cost for first unit for each stage
first_prod_stage_1 = infl*0.3024*(mi1)^0.662;
first_prod_stage_2 = infl*0.3024*(mi2)^0.662;
% now sum all stage 1 and stage 2 costs individually for each unit
% required for each flight
stage1_costs = NRcost1;
stage2_costs = NRcost2;
for ind=1:flights
stage1_costs = stage1_costs + first_prod_stage_1*ind^p;
stage2_costs = stage2_costs + first_prod_stage_2*ind^p;
end
% sum both stage 1 and stage 2 costs for total cost
cost = stage1_costs + stage2_costs;
end

% gross_mass_min: finds the minimum gross mass given the two stage
% propellant types. all masses returned are in kg.
% Parameters:
% ind1: index of the stage 1 propellant
% ind2: index of the stage 2 propellant
% prop_matrix: relevant prop matrix defined above
% Returns:
% gross_mass_min: minimized gross mass
% dv_frac: relevant dv fraction (x-frac) that returns minimal mass solution
% m1: mass of stage 1 (inert + propellant) " "
% m2: mass of stage 2 " "
% mpr1: mass of stage 1 propellant " "
% mpr2: mass of stage 2 propellant " "
function [gross_mass_min, dv_frac, m1, m2, mpr1, mpr2] = gross_mass_min(ind1, ind2, prop_matrix)
% x-fraction grid to use during search for minimum mass
dvgrid = 0:0.001:1;
% results of vary_dv_isp: shape = (dv fraction variation, res params)
res = zeros(length(dvgrid), 6);
% relevant Isps based on prop_matrix
isp1 = prop_matrix(2,ind1);
isp2 = prop_matrix(2,ind2);
% populate res matrix to gather all relevant mass data while varying
% dv_frac
for x_ind=1:length(dvgrid)
x = dvgrid(x_ind);
[res(x_ind, 1), res(x_ind, 2), res(x_ind, 3), res(x_ind, 4)] = vary_dv_isp(x, isp1, isp2);
end
mpl = 63000;
% search for the minimum mass solution within boundarys of 0.3 < xfrac < 0.7
start_ind = find(dvgrid==0.3);
end_ind = find(dvgrid==0.7);
% find min gross mass solution and relevant index
[gross_mass_min, min_ind] = min(res(start_ind:end_ind, 1) + res(start_ind:end_ind, 2) + mpl);
min_ind = min_ind + start_ind;
% calculate and return relevant masses and dv_frac
dv_frac = dvgrid(min_ind);
m1 = res(min_ind, 1);
m2 = res(min_ind, 2);
mpr1 = res(min_ind, 3);
mpr2 = res(min_ind, 4);
end
% vary_dv_isp: calculates estimated mass values for rocket that satisfies
% design parameters given required xfrac and isps. some results may be
% non-physical, but analysis of this function will be only on data
% generated with reasonable x-fractions, namely, with minimums between
% an x-frac of 0.3 to 0.7. all masses returned are in kg.
% Parameters:
% xfrac: fraction of delta V that stage 1 must provide
% isp1: isp of the first stage
% isp2: isp of the second stage
% Returns:
% m1: mass of stage 1 (inert and propellant)
% m2: mass of stage 2
% mpr1: mass of propellant of stage 1
% mpr2: mass of propellant of stage 2
function [m1, m2, mpr1, mpr2] = vary_dv_isp(xfrac, isp1, isp2)
% design parameters as defined in the problem statement
delt1 = 0.05;
delt2 = 0.08;
mpl = 63000;
vtot = 9.8e3; % m/s
g = 9.81;
% solving the system of equations based on definitions of stage inert mass
% fraction and dv
options = optimset('Display','off');
res = fsolve(@(x)eqs(x, xfrac, delt1, delt2, isp1, isp2, mpl, vtot, g), [1e5; 1e5;
1e5; 1e5], options);
% calculate and return the relevant mass values
min1 = res(1);
min2 = res(2);
mpr1 = res(3);
mpr2 = res(4);
m2 = min2 + mpr2;
m1 = min1 + mpr1;
end

% equations used for the fsolve function in vary_dv_isp
function F = eqs(x, xfrac, delt1, delt2, isp1, isp2, mpl, vtot, g)
% x(1) = min1
% x(2) = min2
% x(3) = mpr1
% x(4) = mpr2
% 2 equations defining stage inert mass fraction
% delt1 = min1 / (min1 + mpr1 + min2 + mpr2 + mpl)
% delt2 = min2 / (min2 + mpl + mpr2)
% 2 equations defining the dv1 and dv2
% dv1 = -isp1*g*ln( (min1 + min1 + min2 + mpl + mpr2)/(min1 + min1 + mpl1 + min2 + mpl + mpr2)) = x * vtot
% dv2 = -isp2*g*ln((min1 + mpl)/(min1 + mpl + mpr2)) = (1-x) * vtot
F(1) = x(1) / (x(1) + x(3) + x(2) + x(4) + mpl) - delt1;
F(2) = x(2) / (x(2) + mpl + x(4)) - delt2;
F(3) = -isp1*g*log( (x(1) + x(2) + mpl + x(4))/(x(1) + x(3) + x(2) + mpl + x(4))) - xfrac * vtot;
F(4) = -isp2*g*log((x(2) + mpl)/(x(2) + mpl + x(4))) - (1-xfrac) * vtot;
end
