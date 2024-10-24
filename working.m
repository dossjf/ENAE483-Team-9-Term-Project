clear 
clc
close all
%% Propellant Data
% OX/FUEL RATIO
% SPECIFIC IMPULSE, SEA LVL (S)
% STAGE 1 THRUST (MN)
% STAGE 2 THRUST (MN)
% STAGE 1 De
% STAGE 2 De
% STAGE 1 P0
% STAGE 2 P0
% STAGE 1 NOZZLE AREA RATIO 
% STAGE 2 NOZZLE AREA RATIO 
prop_matrix = [3.6      6.03    2.72    NaN     2.67;
              327      366     311     269     285;
              2.26     1.86    1.92    4.5     1.75;
              0.745    0.099   0.061   2.94    0.067;
              2.4      2.4     3.7     6.6     1.5;
              1.5      2.15    0.92    2.34    1.13;
              35.16    20.64   25.8    10.5    15.7;
              10.1     4.2     6.77    5       14.7;
              34.34    78      37      16      26.2;  
              45       84      14.5    56      81.3];

LOX_CH4=prop_matrix(2,1);
LOX_LH2=prop_matrix(2,2);
LOX_RP1=prop_matrix(2,3);
Solid_prop=prop_matrix(2,4);
Storables_prop=prop_matrix(2,5);


propellant_combinations = [
    LOX_CH4, LOX_CH4;
    LOX_CH4, LOX_LH2;
    LOX_CH4, LOX_RP1;
    LOX_CH4, Solid_prop;
    LOX_CH4, Storables_prop;
    
    LOX_LH2, LOX_CH4;
    LOX_LH2, LOX_LH2;
    LOX_LH2, LOX_RP1;
    LOX_LH2, Solid_prop;
    LOX_LH2, Storables_prop;
    
    LOX_RP1, LOX_CH4;
    LOX_RP1, LOX_LH2;
    LOX_RP1, LOX_RP1;
    LOX_RP1, Solid_prop;
    LOX_RP1, Storables_prop;
    
    Solid_prop, LOX_CH4;
    Solid_prop, LOX_LH2;
    Solid_prop, LOX_RP1;
    Solid_prop, Solid_prop;
    Solid_prop, Storables_prop;
    
    Storables_prop, LOX_CH4;
    Storables_prop, LOX_LH2;
    Storables_prop, LOX_RP1;
    Storables_prop, Solid_prop;
    Storables_prop, Storables_prop
];
%% Data generation
results=data_gen(propellant_combinations);
% results(iteration step{1:71}, propellant number{1:25}, {1,2,3,4,5,6,7,8})
%  {xfrac->1,m01_val->2, m02_val->2, mpr1_val->3, mpr2_val->4,mtot->5,mi1->6,mi2->7}

%Data compliation
[design_matrix,total_mass,cost,v1frac]=design(results); %total mass ,cost and v1frac for debugging
% and to examine one propellant in particular 
% design_matrix(propellant number{1:25}, {1,2,3,4,5,6}) 
% design_matrix(j,:)={min_mass/1000,v1_frac_mass,min_mass_cost,min_cost,v1_frac_cost,mass_min_cost/1000};
%% Sub system analysis 
[num_iterations, num_combinations, ~] = size(results);

dv1_frac=zeros(num_combinations,1);
m01=zeros(num_combinations,1); % Stores total mass of stg1 b4 mass estimation
% for comparsion
m02=zeros(num_combinations,1);% Stores total mass of stg2 b4 mass estimation
% for comparsion
isp1=zeros(num_combinations,1);
isp2=zeros(num_combinations,1);
sub_sys1=zeros(num_combinations,2); % Stores stage 1 values in this
sub_sys2=zeros(num_combinations,2); % stores stage 2 values

mpr1=zeros(num_combinations,1); %Need for mass estimation(i think)
mpr2=zeros(num_combinations,1);

%densities
rho_LH2=71;
rho_RP=820;
rho_CH4=423;
rho_solids=1680;
rho_storables=791;



for i=1:num_combinations
    dv1_frac(i,1)=design_matrix(i,2);
    for j=1:num_iterations
        if dv1_frac(i,1)==results(j,i,1)
            m01(i)=results(j,i,2);
            m02(i)=results(j,i,3);
            mpr1(i)=results(j,i,4);
            mpr2(i)=results(j,i,5);
        end
    end

    isp1(i)=propellant_combinations(i,1);
    isp2(i)=propellant_combinations(i,2);

    % Stage 1 mass estimation
    if isp1(i)==LOX_LH2
        %Change these 
        mpr=mpr1;
        Ae_At=prop_matrix(2,8); %Nozzle area ratio
        Po1=prop_matrix(2,6)*1e6; % Chamber pressure (Pa)
        T=prop_matrix(2,2)*1e6; %Thrust (N)
        mfuel=mpr/6.03; %mass of fuel
        mLOX=mpr-(mpr/6.03);%mass of oxidizer
        
        % Change insulation 
        mass_avions=10*(mpr^0.361);
        V=(1/rho_LH2)*mfuel+(1/rho_LOX)*mLOX % Need to check this
        r=(V/(4*pi/3))^(1/3); % radius
        A=4*pi*r^2; % area
        mass_tank=9.09*V; %mass of tank
        insulation_mass=2.88*A; % tank insulation mass
        mass_engine=(7.81e-4)*T+(3.37e-5)*T*sqrt(Ae_At);
        mass_motor_casing=0.135*fuel;
        mass_thrust_struct=(2.55e-4)*T;
        mass_gimbals=237.8*(T/Po1)^.9375;
    elseif(isp1(i)==LOX_CH4)
          
    elseif(isp1(i)==LOX_RP1)
          
    elseif(isp1(i)==Storables_prop)
           
    elseif(isp1(i)==Solid_prop)
         
    end
    % Stage 2 mass estimation
     if isp2(i)==LOX_LH2
        sub_sys(j,2)=(1/rho_LH2)*mpr2;
    elseif(isp2(i)==LOX_CH4)

    elseif(isp2(i)==LOX_RP1)

    elseif(isp2(i)==Storables_prop)

    elseif(isp2(i)==Solid_prop)

     end

     % Sum of values 
     %sum(...)
end








%% Functions support
%creates mass results by varying stage1 (unfiltered)
function results=data_gen(propellant_combinations)
num_iterations = 0.3:0.01:1; 
results = zeros(length(num_iterations), size(propellant_combinations, 1), 8);
m=1;
% Runs through propllant values and stores them in results(71,6,5)
for i = num_iterations
    % Loop through the 6 combinations of isp1 and isp2
    for j = 1:size(propellant_combinations, 1)
        isp1 = propellant_combinations(j, 2);
        isp2 = propellant_combinations(j, 1);
 
        
        % Call the function with the current delta-v and specific impulses
        [xfrac,m01_val, m02_val, mpr1_val, mpr2_val,min1,min2] = vary_dv_isp(i, isp1, isp2);
        % Stores results in array results;
        mtot=m01_val+m02_val;
        mi1=min1;
        mi2=min2;
        results(m, j, :) = [xfrac,m01_val, m02_val, mpr1_val, mpr2_val,mtot,mi1,mi2];
    end
    m=m+1;
    
end

end
%Calculates costs associated with all stages
function [design_matrix,total_mass,cost,xfrac]=design(results)

    flights=10;
    [num_iterations, num_combinations, ~] = size(results);
    design_matrix=zeros(num_combinations,6);
    % Loop over each propellant combination
    for j = 1:num_combinations
        % Initialize arrays to hold dv1 and masses for this combination
        %dv1_vals = zeros(num_iterations, 1);
        mpr1=zeros(num_iterations,1);
        mpr2=zeros(num_iterations,1);
        mi1=zeros(num_iterations,1);
        mi2=zeros(num_iterations,1);
        mass_stage1 = zeros(num_iterations, 1);
        mass_stage2 = zeros(num_iterations, 1);
        total_mass = zeros(num_iterations, 1);
        cost=zeros(num_iterations,1);
        xfrac=zeros(num_iterations,1);
        for i = 1:num_iterations
            % Extracting values from the results array
            v1frac = results(i, j, 1);
            m01_val = results(i, j, 2);
            m02_val = results(i, j, 3);
            mpr1_val = results(i, j, 4);
            mpr2_val = results(i, j, 5);
            total_mass_val=results(i,j,6);
            mi1_val=results(i,j,7);
            mi2_val=results(i,j,8);
         
            % Calculate the total mass of each stage and the total mass
            mpr1(i)=mpr1_val;
            mpr2(i)=mpr2_val;
            mass_stage1(i) = m01_val;
            mass_stage2(i) = m02_val;
            mi1(i)=mi1_val;
            mi2(i)=mi2_val;
            total_mass(i) = total_mass_val;
            xfrac(i)=v1frac;
            [~, ~, cost(i)] = total_cost_estimate(mi1(i), mi2(i), flights);
        end % For loop j iteration end (combinations)
        for i=1:num_iterations
            [~, ~, cost(i)] = total_cost_estimate(mi1(i), mi2(i), flights);
   
        end
    [min_mass, idx_mass] = min(total_mass); % min mass and its index 
    v1_frac_mass=results(idx_mass,j,1); %v1 fraction
    min_mass_cost=cost(idx_mass);

    [min_cost,idx_cost]=min(cost); % min cost and its index
    v1_frac_cost=results(idx_cost);
    mass_min_cost=total_mass(idx_cost);

    design_matrix(j,:)=[min_mass/1000,v1_frac_mass,min_mass_cost,min_cost,v1_frac_cost,mass_min_cost/1000];
    end % END STORING variables

end

function [stage1_costs, stage2_costs, cost] = total_cost_estimate(mi1, mi2, flights)
    curve = 0.85;
    p = log(curve)/log(2);
    
    
    NRcost1 = (12.73*(mi1)^0.55);
    NRcost2 = (12.73*(mi2)^0.55);
    first_prod_stage_1 = 0.3024*(mi1)^0.662;
    first_prod_stage_2 = 0.3024*(mi2)^0.662;
    
    stage1_costs = 0;
    stage2_costs = 0;
    for ind=1:flights
        stage1_costs = stage1_costs + first_prod_stage_1*ind^p;
        stage2_costs = stage2_costs + first_prod_stage_2*ind^p;
    end
    
    cost = (stage1_costs+NRcost2+NRcost1 + stage2_costs)*1.025;

end


%Establishes eqns used, passed into vary_dv_isp
function F = eqs(x, xfrac, delt1, delt2, isp1, isp2, mpl, vtot, g)
% x(1) = min1
% x(2) = min2
% x(3) = mpr1
% x(4) = mpr2
F(1) = x(1) / (x(1) + x(3) + x(2) + x(4) + mpl) - delt1;
F(2) = x(2) / (x(2) + mpl + x(4)) - delt2;
F(3) = -isp1*g*log( (x(1) + x(2) + mpl + x(4))/(x(1) + x(3) + x(2) + mpl + x(4))) - xfrac * vtot;
F(4) = -isp2*g*log((x(2) + mpl)/(x(2) + mpl + x(4))) - (1-xfrac) * vtot;

end

%Solves eqns defined by eqs function, calls eqs
function [xfrac,m1, m2, mpr1, mpr2,min1,min2] = vary_dv_isp(xfrac, isp1, isp2)

delt1 = 0.05;
delt2 = 0.08;
mpl = 63000;
vtot = 9.8e3; % m/s
g = 9.81; 

% delt1 = min1 / (min1 + mpr1 + min2 + mpr2 + mpl)
% delt2 = min2 / (min2 + mpl + mpr2)
% dv1 = -isp1*g*ln( (min1 + min1 + min2 + mpl + mpr2)/(min1 + min1 + mpl1 + min2 + mpl + mpr2)) = x * vtot
% dv2 = -isp2*g*ln((min1 + mpl)/(min1 + mpl + mpr2)) = (1-x) * vtot

% syms min1 min2 mpr1 mpr2
% eq1 = delt1 == min1 / (min1 + mpr1 + min2 + mpr2 + mpl);
% eq2 = delt2 == min2 / (min2 + mpl + mpr2);
% eq3 = -isp1*g*log( (min1 + min2 + mpl + mpr2)/(min1 + mpr1 + min2 + mpl + mpr2)) == x * vtot;
% eq4 = -isp2*g*log((min2 + mpl)/(min2 + mpl + mpr2)) == (1-x) * vtot;

res = fsolve(@(x)eqs(x, xfrac, delt1, delt2, isp1, isp2, mpl, vtot, g), [1e5; 1e5; 1e5; 1e5]);
if any(res)
    min1 = res(1);
    min2 = res(2);
    mpr1 = res(3);
    mpr2 = res(4);

    m2 = min2 + mpr2+mpl;
    m1 = min1 + mpr1;
else
    min1=NaN;
    min2=NaN;
    m1 = NaN;
    m2 = NaN;
    mpr1 = NaN;
    mpr2 = NaN;
end

end
