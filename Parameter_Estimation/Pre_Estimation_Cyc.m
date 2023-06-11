clear;
clc;

%% Constants
R_g = 8.314; % J/mol-K (Universal gas Constant)
T_ref = 298.15; % K (Reference Temperature)
F = 96485; % C/mol (Faradays Constant)
x_SOC_0 = 0.0085; % Degree of Lithiation at 0 % SOC
x_SOC_100 = 0.78; % Degree of Lithiation at 100 % SOC
U_a_Ref = 0.1233; % V (reverence OCV of anode)

%% Previous results of calendar aging
k_Cal_Ref = 3.6940e-04;
E_a_Cal = 2.0493e+04;
alpha = 0.3840;
k0 = 0.142;

%% High Temperature
Cyc_HT_Act_Data = readmatrix('Cycle_Data.xlsx','Sheet','Pre_Est_HT');
Temperature_Act_Total = Cyc_HT_Act_Data(:,1); % temperature in kelvin
Cycles_Act_Total = Cyc_HT_Act_Data(:,2); % number of cycles
Fade_Act_Total = Cyc_HT_Act_Data(:,3)/100; % fractional capacity fade
Qt = 3*2*Cycles_Act_Total; % Total charge throughput, 3 because C0 = 3 Ah
Number_of_Rows = numel(Cycles_Act_Total);
j = 1; % Counting index
for i = 1:(Number_of_Rows-1)
    % Following if condition notes where temperature in changing
    if Temperature_Act_Total(i) ~= Temperature_Act_Total(i+1)
        Break_Points(j) = i; % This arrays keeps track of temperature changes
        j = j + 1;
    end
end
Break_Points(j) = i+1; % This array contains where a particular temperature set ends and new one starts. The last point is added manually.
n = numel(Break_Points); % Number of sets of data to be analysed.
Temperature_Act = zeros(1,n);
j = 1; % Index keeping track of number of stress factors (or temperature sets)
i = 1; % Index keeping track of cycles and fade data in each temperature set
% The following while loop divides data set according to temperature and
% calculates stress factor at each temperature data set.
while j <= numel(Break_Points)
    Temperature_Act(j) = Temperature_Act_Total(i); % Selecting temperature
    Cycles_Act = Cycles_Act_Total(i:Break_Points(j)); % Cycles array for above temperature
    Fade_Act = Fade_Act_Total(i:Break_Points(j)); % Capacity fade for above temperature
    wait_time = 30*60; % seconds
    C_rate = 1; % C
    total_time = Cycles_Act*2*(60*60/C_rate + wait_time); % seconds
    SOC = 0.5; % For calendar aging effects, 50% SOC is assumed
    x_a_curr = Degree_of_Lithiation(SOC,x_SOC_0,x_SOC_100); % Current degree of lithiation
    U_a_curr = OCV_Anode(x_a_curr); % Current OCV of anode
    k_Cal_Cyc_HT_E_a = Stress_Factor_Calendar_Aging(Temperature_Act(j),U_a_curr,k_Cal_Ref,E_a_Cal,R_g,T_ref,alpha,F,U_a_Ref,k0); % Calendar aging stress factor
    Charge_Data = 3*2*Cycles_Act; % 3Ah*Two times charge throughput* No of Cycles
    options_Act = optimoptions('lsqcurvefit','OptimalityTolerance',1e-20,'FunctionTolerance',1e-20);
    Cal_HT_Fade = k_Cal_Cyc_HT_E_a*sqrt(total_time/3600); % Calendar aging capacity fade at j th temperature set
    Cyc_HT = Fade_Act - Cal_HT_Fade; % Cycle aging capacity fade
    l = 1; % Counting index
    for k = i:j*numel(Cyc_HT)
        Q1(k) = Cyc_HT(l); % Q1 represents pure cycle aging capacity fade. Since it is obtained for one temperature at a time
        % in order to consolidate it for all temperatures, this variable is
        % created
        l = l + 1;
    end
    i = Break_Points(j) + 1;
    j = j + 1;
end
T1 = (1./Temperature_Act_Total - 1/T_ref)/R_g; % Linear equation after taking log
Qt1 = log(Qt)/2; % Linear equation after taking log
m = [T1,Qt1]; % inputs
mdl = fitlm(m,log(Q1)); % model
High_T_Ref_Est = exp(mdl.Coefficients.Estimate(1)); % Reference stress factor
High_T_Ea_Est = -1*mdl.Coefficients.Estimate(2); % Activation energy

%% Low Temperature
Cyc_LT_Act_Data = readmatrix('Cycle_Data.xlsx','Sheet','Pre_Est_LT');
Temperature_Act_Total_LT = Cyc_LT_Act_Data(:,1); % temperature in kelvin
Charging_Rate_Total_LT = Cyc_LT_Act_Data(:,2); % A
Cycles_Act_Total_LT = Cyc_LT_Act_Data(:,3); % number of cycles
Fade_Act_Total_LT = Cyc_LT_Act_Data(:,4)/100; % fractional capacity fade
Qch = 3*Cycles_Act_Total_LT; % Charge throughput, 3 because C0 = 3 Ah
Number_of_Rows = numel(Cycles_Act_Total_LT);
j = 1; % Counting index
for i = 1:(Number_of_Rows-1)
    % Following if condition notes where temperature in changing
    if Temperature_Act_Total_LT(i) ~= Temperature_Act_Total_LT(i+1)
        Break_Points_LT(j) = i; % This arrays keeps track of temperature changes
        j = j + 1;
    end
end
Break_Points_LT(j) = i+1; % This array contains where a particular temperature set ends and new one starts. The last point is added manually.
n = numel(Break_Points_LT); % Number of sets of data to be analysed.
Temperature_Act_LT = zeros(1,n);
j = 1; % Index keeping track of number of stress factors (or temperature sets)
i = 1; % Index keeping track of cycles and fade data in each temperature set
% The following while loop divides data set according to temperature and
% calculates stress factor at each temperature data set.
while j <= numel(Break_Points_LT)
    Temperature_Act_LT(j) = Temperature_Act_Total_LT(i);
    Cycles_Act_LT = Cycles_Act_Total_LT(i:Break_Points_LT(j));
    Fade_Act_LT = Fade_Act_Total_LT(i:Break_Points_LT(j));
    C_Rate = Charging_Rate_Total_LT(i:Break_Points_LT(j))/3; % Dividing by 3 because C0 = 3 Ah
    wait_time = 30*60; % seconds
    C_rate = C_Rate(1); % C
    total_time_LT = Cycles_Act_LT*2*(60*60/C_rate + wait_time); % seconds
    SOC = 0.5; % For calendar aging effects, 50% SOC is assumed
    x_a_curr_LT = Degree_of_Lithiation(SOC,x_SOC_0,x_SOC_100); % Current degree of lithiation
    U_a_curr_LT = OCV_Anode(x_a_curr_LT); % Current OCV of anode
    k_Cal_Cyc_LT_E_a = Stress_Factor_Calendar_Aging(Temperature_Act_LT(j),U_a_curr_LT,k_Cal_Ref,E_a_Cal,R_g,T_ref,alpha,F,U_a_Ref,k0); % Calendar aging stress factor
    Charge_Data_LT = 3*Cycles_Act_LT; % 3Ah* No of Cycles
    options_Act = optimoptions('lsqcurvefit','OptimalityTolerance',1e-20,'FunctionTolerance',1e-20);
    Cal_HT_Fade_LT = k_Cal_Cyc_LT_E_a*sqrt(total_time_LT/3600);
    k_Cyc_HT = Stress_Factor_Cycle_Aging_High_Temperature(Temperature_Act_LT(j),R_g,T_ref,High_T_Ref_Est,High_T_Ea_Est);
    Cyc_HT_LT = k_Cyc_HT*sqrt(3*2*Cycles_Act_LT); % Capacity fade due to cycle aging at high temperature
    Cyc_LT = Fade_Act_LT - Cal_HT_Fade_LT - Cyc_HT_LT; % Capacity fade due to lcycle aging at low temperature
    l = 1; % Counting index
    for k = i:j*numel(Cyc_LT)
        Q1_LT(k) = Cyc_LT(l); % Pure low temperature cycle aging capacity fade
        l = l + 1;
    end
    i = Break_Points_LT(j) + 1;
    j = j + 1;
end
T1_LT = (1./Temperature_Act_Total_LT - 1/T_ref)/R_g; % Linear equation after taking log
Qc1_LT = log(Qch)/2; % Linear equation after taking log
Ic1_LT = (Charging_Rate_Total_LT - 3)/3; % Linear equation after taking log
m = [T1_LT,Ic1_LT,Qc1_LT]; % Inputs
mdl_LT = fitlm(m,log(Q1_LT)); % Model
Low_T_Ref_Est = exp(mdl_LT.Coefficients.Estimate(1)); % Reference stress factor
Low_T_Ea_Est = mdl_LT.Coefficients.Estimate(2); % Activation Energy
Low_T_beta_Est = mdl_LT.Coefficients.Estimate(3); % beta

%% Results
High_T_Ref_Est
High_T_Ea_Est
Low_T_Ref_Est
Low_T_Ea_Est
Low_T_beta_Est

%% Functions

% Degree of lithiation (x_a)
function x_a = Degree_of_Lithiation(SOC,x_a_SOC_0,x_a_SOC_100)
    x_a = x_a_SOC_0 + (SOC*(x_a_SOC_100 - x_a_SOC_0)); 
end
% OCV of anode
function U_a = OCV_Anode(x_a)
    U_a = 0.6379+(0.5416*exp(-305.5309*x_a))+(0.044*tanh((0.1958-x_a)/0.1088))-(0.1978*tanh((x_a-1.0571)/0.0854))-(0.6875*tanh((x_a+0.0117)/0.0529))-(0.0175*tanh((x_a-0.5692)/0.0875));
end
% Stress Factor - Calendar Aging
function SF_CA = Stress_Factor_Calendar_Aging(Temperature,U_a,k_Cal_Ref,E_a_Cal,R_g,T_ref,alpha,F,U_a_Ref,k0)
    SF_CA = k_Cal_Ref*(exp(-E_a_Cal*((1/Temperature)-(1/T_ref))/R_g))*(exp(alpha*F*(U_a_Ref-U_a)/(R_g*T_ref)) + k0);
end
% Stress Factor - Cycle Aging - High temperature
function SF_CA_HT = Stress_Factor_Cycle_Aging_High_Temperature(Temperature,R_g,T_ref,k_cyc_h_T_ref,E_a_cyc_h_T)
    SF_CA_HT = k_cyc_h_T_ref*(exp(-E_a_cyc_h_T*((1/Temperature)-(1/T_ref))/R_g));
end