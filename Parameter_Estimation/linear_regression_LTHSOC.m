clear;
clc;

%% Constants
R_g = 8.314; % J/mol-K (Universal gas Constant)
T_ref = 298.15; % K (Reference Temperature)
F = 96485; % C/mol (Faradays Constant)
x_SOC_0 = 0.0085; % Degree of Lithiation at 0 % SOC
x_SOC_100 = 0.78; % Degree of Lithiation at 100 % SOC
U_a_Ref = 0.1233; % V (reverence OCV of anode)

%% Previous results
k_Cal_Ref = 3.6940e-04;
E_a_Cal = 2.0493e+04;
alpha = 0.3840;
k0 = 0.142;
k_cyc_HT_ref = 1.9447e-04;
E_a_HT = 1.4958e+04;
k_cyc_LT_ref = 3.2208e-04;
E_a_LT = 6.2826e+04;
beta_LT = 3.1831;

%% Low temperature and high SOC
Cyc_LTHSOC_Act_Data = readmatrix('Cycle_Data.xlsx','Sheet','Pre_Est_LT_HSOC');
Temperature_Act_Total_LTHSOC = Cyc_LTHSOC_Act_Data(:,1); % temperature in kelvin
Charging_Rate_Total_LTHSOC = Cyc_LTHSOC_Act_Data(:,2);
Cycles_Act_Total_LTHSOC = Cyc_LTHSOC_Act_Data(:,3); % number of cycles
Fade_Act_Total_LTHSOC = Cyc_LTHSOC_Act_Data(:,4)/100; % fractional capacity fade
Qch_HSOC = 3*Cycles_Act_Total_LTHSOC;
Number_of_Rows = numel(Cycles_Act_Total_LTHSOC);
j = 1; % Counting index
for i = 1:(Number_of_Rows-1)
    % Following if condition notes where temperature in changing
    if Temperature_Act_Total_LTHSOC(i) ~= Temperature_Act_Total_LTHSOC(i+1)
        Break_Points_LTHSOC(j) = i; % This arrays keeps track of temperature changes
        j = j + 1;
    end
end
Break_Points_LTHSOC(j) = i+1; % This array contains where a particular temperature set ends and new one starts. The last point is added manually.
n = numel(Break_Points_LTHSOC); % Number of sets of data to be analysed.
Temperature_Act_LTHSOC = zeros(1,n);
k_Cyc_LTHSOC = zeros(1,n);
j = 1; % Index keeping track of number of stress factors (or temperature sets)
i = 1; % Index keeping track of cycles and fade data in each temperature set
% The following while loop divides data set according to temperature and
% calculates stress factor at each temperature data set.
k_Cyc_0 = 0; % Initial guess for k_Cal
count = 0;
while j <= n
    Temperature_Act_LTHSOC(j) = Temperature_Act_Total_LTHSOC(i);
    Cycles_Act_LTHSOC = Cycles_Act_Total_LTHSOC(i:Break_Points_LTHSOC(j));
    Fade_Act_LTHSOC = Fade_Act_Total_LTHSOC(i:Break_Points_LTHSOC(j));
    wait_time = 30*60; % seconds
    C_Rate = Charging_Rate_Total_LTHSOC(i:Break_Points_LTHSOC(j))/3;
    C_rate = C_Rate(1); % C
    total_time_LTHSOC = Cycles_Act_LTHSOC*2*(60*60/C_rate + wait_time); % seconds
    SOC = 0.5; % For calendar aging effects, 50% SOC is assumed
    x_a_curr_LT = Degree_of_Lithiation(SOC,x_SOC_0,x_SOC_100); % Current degree of lithiation
    U_a_curr_LT = OCV_Anode(x_a_curr_LT); % Current OCV of anode
    k_Cal = Stress_Factor_Calendar_Aging(Temperature_Act_LTHSOC(j),U_a_curr_LT,k_Cal_Ref,E_a_Cal,R_g,T_ref,alpha,F,U_a_Ref,k0); % Calendar aging stress factor
    Charge_Data_LTHSOC = 3*Cycles_Act_LTHSOC; % 3Ah* No of Cycles
    options_Act = optimoptions('lsqcurvefit','OptimalityTolerance',1e-20,'FunctionTolerance',1e-20);
    Cal_Fade = k_Cal*sqrt(total_time_LTHSOC/3600);
    k_Cyc_HT_HSOC = Stress_Factor_Cycle_Aging_High_Temperature(Temperature_Act_LTHSOC(j),R_g,T_ref,k_cyc_HT_ref,E_a_HT);
    Cyc_HT_LT_HSOC = k_Cyc_HT_HSOC*sqrt(3*2*Cycles_Act_LTHSOC);
    k_Cyc_LT_HSOC = Stress_Factor_Cycle_Aging_Low_Temperature(Temperature_Act_LTHSOC(j),Charging_Rate_Total_LTHSOC(i),R_g,T_ref,k_cyc_LT_ref,E_a_LT,3,3,beta_LT);
    Cyc_LT_HSOC = k_Cyc_LT_HSOC*sqrt(3*Cycles_Act_LTHSOC);
    Cyc_LTHSOC = Fade_Act_LTHSOC - Cal_Fade - Cyc_HT_LT_HSOC - Cyc_LT_HSOC;
    l = 1;
    for k = i:j*numel(Cyc_LTHSOC)
        Q1_LTHSOC(k) = Cyc_LTHSOC(l);
        l = l + 1;
    end
    i = Break_Points_LTHSOC(j) + 1;
    j = j + 1;
end
n1 = numel(Q1_LTHSOC);
j = 1;
for i = 1:n1
    if Q1_LTHSOC(i) > 0
        Qnew(j) = Q1_LTHSOC(i);
        Tnew(j) = (1/Temperature_Act_Total_LTHSOC(i) - 1/T_ref)/R_g;
        Qchnew(j) = log(Qch_HSOC(i));
        Inew(j) = (Charging_Rate_Total_LTHSOC(i) - 3)/3;
        j = j + 1;
    end
end
mHSOC = [Tnew;Inew;Qchnew]';
q = log(Qnew);
mdl_LTHSOC = fitlm(mHSOC,q);
Low_T_HSOC_Ref_Est = exp(mdl_LTHSOC.Coefficients.Estimate(1))
Low_T_HSOC_Ea_Est = mdl_LTHSOC.Coefficients.Estimate(2)
Low_T_HSOC_beta_Est = mdl_LTHSOC.Coefficients.Estimate(3)

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
% Stress Factor - Cycle Aging - Low Temperature
function SF_CA_LT = Stress_Factor_Cycle_Aging_Low_Temperature(Temperature,I_ch,R_g,T_ref,k_cyc_l_T_ref,E_a_cyc_l_T,C0,I_ch_ref,beta_l_T)
    SF_CA_LT = k_cyc_l_T_ref*(exp(E_a_cyc_l_T*((1/Temperature)-(1/T_ref))/R_g))*(exp(beta_l_T*(I_ch-I_ch_ref)/C0));
end