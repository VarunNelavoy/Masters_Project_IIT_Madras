clear;
clc;

%% Initial guesses taken from linear regression in pre-estimation
p0_Main_Est = [-8.4729 10.1256 -8.1494 11.0368 2.9324]; % Please note that first four parameters should be after applying natural log. Last parameter does not require any change.

%% Constants
R_g = 8.314; % J/mol-K (Universal gas Constant)
T_ref = 298.15; % K (Reference Temperature)
F = 96485; % C/mol (Faradays Constant)
x_SOC_0 = 0.0085; % Degree of Lithiation at 0 % SOC
x_SOC_100 = 0.78; % Degree of Lithiation at 100 % SOC
U_a_Ref = 0.1233; % V (reverence OCV of anode)

%% Previous results - Calendar Aging
k_Cal_Ref = 3.6940e-04;
E_a_Cal = 2.0493e+04;
alpha = 0.3840;
k0 = 0.142;

%% Calculations
Main_Est_Data = readmatrix('Cycle_Data.xlsx','Sheet','Main_Est');
Temperature_Main_Est = Main_Est_Data(:,1); % temperature in kelvin
Charging_Rate_Main_Est = Main_Est_Data(:,2);
Cycles_Main_Est = Main_Est_Data(:,3); % number of cycles
Fade_Main_Est = Main_Est_Data(:,4)/100; % fractional capacity fade
Qt = 3*2*Cycles_Main_Est;
Qch = 3*Cycles_Main_Est;
Number_of_Rows = numel(Cycles_Main_Est);
j = 1; % Counting index
for i = 1:(Number_of_Rows-1)
    % Following if condition notes where temperature in changing
    if Temperature_Main_Est(i) ~= Temperature_Main_Est(i+1)
        Break_Points(j) = i; % This arrays keeps track of temperature changes
        j = j + 1;
    end
end
Break_Points(j) = i+1; % This array contains where a particular temperature set ends and new one starts. The last point is added manually.
n = numel(Break_Points); % Number of sets of data to be analysed.
Temperature_Main_Est_Sub = zeros(1,n);
k_Cyc = zeros(1,n);
j = 1; % Index keeping track of number of stress factors (or temperature sets)
i = 1; % Index keeping track of cycles and fade data in each temperature set
% The following while loop divides data set according to temperature and
% calculates stress factor at each temperature data set.
k_Cyc_0 = 0; % Initial guess for k_Cal
l = 1;
while j <= numel(Break_Points)
    Temperature_Main_Est_Sub(j) = Temperature_Main_Est(i);
    Cycles_Main_Est_Sub = Cycles_Main_Est(i:Break_Points(j));
    Fade_Main_Est_Sub = Fade_Main_Est(i:Break_Points(j));
    wait_time = 30*60; % seconds
    I_ch = Charging_Rate_Main_Est(i); % A
    C0 = 3; % Ah
    C_rate = I_ch/C0; % C
    total_time_Main_Est = Cycles_Main_Est_Sub*2*(60*60/C_rate + wait_time); % seconds
    SOC = 0.5; % For calendar aging effects, 50% SOC is assumed
    x_a_curr_Main_Est = Degree_of_Lithiation(SOC,x_SOC_0,x_SOC_100); % Current degree of lithiation
    U_a_curr_Main_Est = OCV_Anode(x_a_curr_Main_Est); % Current OCV of anode
    k_Cal_Cyc_Main_Est = Stress_Factor_Calendar_Aging(Temperature_Main_Est_Sub(j),U_a_curr_Main_Est,k_Cal_Ref,E_a_Cal,R_g,T_ref,alpha,F,U_a_Ref,k0); % Calendar aging stress factor
    Cal_Main_Est_Fade = k_Cal_Cyc_Main_Est*sqrt(total_time_Main_Est/3600);
    Cyc_Main_Est = Fade_Main_Est_Sub - Cal_Main_Est_Fade;
    for k = 1:numel(Cyc_Main_Est)
        Q1(l) = Cyc_Main_Est(k);
        l = l + 1;
    end
    i = Break_Points(j) + 1;
    j = j + 1;
end
x_Main_Est_Data = [Temperature_Main_Est,Qt,Charging_Rate_Main_Est,Qch]; % Inputs
p_Main_Est = lsqcurvefit(@Cyc,p0_Main_Est,x_Main_Est_Data,Q1',[-10 1 -10 1 0],[]); % Estimation

%% Results
k_cyc_h_T_ref = exp(p_Main_Est(1))
E_a_cyc_h_T = exp(p_Main_Est(2))
k_cyc_l_T_ref = exp(p_Main_Est(3))
E_a_cyc_l_T = exp(p_Main_Est(4))
beta_l_T = p_Main_Est(5)

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
% Cycle aging fade function
function f = Cyc(p,x)
    f = (exp(p(1)).*exp(-exp(p(2)).*(1./x(:,1) - 1/298.15)/8.314).*sqrt(x(:,2))) + (exp(p(3)).*exp(exp(p(4)).*(1./x(:,1) - 1/298.15)/8.314).*exp(p(5).*(x(:,3) - 3)/3).*sqrt(x(:,4)));
end