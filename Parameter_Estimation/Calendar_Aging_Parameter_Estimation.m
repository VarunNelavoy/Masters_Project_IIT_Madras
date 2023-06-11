clear;
clc;

%% Constants

R_g = 8.314; % J/mol-K (Universal gas Constant)
T_ref = 298.15; % K (Reference Temperature)
F = 96485; % C/mol (Faradays Constant)
x_SOC_0 = 0.0085; % Degree of Lithiation at 0 % SOC
x_SOC_100 = 0.78; % Degree of Lithiation at 100 % SOC
U_a_Ref = 0.1233; % V (reverence OCV of anode)

%% Calendar Aging Reference Stress Factor Estimation

Cal_Ref_Data = readmatrix('Calendar_Data.xlsx','Sheet','Reference'); % This sheet contains time (days) vs capacity fade (%)
Cal_Ref_Time = Cal_Ref_Data(:,1)*24; % time in hours
Cal_Ref_Fade = Cal_Ref_Data(:,2)/100; % fractional capacity fade
Cal_Ref_Func = @(k_Cal_Ref,Time_Data)k_Cal_Ref*sqrt(Time_Data); % Function relating capacity fade and time for calendar aging
k_Cal_Ref0 = 1; % Initial guess for k_Cal_Ref
k_Cal_Ref = lsqcurvefit(Cal_Ref_Func,k_Cal_Ref0,Cal_Ref_Time,Cal_Ref_Fade); % h^(-0.5)

%% Calender Aging Activation Energy Estimation

% At each temperature, we have time vs fade data each of which set has its
% own stress factor. First data (time and fade) has to be split according
% to temperature and stress factor has to be calculated for each set. These
% stress factors are then used to calculate activation energy by fitting it
% to arrhenius model.
Cal_Act_Data = readmatrix('Calendar_Data.xlsx','Sheet','ActivationEnergy'); % This sheet contains temperature (K), time (days) and capacity fade (%).
Temperature_Act_Total = Cal_Act_Data(:,1); % temperature in kelvin
Time_Act_Total = Cal_Act_Data(:,2)*24; % time in hours
Fade_Act_Total = Cal_Act_Data(:,3)/100; % fractional capacity fade
Number_of_Rows = numel(Time_Act_Total);
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
k_cal = zeros(1,n);
j = 1; % Index keeping track of number of stress factors (or temperature sets)
i = 1; % Index keeping track of time and fade data in each temperature set
% The following while loop divides data set according to temperature and
% calculates stress factor at each temperature data set.
k_Cal_0 = 0; % Initial guess for k_Cal
while j <= numel(Break_Points)
    Temperature_Act(j) = Temperature_Act_Total(i);
    Time_Act = Time_Act_Total(i:Break_Points(j));
    Fade_Act = Fade_Act_Total(i:Break_Points(j));
    Cal_Func = @(k_Cal,Time_Data)k_Cal*sqrt(Time_Data);
    k_Cal(j) = lsqcurvefit(Cal_Func,k_Cal_0,Time_Act,Fade_Act); % h^(-0.5)
    i = Break_Points(j) + 1;
    j = j + 1;
end
Act_Func = @(E_a_Cal,T_data)k_Cal_Ref*exp(-E_a_Cal*(1./T_data - 1/T_ref)/R_g);
options_Act = optimoptions('lsqcurvefit','OptimalityTolerance',1e-12,'FunctionTolerance',1e-12);
LB_Act = []; % Lower bound for activation energy
UB_Act = []; % Upper bound for activation energy
E_a_0 = 10000; % Initial guess for activation energy
E_a_Cal = lsqcurvefit(Act_Func,E_a_0,Temperature_Act,k_Cal,LB_Act,UB_Act,options_Act);

%% Calendar Aging alpha and k0 estimation

% At each SOC, we have time vs fade data each of which set has its
% own stress factor. First data (time and fade) has to be split according
% to SOC and stress factor has to be calculated for each set. These
% stress factors are then used to calculate alpha and k0 by fitting it
% to modified tafel equation.
Cal_Other_Data = readmatrix('Calendar_Data.xlsx','Sheet','OtherParameters'); % This sheet contains SOC (%), time (days) and fade (%).
SOC_Other_Total = Cal_Other_Data(:,1)/100; % Fractional SOC
Time_Other_Total = Cal_Other_Data(:,2)*24; % time in hours
Fade_Other_Total = Cal_Other_Data(:,3)/100; % Fractional fade
Number_of_Rows = numel(SOC_Other_Total);
j = 1; % Index to count number of distinct SOC data sets
for i = 1:(Number_of_Rows-1)
    if SOC_Other_Total(i) ~= SOC_Other_Total(i+1)
        Break_Points_Other(j) = i;
        j = j + 1;
    end
end
Break_Points_Other(j) = i+1; % Similar to estimation of activation energy, adding last set index manually
n = numel(Break_Points_Other);
SOC_Other = zeros(1,n);
k_cal_other = zeros(1,n);
j = 1; % Index of distinct SOC
i = 1; % Index in original data set
k_Cal_Other0 = 0; % Initial guess for k_Cal
while j <= n
    SOC_Other(j) = SOC_Other_Total(i);
    Time_Other = Time_Other_Total(i:Break_Points_Other(j));
    Fade_Other = Fade_Other_Total(i:Break_Points_Other(j));
    Cal_Func = @(k_Cal,Time_Data)k_Cal*sqrt(Time_Data);    
    k_cal_other(j) = lsqcurvefit(Cal_Func,k_Cal_Other0,Time_Other,Fade_Other); % h^(-0.5)
    i = Break_Points_Other(j) + 1;
    j = j + 1;
end
x_a = Degree_of_Lithiation(SOC_Other,x_SOC_0,x_SOC_100);
U_a = OCV_Anode(x_a);
p0 = [1 1]; % Initial guess for alpha and k0
LB_other = [];
UB_other = [];
options = optimoptions('lsqcurvefit','OptimalityTolerance',1e-12,'FunctionTolerance',1e-12);
k_Func = @(p,x_data)k_Cal_Ref*(exp(p(1)*F*(U_a_Ref - x_data)/(R_g*T_ref)) + p(2)); % Fitting modified tafel equation
p = lsqcurvefit(k_Func,p0,U_a,k_cal_other,LB_other,UB_other,options); % First item in this array is alpha and second is k0

%% Results
k_Cal_Ref
E_a_Cal
alpha = p(1)
k0 = p(2)

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