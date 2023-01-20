% -----------------------
% Box Model to Check Calculations of Changes in Oxygen under Warming
% Written by: Andreas Andersson and Ariel Pezner
% Last updated: 11 January 2023
% -----------------------

clear all
close all
% global variables
global oceanO2 reefvol tau resp
  
meanO2=[]; dresp=[]; O2data=[]; n=1;
for i=1:4:5                 %sensitivity runs changing residence time 1 and 5 hours
    tau=i;                  
    m=1;
    for j=0:3:6             %sensitivity runs changing respiration rate based on Q10 at 25 C, 28 C, 31 C
        %Set up initial variables
        t0 = 0;                     % time (start)
        tfinal = 168;               % time hours (end)
        tspan=t0:0.5:tfinal;        % solutions for every 0.5 h
        temp = 25+j;                  % Temperature (°C)
        sal = 35;                   % Salinity
        o2_sol = oxy_sol(temp,sal); % O2 Equilibrium concentration (mumol kg-1)

        %Reef/mesocosm dimensions 1x1x1m
        reefvol= 1000;              % Volume L
        reefo2=o2_sol;              % reef assuming equilibrium; 
        oceanO2=o2_sol;             % ocean assuming equilibrium; 

        F0=reefo2*1e-06 * reefvol;  %Total O2 mol and initial reservoir value
        
        % Q10 calculations
        Q10 = 2; % Q10
        R1=12*1e-03;        %respiration rate in base scenario 12 mmol/m2/h
        resp = R1*Q10.^((j)/10);
        Rper=(resp-R1)/R1*100;   %percent change relative to R1
        
        % solve differential equations.
        % use matlab routine ode45 (Runge-Kutta)
        % with function 'reefo2dif' containing the difeq.
        [tv,F] = ode45('reefo2dif',[tspan],F0);
        O2=F./reefvol.*1e+06;
        
        %Plot data
        subplot(1,2,1)
        plot(tv, O2); hold on 
        legend('1 hr | 25ºC','1 hr | 28ºC','1 hr | 31ºC','5 hr | 25ºC','5 hr | 28ºC','5 hr | 31ºC')
        %set(gca,'XTick',0:24:168)
        ylabel('Dissolved Oxygen (µmol kg^{-1})')
        xlabel('Hour')
       
        %Calculate mean O2 (meanO2), %change in R, and save each run (O2data)
        meanO2(m,n)=mean(O2);
        dresp(m,n)=Rper;
        x=size(O2data,2); O2data(:,x+1)=O2;
        m=m+1;
    end
    n=n+1;
    
    
end

%% Compare box model output to calculation approach

%tmp=(reefo2-min(O2data(:,1)))*(dresp(2,1)/100+1);
%tmp1=(reefo2-min(O2data(:,1)))-tmp;

%Residence time scenario 1(1h) +3C
tmp2=(mean(O2data(:,1))-min(O2data(:,1)))*(1-(dresp(2,1)/100+1)); %combine tmp and tmp1 eqs
tmp3=O2data(:,1)+tmp2-(oxy_sol(25,sal)-oxy_sol(28,sal)); % also subtract decrease from solubility change at temp+3
tmp4=O2data(:,2)-tmp3;

%Residence time scenario 1(1h) +6C
tmp2b=(mean(O2data(:,1))-min(O2data(:,1)))*(1-(dresp(3,1)/100+1)); 
tmp3b=O2data(:,1)+tmp2b-(oxy_sol(25,sal)-oxy_sol(31,sal));
tmp4b=O2data(:,3)-tmp3b;


%Residence time scenario 2(5h) +3C
tmp2c=(mean(O2data(:,4))-min(O2data(:,4)))*(1-(dresp(2,2)/100+1)); 
tmp3c=O2data(:,4)+tmp2c-(oxy_sol(25,sal)-oxy_sol(28,sal));
tmp4c=O2data(:,5)-tmp3c;

%Residence time scenario 2(5h) +6C
tmp2d=(mean(O2data(:,4))-min(O2data(:,4)))*(1-(dresp(3,2)/100+1)); 
tmp3d=O2data(:,4)+tmp2d-(oxy_sol(25,sal)-oxy_sol(31,sal));
tmp4d=O2data(:,6)-tmp3d;
    
% Make a plot with each scenario
subplot(1,2,2)
plot(tv,tmp4,'linewidth',2)
hold on
plot(tv,tmp4b,'linewidth',2)
plot(tv,tmp4c,'linewidth',2)
plot(tv,tmp4d,'linewidth',2)
legend('1hr,+3ºC','1hr,+6ºC','5hr,+3ºC','5hr,+6ºC')
ylabel('Deviation from Calculations (µmol kg^{-1})')
ylim([-10,30])

