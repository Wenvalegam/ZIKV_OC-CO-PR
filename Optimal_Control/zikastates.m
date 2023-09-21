function [ dx ] = zikastates(t,x,tvec,u1,u2,phi,pars,fixpar,IC)
%% Parameter definitions 
%%%% Parameteres based on Dr. Rios student%%%%% 
% mu=1/(52*80); %Birth and death rate
% %b1=2.85; %mosquito human contact rate
% %a1=0.63; %sexual transmission rate
% %k =1.19 ; %incubation period 
% eff=1; %Vaccine efficacy
% % a=0;b=0; % upper and lower bounds for u1 control (condom)
% % c=0;d=0;  %upper and lower bounds for u2 control (repellent)
% % e=0; f=0;  %upper and lower bounds for phi control (vaccine)
% %phi=0.4; %vaccination rate
% omega=0.; % wane of immunity 
% a2=1.40; %recovery rate
% %b2=1.45; %human mosquito contact rate
%c1= 0.90; %birth and death mosquito rate
%c2=0.67; %rate at which exposed mosquito become infectious
% A1=1; %Infectious weight parameter
% B1=0.00001; %u1 nonlinear weight parameter
% B2=0.001; %u2 nolinear weight parameter
% B3=0.1; %phi nonlinear weight parameter
% B4=0; %u1 linear weight parameter
% B5=0; %u2 linear weight parameter
% B6=0; %phi linear weight parameter
%% fixed parmateres
mu=fixpar(1);
omega=fixpar(2);
a2=fixpar(3);
c1=fixpar(4);
%% parameters in the state system
b1=pars(1); %mosquito human contact rate
b2=pars(2); %human mosquito contact rate
a1=pars(3); %sexual transmission rate
k=pars(4);  %incubation period 
c2=pars(5); %rate at which exposed mosquito become infectious
%pars=[ b1  b2   a1   k    c2];
%% Initial conditions 
S0=IC(1);
V0=IC(2);
E0=IC(3);
I0=IC(4);
Sv0=IC(5); 
Ev0=IC(6);
Iv0=IC(7);
R0=0;
%% Total human and mosquito pop are constant
N=S0+V0+E0+I0+R0;
Nv=Sv0+Ev0+Iv0;
%% Interpolating controls
u1=pchip(tvec,u1,t);
u2=pchip(tvec,u2,t);
phi=pchip(tvec,phi,t);
%% State ODE system 
dx=zeros(7,1);
    dx(1)= mu*N - b1*(x(1)*x(7)*(1-u2))/Nv - a1*(x(1)*x(4)*(1-u1))/N + omega*x(2) - (mu+phi)*x(1);
    dx(2) = phi*x(1)-(omega+mu)*x(2);
    dx(3)=  b1*x(1)*x(7)*(1-u2)/Nv  + a1*x(1)*x(4)*(1-u1)/N - k*x(3)-mu*x(3);
    dx(4) = k*x(3) - a2*x(4)-mu*x(4);
    dx(5) = c1*(Nv-x(5)) - b2*x(5)*x(4)/N;
    dx(6) = b2*x(5)*x(4)/N - (c1+c2)*x(6);
    dx(7)= c2*x(6) - c1*x(7);

end
