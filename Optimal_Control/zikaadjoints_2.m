function [ dlambda ] = zikaadjoints_2(t,lambda,tvec,x,u1,u2,phi,pars,fixpar,weights,IC)
%% Parameter definitions
%%%% Parameteres based on Dr. Rios student%%%%% 
% mu=1/(52*80); %Birth and death rate
% %b1=2.85; %mosquito human contact rate
% %a1=0.63; %sexual transmission rate
% %k =1.19 ; %incubation period 
% % a=0;b=0; % upper and lower bounds for u1 control (condom)
% % c=0;d=0;  %upper and lower bounds for u2 control (repellent)
% % e=0; f=0;  %upper and lower bounds for phi control (vaccine)
% %phi=0.4; %vaccination rate
% omega=0.; % wane of immunity 
% a2=1.40; %recovery rate
% %b2=1.45; %human mosquito contact rate
%c1= 0.90; %birth and death mosquito rate
%c2=0.67; %rate at which exposed mosquito become infectious
A1=weights(1); %Infectious weight parameter
A2=weights(8); %New cases weight parameter
B1=weights(2); %u1 nonlinear weight parameter
B2=weights(3); %u2 nolinear weight parameter
B3=weights(4); %phi nonlinear weight parameter
B4=weights(5); %u1 linear weight parameter
B5=weights(6); %u2 linear weight parameter
B6=weights(7); %phi linear weight parameter
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
%pars=[ 2.85  1.45  0.63   1.19    0.67];
%pars=[ b1  b2   a1   k    c2];
%% Initial conditions (note that transversality cond. are input in the main code ode45zikaOC)
S0=IC(1);
V0=IC(2);
E0=IC(3);
I0=IC(4);
Sv0=IC(5); 
Ev0=IC(6);
Iv0=IC(7);
R0=0;
%% Total human and mosquito population are constant
N=S0+V0+E0+I0+R0;
Nv=Sv0+Ev0+Iv0;
%% Interpolating controls and states (not sure how it works yet)
x=interp1(tvec,x,t);
u1=pchip(tvec,u1,t);
u2=pchip(tvec,u2,t);
phi=pchip(tvec,phi,t);
%% Adjoint ODE system
dlambda=zeros(7,1);
    dlambda(1) = (lambda(1)-lambda(3)-A2)*(b1*x(7)*(1-u2)/Nv + a1*x(4)*(1-u1)/N) +...
                 (lambda(1)-lambda(2))*phi + mu*lambda(1)-(B6*phi+B4*u1+B5*u2);
    dlambda(2) = (lambda(2)-lambda(1))*omega + lambda(2)*mu;
    dlambda(3) = (lambda(3)-lambda(4))*k + mu*lambda(3);
    dlambda(4) = a1*x(1)*(1-u1)*(lambda(1)-lambda(3)-A2)/N + lambda(4)*(a2+mu) + ...
                (lambda(5)-lambda(6))*b2*x(5)/N ;
    dlambda(5) = c1*lambda(5) + ((lambda(5)-lambda(6))*b2*x(4))/N;
    dlambda(6) = c2*(lambda(6)-lambda(7)) + c1*lambda(6);
    dlambda(7) = (lambda(1)-lambda(3)-A2)*b1*x(1)*(1-u2)/Nv + c1*lambda(7);
end
