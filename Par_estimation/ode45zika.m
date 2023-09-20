function [cumnewcases, newcasesbyweek, x, t] = ode45zika(pars,country);
%% Parameter values and definitions
mu=1/(52*80); %Birth and death rate
%b1=2.85; %mosquito human contact rate
%a1=0.63; %sexual transmission rate
%k =1.19 ; %incubation period 
% a=0.0;b=0.; % upper and lower bounds for u1 control (condom)
% c=0.0;d=0.;  %upper and lower bounds for u2 control (repellent)
% e=0.0; f=0.0;  %upper and lower bounds for phi control (vaccine)
%phi=0.4; %vaccination rate
omega=0.; % wane of immunity 
a2=1.40; %recovery rate
%b2=1.45; %human mosquito contact rate
c1= 0.4; %birth and death mosquito rate
%c2=0.67; %rate at which exposed mosquito become infectious
% A1=1; %Infectious weight parameter
% B1=0.00001; %u1 nonlinear weight parameter
% B2=0.001; %u2 nolinear weight parameter
% B3=0.1; %phi nonlinear weight parameter
% B4=0; %u1 linear weight parameter
% B5=0; %u2 linear weight parameter
% B6=0; %phi linear weight parameter
%%
fixpar=[mu omega a2 c1];
%% Parameters that are being estimated
b1=pars(1); %mosquito human contact rate
b2=pars(2); %human mosquito contact rate
a1=pars(3); %sexual transmission rate
k=pars(4);  %incubation period 
c2=pars(5); %rate at which exposed mosquito become infectious
r=pars(6); %reporting rate
%S0=pars(7); %Initial susceptible pop.
%% Time Domain and IC for each country
% Time horizon and Initial conditions for PR
if country == 'PR'
    T=104; %weeks
    IC= [2800000 0 ceil(46/r)/3 ceil(46/r) 3000000*2 0 ceil((46*2)/r)];
end
% Time horizon and Initial conditions for Colombia
if country == 'CO'
    T=63; %Weeks for CO
    IC= [35000000 0 ceil(ceil(389/r)/3) ceil(389/r) 90000000 0 ceil(780/r)]; %All CO is susceptible
end
%%
S0=IC(1);
V0=IC(2);
E0=IC(3);
I0=IC(4);
Sv0=IC(5); 
Ev0=IC(6);
Iv0=IC(7);
R0=0;
%% Human and Mosquito population are constant
 N=S0+V0+E0+I0+R0;
 Nv=Sv0+Ev0+Iv0;
%% Setting up the code with ODE45 %%%%%%%%%
 
M = T-1;
tvec=linspace(1,T,M+1)';

x=zeros(M+1,7);
u1=zeros(M+1,1);
u2=zeros(M+1,1);
phi=zeros(M+1,1);

    [t,x] = ode45(@(t,x) zikastates(t,x,tvec,u1,u2,phi,pars,fixpar,IC),tvec,IC);
    
t=tvec;
S=x(:,1);
V=x(:,2);
E=x(:,3);
I=x(:,4);
Sv=x(:,5) ;
Ev=x(:,6) ;
Iv=x(:,7);
R=N-(S+V+E+I);
%% Total number of new infected cases
ratenewcases = b1*S.*Iv.*(1-u2)/Nv  + a1*S.*I.*(1-u1)/N ;
totalnumberinfectedcases=trapz(t,ratenewcases);
cumnewcases=zeros(M+1,1);
newcasesbyweek=zeros(M+1,1);
cumnewcases(1)=I0;
for i=2:M+1
    cumnewcases(i)=trapz(t(1:i),ratenewcases(1:i));
end
x(:,8)=cumnewcases;
for i=1:M
    newcasesbyweek(1)=cumnewcases(1);
    newcasesbyweek(i+1)=cumnewcases(i+1)-cumnewcases(i);
end
x(:,9)=newcasesbyweek;
Table=table(transpose(t)',transpose(ratenewcases)',transpose(cumnewcases)',...
    transpose(newcasesbyweek)','VariableNames',{'time','ratenewcases',...
    'cumulativenewcases','newcasesbyweek'});
end