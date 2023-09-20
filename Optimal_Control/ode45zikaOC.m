function [I,tvec,newcasesbyday,totalnumberinfectedcases,u1,u2,phi,J,J1,...
    J2,totalvac] = ode45zikaOC(pars,b,d,f,delay,country,objective) 
%% Parameter values and definitions
mu=1/(52*80); %Birth and death rate
%b1 ; %mosquito human contact rate
%a1 ; %sexual transmission rate
%k ; %incubation period 
eff=0.9; %Vaccine efficacy
a=0.0; %b=; % upper and lower bounds for u1 control (condom)
c=0.0; %d=;  %upper and lower bounds for u2 control (repellent)
e=0.0; %f=;  %upper and lower bounds for phi control (vaccine)
%phi=0.4; %vaccination rate
omega=0; % wane of immunity 
a2=1.40; %recovery rate
%b2; %human mosquito contact rate
c1= 0.4; %birth and death mosquito rate
%c2; %rate at which exposed mosquito become infectious
A1=10; %Infectious weight parameter
A2=10; %weight parameter associated with newcases
% B1; %u1 nonlinear weight parameter
% B2; %u2 nolinear weight parameter
% B3; %phi nonlinear weight parameter
% B4; %u1 linear weight parameter
% B5; %u2 linear weight parameter
% B6; %phi linear weight parameter
%%
fixpar=[mu omega a2 c1];
%% Parameters that are being estimated
b1=pars(1); %mosquito human contact rate
b2=pars(2); %human mosquito contact rate
a1=pars(3); %sexual transmission rate
k=pars(4);  %incubation period 
c2=pars(5); %rate at which exposed mosquito become infectious
r=pars(6); %reporting rate
%% Time Domain and IC for each country
% Time horizon and Initial conditions for PR
if country == 'PR'
    T=104; %weeks
    IC= [2800000 0 ceil(46/r)/3 ceil(46/r) 3000000*2 0 ceil((46*2)/r)];
    B1=1; %u1 nonlinear weight parameter
    B2=2; %u2 nolinear weight parameter
    B3=5; %phi nonlinear weight parameter
    B4=0.1; %u1 linear weight parameter
    B5=0.1; %u2 linear weight parameter
    B6=0.1; %phi linear weight parameter
end
% Time horizon and Initial conditions for Colombia
if country == 'CO'
    T=63; %Weeks for CO
    IC= [35000000 0 ceil(ceil(389/r)/3) ceil(389/r) 90000000 0 ceil(780/r)]; %All CO is susceptible
    B1=3; %u1 nonlinear weight parameter
    B2=5; %u2 nolinear weight parameter
    B3=1; %phi nonlinear weight parameter
    B4=0.01; %u1 linear weight parameter
    B5=0.01; %u2 linear weight parameter
    B6=0.1; %phi linear weight parameter
end
weights=[A1 B1 B2 B3 B4 B5 B6 A2];
%%
S0=IC(1);
V0=IC(2);
E0=IC(3);
I0=IC(4);
Sv0=IC(5); 
Ev0=IC(6);
Iv0=IC(7);
R0=0;

%% Forward-Backward sweep method using ODE45
%%% Setting up the code %%

N=S0+V0+E0+I0+R0; % Total human pop is constant
Nv=Sv0+Ev0+Iv0;   % Total mosquito pop is constant
 
count=0;
test = -1;


delta = 0.01;
M = T;
tvec=linspace(0,T,T+1)'; %time vector

x=zeros(M+1,7); %zero vector to save state solutions 

%%% zero vectors to save adjoint and control solutions
lambda=zeros(M+1,7);
u1=zeros(M+1,1);
u2=zeros(M+1,1);
phi=zeros(M+1,1);
R=zeros(M+1,1);

while(test < 0)
   
 %%%count number of iteration the code need to converge   
 count = count+1;

 %%% updating solutions
    oldu1 = u1;
    oldu2 = u2;
    oldphi = phi;
    oldx = x;
    oldlambda = lambda;

%%%Solving states system forward and adjoint system backward

    %%% Need zikastates and zikaadjoint/zikaadjoint_2 m-files%%% 
    solx = ode45(@(t,x) zikastates(t,x,tvec,u1,u2,phi,pars,fixpar,IC),tvec,IC);
    x = deval(solx,tvec)';
if objective == 1

    sollamb = ode45(@(t,lambda) zikaadjoints(t,lambda,tvec,x,u1,u2,phi,pars,fixpar,weights,IC),...
                                                [T 0],[0 0 0 0 0 0 0]);
    lambda = deval(sollamb,tvec)';
    
    S=x(:,1);
    I=x(:,4);
    Iv=x(:,7);
    lambdaS=lambda(:,1);
    lambdaV=lambda(:,2);
    lambdaE=lambda(:,3);

%%% Control characterization     
    tempu1 = ((lambdaE-lambdaS).*a1.*S.*I - B4.*S*N)./(2*B1*N);
    u11=min(b,max(a,tempu1));
    u1=0.5*(u11+oldu1);
    
    tempu2 = ((lambdaE-lambdaS).*b1.*S.*Iv - B5.*S*Nv)./(2*B2*Nv);
    u22=min(d,max(c,tempu2));
    u2=0.5*(u22+oldu2);
    
    tempphi =((lambdaS-lambdaV).*S- B6.*S)/(2*B3); 
    phi1=min(f,max(e,tempphi));
    phi=0.5*(phi1+oldphi);
end
if objective == 2
    sollamb = ode45(@(t,lambda) zikaadjoints_2(t,lambda,tvec,x,u1,u2,phi,pars,fixpar,weights,IC),...
                                            [T 0],[0 0 0 0 0 0 0]);
    lambda = deval(sollamb,tvec)';
    
    S=x(:,1);
    I=x(:,4);
    Iv=x(:,7);
    lambdaS=lambda(:,1);
    lambdaV=lambda(:,2);
    lambdaE=lambda(:,3);

%%% Control characterization     
    tempu1 = ((A2+lambdaE-lambdaS).*a1.*S.*I - B4.*S*N)./(2*B1*N);
    u11=min(b,max(a,tempu1));
    u1=0.5*(u11+oldu1);
    
    tempu2 = ((A2+lambdaE-lambdaS).*b1.*S.*Iv - B5.*S*Nv)./(2*B2*Nv);
    u22=min(d,max(c,tempu2));
    u2=0.5*(u22+oldu2);
    
    tempphi =((lambdaS-lambdaV).*S- B6.*S)/(2*B3); 
    phi1=min(f,max(e,tempphi));
    phi=0.5*(phi1+oldphi);
end    

%%% Comment the control you don't wish to delay
for i=1:delay
    u1(i)=0;
    u2(i)=0;
    phi(i)=0;
end
%%% Once test is bigger or equal than zero the code stops
    test=min([delta*norm(u1,1)-norm(oldu1-u1,1) delta*norm(u2,1)-norm(oldu2-u2,1) delta*norm(phi,1)-norm(oldphi-phi,1) delta*norm(x,1)-norm(oldx-x,1) delta*norm(lambda,1)-norm(oldlambda-lambda,1)]);
count, test
end
%% Saving vector of solution with states, adjoints and controls
t=tvec;
S=x(:,1);
V=x(:,2);
E=x(:,3);
I=x(:,4);
Sv=x(:,5) ;
Ev=x(:,6) ;
Iv=x(:,7);
lambdaS=lambda(:,1) ;
lambdaV=lambda(:,2) ;
states=[tvec S V E I];
R=N-(S+V+E+I);
%% Total number of new infected cases
ratenewcases = b1*S.*Iv.*(1-u2)/Nv  + a1*S.*I.*(1-u1)/N ;
totalnumberinfectedcases=trapz(t,ratenewcases);
cumnewcases=zeros(M+1,1);
newcasesbyday=zeros(M+1,1);
for i=2:M+1
    cumnewcases(i)=trapz(t(1:i),ratenewcases(1:i));
end
for i=1:M
    newcasesbyday(1)=cumnewcases(1);
    newcasesbyday(i+1)=cumnewcases(i+1)-cumnewcases(i);
end
newcasesbyday;
%% Total Vaccinated Individuals
rateofvac=phi.*V;
totalvac=trapz(t,rateofvac);
%% table with rate of new cases, cumulative cases and cases by day 
T=table(transpose(t)',transpose(ratenewcases)',transpose(cumnewcases)',transpose(newcasesbyday)','VariableNames',{'time','ratenewcases','cumulativenewcases','newcasesbyday'});
%% Objective functional
integrand1=A1*I;
integrand2=B4*u1.*S+B1*u1.^2+B5*u2.*S+B2*u2.^2+B6*phi.*S+B3*phi.^2;
J1=trapz(t,integrand1);
J2=trapz(t,integrand2);
J=J1+J2;
%% Plotting the mosquito and human population
if objective == 0
figure
subplot(2,1,1);plot(t,S,'b--',t,V,'b',t,E,'.g',t,I,'r',t,R,'m','LineWidth',1.5)
if country == 'CO'
    subplot(2,1,1);title('Population Dynamics in Colombia')
else
subplot(2,1,1);title('Population Dynamics in Puerto Rico')
end
subplot(2,1,1);xlabel('Weeks');
subplot(2,1,1);ylabel('Human Classes');
subplot(2,1,1);axis tight
subplot(2,1,1);legend('$S$','$V$','$E$','$I$','$R$','interpreter','latex')
set(gca,'Fontsize',13)

subplot(2,1,2);plot(t,Sv,'b--',t,Ev,'.g',t,Iv,'r','LineWidth',1.5)
if country == 'CO'
    subplot(2,1,2);title('Mosquito Population in Colombia')
else
    subplot(2,1,2);title('Mosquito Population in Puerto Rico')
end
subplot(2,1,2);xlabel('Weeks');
subplot(2,1,2);ylabel('Mosquito Classes');
subplot(2,1,2);axis tight
subplot(2,1,2);legend('$S_v$','$E_v$','$I_v$','interpreter','latex')
set(gca,'Fontsize',13)

figure
plot(t,E,'.k',t,I,'r','LineWidth',1.5)
title('Population Dynamics of $E$ and $I$')
xlabel('Weeks');
ylabel('Human Classes');
axis tight
legend('$E$','$I$','interpreter','latex')
set(gca,'Fontsize',15)
end
end 
