%% 'Parameters used for CO'
pars_CO=[0.7073    1.5138    0.8745    3.3951    0.5463    0.0042];
pars = pars_CO;
country = 'CO';
set(0,'defaultTextInterpreter','latex');
%% 'Parameters used for PR'
pars_PR=[0.7015    2.8235    0.0073    1.0005    0.6347    0.0201];
pars = pars_PR;
country = 'PR';
%% Setting up late controls
if country == 'CO'
    delay1 = 10;
    delay2 = 20;
else
    delay1 = 20;
    delay2 = 40;
end
%% Population dynamics without controls 
ode45zikaOC(pars,0,0,0,0,country,0);
%% Solving and storing optimization solutions for J_1 or J_2
J = 1; %paper graphs are with J_1
[I0, tvec1,newcasesbyday0,total0,u10,u20,phi10,J0,J01,J02,vac0]=ode45zikaOC(pars,0,0,0,0,country,J);
%
[I1, tvec1,newcasesbyday1,total1,u11,u21,phi1,J1,J11,J12,vac1]=ode45zikaOC(pars,0.2,0,0.01,0,country,J);
% Only u1 and u2
[I2,tvec2,newcasesbyday2,total2,u12,u22,phi2,J2,J21,J22,vac2]=ode45zikaOC(pars,0,0.2,0.01,0,country,J);
% u1, u2 and phi
[I3,tvec3,newcasesbyday3,total3,u13,u23,phi3,J3,J31,J32,vac3]=ode45zikaOC(pars,0.2,0.2,0.01,0,country,J);
%% Solving and storing optimization solutions for J_1 implementing late control
[I0, tvec1,newcasesbyday0,total0,u10,u20,phi10,J0,J01,J02,vac0]=ode45zikaOC(pars,0,0,0,0,country);
[I3,tvec1,newcasesbyday3,total3,u13,u23,phi3,J3,J31,J32,vac3]=ode45zikaOC(pars,0.2,0.2,0.01,0,country);
[I1,tvec1,newcasesbyday1,total1,u11,u21,phi1,J1,J11,J12,vac1]=ode45zikaOC(pars,0.2,0.2,0.01,delay1,country);
[I2,tvec2,newcasesbyday2,total2,u12,u22,phi2,J2,J21,J22,vac2]=ode45zikaOC(pars,0.2,0.2,0.01,delay2,country);
% for the graph with only late vaccination phi, go to ode45zikaOC and
% coment the u1 and u2 controls on the delay part, then run this cell
%% Control Figure: only u1 and phi, u2 and phi, and u1,u2 and phi
t=tvec1;
figure
subplot(3,3,1);plot(t,u11,'b','LineWidth',1.5);title('$u_1$ and $\phi$ ','interpreter','latex');...
    ylabel('$u_1$');set(gca,'Fontsize',14);
subplot(3,3,4);plot(t,phi1,'b','LineWidth',1.5);title('$u_1$ and $\phi$ ','interpreter','latex');...
    ylabel('$\phi$');xlabel('Weeks');set(gca,'Fontsize',14);
subplot(3,3,2);plot(t,u22,'color',"#EDB120",'LineWidth',1.5);title('$u_2$ and $\phi$','interpreter','latex');...
    ylabel('$u_2$');set(gca,'Fontsize',14);
subplot(3,3,5);plot(t,phi2,'color',"#EDB120",'LineWidth',1.5);title('$u_2$ and $\phi$','interpreter','latex');...
    ylabel('$\phi$');xlabel('Weeks');set(gca,'Fontsize',14);
subplot(3,3,3);plot(t,u13,'m','LineWidth',1.5);title('$u_1$, $u_2$ and $\phi$','interpreter','latex');...
    ylabel('$u_1$');set(gca,'Fontsize',14);
subplot(3,3,6);plot(t,u23,'m','LineWidth',1.5);...
    ylabel('$u_2$');set(gca,'Fontsize',14);
subplot(3,3,9);plot(t,phi3,'m','LineWidth',1.5);...
    ylabel('$\phi$');xlabel('Weeks');set(gca,'Fontsize',14);
%legend('without control', 'with u1 and u2','with u1, u2 and u3')
axis tight
%% Figure: Zika Virus Incidence for Colombia/PR Outbreak with three controls
t=tvec1;
h=figure; 
plot(t,newcasesbyday0,'r','LineWidth',1.5)
hold on
plot(t,newcasesbyday1,'-*b','LineWidth',0.3)
hold on
plot(t,newcasesbyday2,'color',"#EDB120","LineStyle","--",'LineWidth',2)
hold on
plot(t,newcasesbyday3,':m','LineWidth',2)
xlabel('Weeks')
ylabel('Incidence')
hold on 
axis tight
legend('No Control','$u_1$ and $\phi$', '$u_2$ and $\phi$','$u_1$, $u_2$ and $\phi$','interpreter','latex')
if country == 'CO'
    title('Zika Virus Incidence for 2015-16 Colombia Outbreak')
else 
    title('Zika Virus Incidence for 2016-17 Puerto Rico Outbreak')
end
set(gca,'Fontsize',14)
%% Control Figure: u1, u2 and phi
t = tvec1;
figure
subplot(1,3,1);plot(t,u11,'-b',t,u12,':k','LineWidth',1.5);xlabel('Weeks');...
    ylabel('$u_1$');set(gca,'Fontsize',14);
subplot(1,3,2);plot(t,u21,'-b',t,u22,':k','LineWidth',1.5);xlabel('Weeks');...
    ylabel('$u_2$');set(gca,'Fontsize',14);
subplot(1,3,3);plot(t,phi1,'-b',t,phi2,':k','LineWidth',1.5);...
    ylabel('$\phi$');xlabel('Weeks');set(gca,'Fontsize',14);
legend(['$t= $',num2str(delay1)],['$t= $',num2str(delay2)],'interpreter','latex')
axis tight
%% Slow Control Strategy Implementation
t=tvec1;
h=figure; 
plot(t,newcasesbyday0,'r','LineWidth',1.5)
legend('without control')
hold on
plot(t,newcasesbyday3,':m','LineWidth',2)
hold on
plot(t,newcasesbyday1,'--b','LineWidth',2)
hold on
plot(t,newcasesbyday2,'-*','LineWidth',0.5)
xlabel('Weeks')
ylabel('Incidence')
hold on 
axis tight
legend('No Control','$u_1$, $u_2$ and $\phi$ at $t=0$',['$u_1$, $u_2$ and $\phi$ at $t= $',num2str(delay1)],...
    ['$u_1$, $u_2$ and $\phi$ at $t= $',num2str(delay2)],'interpreter','latex')
title('Slow Control Strategy Implementation')
set(gca,'Fontsize',14)
%% u1 and u2 with late vaccination
t=tvec1;
h=figure; 
plot(t,newcasesbyday0,'r','LineWidth',1.5)
legend('without control')
hold on
plot(t,newcasesbyday3,':m','LineWidth',2)
hold on
plot(t,newcasesbyday1,'--b','LineWidth',2)
hold on
plot(t,newcasesbyday2,'-*','LineWidth',0.5)
xlabel('Weeks')
ylabel('Incidence')
hold on 
axis tight
legend('No Control','$\phi$ at $t=0$',['$\phi$ at $t= $',num2str(delay1)],...
    ['$\phi$ at $t= $',num2str(delay2)],'interpreter','latex')
title('$u_1$ and $u_2$ with Late Vaccination')
set(gca,'Fontsize',14)
%% Control Figure: only u1, u2 and both
figure
subplot(2,2,1);plot(t,u11,'b','LineWidth',1.5);title('Only $u_1$','interpreter','latex');...
    ylabel('$u_1$');set(gca,'Fontsize',14);
subplot(2,2,3);plot(t,u22,'color',"#EDB120",'LineWidth',1.5);title('Only $u_2$','interpreter','latex');...
    ylabel('$u_2$');xlabel('Weeks');set(gca,'Fontsize',14);
subplot(2,2,2);plot(t,u13,'m','LineWidth',1.5);title('Both $u_1$ and $u_2$','interpreter','latex');...
    ylabel('$u_1$');set(gca,'Fontsize',14);
subplot(2,2,4);plot(t,u23,'m','LineWidth',1.5);xlabel('Weeks');ylabel('$u_2$');...
   set(gca,'Fontsize',14)
axis tight
%%
rat([J0 J01  J02])
%%
rat([J1 J11  J12]) 
%%
rat([J2 J21  J22]) 
%%
rat([J3 J31  J32])
%%
disp(rat([total0 total3 total1 total2]))
%%
rat([vac0 vac3  vac1 vac2])
%%
disp(rat([J0 J3 J1 J2]))
%%
disp(rat([total11 total22 total33]))
