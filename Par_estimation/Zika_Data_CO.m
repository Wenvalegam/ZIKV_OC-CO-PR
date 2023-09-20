%% CO Data 2015-2016 Zika Outbreak
%% Confirmed Cases of 2015-2016 ZIKV (cumulative) by weeks
%% 2015 from week 42-52 (11 weeks)
cumulative2015=[389
827
1241
2077
3213
4278
5678
7516
8992
9868
10026]; 
%% 2016-52 weeks 
cumulative2016=[13531
16419
20297
25645
31555
37011
42706
47789
51473
55724
58838
61393
64839
68630
71952
75187
78085
80793
83889
87355
88945
91156
93242
94946
96494
97628
98788
99721
100466
101145
101668
102107
102341
102654
102938
103550
103955
104238
104465
104619
104724
104755
105085
105247
105372
105518
105555
105864
105962
106558
106562
106659];
%% Define parameters 
% pars=[zikaparameters(1,:)];
 pars=[0.7073    1.5138    0.8745    3.3951    0.5463    0.0042];
 %pars=[0.7073    1.5138    0.8745    3.3951    0.5463    0.015];
[cumulative_model newcases_model]= ode45zika(pars,'CO');
%% Concatenating cummulative data 
%r=0.02; %reporting rate from 0 to 1
r=pars(6);
cumulative=[cumulative2015; cumulative2016]/r;
%%
a=length(cumulative);
t=linspace(1,a,a);
%% adjusting data such that M(i)<= M(i+1)
M=[];
for i=2:a
    if cumulative(i)-cumulative(i-1)<0
       M(i)=cumulative(i-1);
    else
       M(i)=cumulative(i);
    end
M(1)=cumulative(1);    
end
M=M';
 %% Incidence or New Cases 
 newcases=[];
 newcases(1)=M(1);
 for i=2:a
     newcases(i)=M(i)-M(i-1);
     if newcases(i)< 0
         newcases(i)=0;
     end
 end
newcases=newcases';
%%
figure
 scatter(t,newcases,'b','LineWidth',1.5)
 legend('Data','Nort')
 hold on
 plot(t,newcases_model,'r','LineWidth',1.5)
 xlabel('Weeks');ylabel('Incidence')
 axis([1 a 0 max(max(newcases_model), max(newcases))+10])
 legend('Data','Approximation','Location','northeast')
 title('2015-16 Zika Outbreak in Colombia')
 set(gca,'FontSize',13)
 grid on