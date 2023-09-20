function value=errorfun_PR(pars)
%%% PR Data 2016 Zika Outbreak
%%% Confirmed Cases of 2016 ZIKV (cumulative) by weeks
%%%Loading the data
r=pars(6);
cumulative2016=[4
9
12
21
54
102
147
192
192
340
362
426
540
616
697
775
915
1098
1160
1342
1491
1716
2152
2377
3091
4427
5572
7286
8766
10680
13176
14324
16527
17861
19957
22348
24117
26691
29074
29965
31454
32730
33445
34060
34552
34815
35126
35638
35860
36316
36365
36921];
cumulative2017=[37478
37889
38297
38461
38733
38940
39339
39529
39559
39636
39815
39839
39984
40067
40134
40139
40140
40274
40301
40330
40398
40374
40357
40365
40383
40392
40421
40427
40460
40545
40569
40570
40576
40582
40588
40562
40562
40562
40562
40562
40562
40562
40564
40573
40573
40601
40610
40630
40630
40630
40630
40630];
cumulative=[cumulative2016;cumulative2017]/r;
a=length(cumulative);
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
%% Loading approximation by the model 
[cumulative_model newcasesbyweek]=ode45zika(pars,'PR');
%% Square error between incidence from ode and data
newcases=[];
 newcases(1)=M(1);
 for i=2:a
     newcases(i)=M(i)-M(i-1);
 end
for i=1:a
    if newcases(i)<0
        newcases(i)=0;
    end
end
newcases=newcases';
%% Incidence Square error
diff_inc=newcasesbyweek-newcases;
value=norm(diff_inc,2)/norm(newcases,2);
end