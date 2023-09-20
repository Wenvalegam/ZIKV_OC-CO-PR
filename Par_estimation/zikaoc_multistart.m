%Parameter vector we are approximating:
% pars = [ b1  b2   a1   k    c2  r];
country = 'PR'; % or 'PR'
if country == 'PR'
    % PR
    LowerBounds = [ 0.7  0.7  0.007   1    7/13   0.01];
    UpperBounds = [ 15   14   0.7   7/2   7/8    0.3];
else
    % CO
    LowerBounds = [ 0.7  1.5  0.007   1    7/13 0.001];
    UpperBounds = [ 6   6     3     7/2     7/8  0.3];
end
%%
xstart=.5*(LowerBounds+UpperBounds); %start in the middle

if country == 'PR'
    problem = createOptimProblem('fmincon','x0',xstart,'objective',@errorfun_PR...
    ,'lb',LowerBounds,'ub',UpperBounds);
else 
    problem = createOptimProblem('fmincon','x0',xstart,'objective',@errorfun_CO...
    ,'lb',LowerBounds,'ub',UpperBounds);
end
problem.options = optimoptions(problem.options,'MaxFunEvals',3000,'MaxIter',3000);%'MaxFunEvals',99999,'MaxIter',99999);%,'TolFun',0,'TolCon',0)

numstartpoints=20; % Test this first to check it runs, then comment and run 200 startpoints 
%numstartpoints=200; 
%ms=MultiStart('Display','iter'); 
ms=MultiStart('UseParallel',true,'Display','iter');    %defines a multistart problem

[parfit_cum2,fval,exitflag,output,manymins]=run(ms,problem,numstartpoints);  %runs the multistart 
%problem  with "k" instances. It uses start values specified by xstart
%as well as a uniform mesh from the lower bounds to upper bounds.

%from "run" literature: MultiStart generates k - 1 start points using the same
%algorithm as list for a RandomStartPointSet object. MultiStart also uses 
%the x0 point from the problem structure.

%"run" refers to "points". From points literature: A k-by-n matrix. 
%The number of rows k is the number of start points that RandSet specifies.
% The number of columns n is the dimension of the start points. n is equal to
% the number of elements in the x0 field in problem. The MultiStart algorithm 
% uses each row of points as an initial point in an optimization.
%i chose to use "manymins" instead of "solutions" because it allows me to
%more easily list solutions in a matrix when there are numerous runs (100
%for example)


% the following takes solutions from manymins and makes a matrix out of them
global zikaparameters

for i=1:length(manymins)
    zikaparameters(i,:)=manymins(i).X;
end

for i=1:length(manymins)
    fval(i)=manymins(i).Fval;
end
fval=fval';

for i=1:length(manymins)
    EF(i)=manymins(i).Exitflag;
end
EF=EF';

beep on;
beep
beep
beep
beep
beep
beep