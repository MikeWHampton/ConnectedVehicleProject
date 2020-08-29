alpha01 = 0.3; beta01 = 0.4; beta02 = 0.3; gamma01 = 0.3; gamma02 = 0.2;
x0 = [alpha01 beta01 beta02 gamma01 gamma02]
rng default % For reproducibility
gs = GlobalSearch;
func = @(x)EConsumpFunc(x);
problem = createOptimProblem('fmincon','x0',x0,...
    'objective',func,'lb',[0.1,0.1,0.1,0.1,0.1],'ub',[5,2,2,2,2]);
x = run(gs,problem)
