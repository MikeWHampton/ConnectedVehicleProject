syms alpha01 kappa01 alpha01 beta01 gamma01 s sigma
a0 = [0 -1; 0 0];
asig = [0 0; alpha01*kappa01 -(alpha01+beta01)];
bhatsig = [0; gamma01];
b0 = [1; 0];
bsig = [0; beta01];
c = [0 1];

TF = simplify(c*inv(s*eye(2)-a0-asig*exp(-s*sigma))*(b0+bsig*exp(-s*sigma)+bsig*s*exp(-s*sigma)))
alpha01 = 0.1; beta01 = 0.25; sigma = 0.6; kappa01 = 0.6;
CharEq1 = alpha01*kappa01 + alpha01*s + beta01*s + s^2*exp(s*sigma) == 0
S1 = solve(CharEq1,s)

alpha01 = 0.001; beta01 = 0; beta02 = 0.01; beta03 = 0.3; gamma01 = 0.03; gamma02 = 0; gamma03 = 0.1;
CharEq2 = s^2*exp(s*sigma)+(alpha01+beta01+beta02+beta03)*s+alpha01*kappa01
S2 = solve(CharEq2,s)