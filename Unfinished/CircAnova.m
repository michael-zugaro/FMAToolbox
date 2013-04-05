% Circular ANOVA with significance calculated from 
% randomization method
% 
% p = CircAnova(Th, Gp)
% or 
% p = CircAnova({Th1, Th2, Th3, ...})
%
% computes mean angles for each group and residual sum as
% sum(cos(th-pred))
%
% significance is calculated by randomizing groups
%
% CircAnova(..., nRands) specifies number of randomizations.

function [p, r0, r] = CircAnova(Arg1, Arg2, Arg3)

nRands = 1000; % default value

% deal with argument form 2
if iscell(Arg1)
    nGps = length(Arg1);
    
    Th = []; Gp = [];
    for g=1:nGps
        Th = [Th ; Arg1{g}(:)];
        Gp = [Gp ; g*ones(length(Arg1{g}),1)];
    end
   if nargin==2
       nRands = Arg2;
   end
else
    Th = Arg1;
    OrigGp = Arg2;
    % relabel groups 1...nGps
    [dummy dummy2 Gp] = unique(OrigGp);
    nGps = max(Gp);
   if nargin==3
       nRands = Arg3;
   end
end

ExpAng = exp(j*Th);

r0 = resid(ExpAng, Gp);
N = length(Th);
for i=1:nRands
    p = randperm(N);
    r(i) = resid(ExpAng, Gp(p));
end

p = 1-Quantile(r, r0);

% hold off; hist(r, 100);
% hold on; plot([r0 r0], ylim, 'r');

function r = resid(ExpAng, Gp)

r = 0;
for g=1:max(Gp)
    My = ExpAng(find(Gp==g));
    if length(My)==0
        continue;
    end
    mu = mean(My); mu0 = mu/abs(mu);
    r = r + sum(real(My./mu0));
end