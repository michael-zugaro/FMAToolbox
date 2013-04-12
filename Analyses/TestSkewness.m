function [h,p,stats] = TestSkewness(control,repeat,test,varargin)

%TestSkewness - Test if firing field skewness changes between two conditions.
%
%  Current implementation only tests 1D environments.
%
%  USAGE
%
%    [h,p,cr,ca] = TestSkewness(control,repeat,test,<options>)
%
%    control        firing fields in control condition (MxN: M fields, N bins)
%    repeat         optional firing fields in repeated control condition
%    test           firing fields in test condition
%    <options>      optional list of property-value pairs (see table below)
%
%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%     'alpha'       significance level (default = 0.05)
%     'iterations'  number of iterations for bootstrap (default = 150)
%     'show'        plot results (default = 'off')
%    =========================================================================
%
%  OUTPUT
%
%    h                  1 if H0 (random remapping) can be rejected
%    p                  p-value of bootstrap test
%    stats.control.m	   median skewness (control)
%    stats.control.s	   standard error of the median skewness (control)
%    stats.repeat.m	   median skewness (repeat)
%    stats.repeat.s	   standard error of the median skewness (repeat)
%    stats.test.m	      median skewness (test)
%    stats.test.s	      standard error of the median skewness (test)

% Copyright (C) 2013 by Michaël Zugaro and Anne Cei
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

% Defaults
show = 'off';
nBootstrap = 150;
alpha = 0.05;

% Check number of parameters
if nargin < 3 | mod(length(varargin),2) ~= 0,
  error('Incorrect number of parameters (type ''help <a href="matlab:help TestSkewness">TestSkewness</a>'' for details).');
end

% Check parameter sizes
if ~isdmatrix(control) || ~isdmatrix(repeat) || ~isdmatrix(test),
	error('All firing fields should be MxN matrices (type ''help <a href="matlab:help TestSkewness">TestSkewness</a>'' for details).');
end

% Parse parameter list
for i = 1:2:length(varargin),
	if ~ischar(varargin{i}),
		error(['Parameter ' num2str(i+2) ' is not a property (type ''help <a href="matlab:help TestSkewness">TestSkewness</a>'' for details).']);
	end
	switch(lower(varargin{i})),
		case 'alpha',
			alpha = varargin{i+1};
			if ~isdscalar(alpha,'>=0','<=1'),
				error('Incorrect value for property ''alpha'' (type ''help <a href="matlab:help TestSkewness">TestSkewness</a>'' for details).');
			end
		case 'iterations',
			nBootstrap = varargin{i+1};
			if ~isiscalar(nBootstrap,'>0'),
				error('Incorrect value for property ''iterations'' (type ''help <a href="matlab:help TestSkewness">TestSkewness</a>'' for details).');
			end
		case 'show',
			show = varargin{i+1};
			if ~isstring(show,'on','off'),
				error('Incorrect value for property ''show'' (type ''help <a href="matlab:help TestSkewness">TestSkewness</a>'' for details).');
			end
		otherwise,
			error(['Unknown property ''' num2str(varargin{i}) ''' (type ''help <a href="matlab:help TestSkewness">TestSkewness</a>'' for details).']);
	end
end

% Matrix sizes
[mc,nc] = size(control);
[mr,nr] = size(repeat);
[mt,nt] = size(test);

% Normalize: transform firing rates into firing probabilities
control = control ./ repmat(sum(control,2),1,nc);
repeat = repeat ./ repmat(sum(repeat,2),1,nr);
test = test ./ repmat(sum(test,2),1,nt);

% Compute skewness
skewnessControl = UnbiasedSkewness(control);
skewnessRepeat = UnbiasedSkewness(repeat);
skewnessTest = UnbiasedSkewness(test);

% Test distribution differences and compute medians
[hc,pc] = kstest2(skewnessControl,skewnessRepeat,alpha);
[h,p] = kstest2(skewnessControl,skewnessTest,alpha);
stats.control.m = median(skewnessControl);
stats.repeat.m = median(skewnessRepeat);
stats.test.m = median(skewnessTest);
stats.control.s = semedian(skewnessControl);
stats.repeat.s = semedian(skewnessRepeat);
stats.test.s = semedian(skewnessTest);

% Plot results
if strcmp(show,'on'),
	figure;hold on;
	% Test results
	if hc,
		sc = ['p < ' num2str(alpha)];
	else
		sc = 'NS';
	end
	if h,
		st = ['p < ' num2str(alpha)];
	else
		st = 'NS';
	end
	% Histograms
	x = linspace(-1,1,20);
	dx = x(2)-x(1);
	hControl = hist(skewnessControl,x)/mc;
	hRepeat = hist(skewnessRepeat,x)/mr;
	hTest = hist(skewnessTest,x)/mt;
	% Plot control and repeat
	subplot(2,1,1);hold on;
	bar(x,hControl,0.5,'b');
	bar(x+dx/2,hRepeat,0.5,'k');
	PlotHVLines(stats.control.m,'v','b');
	PlotHVLines(stats.repeat.m,'v','k');
	PlotHVLines(0,'v','k--');
	legend(['control ' mean2str(stats.control.m,stats.control.s,mc)],['repeat ' mean2str(stats.repeat.m,stats.repeat.s,mr)],'location','northwest');
	title(['Control vs Repeat (KS test: ' sc ')']);
	xlabel('Skewness');
	ylabel('Probability');
	% Plot control and test
	subplot(2,1,2);hold on;
	bar(x,hControl,0.5,'b');
	bar(x+dx/2,hTest,0.5,'r');
	PlotHVLines(stats.control.m,'v','b');
	PlotHVLines(stats.test.m,'v','r');
	PlotHVLines(0,'v','k--');
	legend(['control ' mean2str(stats.control.m,stats.control.s,mc)],['repeat ' mean2str(stats.test.m,stats.test.s,mt)],'location','northwest');
	title(['Control vs Test (KS test: ' st ')']);
	xlabel('Skewness');
	ylabel('Probability');
end

function s = UnbiasedSkewness(distributions)

[m,n] = size(distributions);
x = repmat(1:n,m,1);

% Mean for each distribution
mu = sum(x.*distributions,2);
% Centered x
x0 = x-repmat(mu,1,n);
% Third central moment
m3 = sum( x0.^3 .* distributions,2 );
% Second central moment
s2 = sum( x0.^2 .* distributions,2 );
% Biased skewness
s = m3./s2.^1.5;
% Unbiased skewness
s = sqrt(n*(n-1))/(n-2) * s;
