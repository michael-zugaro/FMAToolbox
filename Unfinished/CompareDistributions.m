%CompareDistributions - Non-parametric comparison of two N-dimensional distributions.
%
%  This test determines in which dimensions the two distributions are
%  significantly different. For example, the two distributions could be
%  leftward vs rightward trajectories in a T maze (each 'dimension' of these
%  distributions corresponds to a spatial bin), and the test would determine
%  where the trajectories are significantly different.
%
%  This test is based on the bootstrap method developed in Fujisawa et al.
%  (2008). This compares the observed difference between the two distribution
%  averages and averages and confidence intervals for surrogate data (where
%  individual observations are shuffled across groups). Both pointwise and
%  global confidence intervals are computed.
%
%  USAGE
%
%    [h,cross,prop] = CompareDistributions(group1,group2,<options>)
%
%    group1,group2  values for the two groups: lines are observations, columns
%                   are dimensions (e.g. spatial bins).
%    <options>      optional list of property-value pairs (see table below)
%
%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%     'shuffles'    number of shuffles to estimate the distribution
%                   (default = 5000)
%     'alpha'       confidence level (default = 0.05)
%     'max'         maximum number of iterations for global confidence
%                   intervals (default = 6)
%     'tolerance'   maximum difference (in %) between target and computed
%                   global alpha (defalut = 0.8)
%     'tail'        one or two-tailed ditributions (default = two)
%                   (to be adapted for one)
%     'show'        show figure (default 'off')
%    =========================================================================
%
%  OUTPUT
%
%    h                  h = 0 if the null hypothesis is accepted
%                       h = 1 if the null hypothesis is rejected
%                       (for at least one bin)
%    cross.sup          cross.sup = 1 if the null hypothesis is rejected,
%                       because the tested value exceed the superior limit of
%                       confidence. If the input are matrices, then cross.sup
%                       is a vector (one value per bin) whose ones correspond
%                       to points where the null hypothesis is rejected.
%    cross.inf          cross.inf = 1 if the null hypothesis is rejected,
%                       because the tested value exceed the inferior limit of
%                       confidence. If the input are matrices, then cross.inf
%                       is a vector (one value per bin) whose ones correspond
%                       to points where the null hypothesis is rejected.
%    prop.original      original value of the difference.
%    prop.mean          mean value of the differences distribution (or vector)
%    prop.localCI       upper and lower limit corresponding to the local
%                       level of confidence alpha (or [2Xn] matrix).
%    prop.globalCI      upper and lower limit corresponding to the global
%                       level of confidence alpha (or [2Xn] matrix).
%    prop.gAlfa         local alpha which lead to the global alpha, given for
%                       each iteration
%    prop.gPerc         global alpha multiplied per 100, given for each
%                       iteration
%

% Copyright (C) 2010-2011 by Erika Cerasti, 2011 by Michaël Zugaro
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

function [h,cross,prop] = CompareDistributions(group1,group2,varargin)

% Default values
nShuffles = 5000;
alpha = 0.05;
tail = 'two';
show = 'off';
maxIterations = 6;
tolerance = 0.8;

% Check number of parameters
if nargin < 2 | mod(length(varargin),2) ~= 0,
    error('Incorrect number of parameters (type ''help <a href="matlab:help CompareDistributions">CompareDistributions</a>'' for details).');
end

if isempty(group1) || isempty(group2), return; end

% Check parameter sizes
if size(group1,2) ~= size(group2,2),
    error('The two groups do not have the same number of columns (type ''help <a href="matlab:help CompareDistributions">CompareDistributions</a>'' for details).');
end

if size(group1,2) > 1,
	nDimensional = 1;       % N dimensions (N>1) => global confidence intervals needed
end

% Parse parameter list
for i = 1:2:length(varargin),
	if ~ischar(varargin{i}),
		error(['Parameter ' num2str(i+1) ' is not a property (type ''help <a href="matlab:help CompareDistributions">CompareDistributions</a>'' for details).']);
	end
	switch(lower(varargin{i})),
		case 'shuffles',
			nShuffles = lower(varargin{i+1});
			if ~isdscalar(nShuffles,'>0'),
				error('Incorrect value for property ''shuffles'' (type ''help <a href="matlab:help CompareDistributions">CompareDistributions</a>'' for details).');
			end
		case 'alpha',
			alpha = varargin{i+1};
			if ~isdscalar(alpha,'>0','<1'),
				error('Incorrect value for property ''alpha'' (type ''help <a href="matlab:help CompareDistributions">CompareDistributions</a>'' for details).');
			end
		case 'max',
			maxIterations = varargin{i+1};
			if ~isdscalar(maxIterations,'>1'),
				error('Incorrect value for property ''max'' (type ''help <a href="matlab:help CompareDistributions">CompareDistributions</a>'' for details).');
			end
		case 'tolerance',
			tolerance = varargin{i+1};
			if ~isdscalar(tolerance,'>0'),
				error('Incorrect value for property ''tolerance'' (type ''help <a href="matlab:help CompareDistributions">CompareDistributions</a>'' for details).');
			end
		case 'tail',
			tail = varargin{i+1};
			if ~isstring(tail,'one','two'),
				error('Incorrect value for property ''tail'' (type ''help <a href="matlab:help CompareDistributions">CompareDistributions</a>'' for details).');
			end
		case 'show',
			show = varargin{i+1};
			if ~isstring(show,'on','off'),
				error('Incorrect value for property ''show'' (type ''help <a href="matlab:help CompareDistributions">CompareDistributions</a>'' for details).');
			end
		otherwise,
			error(['Unknown property ''' num2str(varargin{i}) ''' (type ''help <a href="matlab:help CompareDistributions">CompareDistributions</a>'' for details).']);
	end
end

disp('  Start SHUFFLING process');
% Shuffle data and calculate surrogate differences in averages
n1 = size(group1,1);
n2 = size(group2,1);
all = [group1;group2];
for i = 1:nShuffles,
	shuffled = randperm(n1+n2);

	g1 = group(group(:,2)==1,:);
	indShuf = randperm(gnum);
	Gtot = Gtot(indShuf,:);


    G1shu = Gtot(1:trial1,:);
    G2shu = Gtot(trial1+1:end, :);

    M1=mean(G1shu);
    M2=mean(G2shu);
    AllDiffDist(i,:)=M1-M2;       %matrix of all differences distribution [nShuffles X n]

    step=i-floor(i/1000)*1000;
    if(step==0)
        disp(['  Shuffle #'  num2str(i) ' done']);
    end

end   %end shuffling

disp(' End SHUFFLE process');
disp(' Confidence Interval Computation...  ');

Nbin = size(group1,2);
deviation=100;
iter=0;

while (deviation>tolerance)

    clear distr
    iter=iter+1;
    if(iter > maxIterations)
        disp(['STOP: More than ' num2str(maxIterations) ' iterations']);
        break
    end

    qvect=[(alpha/2), 1-(alpha/2)];
    if strcmp(tail,'one')
        qvect=[1-alpha];
    end

    %---------STATISTICAL POINTWISE ANALYSIS------------
    for n=1:Nbin

        Distr= AllDiffDist(:,n);
        if(Distr==0)
            continue
        end
        CInt=quantile(Distr, qvect);
        CIntBand(:,n)=flipdim(CInt',1);

        clear CInt
    end

    if(nDimensional==0)
        CIntPB = CIntBand;
        break
    end
    %--------STATISTICAL GLOBALWISEANALYSIS--------------
    for sh=1:nShuffles
        cont=0;

        ShufdLine= AllDiffDist(sh,:);
        indsup = find(ShufdLine > CIntBand(1,:));
        indinf = find(ShufdLine < CIntBand(2,:));
        if(~isempty(indsup) || ~isempty(indinf))
            cont=cont+1;
        end
    end

    perc=(cont/nShuffles)*100;
    alfaprc=alpha*100;
    deviation = abs( perc - alfaprc );

    GlobPerc(iter)=perc;
    GlobAlfa(iter)=alpha;

    if(iter==1)
        CIntPB = CIntBand;
        alpha = 0.003;
    else
        ss=sign(perc - alfaprc);
        if(deviation>5)
            alpha = alpha - ss*deviation*0.0000075;
        else
            alpha = alpha - ss*deviation*0.00033;
        end
        if( alpha < 0.001 )
            disp('STOP: new alpha is too small');
            break
        end
    end


end

M1ori=mean(group1);
M2ori=mean(group2);
OriginalDiff = M1ori-M2ori;

prop.original = OriginalDiff;
prop.mean = mean(AllDiffDist,1);
prop.localCI = CIntPB;
prop.globalCI = [];
prop.galfa = [];
prop.gperc = [];

if(iter > 1)
    CIntGB = CIntBand;
    prop.globalCI = CIntBand;
    prop.galfa = GlobAlfa;
    prop.gperc = GlobPerc;
end



%--------Find significant segments--------------
if(nDimensional==0)

    if( OriginalDiff > CIntPB(1) || OriginalDiff < CIntPB(2) )
        h=1;
    else
        h=0;
    end

    cross.sup = [];
    cross.inf = [];

else
    cxsup = OriginalDiff > CIntPB(1,:);       %cxsup=crossing superior band
    cxinf = OriginalDiff < CIntPB(2,:);       %cxinf=crossing inferior band
    cxGsup = OriginalDiff > CIntGB(1,:);
    cxGinf = OriginalDiff < CIntGB(2,:);

    Psup = (cxsup & cxGsup);
    PsupInd = find(cxsup & cxGsup);

    if(~isempty(PsupInd))
        induni = logical([1, (diff(PsupInd)~=1)]);
        PsupStart = PsupInd(induni);

        for t=1: length(PsupStart)
            i=PsupStart(t);
            step=1;
            while(i-step)~=0 && (cxsup(i-step)==1)
                Psup(i-step)=1;
                step=step+1;
            end
            step=1;
            while(i+step)<=length(cxsup) && (cxsup(i+step)==1)
                Psup(i+step)=1;
                step=step+1;
            end
        end
    end


    Pinf = (cxinf & cxGinf);
    PinfInd = find(cxinf & cxGinf);
    if(~isempty(PinfInd))
        induni = logical([1, (diff(PinfInd)~=1)]);
        PinfStart = PinfInd(induni);

        for t=1: length(PinfStart)
            i=PinfStart(t);
            step=1;
            while(i-step)~=0 && (cxinf(i-step)==1)
                Pinf(i-step)=1;
                step=step+1;
            end
            step=1;
            while(i+step)<=length(cxinf) && (cxinf(i+step)==1)
                Pinf(i+step)=1;
                step=step+1;
            end
        end
    end

    cross.sup = Psup;
    cross.inf = Pinf;

    if( all(Psup==0) && all(Pinf==0) )
        h=0;
    else
        h=1;
    end
end

if (nDimensional) && strcmp(show,'on'),

    orange = [1 0.4 0];
    brigthblue = [0.3 0.5 0.9];
    green=[0.6 0.9 0.2];
    darkgreen = [0.15 0.35 0.15];

    Xaxes = [1:length(OriginalDiff)];
    Yaxes = zeros(1,length(OriginalDiff));
    figure
    hold on
    plot(Xaxes, prop.globalCI, 'o-', 'MarkerSize', 4, 'Color', darkgreen, 'MarkerFaceColor', darkgreen, 'MarkerEdgeColor', darkgreen)
    plot(Xaxes, prop.localCI, 'o-', 'MarkerSize', 4, 'Color', green, 'MarkerFaceColor', green, 'MarkerEdgeColor', green)
    plot(Xaxes, OriginalDiff, 'Color', orange,'linewidth',2)
    plot(Xaxes, Yaxes, '--','Color', 'k')
    ylabel('\bf Mean Differences','FontSize',14)
    xlabel('\bf # bin','FontSize',14)

%     ylim([-20 20]);
%     info = ['\bf Global alpha p<' num2str(alpha) 'Local alpha p<' num2str(alpha)];
%     text(x, y, info,'FontSize',12, 'Color', 'k')

    supInt = diff(cross.sup);
    infInt = diff(cross.inf);
    StartSup = find(supInt==1);
    StartInf = find(infInt==1);
    EndSup = find(supInt==-1);
    EndInf = find(infInt==-1);

    if(cross.sup(1) == 1)
        StartSup = [1, StartSup];
    elseif(cross.sup(end) == 1)
        EndSup = [EndSup, length(cross.sup)];
    end

    if(cross.inf(1) == 1)
        StartInf = [1, StartInf];
    elseif(cross.inf(end) == 1)
        EndInf = [EndInf, length(cross.inf)];
    end


    PairsInf = [StartInf' , EndInf'];
    PairsSup = [StartSup' , EndSup'];
    IntervalsPairs = sort([PairsSup; PairsInf]);

    PlotIntervals(IntervalsPairs, 'rectangles')
end

%    vedi
%    if (exist('supGB'))
%    else
%        supGB=supPB;
%        infGB=infPB;
%        GlobAlfa((ind+1),sig)=alfa_new;
%    end

end



