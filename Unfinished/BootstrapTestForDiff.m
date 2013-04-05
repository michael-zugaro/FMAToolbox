%BootstrapTestForDiff - Non parametric test for significant differences between two distributions.
%
%  Compute the boostrap distribution of all possible differences between
%  two means, obtained by grouping randomly the variable values in two
%  different groups several times.
%  The function then computes the confidence intervals for the acceptance of
%  the null hypothesis of no differences between the two groups means.
%  If values corresponding to several bins are provided as input, the function
%  computes for each bin one differences distribution with confidence intervals.
%
%
%  USAGE
%
%    [h, cross, prop] = BootstrapTestForDiff(G1, G2, <options>)
%
%    G1, G2         vectors containing the variable values for the two groups
%                   If G1 and G2 are matrices the test is computed for each
%                   column.
%                   The number of rows of G1 and G2 can be different, but the
%                   number of columns have to be the same.
%                   Example: G1 matrix [tXn] containing the firing rates for
%                   each trial t and for each binned position n, corresponding
%                   to leftwards trajecoties (first group); same for matrix G2
%                   corresponding to rightwards trajectories (second group).
%                   the number for trials
%
%    <options>      optional list of property-value pairs (see table below)
%
%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%     'Nshuff'      number of shuffles used to create the distribution
%                   (default = 5000)
%     'alfa'        confidence level (default = 0.05)
%     'IterLimit'   maximum number of iterations for the global analysis
%                   loop (default = 6)
%     'Tolerance'   maximum accepted difference between the desired global
%                   alfa (in %) and the obtained one (defalut = 0.8)
%     'tail'        distribution type: two-tailed ditribution.
%                   (default = two) (To be adapted for one)
%     'fig'         on off (default 'off')
%    =========================================================================
%
%  OUTPUT
%
%    h                  h = 0 if the null hypothesis is accepted
%                       h = 1 if the null hypothesis is rejected
%                       If the input are matrices, then h = 1 indicates that
%                       the null hypothesis is rejected at least for one bin.
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
%                       level of confidence alfa (or [2Xn] matrix).
%    prop.globalCI      upper and lower limit corresponding to the global
%                       level of confidence alfa (or [2Xn] matrix).
%    prop.gAlfa         local alfa which lead to the global alfa, given for
%                       each iteration
%    prop.gPerc         global alfa multiplied per 100, given for each
%                       iteration
%


% Copyright (C) 2010-2011 by Erika Cerasti
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

function [h, cross, prop] = BootstrapTestForDiff(G1, G2, varargin)


% Default values
Nshuff = 5000;
alfa = 0.05;
tail = 'two';
fig_display = 1;
IterLimit = 6;
Tol =0.8;

% Check number of parameters
if nargin < 2 | mod(length(varargin),2) ~= 0,
    error('Incorrect number of parameters (type ''help <a href="matlab:help BootstrapTestForDiff">BootstrapTestForDiff</a>'' for details).');
end

if isempty(G1) || isempty(G2), return; end

% Check parameter sizes
if size(G1,2) ~= size(G2,2),
    error('"G1" and "G2" have not the same number of columns (type ''help <a href="matlab:help BootstrapTestForDiff">BootstrapTestForDiff</a>'' for details).');
end

if size(G1,2) > 1,
    multipleComp=1;       %multiple comparison = 1 --> global analysis needed
end

% Parse parameter list
for i = 1:2:length(varargin),
    if ~ischar(varargin{i}),
        error(['Parameter ' num2str(i+1) ' is not a property (type ''help <a href="matlab:help BootstrapTestForDiff">BootstrapTestForDiff</a>'' for details).']);
    end
    switch(lower(varargin{i})),
        case 'Nshuff',
            Nshuff = lower(varargin{i+1});
            if ~isdscalar(Nshuff,'>0'),
                error('Incorrect value for property ''Nshuff'' (type ''help <a href="matlab:help BootstrapTestForDiff">BootstrapTestForDiff</a>'' for details).');
            end

        case 'alfa',
            alfa = varargin{i+1};
            if ~isdscalar(alfa,'>0','<1'),
                error('Incorrect value for property ''alfa'' (type ''help <a href="matlab:help BootstrapTestForDiff">BootstrapTestForDiff</a>'' for details).');
            end

        case 'IterLimit',
            IterLimit = varargin{i+1};
            if ~isdscalar(IterLimit,'>1'),
                error('Incorrect value for property ''IterLimit'' (type ''help <a href="matlab:help BootstrapTestForDiff">BootstrapTestForDiff</a>'' for details).');
            end

        case 'Tolerance',
            Tol = varargin{i+1};
            if ~isdscalar(Tol,'>0'),
                error('Incorrect value for property ''Tolerance'' (type ''help <a href="matlab:help BootstrapTestForDiff">BootstrapTestForDiff</a>'' for details).');
            end
        case 'tail',
            tail = varargin{i+1};
            if ~isstring(tail,'one','two'),
                error('Incorrect value for property ''tail'' (type ''help <a href="matlab:help BootstrapTestForDiff">BootstrapTestForDiff</a>'' for details).');
            end

        case 'fig',
            fig_display = varargin{i+1};
            if ~isstring(fig_display,'on','off'),
                error('Incorrect value for property ''fig'' (type ''help <a href="matlab:help BootstrapTestForDiff">BootstrapTestForDiff</a>'' for details).');
            end
            if strcmp(fig_display,'on')
                fig_display=1;
            end
        otherwise,
            error(['Unknown property ''' num2str(varargin{i}) ''' (type ''help <a href="matlab:help BootstrapTestForDiff">BootstrapTestForDiff</a>'' for details).']);
    end
end

disp('  Start SHUFFLING process');
%Create shuffled data and calculate the shuffled LR difference
for sh=1:Nshuff

    %function for shuffling "ShuffleSpike"
    trial1 = size(G1,1);
    trial2 = size(G2,1);
    gnum = trial1+trial2;
    Gtot= [G1; G2];

    indShuf=randperm(gnum);
    Gtot = Gtot(indShuf,:);


    G1shu = Gtot(1:trial1,:);
    G2shu = Gtot(trial1+1:end, :);

    M1=mean(G1shu);
    M2=mean(G2shu);
    AllDiffDist(sh,:)=M1-M2;       %matrix of all differences distribution [Nshuff X n]

    step=sh-floor(sh/1000)*1000;
    if(step==0)
        disp(['  Shuffle #'  num2str(sh) ' done']);
    end

end   %end shuffling

disp(' End SHUFFLE process');
disp(' Confidence Interval Computation...  ');

Nbin = size(G1,2);
deviation=100;
iter=0;

while (deviation>Tol)

    clear distr
    iter=iter+1;
    if(iter > IterLimit)
        disp(['STOP: More than ' num2str(IterLimit) ' iterations']);
        break
    end

    qvect=[(alfa/2), 1-(alfa/2)];
    if strcmp(tail,'one')
        qvect=[1-alfa];
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

    if(multipleComp==0)
        CIntPB = CIntBand;
        break
    end
    %--------STATISTICAL GLOBALWISEANALYSIS--------------
    for sh=1:Nshuff
        cont=0;

        ShufdLine= AllDiffDist(sh,:);
        indsup = find(ShufdLine > CIntBand(1,:));
        indinf = find(ShufdLine < CIntBand(2,:));
        if(~isempty(indsup) || ~isempty(indinf))
            cont=cont+1;
        end
    end

    perc=(cont/Nshuff)*100;
    alfaprc=alfa*100;
    deviation = abs( perc - alfaprc );

    GlobPerc(iter)=perc;
    GlobAlfa(iter)=alfa;

    if(iter==1)
        CIntPB = CIntBand;
        alfa = 0.003;
    else
        ss=sign(perc - alfaprc);
        if(deviation>5)
            alfa = alfa - ss*deviation*0.0000075;
        else
            alfa = alfa - ss*deviation*0.00033;
        end
        if( alfa < 0.001 )
            disp('STOP: new alfa is too small');
            break
        end
    end


end

M1ori=mean(G1);
M2ori=mean(G2);
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
if(multipleComp==0)

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

if (multipleComp) && (fig_display)

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
%     info = ['\bf Global alfa p<' num2str(alfa) 'Local alfa p<' num2str(alfa)];
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



