%	WorkingMemory


function [counts,indices] = WMCounts(data,baitedArms)


RAT = 1;
DAY = 2;
TRIAL = 3;
CONF = 4;
GROUP = 5;
ARM = 6;
WORKB = 6;
WORKU = 7;


% Constants
nDays = max(data(:,DAY));
nTrials = 3;
nArms = 8;
nGroups = max(data(:,GROUP));
nConfigurations = max(data(:,CONF));
nRats = length(unique(data(:,RAT)));


% Let's go...
line = 1;
for rat = 1:nRats,

	% Which group does this rat belong to?
	group = data(data(:,RAT)==rat,GROUP);
	group = group(1);

	for configuration = 1:nConfigurations,

		% Determine baited/unbaited arms for this rat and this configuration
		b = baitedArms(rat,configuration,:);
		baitedArm1 = b(1);
		baitedArm2 = b(2);
		baitedArm3 = b(3);
		unbaitedArms = setdiff(1:8,[baitedArm1 baitedArm2 baitedArm3]);
		unbaitedArm1 = unbaitedArms(1);
		unbaitedArm2 = unbaitedArms(2);
		unbaitedArm3 = unbaitedArms(3);
		unbaitedArm4 = unbaitedArms(4);
		unbaitedArm5 = unbaitedArms(5);

		for day = 1:nDays,
			for trial = 1:nTrials,

				% Select the block of data corresponding the current rat, config, trial and day
				block = data(data(:,RAT)==rat&data(:,CONF)==configuration&data(:,DAY)==day&data(:,TRIAL)==trial,:);

				% Total number of visits
				nVisits = size(block,1);

				% Missing data: skip
				if nVisits == 0,
					continue;
				end

				% Set values of the variables for this line in 'counts'
				counts(line,1:5) = [rat,day,trial,configuration,group];

				% Initially set all measures to zero
				counts(line,WORKB:WORKU) = zeros(1,2);

				% No arm visited yet
				visitedBaitedArm1 = false;
				visitedBaitedArm2 = false;
				visitedBaitedArm3 = false;
				visitedUnbaitedArm1 = false;
				visitedUnbaitedArm2 = false;
				visitedUnbaitedArm3 = false;
				visitedUnbaitedArm4 = false;
				visitedUnbaitedArm5 = false;

				for visit = 1:nVisits,
					% Test successively visited arms to count the number of working memory errors, reference memory errors, etc.
					% In particular, check for repeated visits in baited arms (returning in an already visited baited arm is considered a working memory error)
					currentArm = block(visit,ARM);
					switch currentArm,
						case baitedArm1,
							if visitedBaitedArm1,
								% Returning in an already visited baited arm is considered as a working memory error
								counts(line,WORKB) = counts(line,WORKB)+1;
							else
								visitedBaitedArm1 = true;
							end
						case baitedArm2,
							if visitedBaitedArm2,
								counts(line,WORKB) = counts(line,WORKB)+1;
							else
								visitedBaitedArm2 = true;
							end
						case baitedArm3,
							if visitedBaitedArm3,
								counts(line,WORKB) = counts(line,WORKB)+1;
							else
								visitedBaitedArm3 = true;
							end
						case unbaitedArm1,
							if visitedUnbaitedArm1,
								counts(line,WORKU) = counts(line,WORKU)+1;
							else
								visitedUnbaitedArm1 = true;
							end
						case unbaitedArm2,
							if visitedUnbaitedArm2,
								counts(line,WORKU) = counts(line,WORKU)+1;
							else
								visitedUnbaitedArm2 = true;
							end
						case unbaitedArm3,
							if visitedUnbaitedArm3,
								counts(line,WORKU) = counts(line,WORKU)+1;
							else
								visitedUnbaitedArm3 = true;
							end
						case unbaitedArm4,
							if visitedUnbaitedArm4,
								counts(line,WORKU) = counts(line,WORKU)+1;
							else
								visitedUnbaitedArm4 = true;
							end
						case unbaitedArm5,
							if visitedUnbaitedArm5,
								counts(line,WORKU) = counts(line,WORKU)+1;
							else
								visitedUnbaitedArm5 = true;
							end
					end
				end
				line = line + 1;
			end
		end
	end
end

indices = ComputeIndices(counts);

[p,t,stats terms]=anovan(indices(:,WORKB),{indices(:,GROUP),indices(:,DAY)},'model','full','varnames',{'group','day'});
[p,t,stats terms]=anovan(indices(:,WORKU),{indices(:,GROUP),indices(:,DAY)},'model','full','varnames',{'group','day'});

nDays=15
colors = 'rbkgcmy';
figure;hold on;
for group = 1:3,
	d = indices(indices(:,GROUP)==group,[DAY WORKB]);
	if isempty(d), continue; end
	[m,s] = Accumulate(d(:,1),d(:,2));
	n = Accumulate(d(:,1),1);
	sem = s./sqrt(n);
	p = errorbar(m,sem,colors(group));
	set(p,'marker','.','markersize',24);
	ylabel('index WME Baited Arms');
	xlim([0 nDays+1]);set(gca,'xtick',1:nDays);
end

figure;hold on;
for group = 1:3,
	d = indices(indices(:,GROUP)==group,[DAY WORKU]);
	if isempty(d), continue; end
	[m,s] = Accumulate(d(:,1),d(:,2));
	n = Accumulate(d(:,1),1);
	sem = s./sqrt(n);
	p = errorbar(m,sem,colors(group));
	set(p,'marker','.','markersize',24);
	ylabel('index WME Unbaited Arms');
	xlim([0 nDays+1]);set(gca,'xtick',1:nDays);
end

function indices = ComputeIndices(counts)

WORKB = 6;
WORKU = 7;

indices = counts;
indices(:,[WORKB WORKU]) = sqrt(counts(:,[WORKB WORKU])+1);
%indices(:,[WORKB WORKU]) = counts(:,[WORKB WORKU]);