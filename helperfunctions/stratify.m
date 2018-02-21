%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% stratify
% by Til Ole Bergmann 2017
% last modified 2017/11/21 by TOB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% stratifies a set of n distributions and outputs the respective indices to keep
%
% input = array with observations (e.g., trials) in rows and parameters (e.g. trialinfo) in columns 
%
% putput = input array without rejected rows and with additonal last column indicating the row index of the values keept after stratification
%
% cfg = [];
% cfg.targetcol = integer indicating  the column of the array according to which stratification shall occour (default = 1)
% cfg.conditioncol = integer indicating  the column of the array according to which conditions shall be defined (default = 2)
% cfg.initialtrimtype = string indicating the type of trimming, 'none' (default),'absolute', 'percentile' 
% cfg.initialtrimval = numerical values indicating the upper and lower boundary for trimming [5 95]
% cfg.method = string indicating  the method for startification: 'removefromtop' (default)
% cfg.graphicoutput = string indicating whether to prodcude graphical putput, 'yes' or 'no' (default)
% cfg.stopcriteriontype = string indicating type of stop criterion for stratification process, 'absdev', 'reldev', 'anova' (default)
% cfg.stopcriterionval = integer indicating value for stop criterion of stratification process, e.g., 0.05 (default)
% 
% [output] = stratify(cfg,input);

function [output] = stratify(gcfg, input)


%% check and adjust config
if ~isfield(gcfg,'method'), gcfg.method = 'removefromtop'; end
if ~isfield(gcfg,'targetcol'), gcfg.targetcol = 1; end
if ~isfield(gcfg,'conditioncol'), gcfg.conditioncol = 2; end
if ~isfield(gcfg,'graphicoutput'), gcfg.graphicoutput = 'no'; end
if ~isfield(gcfg,'initialtrimtype'), gcfg.initialtrimtype = 'none'; end
if ~isfield(gcfg,'initialtrimval'), gcfg.initialtrimval = [1 5]; end
if ~isfield(gcfg,'stopcriteriontype'), gcfg.stopcriteriontype = 'anova'; end
if ~isfield(gcfg,'stopcriterionval'), gcfg.stopcriterionval = 0.2; end

output = [];
input = [input, (1:size(input,1))'];
conditions = unique(input(:,gcfg.conditioncol)); 
N = length(conditions);
for iCond = 1:N
    X{iCond} = input(input(:,gcfg.conditioncol)==conditions(iCond),:);
end
numtrial = cellfun(@length,X);

%% sort data and calculate indices
pX = [];
for iCond = 1:N
    [Xs{iCond},Idx{iCond}] = sortrows(X{iCond},gcfg.targetcol);
    pX = [pX;Xs{iCond}];
    Xsmean(iCond) = mean(Xs{iCond}(:,gcfg.targetcol),1);
    Xsmedian(iCond) = median(Xs{iCond}(:,gcfg.targetcol),1);
    Xsstd(iCond) = std(Xs{iCond}(:,gcfg.targetcol),1,1);
end

%% remove Inf and NaN values
for iCond = 1:N
    remIdx = find(isnan(Xs{iCond}(:,gcfg.targetcol)) | isinf(Xs{iCond}(:,gcfg.targetcol)));
    Xs{iCond}(remIdx,:) = [];
    Idx{iCond}(ismember(Idx{iCond},remIdx)) = [];
end

%% calcualte indices 
for iCond = 1:N    
    Xsmean(iCond) = mean(Xs{iCond}(:,gcfg.targetcol),1);
    Xsmedian(iCond) = median(Xs{iCond}(:,gcfg.targetcol),1);
    Xsstd(iCond) = std(Xs{iCond}(:,gcfg.targetcol),1,1);
end
gsmean = mean(Xsmean);
pstd = std(pX(:,gcfg.targetcol),1,1);
dev = mean(abs(Xsmean-gsmean));

% backup old values
oXs = Xs;
ogsmean = gsmean; 
oXsmean = Xsmean;
oXsmedian = Xsmedian;
opX = pX;
onumtrial = numtrial;

% run initial anova
[op, otbl] = anova1(opX(:,gcfg.targetcol),opX(:,gcfg.conditioncol),'off')

%% initial trimming
switch gcfg.initialtrimtype
    case 'none'
        gcfg.initialtrimvalabs = [0 Inf];
	case 'absolute'
        gcfg.initialtrimvalabs = gcfg.initialtrimval;
    case 'percentile'
        pXs = sortrows(pX,gcfg.targetcol);
        gcfg.initialtrimvalabs = [pXs(ceil(gcfg.initialtrimval(1)/100 * length(pXs))), pXs(floor(gcfg.initialtrimval(2)/100 * length(pXs)))];        
end
for iCond = 1:N
    remIdx = find(Xs{iCond}(:,gcfg.targetcol) <= gcfg.initialtrimvalabs(1) | Xs{iCond}(:,gcfg.targetcol) >= gcfg.initialtrimvalabs(2));
    Xs{iCond}(remIdx,:) = [];
    Idx{iCond}(ismember(Idx{iCond},remIdx)) = [];
end

% recalculate mean etc.
pX = [];
for iCond = 1:N    
    pX = [pX;Xs{iCond}];
    Xsmean(iCond) = mean(Xs{iCond}(:,gcfg.targetcol),1);
    Xsmedian(iCond) = median(Xs{iCond}(:,gcfg.targetcol),1);
    Xsstd(iCond) = std(Xs{iCond}(:,gcfg.targetcol),1,1);
end
gsmean = mean(Xsmean);
pstd = std(pX(:,gcfg.targetcol),1,1);
dev = mean(abs(Xsmean-gsmean));

%% stratification
iCount = 0;
p = []; testval = 1;
display(['Iteration ' num2str(iCount) ', deviation = ' num2str(dev) ', stop criterion = ' num2str(gcfg.stopcriterionval)]);
while testval > gcfg.stopcriterionval
    iCount = iCount + 1;
    
    switch gcfg.method
        case 'removefromtop'
            % remove trial with largest value from condition with larges deviation value
            [condVal condPos] = max((Xsmean-gsmean));
            remIdx = size(Xs{condPos},1);
            Xs{condPos}(remIdx,:) = [];
            Idx{condPos}(Idx{condPos} == remIdx) = [];
    end            
    
            % recalculate mean etc.
            Xsmean(condPos) = mean(Xs{condPos}(:,gcfg.targetcol),1);
            Xsmedian(condPos) = median(Xs{condPos}(:,gcfg.targetcol),1);
            Xsstd(condPos) = std(Xs{condPos}(:,gcfg.targetcol),1,1);
            gsmean = mean(Xsmean);
            dev = mean(abs(Xsmean-gsmean));
            
            % calcluate criterion
            switch gcfg.stopcriteriontype
                case 'reldev'                    
                    testval = dev;
                case 'anova'
                    pX = [];
                    for iCond = 1:N
                        pX = [pX;Xs{iCond}];
                    end                    
                    [p, tbl] = anova1(pX(:,gcfg.targetcol),pX(:,gcfg.conditioncol),'off');
                    testval = 1-p;
            end
            
            % display info
            display(['Iteration ' num2str(iCount) ', removed trial ' num2str(remIdx) ' from condition ' num2str(condPos) ' ,  deviation = ' num2str(dev) ', stop criterion = ' num2str(gcfg.stopcriterionval)]);

end
numtrial = cellfun(@length,Xs);

%% produce graphical output
if strcmp(gcfg.graphicoutput,'yes')
    h_stratresults = figure;
    subplot(2,1,1);
    title('before stratification');
    boxplot(opX(:,gcfg.targetcol),opX(:,gcfg.conditioncol));
    set(gca,'XTicklabel',num2cell(onumtrial));
    subplot(2,1,2);
    title('after stratification');
    boxplot(pX(:,gcfg.targetcol),pX(:,gcfg.conditioncol));
    set(gca,'XTicklabel',{1,2,3,4,5,6,7,8,9,10});
    set(gca,'XTicklabel',num2cell(numtrial));
end

output = sortrows(pX,size(pX,2));

end % of function

