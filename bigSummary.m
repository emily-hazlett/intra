%% Description
% This script is for automatically dropping reps that occur before
% penetration of a cell or after seal degrades, and then saving a cleaned
% .txt file without bad reps. Data should be saved as .txt files with each column
% being adifferent presentation of the stimulus and each row being a sample.
% .Txt file can contain header data before the traces begin, but the number
% of numeric headers must be specified.  This script will batch through all
% files matching a specified naming structure in a specifies folder and
% create and save a .tiff image of the summary figure.

% Created by EHazlett, last edited 2018-06-07
close all
clear all

folderold = cd;
%% User editted info
cd('C:\Data Processing\Processing\1239\'); % Look for files in this folder
Files = dir('*_All_trace.txt'); % Find txt files containing this phrase to batch through

badrmp = -20; % Drop reps above this threshold
driftrmpPOS = 35; % amount rmp is allowed to drift negative
driftrmpNEG = 15; % amount rmp is allowed to drift negative
minspike = 40; % minimum spike height allowed
smoother = 3; % How many reps to smooth over
headers = 3; % number of rows containing numeric data in ascii file before the traces start

%% Batch through all files in the folder
for ii = 1:length(Files)
    %% Import data
    filename = Files(ii).name;
    traces = importdata(filename);
    Reps.stim = strrep(traces.textdata(1,2:end),' ','');
    Reps.pulsepolarity = traces.data(1,:);
    Reps.pulsevoltage = traces.data(2,:);
    Reps.timestamp = traces.data(3,:);
    Reps.trace = (traces.data(headers+1:end,:))/10;
    clear traces
    [samples, reps] = size(Reps.trace);
    if reps < 2
        continue
    end
    Reps.todrop = false(1,reps);
    
    %% Drop reps if not sufficiently negative
    Reps.rmp = mode(round(Reps.trace,1));
    Reps.todrop(Reps.rmp>badrmp) = true;
    
    %% Drop reps if spike height isn't good enough
    threshold = max([min(mode(Reps.rmp(~Reps.todrop)) + 20),-35]); % Set spike threshold
    Reps.meanspikeheight = nan(1,reps);
    for i = 1:reps
        if any(Reps.trace(:,i)>threshold)
            [spikeheight, ~] = findpeaks(Reps.trace(:,i), ...
                'MinPeakHeight',threshold, ...
                'MinPeakDistance', ceil(1.5/((1000)/samples)) ...
                ); %Find spike peaks that break threshold and with a hold time of ~1s
            Reps.meanspikeheight (i) = mean(spikeheight - Reps.rmp(i));
            if Reps.meanspikeheight(i) < minspike
                Reps.todrop(i) = true;
            end
            clear spikeheight
        end
    end
    
    % Drop reps if there's too much of a decrease in spike height (don't
    % drop if overshoot)
    spikeheighter = find(movmean(Reps.meanspikeheight,smoother)<(max(Reps.meanspikeheight)-20) ...
        & (movmean(Reps.meanspikeheight + Reps.rmp,smoother)) < 0);
    if ~isempty(spikeheighter)
        if spikeheighter(1) > find(max(Reps.meanspikeheight))
            Reps.todrop(spikeheighter(1):end) = true;
        end
        if spikeheighter(1) < find(max(Reps.meanspikeheight))
            Reps.todrop(1:spikeheighter(end)) = true;
        end
    end
    clear spikeheighter
    
    %% Drop reps if the rmp drifts
    % Don't include reps already marked to drop in overall measures of rmp
    p = histogram(Reps.trace(:,Reps.pulsepolarity==-1 & ~Reps.todrop));
    [~, i] = max(movmean(p.Values,20));
    valuers = p.BinLimits(1):p.BinWidth:p.BinLimits(2)-p.BinWidth;
    rmp(1) = valuers(i);
    clear p valuers i
    close gcf
    
    p = histogram(Reps.trace(:,Reps.pulsepolarity==0 & ~Reps.todrop));
    [~, i] = max(movmean(p.Values,20));
    valuers = p.BinLimits(1):p.BinWidth:p.BinLimits(2)-p.BinWidth;
    rmp(2) = valuers(i);
    clear p valuers i
    close gcf
    
    p = histogram(Reps.trace(:,Reps.pulsepolarity==1 & ~Reps.todrop));
    [~, i] = max(movmean(p.Values,20));
    valuers = p.BinLimits(1):p.BinWidth:p.BinLimits(2)-p.BinWidth;
    rmp(3) = valuers(i);
    clear p valuers i
    close gcf
    rmp(rmp>badrmp) = NaN;
    
    Reps.todrop(Reps.pulsepolarity==-1 & (Reps.rmp < (rmp(1)-driftrmpNEG) | Reps.rmp > (rmp(1)+driftrmpPOS))) = true; % hyperpolarizing current reps
    Reps.todrop(Reps.pulsepolarity==0 & (Reps.rmp < (rmp(2)-driftrmpNEG) | Reps.rmp > (rmp(2)+driftrmpPOS))) = true; % no current reps
    Reps.todrop(Reps.pulsepolarity==1 & (Reps.rmp < (rmp(3)-driftrmpNEG) | Reps.rmp > (rmp(3)+driftrmpPOS))) = true; % depolarizing current reps
    
    %% Clean droprep index
    % At least two consecutive reps needed to really drop reps
    Reps.todrop(strfind(Reps.todrop, [0 1 0])+1) = false;
    
    % If reps are marked not to drop but are sufficiently within a string
    % of dropped reps, then drop them anyway
    Reps.todrop(strfind(Reps.todrop, [1 0 0 0 1 1 1 1 1])+1) = true;
    Reps.todrop(strfind(Reps.todrop, [1 1 0 0 0 1 1 1 1])+2) = true;
    Reps.todrop(strfind(Reps.todrop, [1 1 1 0 0 0 1 1 1])+3) = true;
    Reps.todrop(strfind(Reps.todrop, [1 1 1 1 0 0 0 1 1])+4) = true;
    Reps.todrop(strfind(Reps.todrop, [1 1 1 1 1 0 0 0 1])+5) = true;
    
    Reps.todrop(strfind(Reps.todrop, [1 0 0 1 1 1])+1) = true;
    Reps.todrop(strfind(Reps.todrop, [1 1 0 0 1 1])+2) = true;
    Reps.todrop(strfind(Reps.todrop, [1 1 1 0 0 1])+3) = true;
    
    Reps.todrop(strfind(Reps.todrop, [1 0 1])+1) = true;
    
    % Only allow the first chunk of workable reps to remain
    cleaner = diff(Reps.todrop);
    mark = strfind(cleaner(cleaner~=0),[1, -1]);
    if ~isempty(mark)
        cleaner = find(diff(Reps.todrop)~=0);
        Reps.todrop(cleaner(mark(1)):end) = true;
        if ~any(Reps.pulsepolarity(~Reps.todrop)==-1)
            rmp(1) = NaN;
        end
        if ~any(Reps.pulsepolarity(~Reps.todrop)==0)
            rmp(2) = NaN;
        end
        if ~any(Reps.pulsepolarity(~Reps.todrop)==1)
            rmp(3) = NaN;
        end
    end
    
    %% Plot and save figure
    figure
    set(gcf,'position', [0, 0, 900, 900])
    ax(1) = subplot(3,1,1);
    histogram(Reps.rmp(~Reps.todrop & Reps.pulsepolarity==-1), 'BinWidth', 1, 'FaceColor', 'b')
    hold on
    histogram(Reps.rmp(~Reps.todrop & Reps.pulsepolarity==1), 'BinWidth', 1, 'FaceColor', 'r')
    histogram(Reps.rmp(~Reps.todrop & Reps.pulsepolarity==0), 'BinWidth', 1, 'FaceColor', 'k')
    title(filename,'Interpreter','none','Fontsize',12)
    xlim([-90 badrmp])
    ylabel('Count')
    xlabel('RMP')
    
    ax(2) = subplot(3,1,2:3);
    constant = ceil((Reps.timestamp(2)-Reps.timestamp(1))/samples);
    for i = 1:reps
        plot(Reps.timestamp(i)/constant:Reps.timestamp(i)/constant+samples-1, Reps.trace(:,i))
        hold on
    end
    ylim([-100 20])
    timeline = (Reps.timestamp/constant)+samples/2;
    plot(timeline, Reps.rmp, 'k', 'Linewidth', 2)
    plot(timeline, repmat(rmp(1) - driftrmpNEG, reps), 'b', 'Linewidth', 1)
    plot(timeline, repmat(rmp(1) , reps), 'b--', 'Linewidth', 1)
    plot(timeline, repmat(rmp(1) + driftrmpPOS, reps), 'b', 'Linewidth', 1)
    plot(timeline, repmat(rmp(2) - driftrmpNEG, reps), 'k', 'Linewidth', 1)
    plot(timeline, repmat(rmp(2) , reps), 'k--', 'Linewidth', 1)
    plot(timeline, repmat(rmp(2) + driftrmpPOS, reps), 'k', 'Linewidth', 1)
    plot(timeline, repmat(rmp(3) - driftrmpNEG, reps), 'r', 'Linewidth', 1)
    plot(timeline, repmat(rmp(3) , reps), 'r--', 'Linewidth', 1)
    plot(timeline, repmat(rmp(3) + driftrmpPOS, reps), 'r', 'Linewidth', 1)
    
    plot(timeline, Reps.todrop * 18, 'g', 'Linewidth', 2)
    plot(timeline, Reps.pulsepolarity * 5 - 90, 'r', 'Linewidth', 2)
    title(['RMP(o) = ', num2str(round(rmp(2),1)), 'mV', ...
        '  RMP(-) = ', num2str(round(rmp(1),1)), 'mV', ...
        '  RMP(+) = ', num2str(round(rmp(3),1)), 'mV'], ...
        'Interpreter','none','Fontsize',9)
    clear timeline
    
    testname = ['x', strrep(filename, '.txt','')];
    print('-dtiff','-r500',[testname,'_dropreps.tif'])
    
    %% Drop Reps and save cleaned .txt file
    if sum(~Reps.todrop)< 3
        continue
    end
    Reps.stim(:,Reps.todrop) = [];
    Reps.trace(:,Reps.todrop) = [];
    Reps.pulsepolarity(:,Reps.todrop) = [];
    Reps.pulsevoltage(:,Reps.todrop) = [];
    Reps.timestamp(:,Reps.todrop) = [];
    
    outputter = strrep(filename, 'trace', 'trace_cleaned');
    fid = fopen(outputter, 'w');
    fprintf(fid,[repmat('%s \t',1,length(Reps.stim)-1), '%s'],Reps.stim{:});
    fprintf(fid,'\n');
    fclose(fid);
    dlmwrite(outputter, [Reps.pulsepolarity; Reps.pulsevoltage; Reps.timestamp; Reps.trace], ...
        'precision','%f', '-append', 'delimiter', '\t','roffset', 0);
    disp([num2str(ii), '/', num2str(length(Files))]);
    clear Reps reps pp ii samples rmp
end
clear headers *stim background
cd(folderold);