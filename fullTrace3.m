% Created by EHazlett, last edited 2018-06-07
% close all
% clear all
warning('off', 'all') % don't tell me that the fit equation is badly conditioned!

cd('C:\Data Processing\Processing\1239\'); % Look for files in this folder
Files = dir('*All_trace.txt'); % Find txt files containing this phrase to batch through

Fs = 500000;
badrmp = -40; % Drop reps above this threshold
driftrmpPOS = 3; % amount rmp is allowed to drift negative
driftrmpNEG = 3; % amount rmp is allowed to drift negative
smoother = 3; % How many reps to smooth over
headers = 3; % number of rows containing numeric data in ascii file before the traces start
loadmat = 1; %1 = load files from mat, 0= load files from txt
saver = 0; % save the txt file output?

%% Batch through all files in the folder
for i = 1:length(Files)
    %% Import data
    filename = Files(i).name;
    if loadmat == 0
        % stim files
        traces = importdata(filename);
        Reps.stim = strrep(traces.textdata(1,2:end),' ','');
        Reps.pulsepolarity = traces.data(1,:);
        Reps.pulsevoltage = traces.data(2,:);
        Reps.timestamp = traces.data(3,:);
        Reps.trace = (traces.data(headers+1:end,:))/10;
        Reps.todrop = false(1,size(Reps.trace, 2));
        stimtimestamp = Reps.timestamp;
        pulsepolarity = Reps.pulsepolarity;
        save('Reps.mat', '-struct', 'Reps')
        clear traces Reps
        
        % long trace
        largetracename = strrep(filename, 'All_trace', 'largetrace');
        largetrace = importdata(largetracename); %'C:\Users\emily\Desktop\importer\test4.txt'); %largetracename); %
        Sweeps.trace = largetrace(2:2:end, :)/10;
        Sweeps.timestamp = largetrace(1:2:end, :);
        Sweeps.rmp = mode(Sweeps.trace, 2);
        Sweeps.rmpSmooth = smoothdata(Sweeps.rmp, 'gaussian', 101);
        save('Sweeps.mat', 'Sweeps')
        clear largetrace
    elseif loadmat == 1
        load('Sweeps.mat');
        Sweeps.rmp = mode(Sweeps.trace, 2);
        Sweeps.rmpSmooth = smoothdata(Sweeps.rmp, 'gaussian', 101);
        s = load('Reps.mat', 'timestamp', 'pulsepolarity');
        stimtimestamp = s.timestamp;
        pulsepolarity = s.pulsepolarity;
        clear s
    end
    
    %% Spikes
    %spike detect
    nSweeps = size(Sweeps.trace, 1);
    spiker = Sweeps.trace;
    spiker(Sweeps.rmpSmooth>badrmp, :) = NaN;% don't spike detect reps with badrmp
    spiker = reshape(spiker', 1, numel(Sweeps.timestamp));
    rmp = mode(spiker);
    
    threshold = -20;
    [spikeheight, spikeindex] = findpeaks(spiker, ...
        'MinPeakHeight',threshold, ...
        'MinPeakDistance', ceil(1.5/((1000)/Fs)) ...
        ); %Find spike peaks that break threshold and with a hold time of ~1s
    spikeTimestamp = spikeindex*20+Sweeps.timestamp(1,1);
    
    %Adjust spike height for rmp
    spikeheightAdj = nan(1, length(spikeheight));
    for ii = 1:length(spikeheight)
        spikeheightAdj(ii) = spikeheight(ii) - ...
            Sweeps.rmp(find(Sweeps.timestamp(:,1)<spikeTimestamp(ii), 1, 'last'));
    end
    
    %% Fit line to segments of test based on current injection
    pulser = [1,find([0, diff(pulsepolarity)]),length(pulsepolarity)];
    chunkerRMP = cell(length(pulser)-1,2);
    outliersRMP = cell(length(pulser-1),1);
    coeffRMP = nan(length(pulser)-1,2);
    for ii = 1:length(pulser)-1
        finder(1) = find(Sweeps.timestamp(:,1)<stimtimestamp(pulser(ii)), 1, 'last');
        finder(2) = find(Sweeps.timestamp(:,1)<stimtimestamp(pulser(ii+1)), 1, 'last');
        chunkerRMP{ii,1} = Sweeps.rmpSmooth(finder(1):finder(2));
        chunkerRMP{ii,2} = Sweeps.timestamp(finder(1):finder(2),1);
        outliersRMP{ii} = excludedata(chunkerRMP{ii,2}, chunkerRMP{ii,1}, 'range', ...
            [mean(chunkerRMP{ii,1})-2*std(chunkerRMP{ii,1}), mean(chunkerRMP{ii,1})+2*std(chunkerRMP{ii,1})]);
        f = fit(chunkerRMP{ii,2},chunkerRMP{ii,1},'poly1', 'exclude', outliersRMP{ii});
        coeffRMP(ii,1:2) = coeffvalues(f);
        clear finder
    end
    
    %trajectory
    chunkerSpike = cell(length(pulser)-1,2);
    outliersSpike = cell(length(pulser-1),1);
    coeffSpike = nan(length(pulser)-1,2);
%     if spikeTimestamp(end) > stimtimestamp(end)
%         for ii = 1:length(pulser)-1
%             finder(1,1) = find(Sweeps.timestamp(:,1)<stimtimestamp(pulser(ii)), 1, 'last');
%             finder(2,1) = find(Sweeps.timestamp(:,1)<stimtimestamp(pulser(ii+1)), 1, 'last');
%             finder(1,2) = find(spikeTimestamp < Sweeps.timestamp(finder(1,1),1), 1, 'last');
%             finder(2,2) = find(spikeTimestamp < Sweeps.timestamp(finder(2,1),1), 1, 'last');
%             chunkerSpike{ii,1} = spikeheight(finder(1,2):finder(2,2))';
%             chunkerSpike{ii,2} = linspace(Sweeps.timestamp(finder(1,1),1),Sweeps.timestamp(finder(2,1),1), numel(chunkerSpike{ii,1}))';
%             outliersSpike{ii} = excludedata(chunkerSpike{ii,2}, chunkerSpike{ii,1}, 'range', ...
%                 [mean(chunkerSpike{ii,1})-2*std(chunkerSpike{ii,1}), mean(chunkerSpike{ii,1})+2*std(chunkerSpike{ii,1})]);
%             if size(chunkerSpike {ii,2}, 1) > 5
%                 f = fit(chunkerSpike{ii,2},chunkerSpike{ii,1},'poly1', 'exclude', outliersSpike{ii});
%                 coeffSpike(ii,1:2) = coeffvalues(f);
%             else
%                 coeffSpike(ii, 1:2) = 0;
%             end
%             clear finder
%         end
%     else
%     end
    
    %% Find center of stablest part of recording
    c = abs(Sweeps.rmpSmooth-mode(Sweeps.rmp(Sweeps.rmp<badrmp)))<2;
    count = 0;
    indexer = [1,length(c)];
    currentstreak = 0;
    maxstreak = 5;
    for ii = 1:length(c)
        if count+c(ii) > count
            currentstreak = currentstreak + 1;
            if currentstreak > maxstreak
                indexer(1,1)=max(ii-currentstreak,1);
                indexer(1,2)=ii;
                maxstreak = currentstreak;
            end
        else
            count = 0;
            currentstreak = 0;
        end
    end
    centerIndex = floor(mean(indexer));
    
    %% RmP
    % penetration
    changer = [false; diff(Sweeps.rmpSmooth) < -1]; % change in mV that counts as change in smoothed rmp
    count = 0;
    penetrationTimestamp = Sweeps.timestamp(1,1);
    penetrationIndex = 1;
    for ii = centerIndex:1
        count = (count+changer(ii))*changer(ii);
        if count > 3 % need 10 reps in a row with increasing rmp
            penetrationTimestamp = Sweeps.timestamp(ii+3,1);
            penetrationIndex = ii+3;
            break
        end
    end
    
    % drift down
    changer = [false; diff(Sweeps.rmpSmooth) < -0.05]; % change in mV that counts as change in smoothed rmp
    count = 0;
    rmpdriftdownTimestamp = NaN;
    rmpdriftdownIndex = NaN;
    for ii = penetrationIndex:centerIndex
        if any(ii-find(diff(pulsepolarity)) <15) && any(ii-find(diff(pulsepolarity)) > -5)
            continue
        end
        chunkadunk = Sweeps.rmpSmooth(i:i+10)-mode(Sweeps.rmpSmooth(i:i+10));
        count = count+any(chunkadunk>(ii))*changer(ii);
        if count > 20 % need 10 reps in a row with increasing rmp
            rmpdriftdownTimestamp = Sweeps.timestamp(ii+9,1);
            rmpdriftdownIndex = ii+9;
            break
        end
    end
    
    % lose cell
    changer = [false; diff(Sweeps.rmpSmooth) > 1]; % change in mV that counts as change in smoothed rmp
    count = 0;
    lostTimestamp = Sweeps.timestamp(end,end);
    lostIndex = length(changer);
    for ii = centerIndex:length(changer)
        count = (count+changer(ii))*changer(ii);
        if count > 3 % need 10 reps in a row with increasing rmp
            lostTimestamp = Sweeps.timestamp(ii-3,1);
            lostIndex = ii-3;
            break
        end
    end
    
    % drift up
    changer = [false; diff(Sweeps.rmpSmooth) > 0.05]; % change in mV that counts as change in smoothed rmp
    count = 0;
    rmpdriftupTimestamp = NaN;
    rmpdriftupIndex = NaN;
    for ii = centerIndex:lostIndex
        if any(ii-find(diff(pulsepolarity)) <15) && any(ii-find(diff(pulsepolarity)) > -5)
            continue
        end
        count = (count+changer(ii))*changer(ii);
        if count > 20 % need 10 reps in a row with increasing rmp
            rmpdriftupTimestamp = Sweeps.timestamp(ii-9,1);
            rmpdriftupIndex = ii-9;
            break
        end
    end
    
    % drifting spike height
    changer = [false, diff(movmean(spikeheightAdj,5)) < -0.1]; % change in mV that counts as change in spike height
    changer(spikeheight>0) = false; % don't count decreases in overshoot spikes
    count = 0;
    spikedriftTimestamp = NaN;
    for ii = 1:length(changer)
        count = (count+changer(ii))*changer(ii);
        if count > 4 % need 5 reps in a row with decreasing spike height
            spikedriftTimestamp = spikeindex(ii-4)*20+Sweeps.timestamp(1,1);
            break
        end
    end
    
    %% plot
    figure
    set(gcf,'position', [0, 0, 1200, 900])
    % trace
    plot(reshape(Sweeps.timestamp', 1, numel(Sweeps.timestamp)), reshape(Sweeps.trace', 1, numel(Sweeps.trace)), ...
        'color', [0.5, 0.5, 0.5])
    hold on
    % rmp as it varies during the recoridng
    plot(Sweeps.timestamp(:,1)+mean(Sweeps.timestamp(1,1)-Sweeps.timestamp(2,1)), ...
        Sweeps.rmp, 'k', 'LineWidth', 2)
    plot(Sweeps.timestamp(:,1)+mean(Sweeps.timestamp(1,1)-Sweeps.timestamp(2,1)), ...
        Sweeps.rmpSmooth, 'r', 'LineWidth', 1.5)
    % flattest part of recording
    line(Sweeps.timestamp(indexer), repmat(Sweeps.rmp(centerIndex)-5, 1, 2), 'color', 'y')
    scatter(Sweeps.timestamp(centerIndex,1),Sweeps.rmp(centerIndex)-5, 100, 'd', 'y', 'filled', 'MarkerEdgeColor', 'k')
    % fit lines for pulses
    for ii = 1:length(pulser)-1
        line([chunkerRMP{ii,2}(1),chunkerRMP{ii,2}(end)], ...
            [chunkerRMP{ii,2}(1)*coeffRMP(ii,1)+coeffRMP(ii,2),chunkerRMP{ii,2}(end)*coeffRMP(ii,1)+coeffRMP(ii,2)], 'color', 'y')
    end
    % mean rmp +/- range
    %     plot(Sweeps.timestamp(:,1)+mean(Sweeps.timestamp(1,1)-Sweeps.timestamp(2,1)), ...
    %         repmat(rmp, 1, length(Sweeps.rmp)), 'k', 'LineWidth', 1)
    %     plot(Sweeps.timestamp(:,1)+mean(Sweeps.timestamp(1,1)-Sweeps.timestamp(2,1)), ...
    %         repmat(rmp+driftrmpPOS, 1, length(Sweeps.rmp)), 'k--', 'LineWidth', 1)
    %     plot(Sweeps.timestamp(:,1)+mean(Sweeps.timestamp(1,1)-Sweeps.timestamp(2,1)), ...
    %         repmat(rmp-driftrmpNEG, 1, length(Sweeps.rmp)), 'k--', 'LineWidth', 1)
    % spike markers
    scatter(spikeTimestamp, spikeheight, '.', 'k')
%     if spikeTimestamp(end) > stimtimestamp(end)
%         for ii = 1:length(pulser)-1
%             line([chunkerSpike{ii,2}(1),chunkerSpike{ii,2}(end)], ...
%                 [chunkerSpike{ii,2}(1)*coeffSpike(ii,1)+coeffSpike(ii,2),chunkerSpike{ii,2}(end)*coeffSpike(ii,1)+coeffSpike(ii,2)], 'color', 'y', 'linewidth', 2)
%         end
%     end
    % stim markers
    scatter(stimtimestamp, repmat(-20, 1, numel(stimtimestamp)), 'filled')
    plot(stimtimestamp,pulsepolarity*5-90)
    % begining of good recording
    line(repmat(max([penetrationTimestamp, rmpdriftdownTimestamp]),1,2), ...
        [-100, 30], 'color','g')
    % end of good recording
    line(repmat(min([lostTimestamp, rmpdriftupTimestamp, spikedriftTimestamp]),1,2), ...
        [-100, 30], 'color','r')
    % other values
    %     scatter(spikedriftTimestamp, 0, 25, 'd', 'r', 'filled')
    %     scatter(rmpdriftdownTimestamp, 0, 25, 'r')
    %     scatter(rmpdriftdownTimestamp, 0, 25, 'r', 'filled')
    %     scatter(penetrationTimestamp, 0, 25, 'y', 'filled')
    %     scatter(lostTimestamp, 0, 25, 'r', 'filled')
    axis tight
    ylim([-100, 30])
    savefig(gcf, num2str(i))
    saveas(gcf, [num2str(i), '.jpg'], 'jpeg')
        close gcf
    
    if saver == 1
        %% Drop Reps and save cleaned .txt file
        dropstimTimestamp = min([spikedriftTimestamp, rmpdriftTimestamp, Sweeps.timestamp(end,end)]);
        Reps.todrop(Reps.timestamp > dropstimTimestamp) = true;
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
    end
    clear chunker* coeff* *Timestamp *Index spike* Sweeps pulser rmp
end


%     % bad rmp (can't do penetration, only losing cell)
%     dropper1 = find(Sweeps.rmp > badrmp);
%     if ~isempty(dropper1)
%         dropIndexLong = find(Long.timestamp == Sweeps.timestamp(dropper1(1), 1));
%         dropper1(1, dropIndexLong:end) = false;
%     else
%         dropIndexLong = length(Long.timestamp);
%     end
%     rmp = mode(Long.trace(1:dropIndexLong));