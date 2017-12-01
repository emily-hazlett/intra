%% Description
% This script is for creating and saving a summary plot for intracellular
% neural data.  Data should be saved as .txt files with each column being a
% different presentation of the stimulus and each row being a sample. .Txt
% file can contain header data before the traces begin, but the number of
% numeric headers must be specified.  This script will batch through all
% files matching a specified naming structure in a specifies folder and
% create and save a .tiff image of the summary figure.

% Created by EHazlett, last edited 2017-11-06
clear all
close all

folderold = cd;
%% User editted info
% cd('D:\Intracellular\data analysis\processing\sub\'); % Look for files in this folder
cd('C:\Data Processing\Processing\sub\1195\'); % Look for files in this folder
Files = dir('1195_112817_3474_1_*_trace_*.txt'); % Find txt files containing this phrase to batch through

prestim = 100; %ms in sweep before stim onset. BACkground discharge calc from prestim window.
poststim = 900; %ms in sweep after stim onset
background = 50; %ms from beginning of sweep to calculate rmp
zone1 = [0, 50]; %Timezone1 = 0-50ms from sound onset
zone2 = [100, 500]; %Timezone2 = 100-500ms from sound onset

headers = 5; % number of rows containing numeric data in ascii file before the traces start
IRTest = 0; %Does this have IR test data?
save = 0; % save the figure as a tiff?

%% Batch through all files in the folder
for ii = 1:length(Files)
    %% Import data
    filename = Files(ii).name;
    testname = ['x', strrep(filename, '.txt','')];
    traces = importdata(filename);
    Reps.stim = regexprep(traces.textdata{1}, '\s+', '');
    if IRTest == 1
        Reps.preIR = traces.data(headers-1,:);
        Reps.postIR =traces.data(headers,:);
    else
        Reps.preIR = NaN(size(traces.data(3,:)),'double'); % Add in NaN for IR data
        Reps.postIR = NaN(size(traces.data(3,:)),'double');
    end
    Reps.trace = traces.data(headers+1:end,:);
    clear traces
    
    %% Descriptives of test
    [samples, reps] = size(Reps.trace);
    samplesPerMS = round(samples/(prestim+poststim));
    
    Reps.background = Reps.trace(1:floor(background*samples/1000),:); % Chunk of trace for calc bg
    for pp = 1:reps
        repper = Reps.background(:,pp);
        Reps.rmp(pp) = round(...
            mean2(repper( ...
            repper < mean2(repper) + std2(repper) )) ...
            ,1); % round to 1 decimal place
        clear repper
    end
    rmp = mean2(Reps.rmp);
    rmpSD = std2(Reps.rmp);
    
    %Mean trace with limited spike contamination
    Reps.traceChop = Reps.trace;
    Reps.traceChop(Reps.traceChop > rmp + 4*rmpSD) = rmp + 4*rmpSD; %Chop spikes at 4 SD above rmp for calculating mean trace
    meanTrace = mean(Reps.traceChop,2);
    
    %% Spike detect
    Reps.spikeIndex = cell(1,reps);
    Reps.spikeHeight = cell(1,reps);
    threshold = min(rmp + 30, -10); % Set spike threshold to the lesser of rmp+30 and =10mV
    for i = 1:reps
        if any(Reps.trace(:,i)>threshold)
            [Reps.spikeHeight{i}, Reps.spikeIndex{i}] = findpeaks(Reps.trace(:,i), ...
                'MinPeakHeight',threshold, ...
                'MinPeakDistance', ceil(1.5/(1000/samples)) ...
                ); %Find spike peaks that are <10 SD above rmp and with a hold time of ~1s (assuming a 1s sweep)
        end
    end
    
    %% Make a raster marking sample with spike peak, count spikes in timezones
    spikeheight = []; %vector of all spike heights for this test
    Reps.raster = zeros(samples, reps);
    Reps.zone1spikes = zeros(1, reps);
    Reps.zone2spikes = zeros(1, reps);
    for i = 1:reps
        if isempty(Reps.spikeIndex{i}) == 0 % Dont run for reps that have no spiking
            for p = 1:length(Reps.spikeIndex{i})
                Reps.raster(Reps.spikeIndex{i}(p),i) = 1;
                spikeheight = cat(1, spikeheight, Reps.spikeHeight{i}(p));
            end
            Reps.zone1spikes(i) = sum(Reps.raster((prestim+zone1(1))*samplesPerMS:(prestim+zone1(2))*samplesPerMS, i));
            Reps.zone2spikes(i) = sum(Reps.raster((prestim+zone2(1))*samplesPerMS:(prestim+zone2(2))*samplesPerMS, i));
        end
    end

    % Calculate background firing rate from last 400 ms of rep
    Reps.background = round( ...
        sum(sum(Reps.raster((poststim-400)*samplesPerMS - 10:end, :))) ...
        / ((poststim-400)*reps/1000)...
        ,1);
    
    %% Sound bar
    stimlength = 1;
    if any(strcmp(Reps.stim,{'None', 'NoneIR', 'Null'})) % how long is the stim duration?
        stimlength = 0;
    elseif any(strcmp(Reps.stim,{'BBN062', 'LFH', 'BBN'}))
        stimlength = 62;
    elseif strcmp(Reps.stim, 'BatAgg')
        stimlength = 48;
    elseif any(strcmp(Reps.stim,{'Chevron', 'ChevronNL', 'FreqStp2', 'Complex'}))
        stimlength = 5;
    elseif strcmp(Reps.stim, 'Noisy')
        stimlength = 52;
    elseif strcmp(Reps.stim, 'MFVh')
        stimlength = 70;
    elseif strcmp(Reps.stim, 'MFVnl')
        stimlength = 90;
    elseif strcmp(Reps.stim, 'Short')
        stimlength = 1;
    elseif any(strcmp(Reps.stim,{'Flat', 'DFM'}))
        stimlength = 2;
    elseif any(strcmp(Reps.stim,{'UFM', 'ChevronRev'}))
        stimlength = 3;
    elseif strcmp(Reps.stim, 'FreqStp1')
        stimlength = 6;
    elseif strcmp(Reps.stim, 'MFVt')
        stimlength = 13;
    end
    barcolor = [0.5 0.5 0.5];
    
    %% Plot figure
    timeaxis = linspace(-prestim, poststim, samples);
    scaleAxis = [rmp - 2*rmpSD, rmp + 2*rmpSD+0.0001];
    
    figure;
    set(gcf, 'Name', testname)
    set(gcf, 'Color', 'w')
    set(gcf, 'GraphicsSmoothing', 'on')
    set(gcf,'position', [0, 0, 500, 900])
    
    % IR
    ax(1) = subplot(20,1,[1,2]);
    plot(Reps.preIR, 'k', 'LineWidth', 1.25)
    hold on
    plot(Reps.postIR, 'r', 'LineWidth', 1.25)
    hold off
    axis tight
    title(testname,'Interpreter','none','Fontsize',12)

    ylabel('IR pre/post stim')
%     ax(1).XTickLabels = [];
    ax(1).TickDir = 'out';
    ax(1).LineWidth = 1.5;
    ax(1).Box = 'off';
    
    % Heatplot
    ax(2) = subplot(20,1,[4,7]);
    imagesc(Reps.trace')
    hold on
    for kk = 1:reps % Mark spikes on top of heatplot
        xpoint = repmat(kk,1,length(Reps.spikeIndex{kk}));
        scatter(Reps.spikeIndex{kk},xpoint,'k', 'LineWidth', 1.5)
        %'.',
    end
    hold off
    colorbar('off')
%     title(testname,'Interpreter','none','Fontsize',12)
    ylabel('Reps')
    colormap(ax(2),'cool')
    ax(2).XTickLabels = [];
    ax(2).XTick = linspace(floor((prestim-100)/(1000/samples)), ...
        floor((poststim+prestim-100)/(1000/samples)), 5);
    ax(2).CLim = scaleAxis;
    ax(2).TickDir = 'out';
    ax(2).LineWidth = 1.5;
    ax(2).Box = 'off';
    
    %Mean trace compared to rmp
    ax(3) = subplot(20,1,[8,10]);
    plot(timeaxis, meanTrace, 'r', 'LineWidth', 0.5)
    hold on
    plot(timeaxis,repmat(rmp,1,samples),'k', 'LineWidth', 1.25)
    ylim([min([rmp - 3*rmpSD, min(meanTrace)]), ...
        max([rmp + 3*rmpSD, max(meanTrace)])])
    xlim([-prestim poststim])
    fill([0, stimlength, stimlength, 0], ...
        [max(ylim), max(ylim), min(ylim), min(ylim)], barcolor, ...
        'EdgeColor', 'none', 'FaceAlpha', 0.1)
    hold off
    ylabel('mV')
    ax(3).XTick = [];
    ax(3).TickDir = 'out';
    ax(3).LineWidth = 1.5;
    ax(3).Box = 'off';
    
    % Individual Traces
    tracePlotter = Reps.trace;
    ax(4) = subplot(20,1,[12,20]);
    for kk = 1:reps
        tracePlotter = tracePlotter -range(tracePlotter(:,kk)) - 5;
        plot(timeaxis,tracePlotter(:,kk),'k')
        hold on
    end
    axis tight
    fill([0, stimlength, stimlength, 0], ...
        [max(ylim), max(ylim), min(ylim), min(ylim)], barcolor, ...
        'EdgeColor', 'none', 'FaceAlpha', 0.1)
    hold off
    title(['RMP = ', num2str(round(rmp,1)), 'mV', ...
        ' BG = ', num2str(Reps.background), 'Hz', ...
        ' Spike height = ',num2str(round(mean(spikeheight),1)), ...
        ' +/- ',num2str(round(std(spikeheight),1)),'mV'], ...
        'Interpreter','none','Fontsize',9)
    ylabel('Individual Traces')
    xlabel('Time around stimulus onset(ms)')
    ax(4).YTick = [];
    ax(4).TickDir = 'out';
    ax(4).LineWidth = 1.5;
    ax(4).Box = 'off';
    clear tracePlotter
    
    % save figure
    if save == 1
        print('-dtiff','-r500',[testname,'_summary.tif'])
    end
    
end
cd(folderold);