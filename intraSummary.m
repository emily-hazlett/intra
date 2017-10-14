%% Description
% This script is for creating and saving a summary plot for intracellular
% neural data.  Data should be saved as .txt files with each column being a
% different presentation of the stimulus and each row being a sample. .Txt
% file can contain header data before the traces begin, but the number of
% numeric headers must be specified.  This script will batch through all
% files matching a specified naming structure in a specifies folder and
% create and save a .tiff image of the summary figure.

% Created by EHazlett, last edited 2017-10-14

folderold = cd;
%% User editted info
% cd('D:\Intracellular\data analysis\processing\sub\'); % Look for files in this folder
cd('C:\Data Processing\Processing\sub\'); % Look for files in this folder
Files = dir('*BBN062_trace_1.txt'); % Find txt files containing this phrase to batch through

prestim = 300; %ms in sweep before stim onset
poststim = 700; %ms in sweep after stim onset
background = 200; %ms starting from end of poststim window to calculate background discharge
headers = 3; % number of rows containing numeric data in ascii file before the traces start

PSTH = 0; % Do you want to make a PSTH of the file?
binsize = 20; % Set binsize of PSTH

%% Batch through all files in the folder
for ii = 1:length(Files)
    %% Import data
    filename = Files(ii).name;
    testname = ['x', strrep(filename, '.txt','')];
    traces = importdata(filename);
    Reps.stim = strrep(traces.textdata{1},' ', '');
    Reps.trace = (traces.data(headers+1:end,:));
    clear traces
    
    %% Descriptives of test
    [samples, reps] = size(Reps.trace);
    rmp = mean2(Reps.trace(Reps.trace>mean2(Reps.trace) - std2(Reps.trace) ...
        & Reps.trace<mean2(Reps.trace) + std2(Reps.trace)));
    rmpSD = std2(Reps.trace);
    
    %Mean trace with limited spike contamination
    Reps.traceChop = Reps.trace;
    Reps.traceChop(Reps.traceChop > rmp + 4*rmpSD) = rmp + 4*rmpSD; %Chop spikes at 4 SD above rmp for calculating mean trace
    meanTrace = mean(Reps.traceChop,2);
    
    %% Spike detect
    Reps.spikeIndex = cell(1,reps);
    Reps.spikeHeight = cell(1,reps);
    threshold = min(rmp + 10*rmpSD, -10); % Set spike threshold to the lesser of rmp+10SD and =10mV
    for i = 1:reps
        if any(Reps.trace(:,i)>threshold)
            [Reps.spikeHeight{i}, Reps.spikeIndex{i}] = findpeaks(Reps.trace(:,i), ...
                'MinPeakHeight',threshold, ...
                'MinPeakDistance', ceil(1.5/(1000/samples)) ...
                ); %Find spike peaks that are <10 SD above rmp and with a hold time of ~1s (assuming a 1s sweep)
        end
    end
    
    %% Make a raster marking sample with spike peak
    spikeheight = []; %vector of all spike heights for this test
    Reps.raster = zeros(samples, reps);
    for i = 1:reps
        if isempty(Reps.spikeIndex{i}) == 0 % Dont run for reps that have no spiking
            for p = 1:length(Reps.spikeIndex{i})
                Reps.raster(Reps.spikeIndex{i}(p),i) = 1;
                spikeheight = cat(1, spikeheight, Reps.spikeHeight{i}(p));
            end
        end
    end

    %% Sound bar
    if any(strcmp(Reps.stim,{'None', 'NoneIR'})) % how long is the stim duration?
        stimlength = 0;
    elseif any(strcmp(Reps.stim,{'BBN062', 'LFH'}))
        stimlength = 62;elseif strcmp(Reps.stim, 'BatAgg')
        stimlength = 48;
    elseif any(strcmp(Reps.stim,{'Chevron', 'ChevronNL', 'FreqStp2', 'Complex'}))
        stimlength = 5;
    elseif strcmp(Reps.stim, 'Noisy')
        stimlength = 52;
    elseif strcmp(Reps.stim, 'MFVh')
        stimlength = 70;
    elseif strcmp(Reps.stim, 'MFVnl')
        stimlength = 90;
    elseif strcmp(Reps.stim, 'short')
        stimlength = 1;
    elseif any(strcmp(Reps.stim,{'Flat', 'DFM'}))
        stimlength = 2;
    elseif any(strcmp(Reps.stim,{'UFM', 'ChevronRev'}))
        stimlength = 3;
    elseif strcmp(Reps.stim, 'FreqStp1')
        stimlength = 6;
    elseif strcmp(Reps.stim, 'MFVt')
        stimlength = 13;
    else
        stimLength = 1;
    end
    barcolor = [0.5 0.5 0.5];
    
    %% Plot figure
    timeaxis = linspace(-prestim, poststim, samples);
    scaleAxis = [rmp - rmpSD, rmp + rmpSD];
    
    figure;
    set(gcf, 'Name', testname)
    set(gcf, 'Color', 'w')
    set(gcf, 'GraphicsSmoothing', 'on')
    set(gcf,'position', [0, 0, 500, 900])
    
    % Heatplot
    ax(1) = subplot(10,1,[1,2]);
    imagesc(Reps.trace')
    hold on
    for kk = 1:reps % Mark spikes on top of heatplot
        xpoint = repmat(kk,1,length(Reps.spikeIndex{kk}));
        scatter(Reps.spikeIndex{kk},xpoint,'k', 'LineWidth', 1.5)
        %'.',
    end
    hold off
    colorbar('off')
    title(testname,'Interpreter','none','Fontsize',12)
    ylabel('Reps')
    colormap(ax(1),'cool')
    ax(1).XTickLabels = [];
    ax(1).XTick = linspace(floor((prestim-100)/(1000/samples)), ...
        floor((poststim+prestim-100)/(1000/samples)), 5);
    ax(1).CLim = scaleAxis;
    ax(1).TickDir = 'out';
    ax(1).LineWidth = 1.5;
    ax(1).Box = 'off';
    
    %Mean trace compared to rmp
    ax(2) = subplot(10,1,[3,4]);
    plot(timeaxis, meanTrace, 'r', 'LineWidth', 0.5)
    hold on
    plot(timeaxis,repmat(rmp,1,samples),'k', 'LineWidth', 1.25)
    ylim([min([rmp - 2*rmpSD, min(meanTrace)]), ...
        max([rmp + 2*rmpSD, max(meanTrace)])])
    xlim([-prestim poststim])
    fill([0, stimlength, stimlength, 0], ...
        [max(ylim), max(ylim), min(ylim), min(ylim)], barcolor, ...
        'EdgeColor', 'none', 'FaceAlpha', 0.1)
    hold off
    ylabel('mV')
    ax(2).XTick = [];
    ax(2).TickDir = 'out';
    ax(2).LineWidth = 1.5;
    ax(2).Box = 'off';
    
    % Individual Traces
    tracePlotter = Reps.trace;
    ax(3) = subplot(10,1,[5,10]);
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
        ' BG = ', num2str(1000 * sum(sum(Reps.raster(samples - floor(background/(1000/samples)):end,:)))/(background*reps)), 'Hz', ...
        ' Spike height = ',num2str(round(mean(spikeheight),1)), ...
        ' +/- ',num2str(round(std(spikeheight),1)),'mV'], ...
        'Interpreter','none','Fontsize',9)
    ylabel('Individual Traces')
    xlabel('Time around stimulus onset(ms)')
    ax(3).YTick = [];
    ax(3).TickDir = 'out';
    ax(3).LineWidth = 1.5;
    ax(3).Box = 'off';
    clear tracePlotter
    
    % save figure
    print('-dtiff','-r500',[testname,'_summary.tif'])
    
    %% PSTH if requested
    if PSTH == 1
        pie = cd;
        cd('C:\Users\emily\OneDrive\Documents\GitHub\intra')
        run('intraPSTH')
        cd(pie)
        clear pie
    end
    clear samples reps Reps rmp* p ii kk meanTrace spikeheight ...
        stimlength timeaxis scaleAxis threshold xpoint testname
end
cd(folderold);