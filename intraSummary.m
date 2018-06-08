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
cd('C:\Data Processing\Processing\'); % Look for files in this folder
Files = dir('1215*BatAgg*trace_*.txt'); % Find txt files containing this phrase to batch through

prestim = 100; %ms in sweep before stim onset
poststim = 900; %ms in sweep after stim onset
background = 100; %ms starting from end of poststim window to calculate background discharge
headers = 3; % number of rows containing numeric data in ascii file before the traces start

binsize = 20; % Set binsize of PSTH
slide = 5; % ms of sliding window for psth

%% Batch through all files in the folder
for ii = 1:length(Files)
    %% Import data
    filename = Files(ii).name;
    testname = ['x', strrep(filename, '.txt','')];
    traces = importdata(filename);
    Reps.stim = regexprep(traces.textdata{1}, '\s+', '');
    Reps.trace = (traces.data(headers+1:end,:));
    clear traces
    [samples, reps] = size(Reps.trace);
   
    %% Spike detect
    Reps.spikeIndex = cell(1,reps);
    Reps.spikeHeight = cell(1,reps);
    threshold = mode(mode(round(Reps.trace,1)))+30; % Set spike threshold to the lesser of rmp+10SD and -20mV
    for i = 1:reps
        if any(Reps.trace(:,i)>threshold)
            [Reps.spikeHeight{i}, Reps.spikeIndex{i}] = findpeaks(Reps.trace(:,i), ...
                'MinPeakHeight',threshold, ...
                'MinPeakDistance', ceil(1.5/((prestim+poststim)/samples)) ...
                ); %Find spike peaks that are <10 SD above rmp and with a hold time of ~1s
        end
    end
    
    %% Make a raster marking sample with spike peak and despike traces
    spikeheight = []; %vector of all spike heights for this test
    Reps.raster = zeros(samples, reps);
    spikelength = ceil(1.5/((prestim+poststim)/samples)); % remove 3ms snippet around spike peak fromm trace
    spikeremover = false(samples, reps);
    
    for i = 1:reps
        if ~isempty(Reps.spikeIndex{i}) % Dont run for reps that have no spiking
            for p = 1:length(Reps.spikeIndex{i})
                Reps.raster(Reps.spikeIndex{i}(p),i) = 1;
                spikeremover(Reps.spikeIndex{i}(p)-spikelength:Reps.spikeIndex{i}(p)+spikelength,i) = true;
                spikeheight = cat(1, spikeheight, Reps.spikeHeight{i}(p));
            end
        end
    end
        
    Reps.traceDespiked = Reps.trace;
    Reps.traceDespiked(spikeremover) = NaN;
    meanTrace = nanmean(Reps.traceDespiked,2);
    clear spikeremover
    
    Reps.rmp = mode(round(Reps.traceDespiked,1));
    
    % Get test RMP from
    p = histogram(Reps.trace);
    [~, i] = max(movmean(p.Values,20));
    valuers = p.BinLimits(1):p.BinWidth:p.BinLimits(2)-p.BinWidth;
    rmp = valuers(i);
    close p
    
    %% PSTH with optional sliding window
    psth = sum(Reps.raster, 2)'; %"raster" created by intraSummary
    timePerSample = (prestim+poststim)/ length(psth); %Assumes 1s sweeplength
    binners = floor(binsize/timePerSample);
    sliders = floor(slide/timePerSample);
    
    bin = 0;
    for p = floor(binners/2):sliders:samples-floor(binners/2)
        bin = bin + 1;
        psthBinSlide (bin) = sum(psth(p-floor(binners/2)+1:p+floor(binners/2)));
    end
    clear p bin
    
    %% Sound bar
    stimlength = 48;
    marks = Reps.stim(1);
    switch marks
        case {'None', 'NoneIR'}
            stimlength = 0;
        case {'BBN062', 'LFH', 'BBN'}
            stimlength = 62;
        case 'BatAgg'
            stimlength = 48;
        case {'Chevron', 'ChevronNL', 'FreqStp2', 'Complex'}
            stimlength = 5;
        case 'Noisy'
            stimlength = 52;
        case 'MFVh'
            stimlength = 70;
        case 'MFVnl'
            stimlength = 90;
        case 'Short'
            stimlength = 1;
        case {'Flat', 'DFM'}
            stimlength = 2;
        case {'UFM', 'ChevronRev'}
            stimlength = 3;
        case 'FreqStp1'
            stimlength = 6;
        case 'MFVt'
            stimlength = 13;
    end
    barcolor = [0.5 0.5 0.5];
    
    %% Plot figure
    timeaxis = linspace(-prestim, poststim, samples);
    scaleAxis = [rmp - 5, rmp + 5];
    
    figure;
    set(gcf, 'Name', testname)
    set(gcf, 'Color', 'w')
    set(gcf, 'GraphicsSmoothing', 'on')
    set(gcf,'position', [0, 0, 500, 900])
    
    
    % PSTH
    ax(1) = subplot(12,1,[1,2]);
    plot(linspace(-prestim, poststim, length(psthBinSlide)),psthBinSlide, 'k')
    hold on
    fill([0, stimlength, stimlength, 0], ...
        [max(ylim), max(ylim), min(ylim), min(ylim)], barcolor, ...
        'EdgeColor', 'none', 'FaceAlpha', 0.1)
    title(testname,'Interpreter','none','Fontsize',12) 
    xlim([-prestim, poststim])
    ylabel(['spikes/ ', num2str(binsize), 'ms bin'])
    ax(1).TickDir = 'out';
    ax(1).XTickLabels = [];
    ax(1).LineWidth = 1.5;
    ax(1).Box = 'off';
    
    %Mean trace compared to rmp
    ax(2) = subplot(12,1,[3,4]);
    hold on
%     for kk = 1:reps
%         plot(timeaxis,Reps.trace(:,kk),'Color',(kk/reps)*[0.5 0.5 0.5])
%         hold on
%     end
    plot(timeaxis,repmat(rmp,1,samples),'Color',[0 0 0], 'LineWidth', 1.25)
    plot(timeaxis, meanTrace, 'r', 'LineWidth', 0.5)
    ylim([rmp-5, rmp+5])
    xlim([-prestim poststim])
    fill([0, stimlength, stimlength, 0], ...
        [max(ylim), max(ylim), min(ylim), min(ylim)], barcolor, ...
        'EdgeColor', 'none', 'FaceAlpha', 0.1)
    hold off
    title(['RMP = ', num2str(round(rmp,1)), 'mV', ...
        ' BG = ', num2str(round(1000 * sum(sum(Reps.raster(samples - floor(background/(1000/samples)):end,:)))/(background*reps)),1), 'Hz', ...
        ' Spike height = ',num2str(round(mean(spikeheight),1)), ...
        ' +/- ',num2str(round(std(spikeheight),1)),'mV'], ...
        'Interpreter','none','Fontsize',9)
    ylabel('mV')
    ax(2).TickDir = 'out';
    ax(2).XTickLabels = [];
    ax(2).LineWidth = 1.5;
    ax(2).Box = 'off';
    
    % Heatplot
    ax(3) = subplot(12,1,[5,6]);
    imagesc(Reps.trace')
    hold on
%     for kk = 1:reps % Mark spikes on top of heatplot
%         xpoint = repmat(kk,1,length(Reps.spikeIndex{kk}));
%         scatter(Reps.spikeIndex{kk},xpoint,'k', 'LineWidth', 1.5)
%         %'.',
%     end
    hold off
    colorbar('off')
    ylabel('Reps')
    colormap(ax(3),'cool')
    ax(3).XTickLabels = [];
    ax(3).XTick = linspace(0, samples, 7);
    ax(3).CLim = scaleAxis;
    ax(3).TickDir = 'out';
    ax(3).LineWidth = 1.5;
    ax(3).Box = 'off';
    
    % Individual Traces
    tracePlotter = Reps.trace;
    ax(4) = subplot(12,1,[7,12]);
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
    
    ylabel('Individual Traces')
    xlabel('Time around stimulus onset(ms)')
    ax(4).YTick = [];
    ax(4).TickDir = 'out';
    ax(4).LineWidth = 1.5;
    ax(4).Box = 'off';
    clear tracePlotter
    
    % save figure
    print('-dtiff','-r500',[testname,'_summary.tif'])
    clear samples reps Reps rmp* p kk meanTrace spikeheight ...
        stimlength timeaxis scaleAxis threshold xpoint testname
    disp([num2str(ii),'/', num2str(length(Files))]);
end
cd(folderold);