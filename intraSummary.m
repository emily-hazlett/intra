%% Description
% This script is for creating and saving a summary plot for intracellular
% neural data.  Data should be saved as .txt files with each column being a
% different presentation of the stimulus and each row being a sample. .Txt
% file can contain header data before the traces begin, but the number of
% numeric headers must be specified.  This script will batch through all
% files matching a specified naming structure in a specifies folder and
% create and save a .tiff image of the summary figure.

% Created by EHazlett, last edited 2017-07-14


folderold = cd;
%% User editted info
cd('D:\Intracellular\data analysis\processing\'); % Look for files in this folder
Files = dir('*trace*.txt'); % Find txt files containing this phrase to batch through

prestim = 100; %ms in sweep before stim onset
poststim = 900; %ms in sweep after stim onset
headers = 3; % number of rows containing numeric data in ascii file before the traces start

%% Batch through all files in the folder
for ii = 1:length(Files)
    %% Import data
    filename = Files(ii).name;
    testname = ['x', strrep(filename, '.txt','')];
    traces = importdata(filename);
    traces = (traces.data(headers+1:end,:))/10;
    
    %% Descriptives of traces
    [samples, reps] = size(traces);
    rmp = mean2(traces(traces>mean2(traces) - std2(traces) & traces<mean2(traces) + std2(traces)));
    tracesSD = std2(traces);
    
    spikeIndex = cell(1,reps);
    spikeHeight = cell(1,reps);
    threshold = min(rmp + 10*tracesSD, -10); % Set spike threshold to the lesser of rmp+10SD and =10mV
    for i = 1:reps
        if any(traces(:,i)>threshold)
            [spikeHeight{i}, spikeIndex{i}] = findpeaks(traces(:,i), ...
                'MinPeakHeight',threshold, ...
                'MinPeakDistance', ceil(samples/1000) ...
                ); %Find spike peaks that are <10 SD above rmp and with a hold time of ~1s (assuming a 1s sweep)
            
        end
    end
    
    %Make a raster marking sample with spike peak
    raster = zeros(samples, reps);
    spikeheight = []; %vector of all spike heights for this test
    for i = 1:reps
        if isempty(spikeIndex{i}) == 0 % Dont run for reps that have no spiking
            for p = 1:length(spikeIndex{i})
                raster(spikeIndex{i}(p),i) = 1;
                spikeheight = cat(1, spikeheight, spikeHeight{i}(p));
            end
        end
    end
    spikeheightM = mean(spikeheight);
    spikeheightSD = std(spikeheight);
    
    %Mean trace with limited spike contamination
    tracesChop = traces;
    tracesChop(tracesChop>rmp + 4*tracesSD) = rmp + 4*tracesSD; %Chop spikes at 4 SD above rmp for calculating mean trace
    meanTrace = mean(tracesChop,2);
    
    %% Plot figure
    timestep = (prestim+poststim)/ samples;
    timeaxis = timestep-prestim:timestep:poststim;
    scaleAxis = [floor(rmp - 1.5*tracesSD), ceil(rmp + 1.5*tracesSD)];
    
    figure;
    set(gcf, 'Name', testname)
    set(gcf, 'Color', 'none')
    set(gcf, 'GraphicsSmoothing', 'on')
    set(gcf,'position', [0, 0, 900, 1200])
    
    % Heatplot
    ax(1) = subplot(10,1,[1,2]);
    imagesc(traces')
    hold on
    for kk = 1:reps % Mark spikes on top of heatplot
        xpoint = repmat(kk,1,length(spikeIndex{kk}));
        scatter(spikeIndex{kk},xpoint,'.','k')
    end
    colorbar('southoutside')
    title(testname,'Interpreter','none','Fontsize',20)
    ylabel('Reps')
    colormap(ax(1),'hot')
    ax(1).XTick = [];
    ax(1).CLim = scaleAxis;
    ax(1).TickDir = 'out';
    ax(1).Box = 'off';
    
    %Mean trace compared to rmp
    ax(2) = subplot(10,1,[3,4]);
    plot(timeaxis,repmat(rmp,1,samples),'k', 'LineWidth', 1.5)
    hold on
    plot(timeaxis, meanTrace, 'r', 'LineWidth', 0.5)
    ylim([rmp - 3*tracesSD, rmp + 3*tracesSD])
    xlim([-100 900])
    ylabel('mV')
    ax(2).XTick = [];
    ax(2).TickDir = 'out';
    ax(2).Box = 'off';
    hold off
    
    % Individual Traces
    ax(3) = subplot(10,1,[5,10]);
    for kk = 1:reps
        traces = traces -range(traces(:,kk)) - 2;
        plot(timeaxis,traces(:,kk),'k')
        hold on
    end
    hold off
    title(['RMP = ', num2str(rmp), ' mV  Spike height = ',num2str(spikeheightM), ' +/- ',num2str(spikeheightSD),' mV'],'Interpreter','none')
    ylabel('Individual Traces')
    xlabel('Time around stimulus onset(ms)')
    axis tight
    ax(3).YTick = [];
    ax(3).TickDir = 'out';
    ax(3).Box = 'off';
    
    % save figure
    %     savefig([testname, '_overlay.fig'])
    print('-dtiff','-r500',[testname,'_overlay.tif'])
    %     clear celly samples reps time* *name
end
cd(folderold);