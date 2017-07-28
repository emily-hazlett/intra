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
    spikeheightM = mean(max(traces));
    spikeheightSD = std(max(traces));
    
    %Mean trace with limited spike contamination
    tracesChop = traces;
    threshold = rmp + 4*(tracesSD); %Chop spikes at 4 SD above rmp for calculating mean trace
    tracesChop(tracesChop>threshold) = threshold;
    meanTrace = mean(tracesChop,2);
    %scaleAxis = [floor(min(min(traces))),floor(spikeheightM-5)];
    
    %Find spikes
    p = [];
    for i = 1:reps
        [m, n] = find(traces(:,i));
        spikeIndex = cat(2, p, m);
    end
    
    %% Plot figure
    timestep = (prestim+poststim)/ samples;
    timeaxis = timestep-prestim:timestep:poststim;
    scaleAxis = [floor(rmp - 2*tracesSD), ceil(rmp + 2*tracesSD)];
    
    figure;
    set(gcf, 'Name', testname)
    set(gcf, 'Color', 'none')
    set(gcf, 'GraphicsSmoothing', 'on')
    set(gcf,'position', [0, 0, 900, 1200])
    
    % Heatplot
    ax(1) = subplot(10,1,[1,2]);
    imagesc(traces')
    colorbar('southoutside')
    title(testname,'Interpreter','none','Fontsize',20)
    ylabel('Reps')
    colormap(ax(1),'hot')
    ax(1).XTick = [];
    ax(1).CLim = scaleAxis;
    ax(1).TickDir = 'out';
    ax(1).Box = 'off';
%     
%     % Individual Traces
%     ax(2) = subplot(10,1,[3,7]);
%     for kk = 1:reps
%         traces = traces -range(traces(:,kk)) - 2;
%         plot(timeaxis,traces(:,kk),'k')
%         hold on
%     end
%     hold off
%     title(['RMP = ', num2str(rmp), ' mV  Max spike height = ',num2str(spikeheightM), ' +/- ',num2str(spikeheightSD),' mV'],'Interpreter','none')
%     ylabel('Individual Traces')
%     axis tight
%     ax(2).XTick = [];
%     ax(2).YTick = [];
%     ax(2).TickDir = 'out';
%     ax(2).Box = 'off';
%     
%     %Mean trace compared to rmp
%     ax(3) = subplot(10,1,[8,10]);
%     plot(timeaxis,repmat(rmp,1,samples),'k', 'LineWidth', 1.5)
%     hold on
%     plot(timeaxis, meanTrace, 'r', 'LineWidth', 0.5)
%     ylim([rmp - 4*tracesSD, rmp + 4*tracesSD])
%     xlim([-100 900])
%     xlabel('Time around stimulus onset(ms)')
%     ylabel('mV')
%     ax(3).TickDir = 'out';
%     ax(3).Box = 'off';
%     hold off
    
    %% save figure
    %     savefig([testname, '_overlay.fig'])
    print('-dtiff','-r500',[testname,'_overlay.tif'])
    % %     clear celly samples reps time* *name
end
cd(folderold);