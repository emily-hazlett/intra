folderold = cd;
cd('D:\Intracellular\data analysis\processing\sub\');
Files = dir('*trace*.txt');

prestim = 100; %ms in sweep before stim onset
poststim = 900; %ms in sweep after stim onset
headers = 3; 

for ii = 1:length(Files)
    filename = Files(ii).name;
    testname = ['x', strrep(filename, '.txt','')];
    traces = importdata(filename);
    traces = traces.data(headers+1:end,:);
    traces = traces/10;
    
    
    tracesChop = traces;
    threshold = mean2(traces) + 5*(std2(traces)); %Chop spikes at 5 SD above mean
    tracesChop(tracesChop>threshold) = threshold;
    meanTrace = mean(tracesChop,2);
    
    [samples, reps] = size(traces);
    timestep = (prestim+poststim)/ samples;
    timeaxis = timestep-prestim:timestep:poststim;
    colors = repmat(linspace(0.1,0.6,reps)',1,3);
        
    figure
    set(gcf, 'Name', testname)
    axes('ColorOrder',colors,'NextPlot','replacechildren')
    plot(timeaxis, traces, 'LineWidth', 1)
    hold on
    plot(timeaxis, meanTrace, 'r', 'LineWidth', 2)
    ylim([-100 20])
    xlim([-100 900])
    title(testname,'Interpreter','none')
    xlabel('Time around stimulus onset(ms)')
    ylabel('mV')
    
    set(gcf,'position', [400, 400, 1500, 400])
    set(gca, 'TickDir', 'out')
    savefig([testname, '_overlay.fig'])
    print('-dtiff','-r500',[testname,'_overlay.tif'])
%     clear celly samples reps time* *name
end
cd(folderold);