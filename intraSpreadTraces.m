folderold = cd;
cd('D:\Intracellular\data analysis\processing\sub\');
Files = dir('*trace.txt');

prestim = 100; %ms in sweep before stim onset
poststim = 900; %ms in sweep after stim onset
headers = 3;

for ii = 1:length(Files)
    filename = Files(ii).name;
    testname = ['x', strrep(filename, '.txt','')];
    traces = importdata(filename);
    traces = traces.data(headers+1:end,:);
    traces = traces/10;
    
    [samples, reps] = size(traces);
    timestep = (prestim+poststim)/ samples;
    timeaxis = timestep-prestim:timestep:poststim;
    
    figure
    set(gcf, 'Name', testname)
    for kk = 1:reps
        traces = traces -range(traces(:,kk)) - 2;
        plot(timeaxis,traces(:,kk))
        hold on
    end
    hold off
    
    set(gca,'YTick',[])
    ylim([rmp-10, rmp+10])
    xlim([-100 900])
    title(testname,'Interpreter','none')
    xlabel ('Time around stimulus onset (ms)')
    ylabel('Individual Traces')
    
    set(gcf,'position', [10, 10, 1500, 1000])
    set(gca, 'TickDir', 'out')
    savefig([testname, '_spread.fig'])
    print('-dtiff','-r500',[testname,'_spread.tif'])
    clear traces samples reps time* *name
end
cd(folderold);