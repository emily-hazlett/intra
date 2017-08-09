folderold = cd;
cd('D:\Intracellular\data analysis\1146_062117_3496_1\');
Files = dir('*full.txt');
Fs = 50000; %sample rate

for ii = 1:length(Files)
    filename = Files(ii).name;
    testname = strrep(filename, '.txt','');
    fileOut = ['x', strrep(filename, '.txt','')];
    traces = importdata(filename);
    traces = reshape(traces.data, 1, []);
    traces = traces/10;
    outie.(fileOut) = traces;
    save(fileOut, '-struct', 'outie');
    
    timeaxis = 0:1/Fs:ceil(length(traces)/Fs);
    timeaxis = timeaxis(1,1:length(traces));
    
    figure;
    set(gcf, 'Name', testname)
    plot(timeaxis,traces,'k');
    axis tight
    title(testname,'Interpreter','none')
    ylim([-1000 200])
    xlabel('Time(ms)')
    ylabel('mV')
    
    set(gcf,'position', [400, 400, 1500, 400])
    set(gca, 'TickDir', 'out')
    savefig([testname, '.fig'])
    print('-dtiff','-r500',[testname,'.tif'])
    
    clear traces
end
cd(folderold);