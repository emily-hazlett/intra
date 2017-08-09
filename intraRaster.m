prestim = 100; %ms before sound onset
poststim = 900; %ms after sound onset

folderold = cd;
cd('D:\Intracellular\data analysis\processing\sub\');
Files = dir('*markertimes.txt');

for ii = 1:length(Files)
    markername = Files(ii).name;
    spikename = strrep(markername,'markertimes','spiketimes');
    testname = strrep(markername,'_markertimes.txt', '');
    
    markertimes = importdata(markername);
    spiketimes = importdata(spikename);
    %spiketimes = spiketimes.data;
    
    reps = length(markertimes);
    if contains(markername,'BBN') == 0
        markertimes = markertimes + 50;
    end
    
    raster = cell(reps,1);
    for mm = 1:reps
        stimtime = markertimes(mm,1);
        spikes = spiketimes(spiketimes<stimtime+poststim);
        spikes = spikes(spikes>stimtime-prestim);
        spikes = spikes - stimtime;
        raster{mm} = spikes;
        clear spikes stimtime
    end
    
    figure
    set(gcf, 'Name', testname)
    for kk = 1:reps
        xpoint = repmat(kk,1,length(spikeIndex{kk}));
        scatter(spikeIndex{kk},xpoint,'k', 'mkr', '+')
        hold on
    end
    %     line([0,0],reps,'k')
    hold off
    xlabel('Time around stimulus onset(ms)')
    ylabel('rep')
    xlim([-prestim-10,poststim+10])
    ylim([0.5,reps+0.5])
    title(testname,'Interpreter','none')
    
    set(gcf,'position', [400, 400, 1500, 150+reps])
    set(gca, 'TickDir', 'out')
    savefig([testname, '_raster.fig'])
    print('-dtiff','-r500',[testname,'_raster.tif'])
end
cd(folderold);