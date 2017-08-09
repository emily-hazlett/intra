% This script turns datawave output into rasters and PSTHs.
% Data should be output as an ASCII file.  First column is spike times,
% Second column is stim names, third column is stim times.
% Timestamp format should be in ms.
% Spikes are detected and sorted in DWave.
%
% Created by EHazlett, last edited 2017-08-08
%

folder = 'C:\Data Processing\Processing\Extracted spikes\';
file = '1175_080217_xxxx_-000010_spikes_spiketimes.txt';
prestim = 100; %ms to examine before stimulus onset
poststim = 900; %ms to examine after stimulus onset
binsize = 20; %ms bins in PSTH

%% Import data
cd(folder)
data = importdata(file);
spiketimes = str2num(cell2mat(data.textdata(2:end, 1)));
stimnames_all = data.textdata(2:end, 2);
stimtimes_all = data.data;
clear data

%% Separate out stim
stimlist = unique(stimnames_all);
nstim = length(stimlist);
for ii = 1:nstim
    index = strcmp(stimnames_all,stimlist(ii));
    stimnames {ii} = stimnames_all(index);
    stimtimes {ii} = stimtimes_all(index);
end
clear stimnames_all stimtimes_all index

%% Turn into 1ms binned raster
for ii = 1:nstim
    stimname = stimnames{ii};
    stimtime = stimtimes{ii};
    reps = length(stimname);
    
    rasterbatch = zeros(reps, 1000);
    spikeIndexbatch = cell(1,reps);
    for mm = 1:reps
        soundonset = stimtime(mm);
        time0 = soundonset - prestim;
        if any(spiketimes > ceil(time0) & spiketimes < floor(soundonset + poststim)) % Only run raster if spikes happened
            spikes = spiketimes(spiketimes>time0 & spiketimes<soundonset+poststim);
            spikes = floor(spikes - time0);
            
            spikeIndexbatch {mm} = spikes;
            rasterbatch(mm,spikes) = 1;
            clear spikes
        end
        clear time0 soundonset
    end
    spikeIndex {ii} = spikeIndexbatch;
    rasters {ii} = rasterbatch;
    clear rasterbatch spikeIndexbatch stimname stimtime
end

%% Bin PSTH
for ii = 1: nstim
    PSTHbatch = sum(rasters {ii}, 2)';
    msPerSample = (prestim + poststim) / length(PSTHbatch);
    binners = floor(binsize/msPerSample);
    
    count = 0;
    for i = binners:binners:length(PSTHbatch)
        count = count + 1;
        binned = sum(PSTHbatch(i-binners+1:i));
        PSTH_binned(count) = binned;
        clear binned
    end
    PSTHs {ii} = PSTH_binned;
    clear PSTHbatch msPerSample binners count PSTHbatch
end

%% Rough PSTH plot

figure;
set(gcf, 'Name', [file, num2str(binsize), 'ms bins'])
set(gcf, 'Color', 'w')
set(gcf, 'GraphicsSmoothing', 'on')
set(gcf,'position', [0, 0, 900, 500])

count = 0;
for ii = 1:nstim
    count = count+1;
    ax(count) = subplot(2, nstim, count);
    plot(linspace(-prestim,poststim,length(PSTHs{ii})), PSTHs{ii})
    ax(count).XLim =[-prestim, poststim];
    ax(count).YLim = [-0.1, max([2, max(PSTHs{ii})+ 0.2])];
    ax(count).TickDir = 'out';
    ax(count).LineWidth = 1.5;
    ax(count).Box = 'off';
    title(stimnames{ii}(1))
end

for ii = 1:nstim
    count = count+1;
    reps = length(stimnames{ii});
    ax(count) = subplot(2, nstim, count);
    hold on
    for kk = 1:reps % Mark spikes on top of heatplot
        xpoint = repmat(kk,1,length(spikeIndex{ii}{kk}));
        scatter(spikeIndex{ii}{kk},xpoint,'.','k', 'LineWidth', 1.1)
        clear xpoint
    end
    hold off
    ax(count).XLim =[-prestim, poststim];
    ax(count).YLim = [0, reps];
    %     ax(count).XTick = linspace(-prestim, poststim, 5);
    ax(count).TickDir = 'out';
    ax(count).LineWidth = 1.5;
    ax(count).Box = 'off';
end


print('-dtiff','-r500',strrep(file, '.txt', '.tif'))% save figure
