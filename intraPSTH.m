% This is to create a binned PSTH from intracellular data after running intraSummary
%
% Created 07/27/17 EHazlett


%Desired bin size in ms
if exist('binsize', 'var') %uses preset binsize if it exists
else
    binsize = 1; %set bin size if doesn't
end

PSTH = sum(raster, 2)'; %"raster" created by intraSummary
msPerSample = 1000/ length(PSTH); %Assumes 1s sweeplength
binners = floor(binsize/msPerSample);

count = 0;
for i = binners:binners:length(PSTH)
    count = count + 1;
    binned = sum(PSTH(i-binners+1:i));
    PSTH_binned(count) = binned;
    clear binned
end

figure;
plot(linspace(-100,900,length(PSTH_binned)), PSTH_binned);
title(["bin size = " num2str(binsize) "ms"])
clear PSTH_binned