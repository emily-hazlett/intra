% prompt = {'Folder containing files','Files to use'};
% dlg_title = 'Input';
% num_lines = [1, 100];
% defAns = {'C:\Users\emily\Desktop\bigone\good\', '*_tracetime.txt'};
% options.Resize='on';
% answer = inputdlg(prompt, dlg_title, num_lines, defAns, options);
% 
% folderOld = cd;
% folderNew = answer{1};
% cd(folderNew);
% Files = dir(answer{2});
% 
% button = questdlg('Start at begining of folder?', ...
%     'Start place', ...
%     'Yes', 'No', 'Yes');
% switch button
%     case 'Yes'
%         ii = 1;
%     case 'No'
%         prompt = {'Which file number to start at?'};
%         dlg_title = 'Input';
%         num_lines = [1, 100];
%         defAns = {['Must be <=', num2str(length(Files))]};
%         options.Resize='on';
%         answer = inputdlg(prompt, dlg_title, num_lines, defAns, options);
%         if isempty(answer)
%             warndlg('Must enter value')
%             return
%         elseif strcmp(answer, ['Must be <=', num2str(length(Files))]) == 1
%             warndlg('Must enter value')
%             return
%         end
%         ii = str2double(answer{1});
% end

folderOld = cd;
folderNew = 'C:\Users\emily\Desktop\bigone\maybe\';
cd(folderNew);
Files = dir('*_tracetime.txt');
ii = 1;

p = figure('units','normalized','outerposition',[0 0 1 1]);

for i = 1:length(Files)
    filename = Files(i).name;
    big = importdata(filename);
    IntraData.trace = reshape(big(2:2:end,:)', 1, numel(big(2:2:end,:)));
    IntraData.timepoints = reshape(big(1:2:end-1, :)', 1, numel(big(1:2:end-1, :)));
    clear big
    
    markername = strrep(filename, 'tracetime', 'marker');
    if exist(markername, 'file') == 0
        continue
    end
    Markers = importdata(markername);
    
    IntraData.stimname = Markers.textdata(1,2:end);
    IntraData.markertime = Markers.data(end, :);
    if size(Markers.data, 1) == 3
        IntraData.pulsedir = Markers.data(1,:);
        IntraData.pulseamp = Markers.data(2,:);
    else
        IntraData.pulsedir = [];
        IntraData.pulseamp = [];
    end
    clear Markers
    
    plot(IntraData.timepoints, IntraData.trace)
    hold on
    scatter(IntraData.markertime, repmat(50, 1, length(IntraData.markertime)), 2)
    plot(IntraData.markertime, (IntraData.pulsedir.*IntraData.pulseamp) - 800)
    hold off
    axis tight
    ylim([-1000 200])
    title(filename, 'interpreter', 'none')
    figure(p)
    keyboard
    
    indicies = get(gca, 'XLim');
    indicies (1) = floor(max(1, indicies(1)));
    indicies (2) = ceil(indicies(2));
    IntraData.traceselect(1) = find(IntraData.timepoints>indicies(1), 1, 'first') - 1;

    d = find(IntraData.timepoints>indicies(2), 1, 'first');
    if isempty(d)
        IntraData.traceselect(2) = length(IntraData.timepoints);
    else
        IntraData.traceselect(2) = d-1;
    end
    
    save(strrep(filename, '_tracetime.txt', ''), 'IntraData')
    clear IntraData indicies d
end
% toc
