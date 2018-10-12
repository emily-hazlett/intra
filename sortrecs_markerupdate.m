prompt = {'Folder containing files','Files to use'};
dlg_title = 'Input';
num_lines = [1, 100];
defAns = {'C:\Users\emily\Desktop\bigone\good\', '*_tracetime.txt'};
options.Resize='on';
answer = inputdlg(prompt, dlg_title, num_lines, defAns, options);

folderOld = cd;
folderNew = answer{1};
cd(folderNew);
Files = dir(answer{2});

button = questdlg('Start at begining of folder?', ...
    'Start place', ...
    'Yes', 'No', 'Yes');
switch button
    case 'Yes'
        ii = 1;
    case 'No'
        prompt = {'Which file number to start at?'};
        dlg_title = 'Input';
        num_lines = [1, 100];
        defAns = {['Must be <=', num2str(length(Files))]};
        options.Resize='on';
        answer = inputdlg(prompt, dlg_title, num_lines, defAns, options);
        if isempty(answer)
            warndlg('Must enter value')
            return
        elseif strcmp(answer, ['Must be <=', num2str(length(Files))]) == 1
            warndlg('Must enter value')
            return
        end
        ii = str2double(answer{1});
end

% folderOld = cd;
% folderNew = 'C:\Users\emily\Desktop\bigone\good\';
% cd(folderNew);
% Files = dir('*3821*_tracetime.txt');

p = figure('units','normalized','outerposition',[0 0 1 1]);

for i = ii:length(Files)
    filename = Files(i).name;
    load(filename);
    markername = strrep(filename, '.mat', '_marker.txt');
    if exist(markername, 'file') == 0
        IntraData.stimname = {};
        IntraData.markertime = [];
        IntraData.pulsedir = [];
        IntraData.pulseamp = [];
        IntraData.markerselect = [];
        save(filename, 'IntraData')
    else
        marker = importdata(markername);
        IntraData.stimname = marker.textdata(1, 2:end);
        IntraData.markertime = marker.data(3,:);
        IntraData.pulsedir = marker.data(1,:);
        IntraData.pulseamp = marker.data(2,:);
        save(filename, 'IntraData')
    end
    clear IntraData marker markername filename
end
