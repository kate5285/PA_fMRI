clc;clear all;
folderPath = 'D:\서울대\5-1\intern\brainhp';  
filePattern = fullfile(folderPath, '*.txt');   
pfiles = dir(filePattern);    
l_allfiles = {};
partici=[];
r_allfiles = {};
for k = 1:length(pfiles)
    baseFileName = pfiles(k).name;                  
    fullFileName = fullfile(pfiles(k).folder, baseFileName);
    fileID = fopen(fullFileName, 'r');
    datatype = {'%s', '%f'};
            fileID = fopen(fullFileName, 'r');
            data = textscan(fileID, [strjoin(datatype, ' ')], 'HeaderLines',0); % put all the data here
            fclose(fileID);
            fileStruct = struct();
            fileStruct.name = data{1};
            fileStruct.vol = data{2};

            if contains(baseFileName, 'lh')
                l_allfiles{end+1, 1} = fileStruct; % "lh"
                partici(end+1, 1) = str2double(regexp(baseFileName, '^\d+', 'match', 'once'));
            elseif contains(baseFileName, 'rh')
                r_allfiles{end+1, 1} = fileStruct; % "rh"
           end
end
clearvars -except r_allfiles l_allfiles pfiles partici

dontlookat=[5,10,11,17,19];
indicesToKeep = setdiff(1:22, dontlookat);
for i=1:numel(r_allfiles)
r_allfiles{i} = arrayfun(@(x) struct('name', {x.name(indicesToKeep)}, 'vol', x.vol(indicesToKeep)), r_allfiles{i});
l_allfiles{i} = arrayfun(@(x) struct('name', {x.name(indicesToKeep)}, 'vol', x.vol(indicesToKeep)), l_allfiles{i});
end
%% combine the head and tail
r_combined = cell(size(r_allfiles));
l_combined = cell(size(l_allfiles));

for i = 1:numel(r_allfiles)
names = r_allfiles{1,1}.name;
vols = r_allfiles{i}.vol;
new_names=[{'Hippocamapal_tail'};{'parasubiculum'};{'Whole_hippocampal_body'};{'Whole_hippocampal_head'};{'Whole_hippocampus'}];
new_vols=[vols(1);vols(8);vols(15);vols(16);vols(17)];
pairs = [2 4; 3 6; 5 7; 9 11; 10 14; 12 13];
    for j = 1:size(pairs, 1)
        dash=strfind(names{pairs(j, 1)},'-');
        combined_name = names{pairs(j, 1)}(1:dash(end) - 1);
        combined_vol = vols(pairs(j, 1)) + vols(pairs(j, 2));
        new_names = [new_names; {combined_name}];
        new_vols = [new_vols; combined_vol];
    end
    r_combined{i} = struct('name', {new_names}, 'vol', new_vols);

names = l_allfiles{1,1}.name;
vols = l_allfiles{i}.vol;
new_names=[{'Hippocamapal_tail'};{'parasubiculum'};{'Whole_hippocampal_body'};{'Whole_hippocampal_head'};{'Whole_hippocampus'}];
new_vols=[vols(1);vols(8);vols(15);vols(16);vols(17)];

    for j = 1:size(pairs, 1)
        dash=strfind(names{pairs(j, 1)},'-');
        combined_name = names{pairs(j, 1)}(1:dash(end) - 1);
        combined_vol = vols(pairs(j, 1)) + vols(pairs(j, 2));
        new_names = [new_names; {combined_name}];
        new_vols = [new_vols; combined_vol];
    end
    l_combined{i} = struct('name', {new_names}, 'vol', new_vols);
end
%%
load("C:\Users\kate5\Downloads\rel.mat"); %rel acc

exfilename= "D:\서울대\5-1\intern\s_updated.xlsx";
behavok= readcell(exfilename, 'Sheet', 1, 'Range', 'C3:C264');
dtalabel= readcell(exfilename, 'Sheet', 1, 'Range', 'E3:E264');
okdrow = find(cellfun(@(x) ismember(x, [2,4,5]), dtalabel)); %mid and old
okrow = find(cellfun(@(x) isequal(x, 'O'), behavok));
allok = intersect(okdrow, okrow);

parti= readcell(exfilename, 'Sheet', 1, 'Range', 'A3:A264');
parti=cellfun(@double, parti);
epartici=zeros(numel(allok),1);
epartici=parti(allok);

[relp, sort_idx] = sort(rel.participant);
rel_filtered = rel(sort_idx, :);
% rel_filtered = rel_filtered(non_zero_indices, :);

Kfilename= "D:\서울대\5-1\intern\survey_cop.xlsx";
datalabel = xlsread(Kfilename, 3, 'A:A');
numRows = numel(datalabel);

sex = readcell(Kfilename, 'Sheet', 3, 'Range',sprintf('C3:C%d', numRows + 2));
age= xlsread(Kfilename, 3, sprintf('D3:D%d', numRows + 2));
PA= xlsread(Kfilename, 3, sprintf('AG3:AG%d', numRows + 2));
education= xlsread(Kfilename, 3, sprintf('E3:E%d', numRows + 2));
subj_nums=xlsread(Kfilename, 3, sprintf('B3:B%d', numRows + 2));
load('D:\서울대\5-1\intern\ICVs.mat')

isMemberArray = ismember(subj_nums,partici);
commonElements = subj_nums(isMemberArray);

sex_non_zero = sex(isMemberArray);
sex_binary = strcmp(sex_non_zero, 'F');
age_non_zero= age(isMemberArray);
PA_non_zero= PA(isMemberArray);
education_non_zero= education(isMemberArray);

load('D:\서울대\5-1\intern\participants.mat')
isMemberArray = ismember(participant,commonElements);
IcommonElements = participant(isMemberArray);
ICV_non_zero=ICVs(isMemberArray);
%% for accuracy
isMemberArray = ismember(rel_filtered.participant,commonElements);
IcommonElements = participant(isMemberArray);

sex_acc = sex_binary(isMemberArray);
age_acc= age_non_zero(isMemberArray);
PA_acc= PA_non_zero(isMemberArray);
education_acc= education_non_zero(isMemberArray);
ICV_acc=ICV_non_zero(isMemberArray);
rel_filtered=rel_filtered(isMemberArray, :);

nm = ["what", "where", "when"];
names = {};
rv = [];
pv = [];
task={};
what=1;%area
where=r_combined{1,1}.name(what);

    for P = 1:3
        for i = 1:numel(r_combined)
            r_roi(i, 1) = sum(r_combined{i, 1}.vol(what(:), 1));
            l_roi(i, 1) = sum(l_combined{i, 1}.vol(what(:), 1));
        roi = (r_roi+l_roi);
        end
        % ROI_adj=adjustROI(roi,ICV_non_zero,non_zero_indices);
        data = table(roi(isMemberArray), age_acc, sex_acc, education_acc, ICV_acc, PA_acc, ...
                     'VariableNames', {'ROI', 'Age', 'Sex', 'Education', 'ICV', 'PA'});
        mdl = fitlm(data, 'linear', 'ResponseVar', 'ROI', 'PredictorVars', {'Sex','Age','ICV','Education'});
        [r, p] = corr(rel_filtered.(['accuracy' num2str(P)]), mdl.Fitted);
        % [r, p] = corr(rel_filtered.(['accuracy' num2str(P)]), ROI_adj);

        names = vertcat(names,where); 
        task = vertcat(task, nm{P}); 
        rv = vertcat(rv, r);
        pv = vertcat(pv, p);

        figure;
        plot(rel_filtered.(['accuracy' num2str(P)]), mdl.Fitted,'.')
        % plot(rel_filtered.(['accuracy' num2str(P)]), ROI_adj,'.')
        hold on;
        x_range = [min(rel_filtered.(['accuracy' num2str(P)])), max(rel_filtered.(['accuracy' num2str(P)]))]; % x 축 범위
        md = fitlm(rel_filtered.(['accuracy' num2str(P)]),mdl.Fitted);
        % md = fitlm(rel_filtered.(['accuracy' num2str(P)]),ROI_adj);
        y_pred = md.Coefficients.Estimate(1) + md.Coefficients.Estimate(2) * x_range; 
        plot(x_range,y_pred,'k-','LineWidth',2);
        hold off;
        xlabel('acc');
        ylabel('vol');
        title(nm{P});
    end

tab = table(names, task,rv, pv, 'VariableNames', {'Name','task', 'r', 'p'});
disp(tab);

nm = ["what", "where", "when"]; 
P=1;%정확도 중 어디를 볼 것인지
rv=[];pv=[];names=[];
for h=1:numel(l_combined{1,1}.name)
 for i = 1:numel(l_combined)
            r_roi(i, 1) = sum(r_combined{i, 1}.vol(h, 1));
            l_roi(i, 1) = sum(l_combined{i, 1}.vol(h, 1));
        roi = (r_roi+l_roi);
 end
        % ROI_adj=adjustROI(roi,ICV_non_zero,non_zero_indices);
        data = table(roi(isMemberArray), age_acc, sex_acc, education_acc, ICV_acc, PA_acc, ...
                     'VariableNames', {'ROI', 'Age', 'Sex', 'Education', 'ICV', 'PA'});
        mdl = fitlm(data, 'linear', 'ResponseVar', 'ROI', 'PredictorVars', {'Sex', 'Education','ICV','Age'});
        [r, p] = corr(rel_filtered.(['accuracy' num2str(P)]), mdl.Fitted);
        % [r, p] = corr(rel_filtered.(['accuracy' num2str(P)]),((data.ROI) ./(data.ICV)));
        pv=vertcat(pv,p);
        rv=vertcat(rv,r);
        name=l_combined{1,1}.name(h,1);
        names = vertcat(names,name); 

        figure;
        plot(rel_filtered.(['accuracy' num2str(P)]),mdl.Fitted,'.')

end
        [~,sorted]=sort(pv,'ascend');
        tableof=table(names(sorted),rv(sorted), pv(sorted), ...
     'VariableNames', {'Name', 'r', 'p'})
%% krbans
Kfilename= "D:\서울대\5-1\intern\KRBANS.xlsx";
Kpartici = xlsread(Kfilename, 1, 'A:A');
numRows = numel(Kpartici);

words = readcell(Kfilename, 1, sprintf('B3:B%d', numRows + 2));
story= xlsread(Kfilename, 3, sprintf('C3:C%d', numRows + 2));
PA= xlsread(Kfilename, 3, sprintf('AG3:AG%d', numRows + 2));
