clear all; clc;
folderpath = 'D:\서울대\5-1\intern\brain_volume\';

c_filename = [dir(fullfile(folderpath, '*_rh.aparc.stats')),dir(fullfile(folderpath, '*_lh.aparc.stats'))];
filename = dir(fullfile(folderpath, '*_aseg.stats'));
allfiles = cell(length(filename), 1);

for files=1:2
g_allfiles = cell(length(c_filename), 1);
volumes = zeros(length(g_allfiles), 1); 
ICVs= zeros(length(g_allfiles), 1);
for k = 1:length(c_filename)
    fileID = fopen(fullfile(folderpath, c_filename(k,files).name), 'r');
    headername = cell(1, 10);%header names go here
    col_header_line = 0;%this is where header line is
    line_num = 0;
    while ~feof(fileID)
        line = fgetl(fileID);
        line_num = line_num + 1;
        if startsWith(line, '# ColHeaders')
            col_header_line = line_num;
            break;
        end
    end
    if col_header_line == 0
        error('# ColHeaders가 나오는 행 못 찾음..!!');%if it can't find the header line
    end

    fileID = fopen(fullfile(folderpath, c_filename(k,files).name), 'r');
    headername = cell(1, 10);
    for i = 1:line_num
        lines = fgetl(fileID);
        if i == line_num
            header_tokens = strsplit(lines, ' ');
            header_tokens = header_tokens(~cellfun('isempty',header_tokens));%sometimes empty space bars get inserted..this is for getting rid of them
            for j = 3:length(header_tokens)
                headername{j-2} = header_tokens{j};
            end
        end
    end
     
    data_line = fgetl(fileID);
    data_tokens = strsplit(data_line, ' ');%this is for data type
    fclose(fileID);
    data_tokens = data_tokens(~cellfun('isempty',data_tokens));%same with the data type..no space
    datatype = cell(1, numel(headername));
    for i = 1:numel(headername)
        if isnan(str2double(data_tokens{i}))%text type
            datatype{i} = '%s';
        elseif isnumeric(str2double(data_tokens{i}))%num
            datatype{i} = '%f';
        end
    end
    fileID = fopen(fullfile(folderpath, c_filename(k,files).name), 'r');
    data = textscan(fileID, [strjoin(datatype, ' ')], 'HeaderLines', col_header_line);%put all the data here
    fclose(fileID);
    
    g_allfiles{k} = struct(); 
    for i = 1:numel(headername)
        nme = headername{i}; 
        g_allfiles{k}.(nme) = data{i}; %for storing all data from all files in that folder!!!:D
    end
%total vol
    fileID = fopen(fullfile(folderpath, c_filename(k,files).name), 'r');
    for i = 1:22
        fgetl(fileID);
    end
    line = fgetl(fileID);
    volume = strsplit(line, ',');
    volume = strtrim(volume{4});
    volumes(k) = str2double(volume);
    
    for i=23:27%ICV
        fgetl(fileID);
    end
    lineI=fgetl(fileID);
    ICV=strsplit(lineI, ',');
    ICV = strtrim(ICV{4});
    ICVs(k) = str2double(ICV);
    fclose(fileID);
end
if files==1
r_allfiles=g_allfiles;
elseif files==2
l_allfiles=g_allfiles;
end
end

% subcortical
for k = 1:length(filename)
    fileID = fopen(fullfile(folderpath, filename(k).name), 'r');
    headername = cell(1, 10);%header names go here
    col_header_line = 0;%this is where header line is
    line_num = 0;
    while ~feof(fileID)
        line = fgetl(fileID);
        line_num = line_num + 1;
        if startsWith(line, '# ColHeaders')
            col_header_line = line_num;
            break;
        end
    end
    if col_header_line == 0
        error('# ColHeaders가 나오는 행 못 찾음..!!');%if it can't find the header line
    end

    fileID = fopen(fullfile(folderpath, filename(k).name), 'r');
    headername = cell(1, 10);
    for i = 1:line_num
        lines = fgetl(fileID);
        if i == line_num
            header_tokens = strsplit(lines, ' ');
            header_tokens = header_tokens(~cellfun('isempty',header_tokens));%sometimes empty space bars get inserted..this is for getting rid of them
            for j = 3:length(header_tokens)
                headername{j-2} = header_tokens{j};
            end
        end
    end
     
    data_line = fgetl(fileID);
    data_tokens = strsplit(data_line, ' ');%this is for data type
    fclose(fileID);
    data_tokens = data_tokens(~cellfun('isempty',data_tokens));%same with the data type..no space
    datatype = cell(1, numel(headername));
    for i = 1:numel(headername)
        if isnan(str2double(data_tokens{i}))%text type
            datatype{i} = '%s';
        elseif isnumeric(str2double(data_tokens{i}))%num
            datatype{i} = '%f';
        end
    end
    fileID = fopen(fullfile(folderpath, filename(k).name), 'r');
    data = textscan(fileID, [strjoin(datatype, ' ')], 'HeaderLines', col_header_line);%put all the data here
    fclose(fileID);
    
    allfiles{k} = struct(); 
    for i = 1:numel(headername)
        nme = headername{i}; 
        allfiles{k}.(nme) = data{i}; %for storing all data from all files in that folder!!!:D
    end
   
end

k=setdiff(who, {'r_allfiles','l_allfiles', 'ICVs','volumes','allfiles'});
clear(k{:});%1

xaxis_limit=1;
cutoff=9000;

PAfilename= "D:\서울대\5-1\intern\survey_cop.xlsx";
age= xlsread(PAfilename, 3, 'D3:D174');
PA= xlsread(PAfilename, 3, 'AG3:AG174');
sex= readcell(PAfilename, 'Sheet', 3, 'Range', 'C3:C174');
education= xlsread(PAfilename, 3, 'E3:E174');
grayvol=zeros(size(PA, 1), 1); 

rh=0;%decide which hemisphere you're going to look at
if rh==1
volumes=r_volumes;
for i=1:size(PA,1)
r_grayvol(i,1)=sum(r_allfiles{i,1}.GrayVol);
end
grayvol=r_grayvol;
ICVs=r_ICVs;
else
volumes=l_volumes;
for i=1:size(PA,1)
l_grayvol(i,1)=sum(l_allfiles{i,1}.GrayVol);
end
grayvol=l_grayvol;
ICVs=l_ICVs;
end

zeroPA_indices = find(isnan(PA));
non_zero_indices = setdiff(1:numel(PA), zeroPA_indices);
PA_non_zero = PA(non_zero_indices);
volumes_non_zero = volumes(non_zero_indices);
sex_non_zero = sex(non_zero_indices);
ICV_non_zero = ICVs(non_zero_indices);
age_non_zero = age(non_zero_indices); 
education_non_zero = education(non_zero_indices); 
grayvol_non_zero=grayvol(non_zero_indices);

    if xaxis_limit==1
    cutoffPA_indices = find(PA_non_zero <= cutoff);
    PA_non_zero = PA_non_zero(cutoffPA_indices);
    volumes_non_zero=volumes_non_zero(cutoffPA_indices);
    ICV_non_zero=ICV_non_zero(cutoffPA_indices);
    sex_non_zero= sex_non_zero(cutoffPA_indices);
    age_non_zero = age_non_zero(cutoffPA_indices);
    education_non_zero=education_non_zero(cutoffPA_indices);
    grayvol_non_zero=grayvol_non_zero(cutoffPA_indices);
    else
    end

sex_non_zero = strrep(sex_non_zero, '남', 'M');
sex_non_zero = strrep(sex_non_zero, '여', 'F');
%% plot
[R,P] = corr(PA_non_zero, grayvol_non_zero);
sprintf("correlation:%f P value:%f",R,P)
mdl = fitlm(PA_non_zero, grayvol_non_zero);
r_squared = mdl.Rsquared.Ordinary;
p_value = mdl.Coefficients.pValue;
plot(PA_non_zero, grayvol_non_zero,'.')
hold on; 
x_range = [min(PA_non_zero), max(PA_non_zero)];
y_pred = mdl.Coefficients.Estimate(1) + mdl.Coefficients.Estimate(2) * x_range; 

plot(x_range, y_pred, 'k-', 'LineWidth', 2);
xlabel('age');
ylabel('grayvol');
hold off
%% looking at different areas for correlation
variant=education;% which variant are you looking at
ratio=0;% are you choosing to look at the ratio of that region?
rh=0;

if rh==1
g_allfiles=r_allfiles;   
else
g_allfiles=l_allfiles;   
end

allvol=cell(14,1);
vol=zeros(size(variant, 1), 1); 

for k=1:numel(g_allfiles{1,1}.StructName)
for i=1:size(variant,1)
    if ratio==1
    vol(i,1)=g_allfiles{i,1}.GrayVol(k)/volumes(i);
    else
    vol(i,1)=g_allfiles{i,1}.ThickAvg(k);
    end
eh=vol(non_zero_indices);
ep=variant(non_zero_indices);
[R,P] = corr(ep(cutoffPA_indices),eh(cutoffPA_indices));
end
allvol{k,1}=[R,P];
end

correlations = cellfun(@(x) abs(x(1)), allvol);
[sortedCorrelations, sortIdx] = sort(correlations, 'descend');

sprintf("Highest correlation area: %s, correlation value: %f, P value: %f", ...
    g_allfiles{1,1}.StructName{sortIdx(1), 1}, ...
    allvol{sortIdx(1), 1}(1), allvol{sortIdx(1), 1}(2))

sprintf("Second highest correlation area: %s, correlation value: %f, P value: %f", ...
     g_allfiles{1,1}.StructName{sortIdx(2), 1}, ...
    allvol{sortIdx(2), 1}(1), allvol{sortIdx(2), 1}(2))

sprintf("Third highest correlation area: %s, correlation value: %f, P value: %f", ...
     g_allfiles{1,1}.StructName{sortIdx(3), 1}, ...
    allvol{sortIdx(3), 1}(1), allvol{sortIdx(3), 1}(2))
%% group (box plot)
groupnum=xlsread(PAfilename, 3, 'A3:A174');
PAgroup=xlsread(PAfilename, 3, 'AP3:AP174');
edugroup= readcell(PAfilename, 'Range', 'AS3:AS174', 'Sheet', 3);

agegroup_non_zero = groupnum(non_zero_indices);
PAgroup_non_zero=PAgroup(non_zero_indices);
edugroup_non_zero=edugroup(non_zero_indices);
    if xaxis_limit==1
        agegroup_non_zero = agegroup_non_zero(cutoffPA_indices);
        PAgroup_non_zero=PAgroup_non_zero(cutoffPA_indices);
        edugroup_non_zero=edugroup_non_zero(cutoffPA_indices);
    end

%for character ones
data = cell(4, 1);
label = cell(4, 1);
for k= 1:2
for i = 1:2
    % for j = 1:2
        agelabel={'mid','old'};
        edulabel={'edu level high','edu level low'};
        PAlabel={'low PA','high PA'};
        sexlabel={'M','F'};
        agegroups = cell(numel(agegroup_non_zero), 1);
        agegroups(agegroup_non_zero == 2 | agegroup_non_zero == 4) = {'mid'};
        agegroups(agegroup_non_zero == 5) = {'old'};

        PAgroups = cell(numel(PAgroup_non_zero), 1);
        PAgroups(PAgroup_non_zero == 1 | PAgroup_non_zero == 2) = {'low PA'};
        PAgroups(PAgroup_non_zero == 3 | PAgroup_non_zero == 4) = {'high PA'};

        idx = strcmp(sex_non_zero,sexlabel{1,i}) & strcmp(agegroups,agelabel{1,k});
        % idx = strcmp(agegroups, agelabel{1,i}) & strcmp(PAgroups,PAlabel{1,j}) & strcmp(edugroup_non_zero,edulabel{1,k});
        % data{(i-1)*4 + (j-1)*2+ k} = grayvol_non_zero(idx);
        data{(i-1)*2 + k} = grayvol_non_zero(idx);

        n=numel(data{(i-1)*2 + k});
        % n=numel(data{(i-1)*4 + (j-1)*2+ k});
        % label{(i-1)*4 + (j-1)*2+ k} = repmat(sprintf('(%s,%s,%s)', agelabel{1,i},PAlabel{1,j},edulabel{1,k}),n, 1);
        label{(i-1)*2 + k} = repmat(sprintf('(%s,%s)', sexlabel{1,i},agelabel{1,k}),n, 1);
    % end
end
end

% only for text ones
all_text = [];
for i = 1:numel(label)
    allt=cellstr(label{i, 1});
    all_text = vertcat(all_text,allt);
end
figure;
boxplot(cell2mat(data),all_text)
ylabel('volume (mm)^{3}');
% xlabel('age,PA');
%% one way ANOVA
[p, tbl, stats] = anova1(cell2mat(data), all_text);
%% pfc regions and ratio/ make sure to load subcortical regions too
load('D:\서울대\5-1\intern\brain volume values.mat')

r_pfc=zeros(numel(r_allfiles), 1); 
l_pfc=zeros(numel(r_allfiles), 1); 

r_hippo=zeros(numel(r_allfiles), 1); 
l_hippo=zeros(numel(r_allfiles), 1); 

for i=1:numel(r_allfiles)
r_pfc(i,1)=sum(r_allfiles{i,1}.GrayVol([3,11,13,26,27,31],1));
l_pfc(i,1)=sum(l_allfiles{i,1}.GrayVol([3,11,13,26,27,31],1));
end

for i=1:numel(allfiles)
r_hippo(i,1)=sum(allfiles{i,1}.Volume_mm3(27,1));
l_hippo(i,1)=sum(allfiles{i,1}.Volume_mm3(12,1));
end

r_hippo_non_zero = r_hippo(non_zero_indices);
l_hippo_non_zero = l_hippo(non_zero_indices);
r_hippo_non_zero=r_hippo_non_zero(cutoffPA_indices);
l_hippo_non_zero=l_hippo_non_zero(cutoffPA_indices);

r_pfc_non_zero = r_pfc(non_zero_indices);
l_pfc_non_zero = l_pfc(non_zero_indices);
r_pfc_non_zero=r_pfc_non_zero(cutoffPA_indices);
l_pfc_non_zero=l_pfc_non_zero(cutoffPA_indices);
pfc_non_zero=l_pfc_non_zero+r_pfc_non_zero;

r_hpfratio_nz=r_pfc_non_zero ./r_hippo_non_zero;
l_hpfratio_nz=l_pfc_non_zero ./l_hippo_non_zero;
hpfratio_nz=pfc_non_zero ./hippo_non_zero;
%% 산점도 for all areas together-pfc
ROI=r_hippo_non_zero;%where are you looking at
variant=age_non_zero;

[R,P] = corr(variant, ROI);
sprintf("correlation:%f P value:%f",R,P)
mdl = fitlm(variant, ROI);
r_squared = mdl.Rsquared.Ordinary;
p_value = mdl.Coefficients.pValue;
plot(variant, ROI,'.')
hold on; 
x_range = [min(variant), max(variant)];
y_pred = mdl.Coefficients.Estimate(1) + mdl.Coefficients.Estimate(2) * x_range; 

plot(x_range, y_pred, 'k-', 'LineWidth', 2);
hold off
 
ylabel('hippo vol');
xlabel('age');
%% for individual areas in the pfc
r_pfc=zeros(numel(r_allfiles), 1); 
l_pfc=zeros(numel(r_allfiles), 1); 
output_dir = 'D:\서울대\5-1\intern\figure\';
allR=[];allP=[];

where=[3,11,13,26,27,31];
for k=1:6
for i=1:numel(r_allfiles)
% r_pfc(i,1)=sum(r_allfiles{i,1}.GrayVol(where(i),1));
r_pfc(i,1)=r_allfiles{i,1}.GrayVol(where(k),1);
l_pfc(i,1)=l_allfiles{i,1}.GrayVol(where(k),1);
end

r_pfc_non_zero = r_pfc(non_zero_indices);
l_pfc_non_zero = l_pfc(non_zero_indices);
r_pfc_non_zero=r_pfc_non_zero(cutoffPA_indices);
l_pfc_non_zero=l_pfc_non_zero(cutoffPA_indices);
pfc_non_zero=l_pfc_non_zero+r_pfc_non_zero;

r_hpfratio_nz=r_pfc_non_zero ./r_hippo_non_zero;
l_hpfratio_nz=l_pfc_non_zero ./l_hippo_non_zero;
hpfratio_nz=pfc_non_zero ./hippo_non_zero;

ROI=r_pfc_non_zero;%where are you looking at
variant=age_non_zero;

[R,P] = corr(variant, ROI);
fprintf("correlation:%f P value:%f\n",R,P)
allR=[allR;R];allP=[allP;P];

mdl = fitlm(variant, ROI);
r_squared = mdl.Rsquared.Ordinary;
p_value = mdl.Coefficients.pValue;

figure;
plot(variant, ROI,'.')
hold on; 
x_range = [min(variant), max(variant)];
y_pred = mdl.Coefficients.Estimate(1) + mdl.Coefficients.Estimate(2) * x_range; 

plot(x_range, y_pred, 'k-', 'LineWidth', 2);
hold off

ylabel(sprintf('%s region gray matter vol', r_allfiles{1, 1}.StructName{where(k), 1}));
xlabel('age');

disp(r_allfiles{1, 1}.StructName{where(k), 1})
end
[sortedCorrelations, sortIdx] = sort(allR, 'descend');[sortedp, sortIdx] = sort(allP, 'descend');

fprintf("Highest correlation area: %s, correlation value: %f, P value: %f \n", ...
    r_allfiles{1,1}.StructName{where(sortIdx(1)), 1}, ...
    sortedCorrelations(1), sortedp(1))

fprintf("Second highest correlation area: %s, correlation value: %f, P value: %f \n", ...
    r_allfiles{1,1}.StructName{where(sortIdx(2)), 1}, ...
    sortedCorrelations(2), sortedp(2))

fprintf("Third highest correlation area: %s, correlation value: %f, P value: %f \n", ...
    r_allfiles{1,1}.StructName{where(sortIdx(3)), 1}, ...
    sortedCorrelations(3), sortedp(3))
% close all
