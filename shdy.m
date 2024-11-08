clear all; 
exfilename= "D:\서울대\5-1\intern\s_updated.xlsx";
behavok= readcell(exfilename, 'Sheet', 1, 'Range', 'C3:C264');
dtalabel= readcell(exfilename, 'Sheet', 1, 'Range', 'E3:E264');
mkdrow = find(cellfun(@(x) ismember(x, [2,4,5]), dtalabel)); %mid and old
% mkrow = find(cellfun(@(x) isequal(x, 'O'), behavok));
% mallok = intersect(mkdrow, mkrow);
mallok=mkdrow;
mparti= readcell(exfilename, 'Sheet', 1, 'Range', 'A3:A264');
mparti=cellfun(@double, mparti);
% allok=[1:172,222,223];
mpartici=zeros(numel(mallok),1);
mpartici=mparti(mallok);

okdrow = find(cellfun(@(x) ismember(x, [5]), dtalabel)); %mid and old
okrow = find(cellfun(@(x) isequal(x, 'O'), behavok));
allok = intersect(okdrow, okrow);
% allok=okdrow;
oparti= readcell(exfilename, 'Sheet', 1, 'Range', 'A3:A264');
oparti=cellfun(@double, oparti);
% allok=[1:172,222,223];
opartici=zeros(numel(allok),1);
opartici=oparti(allok);
%%
folderpath = 'D:\서울대\5-1\intern\brain_volume\';

c_filenames = [dir(fullfile(folderpath, '*_rh.aparc.stats')), dir(fullfile(folderpath, '*_lh.aparc.stats'))];
filenames = dir(fullfile(folderpath, '*_aseg.stats'));

g_allfiles = cell(length(opartici), 1);
volumes = zeros(length(opartici), 1); 
ICVs = zeros(length(opartici), 1);

for files = 1:2
    for k = 1:length(c_filenames)
        filename = c_filenames(k, files).name;
        % extract sbj number from file name
        subj_num_str = regexp(filename, '^\d+', 'match', 'once');
        subj_num = str2double(subj_num_str);
        % check if the sbj number is on the participant list
        if ismember(subj_num, opartici)
            fileID = fopen(fullfile(folderpath, filename), 'r');
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

            fileID = fopen(fullfile(folderpath, filename), 'r');
            headername = cell(1, 10);
            for i = 1:line_num
                lines = fgetl(fileID);
                if i == line_num
                    header_tokens = strsplit(lines, ' ');
                    header_tokens = header_tokens(~cellfun('isempty', header_tokens)); % sometimes empty space bars get inserted..this is for getting rid of them
                    for j = 3:length(header_tokens)
                        headername{j-2} = header_tokens{j};
                    end
                end
            end

            data_line = fgetl(fileID);
            data_tokens = strsplit(data_line, ' '); % this is for data type
            fclose(fileID);
            data_tokens = data_tokens(~cellfun('isempty', data_tokens)); % same with the data type..no space
            datatype = cell(1, numel(headername));
            for i = 1:numel(headername)
                if isnan(str2double(data_tokens{i})) % text type
                    datatype{i} = '%s';
                elseif isnumeric(str2double(data_tokens{i})) % num
                    datatype{i} = '%f';
                end
            end
            fileID = fopen(fullfile(folderpath, filename), 'r');
            data = textscan(fileID, [strjoin(datatype, ' ')], 'HeaderLines', col_header_line); % put all the data here
            fclose(fileID);

            g_allfiles{k} = struct(); 
            for i = 1:numel(headername)
                nme = headername{i}; 
                g_allfiles{k}.(nme) = data{i}; % for storing all data from all files in that folder!!!:D
            end

            % Total vol
            fileID = fopen(fullfile(folderpath, filename), 'r');
            for i = 1:22
                fgetl(fileID);
            end
            % extract brain vol
            line = fgetl(fileID);
            volume = strsplit(line, ',');
            volume = strtrim(volume{4});
            % mm^3 
            volumes(k) = str2double(volume);

            for i = 23:27 % ICV
                fgetl(fileID);
            end
            lineI = fgetl(fileID);
            ICV = strsplit(lineI, ',');
            ICV = strtrim(ICV{4});
            ICVs(k) = str2double(ICV);
            fclose(fileID);
        end
    end
    if files == 1
        r_allfiles = g_allfiles;
    elseif files == 2
        l_allfiles = g_allfiles;
    end
end

% subcortical
for k = 1:length(filenames)
    filename = filenames(k).name;
    subj_num_str = regexp(filename, '^\d+', 'match', 'once');
    subj_num = str2double(subj_num_str);
    if ismember(subj_num, opartici)
        fileID = fopen(fullfile(folderpath, filename), 'r');
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

        fileID = fopen(fullfile(folderpath, filename), 'r');
        headername = cell(1, 10);
        for i = 1:line_num
            lines = fgetl(fileID);
            if i == line_num
                header_tokens = strsplit(lines, ' ');
                header_tokens = header_tokens(~cellfun('isempty', header_tokens)); % sometimes empty space bars get inserted..this is for getting rid of them
                for j = 3:length(header_tokens)
                    headername{j-2} = header_tokens{j};
                end
            end
        end

        data_line = fgetl(fileID);
        data_tokens = strsplit(data_line, ' '); % this is for data type
        fclose(fileID);
        data_tokens = data_tokens(~cellfun('isempty', data_tokens)); % same with the data type..no space
        datatype = cell(1, numel(headername));
        for i = 1:numel(headername)
            if isnan(str2double(data_tokens{i})) % text type
                datatype{i} = '%s';
            elseif isnumeric(str2double(data_tokens{i})) % num
                datatype{i} = '%f';
            end
        end
        fileID = fopen(fullfile(folderpath, filename), 'r');
        data = textscan(fileID, [strjoin(datatype, ' ')], 'HeaderLines', col_header_line); % put all the data here
        fclose(fileID);

        allfiles{k} = struct(); 
        for i = 1:numel(headername)
            nme = headername{i}; 
            allfiles{k}.(nme) = data{i}; % for storing all data from all files in that folder!!!:D
        end
    end
end
ICVs=ICVs(ICVs ~= 0);
volumes=volumes(volumes ~=0);
l_allfiles=l_allfiles(~cellfun('isempty', l_allfiles));
r_allfiles=r_allfiles(~cellfun('isempty', r_allfiles));
allfiles=allfiles(~cellfun('isempty', allfiles));
allfiles=allfiles.';
% k = setdiff(who, {'r_allfiles', 'l_allfiles', 'ICVs', 'volumes', 'allfiles', 'allok', 'partici', 'parti', 'dtalabel', 'exfilename'});
% clear(k{:});

%for other- mid
c_filenames = [dir(fullfile(folderpath, '*_rh.aparc.stats')), dir(fullfile(folderpath, '*_lh.aparc.stats'))];
filenames = dir(fullfile(folderpath, '*_aseg.stats'));

g_allfiles = cell(length(mpartici), 1);
mvolumes = zeros(length(mpartici), 1); 
mICVs = zeros(length(mpartici), 1);

for files = 1:2
    for k = 1:length(c_filenames)
        filename = c_filenames(k, files).name;
        subj_num_str = regexp(filename, '^\d+', 'match', 'once');
        subj_num = str2double(subj_num_str);
        if ismember(subj_num, mpartici)
            fileID = fopen(fullfile(folderpath, filename), 'r');
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

            fileID = fopen(fullfile(folderpath, filename), 'r');
            headername = cell(1, 10);
            for i = 1:line_num
                lines = fgetl(fileID);
                if i == line_num
                    header_tokens = strsplit(lines, ' ');
                    header_tokens = header_tokens(~cellfun('isempty', header_tokens)); % sometimes empty space bars get inserted..this is for getting rid of them
                    for j = 3:length(header_tokens)
                        headername{j-2} = header_tokens{j};
                    end
                end
            end

            data_line = fgetl(fileID);
            data_tokens = strsplit(data_line, ' '); % this is for data type
            fclose(fileID);
            data_tokens = data_tokens(~cellfun('isempty', data_tokens)); % same with the data type..no space
            datatype = cell(1, numel(headername));
            for i = 1:numel(headername)
                if isnan(str2double(data_tokens{i})) % text type
                    datatype{i} = '%s';
                elseif isnumeric(str2double(data_tokens{i})) % num
                    datatype{i} = '%f';
                end
            end
            fileID = fopen(fullfile(folderpath, filename), 'r');
            data = textscan(fileID, [strjoin(datatype, ' ')], 'HeaderLines', col_header_line); % put all the data here
            fclose(fileID);

            g_allfiles{k} = struct(); 
            for i = 1:numel(headername)
                nme = headername{i}; 
                g_allfiles{k}.(nme) = data{i}; % for storing all data from all files in that folder!!!:D
            end

            % Total vol
            fileID = fopen(fullfile(folderpath, filename), 'r');
            for i = 1:22
                fgetl(fileID);
            end
            line = fgetl(fileID);
            mvolume = strsplit(line, ',');
            mvolume = strtrim(mvolume{4});
            mvolumes(k) = str2double(mvolume);

            for i = 23:27 % ICV
                fgetl(fileID);
            end
            lineI = fgetl(fileID);
            mICV = strsplit(lineI, ',');
            mICV = strtrim(mICV{4});
            mICVs(k) = str2double(mICV);
            fclose(fileID);
        end
    end
    if files == 1
        mr_allfiles = g_allfiles;
    elseif files == 2
        ml_allfiles = g_allfiles;
    end
end

% subcortical
for k = 1:length(filenames)
    filename = filenames(k).name;
    subj_num_str = regexp(filename, '^\d+', 'match', 'once');
    subj_num = str2double(subj_num_str);
    if ismember(subj_num, mpartici)
        fileID = fopen(fullfile(folderpath, filename), 'r');
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

        fileID = fopen(fullfile(folderpath, filename), 'r');
        headername = cell(1, 10);
        for i = 1:line_num
            lines = fgetl(fileID);
            if i == line_num
                header_tokens = strsplit(lines, ' ');
                header_tokens = header_tokens(~cellfun('isempty', header_tokens)); % sometimes empty space bars get inserted..this is for getting rid of them
                for j = 3:length(header_tokens)
                    headername{j-2} = header_tokens{j};
                end
            end
        end

        data_line = fgetl(fileID);
        data_tokens = strsplit(data_line, ' '); % this is for data type
        fclose(fileID);
        data_tokens = data_tokens(~cellfun('isempty', data_tokens)); % same with the data type..no space
        datatype = cell(1, numel(headername));
        for i = 1:numel(headername)
            if isnan(str2double(data_tokens{i})) % text type
                datatype{i} = '%s';
            elseif isnumeric(str2double(data_tokens{i})) % num
                datatype{i} = '%f';
            end
        end
        fileID = fopen(fullfile(folderpath, filename), 'r');
        data = textscan(fileID, [strjoin(datatype, ' ')], 'HeaderLines', col_header_line); % put all the data here
        fclose(fileID);

        mallfiles{k} = struct(); 
        for i = 1:numel(headername)
            nme = headername{i}; 
            mallfiles{k}.(nme) = data{i}; % for storing all data from all files in that folder!!!:D
        end
    end
end

mICVs=mICVs(mICVs ~= 0);
mvolumes=mvolumes(mvolumes ~=0);
ml_allfiles=ml_allfiles(~cellfun('isempty', ml_allfiles));
mr_allfiles=mr_allfiles(~cellfun('isempty', mr_allfiles));
mallfiles=mallfiles(~cellfun('isempty', mallfiles));
mallfiles=mallfiles.';
%% 1
r_pfc=zeros(numel(r_allfiles), 1); 
l_pfc=zeros(numel(r_allfiles), 1); 
mr_pfc=zeros(numel(mr_allfiles), 1); 
ml_pfc=zeros(numel(mr_allfiles), 1); 

r_hippo=zeros(numel(r_allfiles), 1); 
l_hippo=zeros(numel(r_allfiles), 1); 
mr_hippo=zeros(numel(mr_allfiles), 1); 
ml_hippo=zeros(numel(mr_allfiles), 1); 
output_dir = 'D:\서울대\5-1\intern\figure\';
where=[3,11,13,26,27,31];%PFC
% where=[3,26,27]; %dlPFC
% where=[15];%para
% for k=1:6
for i=1:numel(r_allfiles)
r_pfc(i,1)=sum(r_allfiles{i,1}.GrayVol(where(:),1));
l_pfc(i,1)=sum(l_allfiles{i,1}.GrayVol(where(:),1));
r_pfc_thick(i,1)=mean(r_allfiles{i,1}.ThickAvg(where(:),1));
l_pfc_thick(i,1)=mean(l_allfiles{i,1}.ThickAvg(where(:),1));
pfc_thick=mean(horzcat(r_pfc_thick,l_pfc_thick),2);
% r_pfc(i,1)=r_allfiles{i,1}.GrayVol(where(k),1);
% l_pfc(i,1)=l_allfiles{i,1}.GrayVol(where(k),1);
end
% end
for i=1:numel(mr_allfiles)
mr_pfc(i,1)=sum(mr_allfiles{i,1}.GrayVol(where(:),1));
ml_pfc(i,1)=sum(ml_allfiles{i,1}.GrayVol(where(:),1));
mr_pfc_thick(i,1)=mean(mr_allfiles{i,1}.ThickAvg(where(:),1));
ml_pfc_thick(i,1)=mean(ml_allfiles{i,1}.ThickAvg(where(:),1));
mpfc_thick=mean(horzcat(mr_pfc_thick,ml_pfc_thick),2);
end

for i=1:numel(r_allfiles)
r_hippo(i,1)=allfiles{i,1}.Volume_mm3(27,1);
l_hippo(i,1)=allfiles{i,1}.Volume_mm3(12,1);
end
for i=1:numel(mr_allfiles)
mr_hippo(i,1)=mallfiles{i,1}.Volume_mm3(27,1);
ml_hippo(i,1)=mallfiles{i,1}.Volume_mm3(12,1);
end

hippo=r_hippo+l_hippo;
pfc=l_pfc+r_pfc;
r_hpfratio=r_pfc ./r_hippo;

l_hpfratio=l_pfc ./l_hippo;

hpf_ratio=pfc ./hippo;

mhippo=mr_hippo+ml_hippo;
mpfc=ml_pfc+mr_pfc;
mr_hpfratio=mr_pfc ./mr_hippo;

ml_hpfratio=ml_pfc ./ml_hippo;

mhpf_ratio=mpfc ./mhippo;

[sort_hpf_ratio, sorted_indices] = sort(hpf_ratio, 'descend');%sort the ratio of pfc/hp
high_to_low_hpf=opartici(sorted_indices);

median_participant= high_to_low_hpf(floor(numel(opartici)/2)+1,1);
high_participants=high_to_low_hpf(1:floor(numel(opartici)/2),1);
low_participants=high_to_low_hpf(floor(numel(opartici)/2)+2:end,1);
% k = setdiff(who, {'high_participants','low_participants','median_participant','high_to_low_hpf'});
% clear(k{:});
 
figure;
 for i=1:numel(opartici)
     plot(i,hpf_ratio(find(opartici==high_to_low_hpf(i)),1),'.')
     hold on;
 end
 hold off;
% sort_hpf_ratio goes from the highest to lowest ratio, high_to_low_hpf is the label for it

% hpfrat = table(high_to_low_hpf, sort_hpf_ratio, 'VariableNames', {'person_idx', 'hpf_ratio'});
% save_path = 'D:\서울대\5-1\intern\old_hpf.mat';
% save(save_path, 'hpfrat_mid_old');
%%
PAfilename= "D:\서울대\5-1\intern\survey_cop.xlsx";
load("C:\Users\kate5\Downloads\rel.mat"); %rel acc
% load("D:\서울대\5-1\intern\all_parti.mat") %raw acc
% rel = all_parti; clear all_parti; 

xaxis_limit=1;

datalabel = xlsread(PAfilename, 3, 'A:A');
numRows = numel(datalabel);

sex = readcell(PAfilename, 'Sheet', 3, 'Range',sprintf('C3:C%d', numRows + 2));
age= xlsread(PAfilename, 3, sprintf('D3:D%d', numRows + 2));
PA= xlsread(PAfilename, 3, sprintf('AG3:AG%d', numRows + 2));
education= xlsread(PAfilename, 3, sprintf('E3:E%d', numRows + 2));
subj_nums=xlsread(PAfilename, 3, sprintf('B3:B%d', numRows + 2));
%% sort out only those with behavioral data

%exclude those with 0 PA
zeroPA_indices = find(isnan(PA));
%if you don't wish to look at PA
% non_zero_indices = 1:numel(allfiles);
non_zero_indices = setdiff(1:numel(PA), zeroPA_indices);
PA_non_zero = PA(non_zero_indices);
subj_nums_non_zero=subj_nums(non_zero_indices);
sex_non_zero = sex(non_zero_indices);
age_non_zero = age(non_zero_indices);
education_non_zero = education(non_zero_indices);

rf = find(ismember(opartici,subj_nums_non_zero));
volumes_non_zero = volumes(rf);
oICV_non_zero = ICVs(rf);

row = find(ismember(subj_nums_non_zero, opartici));
osex_non_zero=sex_non_zero(row);
oage_non_zero=age_non_zero(row);
oPA_non_zero=PA_non_zero(row);
oeducation_non_zero=education_non_zero(row);
osubj_nums_non_zero=subj_nums_non_zero(row);

% ratio = hippo ./ volumes;% set the ratio of hipp vol to whole brain vol
% ratio_non_zero= ratio(non_zero_indices);
% hpf_ratio_non_zero=hpf_ratio(non_zero_indices);
% pfc_non_zero=pfc(non_zero_indices);
% hippo_non_zero=hippo(non_zero_indices);

osex_non_zero = strrep(osex_non_zero, '남', 'M');
osex_non_zero = strrep(osex_non_zero, '여', 'F');
osex_binary = strcmp(osex_non_zero, 'F');

[relp, sort_idx] = sort(rel.participant);
rel_filtered = rel(sort_idx, :);
[com,non]=ismember(rel_filtered.participant, opartici);
orel_filtered=rel_filtered(com, :);
row = find(ismember(orel_filtered.participant, opartici));
orel_filtered=orel_filtered(row,:);

%for mid
% mnon_zero_indices = 1:numel(mallfiles);
mrow = find(ismember(subj_nums_non_zero, mpartici));
msex_non_zero=sex_non_zero(mrow);
mage_non_zero=age_non_zero(mrow);
mPA_non_zero=PA_non_zero(mrow);
meducation_non_zero=education_non_zero(mrow);
msubj_nums_non_zero=subj_nums_non_zero(mrow);

mrow = find(ismember(mpartici,subj_nums_non_zero));
mvolumes_non_zero = mvolumes(mrow);
mICV_non_zero = mICVs(mrow);

msex_non_zero = strrep(msex_non_zero, '남', 'M');
msex_non_zero = strrep(msex_non_zero, '여', 'F');
msex_binary = strcmp(msex_non_zero, 'F');

[com,non]=ismember(rel_filtered.participant, mpartici);
mrel_filtered=rel_filtered(com, :);
mrow = find(ismember(mrel_filtered.participant, mpartici));
mrel_filtered=mrel_filtered(mrow,:);
%% ratio reg out
data = table(hippo,age, sex_binary, education, ICVs, ...
             'VariableNames', {'Hippo','Age', 'Sex', 'Education', 'ICV'});
mdl = fitlm(data, 'linear', 'ResponseVar', 'Hippo', 'PredictorVars', {'Sex', 'Education','Age'});
hp_Fitted=mdl.Fitted;

data = table(pfc,age, sex_binary, education, ICVs, ...
             'VariableNames', {'PFC','Age', 'Sex', 'Education', 'ICV'});
mdl = fitlm(data, 'linear', 'ResponseVar', 'PFC', 'PredictorVars', {'Sex', 'Education','Age'});
PFC_Fitted=mdl.Fitted;

fitted_ratio=PFC_Fitted./hp_Fitted;
corr(fitted_ratio,age)
%% for ratio and age, all ages: data label 1,2,4,5
age1=xlsread(exfilename, 3, sprintf('C1:C39'));
idx=[];
age1_partici=xlsread(exfilename, 3, sprintf('A1:A39'));
AGE = age;
AGE=[AGE;age1];
SUBJ=[subj_nums;age1_partici];
AGEs = table(SUBJ, AGE,zeros(numel(AGE),1), 'VariableNames', {'SUBJ', 'AGE','PH'});
for i = 1:length(opartici)
    idx = find(AGEs.SUBJ(i) == opartici);
        AGEs.PH(i) = (hpf_ratio(idx));
end

%% age
% u shape -parahippocampal gyrus

[b,PAoutliers]=rmoutliers(PA_non_zero);
nooutlierPArow=find(ismember(PA_non_zero, b));
ohpf_ratio_non_zero=hpf_ratio_non_zero(nooutlierPArow);

mdl = fitlm(AGEs.AGE,AGEs.PH,'quadratic')
figure('color','white')
figure1=plot(mdl);
xlabel('Age');
ylabel('para/hp');
figure1(1).Marker="o";
figure1(1).MarkerFaceColor="k";
figure1(1).MarkerFaceAlpha='flat';
figure1(1).MarkerEdgeColor='none';
figure1(1).Marker

box off
legend('off');
%%
[r,w]=rmoutliers(PA_non_zero);
rows=find(w==0);
rel_filtered=rel_filtered(rows,:);
PA_non_zero=PA_non_zero(rows,:);
age_non_zero=age_non_zero(rows,:);

variant=rel_filtered.accuracy1;%where are you looking at
ROI=PA_non_zero;
% age=age([1:159,162:165,167:174],:);
[R,P] = corr(variant, ROI);
sprintf("correlation:%f P value:%f",R,P)
figure;
mdl = fitlm(variant, ROI);
r_squared = mdl.Rsquared.Ordinary;
p_value = mdl.Coefficients.pValue;
plot(variant, ROI,'.')
hold on; 
x_range = [min(variant), max(variant)]; % x 축 범위
y_pred = mdl.Coefficients.Estimate(1) + mdl.Coefficients.Estimate(2) * x_range; 

plot(x_range, y_pred, 'k-', 'LineWidth', 2);
 xlabel("age")
 ylabel("dlPFC")
 % title("What")
hold off
%%
PA_mid_old= table(age,rel_filtered.accuracy1,rel_filtered.accuracy2,rel_filtered.accuracy3,'VariableNames', {'age','rel_accuracy1','rel_accuracy2','rel_accuracy3'});
www=["what","where","when"];
pv=[];
for i = 2:4   % www
    figure;
    [r, p] = corr(PA_mid_old{:,1}, PA_mid_old{:,i});
    pv=[pv;p];
    fprintf("%s: p = %f, r =  %f\n",www(i-1),p,r)
    mdl = fitlm(PA_mid_old{:,1}, PA_mid_old{:,i});
    plot(PA_mid_old{:,1}, PA_mid_old{:,i}, '.')
    hold on
    x_range = [min(PA_mid_old{:,1}), max(PA_mid_old{:,1})];
    y_pred = mdl.Coefficients.Estimate(1) + mdl.Coefficients.Estimate(2) * x_range; 
    plot(x_range, y_pred, 'k-', 'LineWidth', 2);
    xlabel('age');
    ylabel('accuracy');
    hold off
    if i == 2
        title('what')
    elseif i == 3
        title('where')
    elseif i==4
        title('when')
    end
end

%% remove PA outliers
[b,PAoutliers]=rmoutliers(PA_non_zero);
nooutlierPArow=find(ismember(PA_non_zero, b));

rel_filtered=rel_filtered(nooutlierPArow,:);

PAnp = table(rel_filtered.participant, subj_nums_non_zero, b, rel_filtered.accuracy1, rel_filtered.accuracy2, rel_filtered.accuracy3, ...
             'VariableNames', {'sh','SUBJ', 'PA', 'What', 'Where', 'When'});

PAnp_sorted = sortrows(PAnp, 'PA');
n = height(PAnp);
group_size = floor(n / 4);
remainder = mod(n, 4);

groups = cell(4, 3);
start_idx = 1;
for i = 1:4
    if i <= remainder
        end_idx = start_idx + group_size;
    else
        end_idx = start_idx + group_size - 1;
    end
    
    groups{i, 1} = PAnp_sorted.SUBJ(start_idx:end_idx);
    groups{i, 2} = PAnp_sorted.PA(start_idx:end_idx);
    groups{i, 3} = PAnp_sorted.When(start_idx:end_idx);
    start_idx = end_idx + 1;
end

means = zeros(4, 1);
ses = zeros(4, 1);

for i = 1:4
    data = groups{i, 3};
    means(i) = mean(data);
    ses(i) = std(data) / sqrt(length(data));
end

figure('Color', 'white');
bar(means);
hold on;
errorbar(1:4, means, ses, 'k', 'LineStyle', 'none');
xlabel('PA Group');
ylabel('Mean Task Accuracy');
title('When, rel Accuracy, PA outliers removed');
hold off;
%%
PA_mid_old= table(subj_nums_non_zero, b,age_non_zero,education_non_zero,hpf_ratio_non_zero,rel_filtered.accuracy1,rel_filtered.accuracy2,rel_filtered.accuracy3,'VariableNames', {'person_idx', 'PA','Age','Edu','ratio','rel_accuracy1','rel_accuracy2','rel_accuracy3'});
% how about I make a quad funct with PA and age
variables = {'PA', 'Age', 'Edu', 'ratio'};
combinations = nchoosek(variables, 2);

save_path = 'D:\서울대\5-1\intern\figure\';

for i = 1:size(combinations, 1)
    predictorVars = combinations(i, :);

    mdl = fitlm(PA_mid_old, 'quadratic', 'ResponseVar', 'rel_accuracy2', 'PredictorVars', predictorVars)

    figure('Color', 'white');
    figure1 = plot(mdl);
    xlabel(predictorVars{1});
    ylabel('relative accuracy');
    figure1(1).Marker = '.';
    figure1(1).MarkerFaceColor = 'k';
    % figure1(1).MarkerFaceAlpha = 1;
    figure1(1).MarkerEdgeColor = 'none';
    
    box off;
    legend('off');

    % save_filename = [save_path, 'Where ', predictorVars{1}, '_', predictorVars{2}, '.png'];
    % saveas(gcf, save_filename);

    close(gcf);
end

%% rel acc
www=["what","where","when"];
for i = 3:5   % www
    figure;
    % 점
    [r, p] = corr(PA_mid_old{:,2}, PA_mid_old{:,i});
    fprintf("%s: p = %f, r =  %f\n",www(i-2),p,r)
    mdl = fitlm(PA_mid_old{:,2}, PA_mid_old{:,i});
    plot(PA_mid_old{:,2}, PA_mid_old{:,i}, '.')
    hold on
    % 선
    x_range = [min(PA_mid_old{:,2}), max(PA_mid_old{:,2})];
    y_pred = mdl.Coefficients.Estimate(1) + mdl.Coefficients.Estimate(2) * x_range; 
    plot(x_range, y_pred, 'k-', 'LineWidth', 2);
    hold off
    if i == 3
        title('what')
    elseif i == 4
        title('where')
    else
        title('when')
    end
end
%% remove rel acc outlier
for i = 3:5

    if i==3
    [b, outlier] = rmoutliers(PA_mid_old(:,i));
    nooutlierPArow=find(ismember(PA_mid_old(:,i), b));
    rel_accuracy1=rel_accuracy1(nooutlierPArow);
    out=rel_accuracy1;
    rPA=PA_mid_old.PA(nooutlierPArow);
    end

    if i==4
    [b, outlier] = rmoutliers(PA_mid_old(:,i));
    nooutlierPArow=find(ismember(PA_mid_old(:,i), b));
    rel_accuracy2=rel_accuracy2(nooutlierPArow);
    out=rel_accuracy2;
    rPA=PA_mid_old.PA(nooutlierPArow);
    end

    if i==5
    [b, outlier] = rmoutliers(PA_mid_old(:,i));
    nooutlierPArow=find(ismember(PA_mid_old(:,i), b));
    rel_accuracy3=rel_accuracy3(nooutlierPArow);
    out=rel_accuracy3;
    rPA=PA_mid_old.PA(nooutlierPArow);
    end

    figure;
    [r, p] = corr(rPA,out);
    fprintf("p = %f, r =  %f\n",p,r)
    mdl = fitlm(rPA, out);
    plot(rPA, out, '.')
    hold on
    x_range = [min(rPA), max(rPA)];
    y_pred = mdl.Coefficients.Estimate(1) + mdl.Coefficients.Estimate(2) * x_range; 
    plot(x_range, y_pred, 'k-', 'LineWidth', 2);
    hold off
    if i == 3
        title('what')
    elseif i == 4
        title('where')
    else
        title('when')
    end
end
%% regressing out age, education, and gender
% sex_binary = strcmp(sex_non_zero, 'F');
% l_hippo_non_zero=l_hippo(non_zero_indices);
% r_hippo_non_zero=r_hippo(non_zero_indices);
% l_pfc_non_zero=l_pfc(non_zero_indices);
% r_pfc_non_zero=r_pfc(non_zero_indices);
% l_hpfratio_non_zero=l_hpfratio(non_zero_indices);
% r_hpfratio_non_zero=r_hpfratio(non_zero_indices);

% data = table(hippo,l_pfc_non_zero,r_pfc_non_zero,age_non_zero, sex_binary, education_non_zero, ICV_non_zero, PA_non_zero,...
%              'VariableNames', {'dl PFC','l dlPFC','r dlPFC','Age', 'Sex', 'Education', 'ICV','PA'});
data = table(pfc,age_non_zero, sex_binary, education_non_zero, ICV_non_zero, PA_non_zero,...
             'VariableNames', {'dl PFC','Age', 'Sex', 'Education', 'ICV','PA'});
nm=['What ','Where ',"When "];
save_path ="D:\서울대\5-1\intern\figure\";

% [k,rw]=rmoutliers(data.PA);
% rows=find(rw==0);
% data = data(rows, :);

mdl = fitlm(data, 'linear', 'ResponseVar', 'dl PFC', 'PredictorVars', {'Sex', 'Education','Age','ICV'})
for i=1:3
[r, p] = corr(rel_filtered.(['accuracy' num2str(i)]),mdl.Fitted);
% [r, p] = corr(data.PA,mdl.Fitted);
fprintf('%s, r: %.4f, p-value: %.4f\n', nm{i}, r, p);
% fprintf('r: %.4f, p-value: %.4f\n', r, p);
figure;
plot(rel_filtered.(['accuracy' num2str(i)]),mdl.Fitted,'.')
% plot(data.PA,mdl.Fitted,'.')
hold on; 
x_range = [min(rel_filtered.(['accuracy' num2str(i)])), max(rel_filtered.(['accuracy' num2str(i)]))];
md = fitlm(rel_filtered.(['accuracy' num2str(i)]),mdl.Fitted);
% x_range = [min(data.PA), max(data.PA)];
% md = fitlm(data.PA,mdl.Fitted);

y_pred = md.Coefficients.Estimate(1) + md.Coefficients.Estimate(2) * x_range; 
plot(x_range, y_pred, 'k-', 'LineWidth', 2);
title(nm{i});
xlabel('accuracy');
ylabel('phg vol');
hold off
    % save_filename = strcat(save_path, nm{i}, 'dl pfc_age.png');
    save_filename = strcat(save_path, 'r hp and PA.png');
    % saveas(gcf, save_filename);;
end
%% other control regions
% what = [10,12,4,20];%for occipital
% what=[28,7,30,21,24];%for parietal 
what=[3,11,13,26,27,31];%pfc
% what=[3,26,27];%dPFC
% what=[27,26,3,17,19,18,11,13,16,23,31];%frontal
% what=[29,8,14,1,6,33,5,32,15];%temporal
% what=[25,2,22,9];%cingulate
% what=[12,27];%hippo
% what=[23,27];%sma
%% with acc
nm = ["what", "where", "when"];
names = {};
rv = [];
pv = [];
task={};
mnames = {};
mrv = [];
mpv = [];
mtask={};
where='hp';

    for P = 1:3
        for i = 1:numel(r_allfiles)
            r_roi(i, 1) = sum(r_allfiles{i, 1}.GrayVol(what(:), 1));
            l_roi(i, 1) = sum(l_allfiles{i, 1}.GrayVol(what(:), 1));
        end

        for i = 1:numel(mr_allfiles)
            mr_roi(i, 1) = sum(mr_allfiles{i, 1}.GrayVol(what(:), 1));
            ml_roi(i, 1) = sum(ml_allfiles{i, 1}.GrayVol(what(:), 1));
        end

        roi = (r_roi+l_roi)./ICVs;
        mroi = (mr_roi+ml_roi)./mICVs;
        % ROI_adj=adjustROI(roi,ICV_non_zero,non_zero_indices);

        row = find(ismember(opartici,subj_nums_non_zero));
        non_zero_roi=roi(row);
        norel_filtered=orel_filtered(row,:);

        [b,PAoutliers]=rmoutliers(oPA_non_zero);
        rrow=find(ismember(oPA_non_zero, b));
        knorel_filtered=norel_filtered(rrow,:);

        data = table(non_zero_roi(rrow), oage_non_zero(rrow), osex_binary(rrow), oeducation_non_zero(rrow), oICV_non_zero(rrow), oPA_non_zero(rrow), ...
                     'VariableNames', {'ROI', 'Age', 'Sex', 'Education', 'ICV', 'PA'});
        mdl = fitlm(data, 'linear', 'ResponseVar', 'ROI', 'PredictorVars', {'Sex','Age','Education'});
        [r, p] = corr(knorel_filtered.(['accuracy' num2str(P)]), oPA_non_zero(rrow));

        mrow = find(ismember(mpartici,subj_nums_non_zero));
        non_zero_mroi=mroi(mrow);
        mnorel_filtered=mrel_filtered(mrow,:);

        [b,PAoutliers]=rmoutliers(mPA_non_zero);
        mrrow=find(ismember(mPA_non_zero, b));
        kmnorel_filtered=mnorel_filtered(mrrow,:);

        mdata = table(non_zero_mroi(mrrow), mage_non_zero(mrrow), msex_binary(mrrow), meducation_non_zero(mrrow), mICV_non_zero(mrrow), mPA_non_zero(mrrow), ...
                     'VariableNames', {'ROI', 'Age', 'Sex', 'Education', 'ICV', 'PA'});
        mmdl = fitlm(mdata, 'linear', 'ResponseVar', 'ROI', 'PredictorVars', {'Sex','Age','Education'});
        [mr, mp] = corr(kmnorel_filtered.(['accuracy' num2str(P)]),mPA_non_zero(mrrow));
        % [r, p] = corr(rel_filtered.(['accuracy' num2str(P)]), ROI_adj);

        names = vertcat(names,where); 
        task = vertcat(task, nm{P}); 
        rv = vertcat(rv, r);
        pv = vertcat(pv, p);

        mnames = vertcat(mnames,where); 
        mtask = vertcat(mtask, nm{P}); 
        mrv = vertcat(mrv, mr);
        mpv = vertcat(mpv, mp);

        figure;
        hold on;
        plot(oPA_non_zero(rrow),knorel_filtered.(['accuracy' num2str(P)]),'b.')
        plot(mPA_non_zero(mrrow),kmnorel_filtered.(['accuracy' num2str(P)]),'r.')
        % plot(rel_filtered.(['accuracy' num2str(P)]), ROI_adj,'.')

        x_range = [min(oPA_non_zero(rrow)), max(oPA_non_zero(rrow))];
        md = fitlm(oPA_non_zero(rrow),knorel_filtered.(['accuracy' num2str(P)]));

        mx_range = [min(mPA_non_zero(mrrow)), max(mPA_non_zero(mrrow))];
        mmd = fitlm(mPA_non_zero(mrrow),kmnorel_filtered.(['accuracy' num2str(P)]));
        % md = fitlm(rel_filtered.(['accuracy' num2str(P)]),ROI_adj);
        y_pred = md.Coefficients.Estimate(1) + md.Coefficients.Estimate(2) * x_range; 
        my_pred = mmd.Coefficients.Estimate(1) + mmd.Coefficients.Estimate(2) * x_range; 

        % if P==1 || P==2
        plot(x_range,y_pred,'b-','LineWidth',2);
        plot(mx_range,my_pred,'r-','LineWidth',2);
        % elseif P==3
        % plot(x_range,y_pred,'b-','LineWidth',2);
        % end
        hold off;

        set(gca, 'FontSize', 14);
        xlabel('PA','FontSize', 18);
        ylabel('acc','FontSize', 18); 
        % title(nm{P},'FontSize', 15);
    end

tab = table(names, task,rv, pv, 'VariableNames', {'Name','task', 'r', 'p'})

disp('mid')
mtab = table(mnames, mtask,mrv, mpv, 'VariableNames', {'Name','task', 'r', 'p'})
PAnp = table(knorel_filtered.participant, oPA_non_zero(rrow), knorel_filtered.accuracy1, knorel_filtered.accuracy2, knorel_filtered.accuracy3, ...
             'VariableNames', {'sh', 'PA', 'What', 'Where', 'When'});

mPAnp = table(kmnorel_filtered.participant, mPA_non_zero(mrrow), kmnorel_filtered.accuracy1, kmnorel_filtered.accuracy2, kmnorel_filtered.accuracy3, ...
             'VariableNames', {'sh', 'PA', 'What', 'Where', 'When'});

PAnp_sorted = sortrows(PAnp, 'PA');
n = height(PAnp);
group_size = floor(n / 4);
remainder = mod(n, 4);

groups = cell(4, 3);
start_idx = 1;
for i = 1:4
    if i <= remainder
        end_idx = start_idx + group_size;
    else
        end_idx = start_idx + group_size - 1;
    end
    

    groups{i, 1} = PAnp_sorted.What(start_idx:end_idx);
    groups{i, 2} = PAnp_sorted.Where(start_idx:end_idx);
    groups{i, 3} = PAnp_sorted.When(start_idx:end_idx);

    start_idx = end_idx + 1;
end

means = zeros(4, 3);
ses = zeros(4, 3);

for k=1:3
    if k==1%what
    for i = 1:4
    data = groups{i, 1};
    means(i,k) = mean(data);
    ses(i,k) = std(data) / sqrt(length(data));
    end
figure('Color', 'white');
bar(means(:,1));
hold on;
errorbar(1:4, means(:,1), ses(:,1), 'k', 'LineStyle', 'none');
xlabel('PA Group');
ylabel('Mean Task Accuracy');
title('What, rel Accuracy, PA outliers removed');
hold off;

    elseif k==2%where
        for i = 1:4
    data = groups{i, 2};
    means(i,k) = mean(data);
    ses(i,k) = std(data) / sqrt(length(data));    
        end
figure('Color', 'white');
bar(means(:,2));
hold on;
errorbar(1:4, means(:,2), ses(:,2), 'k', 'LineStyle', 'none');
xlabel('PA Group');
ylabel('Mean Task Accuracy');
title('Where, rel Accuracy, PA outliers removed');
hold off;

    elseif k==3%where
        for i = 1:4
    data = groups{i, 3};
    means(i,k) = mean(data);
    ses(i,k) = std(data) / sqrt(length(data));    
        end
figure('Color', 'white');
bar(means(:,3));
hold on;
errorbar(1:4, means(:,3), ses(:,3), 'k', 'LineStyle', 'none');
xlabel('PA Group');
ylabel('Mean Task Accuracy');
title('When, rel Accuracy, PA outliers removed');
hold off;
end
end

%for mid

PAnp_sorted = sortrows(mPAnp, 'PA');
n = height(mPAnp);
group_size = floor(n / 4);
remainder = mod(n, 4);

mgroups = cell(4, 3);
start_idx = 1;
for i = 1:4
    if i <= remainder
        end_idx = start_idx + group_size;
    else
        end_idx = start_idx + group_size - 1;
    end
    
    mgroups{i, 1} = PAnp_sorted.What(start_idx:end_idx);
    mgroups{i, 2} = PAnp_sorted.Where(start_idx:end_idx);
    mgroups{i, 3} = PAnp_sorted.When(start_idx:end_idx);

    start_idx = end_idx + 1;
end

mmeans = zeros(4, 3);
mses = zeros(4, 3);

for k=1:3
    if k==1%what
    for i = 1:4
    data = mgroups{i, 1};
    mmeans(i,k) = mean(data);
    mses(i,k) = std(data) / sqrt(length(data));
    end
figure('Color', 'white');
bar(mmeans(:,1));
hold on;
errorbar(1:4, mmeans(:,1), mses(:,1), 'k', 'LineStyle', 'none');
xlabel('PA Group');
ylabel('Mean Task Accuracy');
title('What, rel Accuracy, PA outliers removed');
hold off;

    elseif k==2%where
        for i = 1:4
    data = mgroups{i, 2};
    mmeans(i,k) = mean(data);
    mses(i,k) = std(data) / sqrt(length(data));    
        end
figure('Color', 'white');
bar(mmeans(:,2));
hold on;
errorbar(1:4, mmeans(:,2), mses(:,2), 'k', 'LineStyle', 'none');
xlabel('PA Group');
ylabel('Mean Task Accuracy');
title('Where, rel Accuracy, PA outliers removed');
hold off;

    elseif k==3%where
        for i = 1:4
    data = mgroups{i, 3};
    mmeans(i,k) = mean(data);
    mses(i,k) = std(data) / sqrt(length(data));    
        end
figure('Color', 'white');
bar(mmeans(:,3));
hold on;
errorbar(1:4, mmeans(:,3), mses(:,3), 'k', 'LineStyle', 'none');
xlabel('PA Group');
ylabel('Mean Task Accuracy');
title('When, rel Accuracy, PA outliers removed');
hold off;
end
end
%%
pgo=1; %for PA, set this value to 1. For age, 0. groups are for old, mgroups are for middle-aged.
P=3;%what is 1, where is 2, when is 3

pjohap=[1 2;1 3;1 4;2 3;2 4;3 4];

if pgo==1
for i=1:size(pjohap,1)
    a =groups{pjohap(i,1),P}; %in a group, 4 is for PA, 3 is for WWWs
    b = groups{pjohap(i,2),P};
    %what
    x_1 = a;
    x_2 = b;
    p = vartest2(x_1,x_2);
    if p == 0
        disp('ttest0')
        [~,p] = ttest2(x_1,x_2)
    else
        disp('ttest1')
        [~,p] = ttest2(x_1,x_2,'Vartype','unequal')
    end
end

elseif pgo==0
for i=1:4
a =groups{i,P};
b = mgroups{i,P};
%what
x_1 = a;
x_2 = b;
p = vartest2(x_1,x_2);
if p == 0
    disp('ttest0')
    [~,p] = ttest2(x_1,x_2)
else
    disp('ttest1')
    [~,p] = ttest2(x_1,x_2,'Vartype','unequal')
end
end

end
%% for looking at specific ROIs
nm = ["what", "where", "when"]; P=3;
rv=[];pv=[];names=[];
for h=1:numel(l_allfiles{1,1}.StructName)
 for i = 1:numel(r_allfiles)
            r_roi(i, 1) = sum(r_allfiles{i, 1}.GrayVol(h, 1));
            l_roi(i, 1) = sum(l_allfiles{i, 1}.GrayVol(h, 1));
        roi = (r_roi+l_roi);
 end
        ROI_adj=adjustROI(roi,ICV_non_zero,non_zero_indices);
        data = table(ROI_adj, age, sex_binary, education_non_zero, ICV_non_zero, PA_non_zero, ...
                     'VariableNames', {'ROI', 'Age', 'Sex', 'Education', 'ICV', 'PA'});
        % mdl = fitlm(data, 'linear', 'ResponseVar', 'ROI', 'PredictorVars', {'Sex', 'Education','ICV','Age'});
        [r, p] = corr(rel_filtered.(['accuracy' num2str(P)]), ROI_adj);
        pv=vertcat(pv,p);
        rv=vertcat(rv,r);
        name=l_allfiles{1,1}.StructName(h,1);
        names = vertcat(names,name); 
        [~,sorted]=sort(pv,'ascend');
        tableof=table(names(sorted),rv(sorted), pv(sorted), ...
     'VariableNames', {'Name', 'r', 'p'});
end
        % figure;
        % plot(rel_filtered.(['accuracy' num2str(P)]), ROI_adj,'.')
        % hold on;
        % x_range = [min(rel_filtered.(['accuracy' num2str(P)])), max(rel_filtered.(['accuracy' num2str(P)]))];
        % md = fitlm(rel_filtered.(['accuracy' num2str(P)]),ROI_adj);
        % y_pred = md.Coefficients.Estimate(1) + md.Coefficients.Estimate(2) * x_range; 
        % plot(x_range,y_pred,'k-','LineWidth',2);
        % hold off;
        % xlabel('acc');
        % ylabel('vol');
        % title(nm{P});
%% with age or PA
factor='Age';
where='HP';
        for i = 1:numel(mr_allfiles)
            r_roi(i, 1) = sum(mr_allfiles{i, 1}.GrayVol(what(:), 1));
            l_roi(i, 1) = sum(ml_allfiles{i, 1}.GrayVol(what(:), 1));
        end
        % roi = (r_roi+l_roi)./mICVs;
        roi = mhippo./mICVs;

        data = table(roi, age(mallok), msex_binary(mallok), education(mallok), mICVs, PA(mallok), ...
                     'VariableNames', {'ROI', 'Age', 'Sex', 'Education', 'ICV', 'PA'});
        mdl = fitlm(data, 'linear', 'ResponseVar', 'ROI', 'PredictorVars', {'Sex', 'Education'});
        [r, p] = corr(data.(factor), mdl.Fitted);
        % [r, p] = corr(data.(factor), mdl.Fitted);
        fprintf('%s r: %.4f, p-value: %.4f\n', where, r, p);

figure;
plot(data.(factor),mdl.Fitted,'.')
% plot(data.(factor),roi,'.')
hold on; 
x_range = [min(data.(factor)), max(data.(factor))];
md = fitlm(data.(factor),mdl.Fitted);

%polynomial fit
% md = fitlm(data.(factor),mdl.Fitted,'poly2')
% figure;
% set(gca,'FontSize',13)
% xlabel('Age','FontSize',25)
% ylabel('PFC/ICV','FontSize',25)
% plot(md)

y_pred = md.Coefficients.Estimate(1) + md.Coefficients.Estimate(2) * x_range; 
plot(x_range, y_pred, 'k-', 'LineWidth', 2);
% title(nm{i});       
set(gca, 'FontSize', 13);
        xlabel(factor,'FontSize', 25);
        ylabel('HP/ICV','FontSize', 25); 
hold off
%%
mdl_hp_age = fitlm(age, hp);
residuals_hp = mdl_hp_age.Residuals.Raw;
pfc_hp_ratio = pfc ./ hp;
mdl_ratio_age = fitlm(age, pfc_hp_ratio);
%% residual per participant
for i=1:numel(allfiles)
dud=(allfiles{i,1}.Volume_mm3(5,1)+allfiles{i,1}.Volume_mm3(23,1));
end
% pfc_adj=adjustROI(pfc,ICV_non_zero,1:83);
% hp_adj=adjustROI(hippo,ICV_non_zero,1:83);

pfc_volume_ratio = pfc ./ ICVs;
hp_volume_ratio = hippo ./ ICVs;

pfcmean=mean(pfc_volume_ratio);
hpmean=mean(hp_volume_ratio);

deltapfc=pfc_volume_ratio   -pfcmean;
deltahp=hp_volume_ratio -hpmean;

data = table(deltapfc, deltahp, age_non_zero, sex_binary, education_non_zero, ICV_non_zero, PA_non_zero, ...
             'VariableNames', {'pfc','hp','Age', 'Sex', 'Education', 'ICV', 'PA'});
% md = fitlm(data, 'linear', 'ResponseVar', 'hp', 'PredictorVars', {'Sex', 'Education','Age'});
% mdl = fitlm(data, 'linear', 'ResponseVar', 'pfc', 'PredictorVars', {'Sex', 'Education','Age'});
  
figure;
% plot(md.Fitted,mdl.Fitted,'.')
plot(data.hp,data.pfc,'.')
hold on; 
% x_range = [min(md.Fitted), max(md.Fitted)];
x_range = [min(data.hp), max(data.hp)];
% Md = fitlm(md.Fitted,mdl.Fitted);
Md = fitlm(data.hp,data.pfc);
y_pred = Md.Coefficients.Estimate(1) + Md.Coefficients.Estimate(2) * x_range; 
plot(x_range, y_pred, 'k-', 'LineWidth', 2);
xlabel('\Delta HP');
ylabel('\Delta pfc');
hold off;
[r,p]=corr(data.hp,data.pfc)

rv=[];pv=[];names=[];

onearea=[9;10;11;14;33;34;37;40];
oneName=[{'3rd-Ventricle'},{'4th-Ventricle'},{'Brain-Stem'},{'CSF'},{'5th-Ventricle'},{'WM-hypointensities'},{'non-WM-hypointensities'},{'Optic-Chiasm'}];
ccarea=41:45;
pairs = [3 21; 4 22; 5 23; 6 24; 7 25; 8 26 ;13 28; 15 29;1 19;2 20;16 30;17 31;18 32;35 36;38 39];
Name= [{'Cerebellum-White-Matter'};{'Cerebellum-Cortex'};{'Thalamus'};{'Caudate'};{'Putamen'};{'Pallidum'};{'Amygdala'};{'Accumbens-area'};{'Lateral-Ventricle'};{'Inf-Lat-Vent'};{'VentralDC'};{'vessel'};{'choroid-plexus'};{'WM-hypointensities'};{'non-WM-hypointensities'}];

for h=1: numel(onearea)   
 for i = 1:numel(l_allfiles)
        r_roi(i, 1) = allfiles{i, 1}.Volume_mm3(onearea(h), 1);
        roi= r_roi;
 end
roi_volume_ratio = roi ./ volumes;
hp_volume_ratio = hippo ./ volumes;

roimean=mean(roi_volume_ratio);
hpmean=mean(hp_volume_ratio);

deltaroi=roi_volume_ratio-roimean;
deltahp=hp_volume_ratio-hpmean;

 [r, p] = corr(deltahp, deltaroi);
        pv=vertcat(pv,p);
        rv=vertcat(rv,r);
        % name=l_allfiles{1,1}.StructName(h,1);
        name=oneName{(h)};
        % names = vertcat(names,name); 
        names =[names;{name}]; 
end

for i = 1:numel(l_allfiles)
        r_roi(i, 1) = sum(allfiles{i, 1}.Volume_mm3(ccarea(:), 1));
        roi= r_roi;
end
roi_volume_ratio = roi ./ volumes;
hp_volume_ratio = hippo ./ volumes;

roimean=mean(roi_volume_ratio);
hpmean=mean(hp_volume_ratio);

deltaroi=roi_volume_ratio-roimean;
% deltahp=hp_volume_ratio-hpmean;

 [r, p] = corr(deltahp, deltaroi);
        pv=vertcat(pv,p);
        rv=vertcat(rv,r);
        % name=l_allfiles{1,1}.StructName(h,1);
        name={'CC'};
        names = vertcat(names,name); 


% for h=1:numel(l_allfiles{1,1}.StructName)
for h=1:numel(pairs(:,1))
 for i = 1:numel(l_allfiles)
        r_roi(i, 1) = (allfiles{i, 1}.Volume_mm3(pairs(h,1), 1)+allfiles{i, 1}.Volume_mm3(pairs(h,2), 1));
        % l_roi(i, 1) = (l_allfiles{i, 1}.GrayVol(h, 1));
        % roi = (r_roi+l_roi);
        roi= r_roi;
 end
        
roi_volume_ratio = roi ./ volumes;
hp_volume_ratio = hippo ./ volumes;

roimean=mean(roi_volume_ratio);
hpmean=mean(hp_volume_ratio);

deltaroi=roi_volume_ratio-roimean;
% deltahp=hp_volume_ratio-hpmean;

 [r, p] = corr(deltahp, deltaroi);
        pv=vertcat(pv,p);
        rv=vertcat(rv,r);
        % name=l_allfiles{1,1}.StructName(h,1);
        name=Name{(h)};
        % names = vertcat(names,name); 
        names =[names;{name}]; 
end
    [~,sorted]=sort(rv,'ascend');
    ctableof=table(names(sorted),rv(sorted), pv(sorted), ...
    'VariableNames', {'Name', 'r', 'p'})
%% function: ICV equation
function ROI_adj = adjustROI(roi, ICV_non_zero, non_zero_indices)
meanICV=mean(ICV_non_zero);        
rawROI=(roi);
rawROI_non_zero=rawROI(non_zero_indices);
mdl = fitlm(ICV_non_zero,rawROI_non_zero);
b=mdl.Coefficients.Estimate(2);
ROI_adj= rawROI_non_zero-b*(ICV_non_zero-meanICV);
end
