clear all; clc;
folderpath = 'D:\서울대\5-1\intern\brain_volume\';
filename = dir(fullfile(folderpath, '*_aseg.stats'));
allfiles = cell(length(filename), 1);
volumes = zeros(length(allfiles), 1); 
ICVs= zeros(length(allfiles), 1);
%1
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
    % 23번째 줄에서 뇌 총 부피 값 추출
    line = fgetl(fileID);
    volume = strsplit(line, ',');
    volume = strtrim(volume{4}); % 값 추출
    % mm^3 단위로 변환하여 저장
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

k=setdiff(who, {'allfiles', 'ICVs','volumes'});
clear(k{:});
%% cutoff and other covariates
xaxis_limit=1;
cutoff=9000;

PAfilename= "D:\서울대\5-1\intern\survey_cop.xlsx";
age= xlsread(PAfilename, 3, 'D3:D174');
PA= xlsread(PAfilename, 3, 'AG3:AG174');
sex= readcell(PAfilename, 'Sheet', 3, 'Range', 'C3:C174');
education= xlsread(PAfilename, 3, 'E3:E174');
hippo=zeros(size(PA, 1), 1); 
for i=1:size(PA,1)
hippo(i,1)=(allfiles{i,1}.Volume_mm3(12,1)+allfiles{i,1}.Volume_mm3(27,1));
end

%PA가 0인 사람들 제외
zeroPA_indices = find(isnan(PA));
non_zero_indices = setdiff(1:numel(PA), zeroPA_indices);
PA_non_zero = PA(non_zero_indices);
volumes_non_zero = volumes(non_zero_indices);
sex_non_zero = sex(non_zero_indices);
ICV_non_zero = ICVs(non_zero_indices);
ratio = hippo ./ volumes;% set the ratio of hipp vol to whole brain vol
ratio_non_zero= ratio(non_zero_indices);
hippo_non_zero=hippo(non_zero_indices);
age_non_zero = age(non_zero_indices); 
education_non_zero = education(non_zero_indices); 
    if xaxis_limit==1
    cutoffPA_indices = find(PA_non_zero <= cutoff);
    hippo_non_zero=hippo_non_zero(cutoffPA_indices);
    PA_non_zero = PA_non_zero(cutoffPA_indices);
    volumes_non_zero=volumes_non_zero(cutoffPA_indices);
    ICV_non_zero=ICV_non_zero(cutoffPA_indices);
    sex_non_zero= sex_non_zero(cutoffPA_indices);
    age_non_zero = age_non_zero(cutoffPA_indices);
    education_non_zero=education_non_zero(cutoffPA_indices);
    ratio_non_zero=ratio_non_zero(cutoffPA_indices);
    else
    end

sex_non_zero = strrep(sex_non_zero, '남', 'M');
sex_non_zero = strrep(sex_non_zero, '여', 'F');
%% volume adjustment using ICV
% b is the slope of the linear regression between Volume(Raw) and ICV->근데 이거 느낌상 b 구할 때 x축에 volume말고 ICV 들어가야할 거 같음.
meanICV=mean(ICV_non_zero);
mdl = fitlm(ICV_non_zero,volumes_non_zero);
b=mdl.Coefficients.Estimate(2);
volumes_non_zero_adj= volumes_non_zero-b*(ICV_non_zero-meanICV);
%% plot
% 산점도
ROI=ratio_non_zero;%where are you looking at
variant=education_non_zero;

[R,P] = corr(variant, ROI);
sprintf("correlation:%f P value:%f",R,P)
mdl = fitlm(variant, ROI);
r_squared = mdl.Rsquared.Ordinary;
p_value = mdl.Coefficients.pValue;
plot(variant, ROI,'.')
hold on; 
x_range = [min(variant), max(variant)]; % x 축 범위
y_pred = mdl.Coefficients.Estimate(1) + mdl.Coefficients.Estimate(2) * x_range; 

% 산점도- 회귀선 그리기
plot(x_range, y_pred, 'k-', 'LineWidth', 2);
hold off
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

% 각 조합에 대한 라벨 생성.. for character ones..
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

        idx = strcmp(edugroup_non_zero,edulabel{1,i}) & strcmp(agegroups,agelabel{1,k});
        % idx = strcmp(agegroups, agelabel{1,i}) & strcmp(PAgroups,PAlabel{1,j}) & strcmp(edugroup_non_zero,edulabel{1,k});
        % data{(i-1)*4 + (j-1)*2+ k} = volumes_non_zero(idx);
        data{(i-1)*2 + k} = hippo_non_zero(idx);

        n=numel(data{(i-1)*2 + k});
        % n=numel(data{(i-1)*4 + (j-1)*2+ k});
        % label{(i-1)*4 + (j-1)*2+ k} = repmat(sprintf('(%s,%s,%s)', agelabel{1,i},PAlabel{1,j},edulabel{1,k}),n, 1);
        label{(i-1)*2 + k} = repmat(sprintf('(%s,%s)', edulabel{1,i},agelabel{1,k}),n, 1);
    % end
end
end

% only for text ones
all_text = [];
for i = 1:numel(label)
    allt=cellstr(label{i, 1});
    all_text = vertcat(all_text,allt); % 각 셀의 문자열에서 양쪽 공백을 제거한 후 공백으로 결합
end
figure;
boxplot(cell2mat(data),all_text)
ylabel('volume (mm)^{3}');
% xlabel('age,PA');

%group 별 크기 확인
for n=1:numel(data)
fprintf('%s: %d/ ',label{n}(1,:),numel(data{n}))
end
%% kmeans?
% 데이터를 오름차순으로 정렬
PA_sorted = sort(education_non_zero);

% 데이터의 총 개수
num_data = numel(PA_sorted);

% 클러스터 수 설정
num_clusters = 2;

% 각 클러스터의 크기
cluster_size = floor(num_data / num_clusters);

% 각 클러스터의 경계값을 찾기
cutoff_values = zeros(num_clusters - 1, 2);
for i = 1:num_clusters-1
    max_val = PA_sorted(i * cluster_size);
    min_val = PA_sorted(i * cluster_size + 1);
    cutoff_values(i, :) = [max_val, min_val];
end

% 각 클러스터의 데이터 포인트 수 확인
for i = 1:num_clusters
    if i < num_clusters
        fprintf('Cluster %d: %d points\n', i, cluster_size);
    else
        fprintf('Cluster %d: %d points\n', i, num_data - (num_clusters-1) * cluster_size);
    end
end

% 각 클러스터의 경계값 출력
fprintf('Cutoff values between the clusters (max of current, min of next):\n');
disp(cutoff_values);


%% one way ANOVA
[p, tbl, stats] = anova1(cell2mat(data), all_text);
%% for getting rid of group 2 (for education might be too clumped)
tlqkf = find(~isnan(PA) & groupnum ~= 2);
PA_tlqkf=PA(tlqkf);
age_tlqkf=age(tlqkf);
education_tlqkf=education(tlqkf);

hippo_tlqkf=hippo(tlqkf);
ICV_tlqkf=ICVs(tlqkf);
ratio_tlqkf=ratio(tlqkf);
volumes_tlqkf=volumes(tlqkf);

cutoffPA_indices = find(PA_tlqkf <= cutoff);
age_tlqkf = age_tlqkf(cutoffPA_indices);
PA_tlqkf = PA_tlqkf(cutoffPA_indices);
education_tlqkf = education_tlqkf(cutoffPA_indices);

hippo_tlqkf=hippo_tlqkf(cutoffPA_indices);
ICV_tlqkf=ICV_tlqkf(cutoffPA_indices);
ratio_tlqkf=ratio_tlqkf(cutoffPA_indices);
volumes_tlqkf=volumes_tlqkf(cutoffPA_indices);

variant=education_tlqkf;
ROI=hippo_tlqkf;

[R,P] = corr(variant, ROI);
sprintf("correlation:%f P value:%f",R,P)
mdl = fitlm(variant, ROI);
r_squared = mdl.Rsquared.Ordinary;
p_value = mdl.Coefficients.pValue;
plot(variant, ROI,'.')
hold on; 
x_range = [min(variant), max(variant)]; % x 축 범위
y_pred = mdl.Coefficients.Estimate(1) + mdl.Coefficients.Estimate(2) * x_range; 

% 산점도- 회귀선 그리기
plot(x_range, y_pred, 'k-', 'LineWidth', 2);
hold off
%% other covariates
    tbl = table(PA_non_zero, volumes_non_zero, hippo_non_zero, age_non_zero, education_non_zero, 'VariableNames', {'PA', 'Volumes', 'Hippo', 'Age','Education'});
    mdl = fitlm(tbl, 'Hippo ~ Age+ Education + PA');  
    disp(mdl.Formula.LinearPredictor)
    plot(mdl);
    xlabel('Adjusted Whole Model');
    title('Added variable plot for whole model');
    legend('data points','fit','95% conf. bounds');
    ylabel('hippocampus vol(mm)^{3}');
    sprintf("R squared Ordinary:%f R squared Adjusted:%f P value:%f",mdl.Rsquared.Ordinary,mdl.Rsquared.Adjusted,mdl.ModelFitVsNullModel.Pvalue)

% for 3D graph and mesh
%     tbl = table(PA_non_zero, volumes_non_zero, age_non_zero, 'VariableNames', {'PA', 'Volumes', 'Age'});
%     mdl = fitlm(tbl, 'Volumes ~ PA + Age');
% 
% figure;
% scatter3(PA_non_zero, age_non_zero,  volumes_non_zero, 'filled');
% hold on;
% [PA_mesh, age_mesh] = meshgrid(linspace(min(PA_non_zero), max(PA_non_zero), 10), linspace(min(age_non_zero), max(age_non_zero), 10));
% volumes_pred = predict(mdl, table(PA_mesh(:), age_mesh(:), 'VariableNames', {'PA', 'Age'}));
% volumes_mesh = reshape(volumes_pred, size(PA_mesh));
% mesh(PA_mesh, age_mesh, volumes_mesh, 'FaceAlpha', 0.5);
% xlabel('Physical Activity');
% ylabel('Age');
% zlabel('Brain Volume');
% title('Brain Volume vs Physical Activity and Age');
% hold off;
%% for quad function regression
tbl = table(PA_non_zero, volumes_non_zero, hippo_non_zero, age_non_zero, education_non_zero, 'VariableNames', {'PA', 'Volumes', 'Hippo', 'Age','Education'});

% Create squared age variable
tbl.Age_squared = tbl.Age.^2;

% Fit the quadratic model
mdl_quad = fitlm(tbl, 'Hippo ~ Age + Age_squared + Age * Education + Volumes + PA');

% Define range for age and squared age
age_range = min(tbl.Age):0.1:max(tbl.Age);
age_range_squared = age_range.^2;

% Create predictor variables for the range
volumes_range = repmat(mean(tbl.Volumes), size(age_range));
pa_range = repmat(mean(tbl.PA), size(age_range));
education_range = repmat(mean(tbl.Education), size(age_range)); 

% Plot the actual data and the predicted curve

predicted_hippo_volume = predict(mdl_quad, table(age_range', age_range_squared', volumes_range', education_range', pa_range', 'VariableNames', {'Age', 'Age_squared', 'Volumes', 'Education', 'PA'}));
hold on;
plot(age_range, predicted_hippo_volume,'r-', 'DisplayName', 'Predicted Hippo Volume ');
plot(tbl.Age, tbl.Hippo, 'bo', 'DisplayName', 'Actual Data');
hold off;
xlabel('Age');
ylabel('Hippocampal Volume');
title('Age vs. Hippocampal Volume (Quadratic Fit)');
legend('show');
%% divide into different education levels

colors = {'r-', 'b-'};
colordot={'ro','bo'};
figure;

for i=1:numel(sexlabel)
edugroupforthis= find(strcmp(sex_non_zero,sexlabel{i}));

% if i==1
%     exclde=[81, 98, 97, 36];
%     edugroupforthis(exclde,:)=[];
% end

tbl = table(PA_non_zero(edugroupforthis), volumes_non_zero(edugroupforthis), hippo_non_zero(edugroupforthis), age_non_zero(edugroupforthis), education_non_zero(edugroupforthis), 'VariableNames', {'PA', 'Volumes', 'Hippo', 'Age','Education'});
tbl.Age_squared = tbl.Age.^2;
tbl.PA_squared = tbl.PA.^2;
mdl_quad = fitlm(tbl, 'Hippo ~ Age + Age_squared + Age * PA + Volumes + PA + Education');

% age_range = 0:0.1:100; %나이를 0에서 100까지 연장해보기
% age_range_squared = age_range.^2;
% volumes_range = repmat(mean(tbl.Volumes), size(age_range));
% pa_range = repmat(mean(tbl.PA), size(age_range));
% education_range = repmat(mean(tbl.Education), size(age_range)); 

age_range = min(age_non_zero):0.1:max(age_non_zero); 
age_range_squared = age_range.^2;
volumes_range = repmat(mean(tbl.Volumes), size(age_range));
PA_range = repmat(mean(tbl.PA), size(age_range));
education_range = repmat(mean(tbl.Education), size(age_range)); 

% predicted_hippo_volume = predict(mdl_quad, table(age_range', age_range_squared', volumes_range', education_range', pa_range', 'VariableNames', {'Age', 'Age_squared', 'Volumes', 'Education', 'PA'}));
predicted_hippo_volume = predict(mdl_quad, table(PA_range', age_range_squared', volumes_range', education_range', age_range', 'VariableNames', {'PA', 'Age_squared', 'Volumes', 'Education', 'Age'}));

hold on;
plot(age_range, predicted_hippo_volume, colors{i}, 'DisplayName', ['Predicted Hippo Volume of ' sexlabel{i}]);
plot(tbl.Age, tbl.Hippo, colordot{i}, 'DisplayName', ['Actual Data of ' sexlabel{i}]);
end
xlabel('PA');
ylabel('Hippocampal Volume');
title('PA vs. Hippocampal Volume (Quadratic Fit)');
legend('show');
hold off;
%% maybe a bar graph?
reg=0;%are you doing a regression?

bar_width = 0.5; 
group_spacing = 1.0; 
colors={'red','blue'};

figure;

for k=1:numel(PAlabel)
for i=1:numel(sexlabel)
edugroupforthis= find(strcmp(sex_non_zero,sexlabel{i}) & strcmp(PAgroups,PAlabel{k}));
tbl = table(PA_non_zero(edugroupforthis), volumes_non_zero(edugroupforthis), hippo_non_zero(edugroupforthis), age_non_zero(edugroupforthis), education_non_zero(edugroupforthis), 'VariableNames', {'PA', 'Volumes', 'Hippo', 'Age','Education'});
tbl.PA_squared = tbl.PA.^2;
tbl.Age_squared = tbl.Age.^2;
mdl_quad = fitlm(tbl, 'Hippo ~ PA + Age_squared + PA * Age + Volumes + Age + Education');

x_pos = (k-1) * (2 * bar_width + group_spacing) + (i-1) * (bar_width);

hold on;
if reg==1
h(i)=bar(x_pos, median(mdl_quad.Fitted), bar_width, 'FaceColor', colors{i});
std_error = std(mdl_quad.Fitted) / sqrt(length(mdl_quad.Fitted));
errorbar(x_pos, median(mdl_quad.Fitted), std_error, 'k', 'LineWidth', 1);
else
h(i)=bar(x_pos, median(hippo_non_zero(edugroupforthis)), bar_width, 'FaceColor', colors{i});
std_error = std(hippo_non_zero(edugroupforthis)) / sqrt(length(edugroupforthis));
errorbar(x_pos, median(hippo_non_zero(edugroupforthis)), std_error, 'k', 'LineWidth', 1);
end

plot(x_pos,hippo_non_zero(edugroupforthis),'.');
fprintf('%s,%s\n',sexlabel{i},PAlabel{k})
% disp(numel(edugroupforthis))

end
end
hold off;
set(gca, 'XTick', (0:numel(PAlabel)-1) * (2 * bar_width + group_spacing) + bar_width / 2);
set(gca, 'XTickLabel', PAlabel);%k를 따라감
xlabel('Education');
ylim([6000 11000])
ylabel('Hippocampal Volume');
legend(h,sexlabel);
hold off;

%% t test
group1=sexlabel;
groupofint=agelabel;

group1d=sex_non_zero;
groupintd=agegroups;

for k=1:numel(group1)
for i=1:numel(groupofint)
    edugroupforthis= find(strcmp(groupintd,groupofint{i}) & strcmp(group1d,group1{k}));
        if i > 1 % 보고 싶은 부분들 위주로 고치기
            prev_edugroupforthis = find(strcmp(groupintd, groupofint{i-1}) & strcmp(group1d, group1{k}));
            [~,p_value] = ttest2(hippo_non_zero(edugroupforthis), hippo_non_zero(prev_edugroupforthis));
            all_p_values = [all_p_values; p_value];
            fprintf('T-test between %s,%s and %s,%s: p-value = %f\n',groupofint{i-1}, group1{k}, groupofint{i}, group1{k}, p_value);
        end
end
end
%% which area has max corr?
variant=education;% which variant are you looking at
ratio=0;% are you choosing to look at the ratio of that region?

allvol=cell(14,1);
vol=zeros(size(variant, 1), 1); 
wheretopick = cell(14, 1);
wheretopick{2} = struct('name','Cerebellum','range',[3, 4, 21, 22]);
wheretopick{3} = struct('name','Thalamus','range',[5, 23]);
wheretopick{4} = struct('name','Caudate','range',[6,24]);
wheretopick{5} =struct('name','Putamen','range',[7,25]);
wheretopick{6} =struct('name','Pallidum','range',[8,26]);
wheretopick{7} =struct('name','Brainstem','range',11);
wheretopick{8} =struct('name','Hippocampus','range',[12,27]);
wheretopick{9} =struct('name','Amygdala','range',[13,28]);
wheretopick{10} =struct('name','Csf','range',14);
wheretopick{11} =struct('name','Accumbens','range',[15,29]);
wheretopick{12} =struct('name','ventral_IDC','range',[16,30]);
wheretopick{13} =struct('name','Optic chiasm', 'range',40);
wheretopick{14} =struct('name','Corpus callosum(CC)','range',[41,42,43,44,45]);

for k=2:14
for i=1:size(variant,1)
    if ratio==1
vol(i,1)=sum(allfiles{i,1}.Volume_mm3(wheretopick{k,1}.range))/volumes(i);
    else
    vol(i,1)=sum(allfiles{i,1}.Volume_mm3(wheretopick{k,1}.range));
    end
eh=vol(non_zero_indices);
ep=variant(non_zero_indices);
[R,P] = corr(ep(cutoffPA_indices),eh(cutoffPA_indices));
end
allvol{k,1}=[R,P];
end

wheretopick{1} = struct('name','total volume');
[R,P] = corr(volumes_non_zero,PA_non_zero);
allvol{1,1}=[R,P];

correlations = cellfun(@(x) abs(x(1)), allvol);
[sortedCorrelations, sortIdx] = sort(correlations, 'descend');

sprintf("Highest correlation area: %s, correlation value: %f, P value: %f", ...
    wheretopick{sortIdx(1), 1}.name, ...
    allvol{sortIdx(1), 1}(1), allvol{sortIdx(1), 1}(2))

sprintf("Second highest correlation area: %s, correlation value: %f, P value: %f", ...
    wheretopick{sortIdx(2), 1}.name, ...
    allvol{sortIdx(2), 1}(1), allvol{sortIdx(2), 1}(2))

sprintf("Third highest correlation area: %s, correlation value: %f, P value: %f", ...
    wheretopick{sortIdx(3), 1}.name, ...
    allvol{sortIdx(3), 1}(1), allvol{sortIdx(3), 1}(2))
