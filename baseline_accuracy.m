T_1 = readtable('\\147.47.66.155\public\Interns\전시현\WWW_survey.xlsx', 'Sheet', 'data_rejection');

%% sbj num/What/Where/When
common_folder = '\\147.47.66.155\public\Interns\전시현\www_behavior';
folders = dir(fullfile(common_folder, '*_*_*'));
sbj = 1;

for i = 1:length(folders)
    current_folder = fullfile(common_folder, folders(i).name);
    
    [~, folder_name, ~] = fileparts(current_folder);
    folder_name_parts = strsplit(folder_name, '_');
    if length(folder_name_parts) == 3
        folder_number_str = folder_name_parts(3);

    elseif length(folder_name_parts) == 4
        folder_number_str = folder_name_parts(4);
    end

    folder_number = str2double(folder_number_str);

    if cell2mat(T_1{T_1{:,1} == folder_number,3}) == 'O'
        files = dir(fullfile(current_folder, '*_results.mat'));
        acc(sbj,1) = folder_number;
    
        if ~isempty(files)
            current_file = fullfile(current_folder, files.name);
        
            data = load(current_file);
            acc(sbj,2) = mean(data.results.response.OX{1});
            acc(sbj,3) = mean(data.results.response.OX{2});
            acc(sbj,4) = mean(data.results.response.OX{3});
        end
        sbj = sbj + 1;
    else
        continue
    end
end

nonzero_idx = find(acc(:, 2) ~= 0);

%% sort by data_label
unique_num = unique(T_1(:,5));
label_number = cell(height(unique_num),1);

for i = 1:height(unique_num)
    label_number{i} = T_1(T_1{:,5} == i,1);
end

%% sort acc by data label
grouped_acc = cell(length(label_number), 1);

for i = 1:length(label_number)
    current_participants = label_number{i};
    [~, idx_acc, idx_participants] = intersect(array2table(acc(:,1)), current_participants);

    current_data = table();
    current_data.participant = acc(intersect(idx_acc, nonzero_idx),1);
    current_data.accuracy1 = acc(intersect(idx_acc, nonzero_idx),2);
    current_data.accuracy2 = acc(intersect(idx_acc, nonzero_idx),3);
    current_data.accuracy3 = acc(intersect(idx_acc, nonzero_idx),4);

    grouped_acc{i} = current_data;
end

%% calculate mean accuracy
mean_acc = cell(length(grouped_acc), 3);

for i = 1:length(grouped_acc)
    acc1 = grouped_acc{i}.accuracy1;
    acc2 = grouped_acc{i}.accuracy2;
    acc3 = grouped_acc{i}.accuracy3;

    mean_acc{i,1} = mean(acc1);
    mean_acc{i,2} = mean(acc2);
    mean_acc{i,3} = mean(acc3);
end

relative_dementia = cell(3,1);
relative_label = table();
relative_label.participant = grouped_acc{2}.participant;
relative_label.accuracy1 = grouped_acc{2}.accuracy1/mean_acc{1,1};
relative_label.accuracy2 = grouped_acc{2}.accuracy2/mean_acc{1,2};
relative_label.accuracy3 = grouped_acc{2}.accuracy3/mean_acc{1,3};

relative_dementia{1} = relative_label;

relative_label = table();
relative_label.participant = grouped_acc{4}.participant;
relative_label.accuracy1 = grouped_acc{4}.accuracy1/mean_acc{3,1};
relative_label.accuracy2 = grouped_acc{4}.accuracy2/mean_acc{3,2};
relative_label.accuracy3 = grouped_acc{4}.accuracy3/mean_acc{3,3};

relative_dementia{2} = relative_label;

relative_label = table();
relative_label.participant = grouped_acc{5}.participant;
relative_label.accuracy1 = grouped_acc{5}.accuracy1/mean_acc{3,1};
relative_label.accuracy2 = grouped_acc{5}.accuracy2/mean_acc{3,2};
relative_label.accuracy3 = grouped_acc{5}.accuracy3/mean_acc{3,3};

relative_dementia{3} = relative_label;
