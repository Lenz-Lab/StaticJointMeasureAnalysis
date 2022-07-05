clear, clc%, close all

% Load in the .mat file
[file,path] = uigetfile('*.*','Select .mat file');
addpath(path);
load(file);

% Get all CPs
fns = fieldnames(Data);

% list_names = {fns{1:19}, 'Population'};
list_names = {fns{1:13}};
list_distcong = {'Distance','Congruency'}; %,'Coverage'};
list_joint = {'Subtalar','Talonavicular','Calcaneocuboid','NavicularMedial','NavicularIntermediate','NavicularLateral'};

[joint_indx,~] = listdlg('PromptString',{'Select which joint you selected'}, 'ListString', list_joint, 'SelectionMode','single');
[names_indx,~] = listdlg('PromptString',{'Select which name you want'}, 'ListString', list_names, 'SelectionMode','single');
[distcong_indx,~] = listdlg('PromptString',{'Select which calculation you want'}, 'ListString', list_distcong, 'SelectionMode','single');

if joint_indx > 3
    bone_name = 'Navicular';
end

if names_indx == 20
    list_meanstd = {'Mean','+1_SD','-1_SD'};
    [meanstd_indx,~] = listdlg('PromptString',{'Select which you would like'}, 'ListString', list_meanstd, 'SelectionMode','single');
end

if names_indx == 20
    %% Population Calculations
    for n = 20:length(fns)
        temp = Data.(fns{n}).MeasureData(:,1);

        all_cp(1:length(temp),n) = temp;
    end

    all_cp = nonzeros(all_cp);

    % sort_all_cp = sort(all_cp(:,1),'ascend');
    sort_all_cp = all_cp;

    C = nonzeros(hist(sort_all_cp,1:max(sort_all_cp)));

    all_cp_unique = [unique(sort_all_cp) C];

    for n = 1:length(all_cp_unique)
        if C(n) <= 2
            all_cp_unique_clear(n,:) = [0 0];
        else
            all_cp_unique_clear(n,:) = all_cp_unique(n,:);
        end
    end

    new_cp = nonzeros(all_cp_unique_clear(:,1));

    for p = 20:length(fns)
        for m = 1:length(new_cp)
            for n = 1:length(Data.(fns{p}).MeasureData)
                if new_cp(m) == Data.(fns{p}).MeasureData(n)
                    all_dist(m,(p-19)) = [Data.(fns{p}).MeasureData(n,3)];
                    all_cong(m,(p-19)) = [Data.(fns{p}).MeasureData(n,4)];
                end
            end
        end
    end

    all_dist(all_dist == 0) = NaN;
    all_cong(all_cong == 0) = NaN;

    for n = 1:length(all_dist)
        all_mean_dist(n,1) = mean(all_dist(n,:),'omitnan');
        all_std_dist(n,1) = std(all_dist(n,:),'omitnan');
    end

    for n = 1:length(all_cong)
        all_mean_cong(n,1) = mean(all_cong(n,:),'omitnan');
        all_std_cong(n,1) = std(all_cong(n,:),'omitnan');
    end

    %     for n = 20:length(fns)
    %         for m = 1:length(Data.(fns{n}).MeasureData(:,1))
    %             for p = 1:length(new_cp)
    %                 if Data.(fns{n}).MeasureData(m,1) == new_cp(p,:)
    %                     temp(p,(n-19)) = Data.(fns{n}).MeasureData(m,2);
    %                 end
    %             end
    %         end
    %     end
    %
    %     k = 1;
    %     for n = 1:length(temp(:,1))
    %         average_across(k,:) = mean(nonzeros(temp(n,:)));
    %         k = k+1;
    %     end
    %
    %     average_across = round(average_across,0);
    %     NodalIndex = average_across;


    %     total_mean_dist = mean(all_mean_dist);
    %     total_mean_cong = mean(all_mean_cong);

    ROI = Data.Mean.(string(bone_name)).CP(new_cp,:);
    NodalIndex = 0;
    %     if distcong_indx == 1
    if meanstd_indx == 1
        AllNodalData = [zeros(length(new_cp),1) zeros(length(new_cp),1) all_mean_dist all_mean_cong];
    elseif meanstd_indx == 2
        AllNodalData = [zeros(length(new_cp),1) zeros(length(new_cp),1) all_mean_dist+all_std_dist all_mean_cong+all_std_cong];
    elseif meanstd_indx == 3
        AllNodalData = [zeros(length(new_cp),1) zeros(length(new_cp),1) all_mean_dist-all_std_dist all_mean_cong-all_std_cong];
    end

    %     elseif distcong_indx == 2
    %         if meanstd_indx == 1
    %             data_values = all_mean_cong;
    %         elseif meanstd_indx == 2
    %             data_values = all_mean_cong + all_std_cong;
    %         elseif meanstd_indx == 3
    %             data_values = all_mean_cong - all_std_cong;
    %         end
    %     end
    Bone.Points = Data.Mean.(string(bone_name)).(string(bone_name)).Points;
    Bone.ConnectivityList = Data.Mean.(string(bone_name)).(string(bone_name)).ConnectivityList;

else
    Bone.Points = Data.(fns{names_indx}).(string(bone_name)).(string(bone_name)).Points;
    Bone.ConnectivityList = Data.(fns{names_indx}).(string(bone_name)).(string(bone_name)).ConnectivityList;
    ROI = [];
    if distcong_indx == 1
        data_values = Data.(fns{names_indx}).MeasureData(:,3);
        range = [0 6];
    elseif distcong_indx == 2
        data_values = Data.(fns{names_indx}).MeasureData(:,4);
        range = [0 1];
    end
    NodalIndex = Data.(fns{names_indx}).MeasureData(:,2);
    AllNodalData = Data.(fns{names_indx}).MeasureData;
end

RainbowFish(Bone,NodalIndex,AllNodalData,range,distcong_indx,ROI)
