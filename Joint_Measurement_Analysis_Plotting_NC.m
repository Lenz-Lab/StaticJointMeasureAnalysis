clear, clc, close all

% Load in the .mat file
[file,path] = uigetfile('*.*','Select .mat file');
addpath(path);
load(file);

% Get all CPs
fns = fieldnames(Data);

list_names = {fns{1:13}, 'Population'};
list_distcong = {'Distance','Congruency'};
list_joint = {'NavicularCuneiform','MedialIntermediate','LateralIntermediate'};

[joint_indx,~] = listdlg('PromptString',{'Select which joint you selected'}, 'ListString', list_joint, 'SelectionMode','single');
[names_indx,~] = listdlg('PromptString',{'Select which name you want'}, 'ListString', list_names, 'SelectionMode','single');
[distcong_indx,~] = listdlg('PromptString',{'Select which calculation you want'}, 'ListString', list_distcong, 'SelectionMode','single');

if distcong_indx == 1
    range = [0 4];
elseif distcong_indx == 2
    range = [0 1];
end

if joint_indx == 3
    bone_name = 'Intermediate';
    views = [-115 -15];
elseif joint_indx == 2
    bone_name = 'Intermediate';
    views = [70 -15];
else
    bone_name = 'Navicular';
    views = [-15 -25];
end

%% Population Calculations
if names_indx == 14
    list_meanstd = {'Mean','+1_SD','-1_SD','+2_SD','-2_SD'};
    [meanstd_indx,~] = listdlg('PromptString',{'Select which you would like'}, 'ListString', list_meanstd, 'SelectionMode','single');

    for n = 1:length(fns)
        if n > 13
            if joint_indx == 1
                temp = [Data(1).(fns{n}).MeasureData(:,1); Data(2).(fns{n}).MeasureData(:,1); Data(3).(fns{n}).MeasureData(:,1)];
                all_cp(1:length(temp),n) = temp;
            else
                temp = Data.(fns{n}).MeasureData(:,1);
                all_cp(1:length(temp),n) = temp;
            end
        end
    end

    all_cp = nonzeros(all_cp);

    % sort_all_cp = sort(all_cp(:,1),'ascend');
    sort_all_cp = all_cp;

    C = nonzeros(hist(sort_all_cp,1:max(sort_all_cp)));

    all_cp_unique = [unique(sort_all_cp) C];

    k = 1;
    for n = 1:length(all_cp_unique)
        if C(n) <= 20
            all_cp_unique_clear(k,:) = [0 0];
        else
            all_cp_unique_clear(k,:) = all_cp_unique(n,:);
        end
        k = k+1;
    end

    new_cp = nonzeros(all_cp_unique_clear(:,1));

    if joint_indx == 1
        tempMeasureData = [Data(1).(fns{1}).MeasureData; Data(2).(fns{1}).MeasureData; Data(3).(fns{1}).MeasureData;...
            Data(1).(fns{2}).MeasureData; Data(2).(fns{2}).MeasureData; Data(3).(fns{2}).MeasureData;...
            Data(1).(fns{3}).MeasureData; Data(2).(fns{3}).MeasureData; Data(3).(fns{3}).MeasureData;...
            Data(1).(fns{4}).MeasureData; Data(2).(fns{4}).MeasureData; Data(3).(fns{4}).MeasureData;...
            Data(1).(fns{5}).MeasureData; Data(2).(fns{5}).MeasureData; Data(3).(fns{5}).MeasureData;...
            Data(1).(fns{6}).MeasureData; Data(2).(fns{6}).MeasureData; Data(3).(fns{6}).MeasureData;...
            Data(1).(fns{7}).MeasureData; Data(2).(fns{7}).MeasureData; Data(3).(fns{7}).MeasureData;...
            Data(1).(fns{8}).MeasureData; Data(2).(fns{8}).MeasureData; Data(3).(fns{8}).MeasureData;...
            Data(1).(fns{9}).MeasureData; Data(2).(fns{9}).MeasureData; Data(3).(fns{9}).MeasureData;...
            Data(1).(fns{10}).MeasureData; Data(2).(fns{10}).MeasureData; Data(3).(fns{10}).MeasureData;...
            Data(1).(fns{11}).MeasureData; Data(2).(fns{11}).MeasureData; Data(3).(fns{11}).MeasureData;...
            Data(1).(fns{12}).MeasureData; Data(2).(fns{12}).MeasureData; Data(3).(fns{12}).MeasureData;...
            Data(1).(fns{13}).MeasureData; Data(2).(fns{13}).MeasureData; Data(3).(fns{13}).MeasureData;...
            Data(1).(fns{27}).MeasureData; Data(2).(fns{27}).MeasureData; Data(3).(fns{27}).MeasureData;...
            Data(1).(fns{28}).MeasureData; Data(2).(fns{28}).MeasureData; Data(3).(fns{28}).MeasureData;...
            Data(1).(fns{29}).MeasureData; Data(2).(fns{29}).MeasureData; Data(3).(fns{29}).MeasureData;...
            Data(1).(fns{30}).MeasureData; Data(2).(fns{30}).MeasureData; Data(3).(fns{30}).MeasureData;...
            Data(1).(fns{31}).MeasureData; Data(2).(fns{31}).MeasureData; Data(3).(fns{31}).MeasureData;...
            Data(1).(fns{32}).MeasureData; Data(2).(fns{32}).MeasureData; Data(3).(fns{32}).MeasureData;...
            Data(1).(fns{33}).MeasureData; Data(2).(fns{33}).MeasureData; Data(3).(fns{33}).MeasureData;...
            Data(1).(fns{34}).MeasureData; Data(2).(fns{34}).MeasureData; Data(3).(fns{34}).MeasureData;...
            Data(1).(fns{35}).MeasureData; Data(2).(fns{35}).MeasureData; Data(3).(fns{35}).MeasureData;...
            Data(1).(fns{36}).MeasureData; Data(2).(fns{36}).MeasureData; Data(3).(fns{36}).MeasureData;...
            Data(1).(fns{37}).MeasureData; Data(2).(fns{37}).MeasureData; Data(3).(fns{37}).MeasureData;...
            Data(1).(fns{38}).MeasureData; Data(2).(fns{38}).MeasureData; Data(3).(fns{38}).MeasureData;...
            Data(1).(fns{39}).MeasureData; Data(2).(fns{39}).MeasureData; Data(3).(fns{39}).MeasureData;...
            Data(1).(fns{40}).MeasureData; Data(2).(fns{40}).MeasureData; Data(3).(fns{40}).MeasureData];

    end

    if joint_indx == 1
        tempMeasureData(tempMeasureData(:,1) == tempMeasureData(:,2), :) = [];
    end


    if joint_indx == 1
        for m = 1:length(new_cp)
            for n = 1:length(tempMeasureData)
                if new_cp(m) == tempMeasureData(n)
                    all_dist(n,:) = [tempMeasureData(n,3), new_cp(m)];
                    all_cong(n,:) = [tempMeasureData(n,4), new_cp(m)];
                end
            end
        end
    else
        for p = 1:length(fns)
            if p > 14 || p < 26
                for m = 1:length(new_cp)
                    for n = 1:length(Data.(fns{p}).MeasureData)
                        if new_cp(m) == Data.(fns{p}).MeasureData(n)
                            all_dist(m,(p)) = Data.(fns{p}).MeasureData(n,3);
                            all_cong(m,(p)) = Data.(fns{p}).MeasureData(n,4);
                        end
                    end
                end
            end
        end
    end

     if joint_indx == 1
         for k = 1:length(all_dist)
             temp = mean(all_dist(find(all_dist(k,2) == all_dist(:,2))));
             temp2 = std(all_dist(find(all_dist(k,2) == all_dist(:,2))));
             indx = find(all_dist(k,2) == new_cp(:,:));
             all_mean_dist(indx,:) = temp;
             all_std_dist(indx,:) = temp2;
         end
     clear temp temp2 indx
     for k = 1:length(all_cong)
         if all_cong(k,2) ~= 0
             temp = mean(all_cong(find(all_cong(k,2) == all_cong(:,2))));
             temp2 = std(all_cong(find(all_cong(k,2) == all_cong(:,2))));
             indx = find(all_cong(k,2) == new_cp(:,:));
             all_mean_cong(indx,:) = temp;
             all_std_cong(indx,:) = temp2;
         end
     end

                  all_mean_dist(all_mean_dist == 0) = NaN;
         all_std_dist(all_std_dist == 0) = NaN;
         all_mean_cong(all_mean_cong == 0) = NaN;
         all_std_dist(all_std_dist == 0) = NaN;

%          for n = 1:length(all_mean_dist)
%              all_mean_dist(n,1) = mean(all_mean_dist(n,:),'omitnan');
%              all_std_dist(n,1) = std(all_mean_dist(n,:),'omitnan');
%          end
% 
%          for n = 1:length(all_mean_cong)
%              all_mean_cong(n,1) = mean(all_mean_cong(n,:),'omitnan');
%              all_std_cong(n,1) = std(all_mean_cong(n,:),'omitnan');
%          end

     else

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

     for n = 1:length(fns)
         if n > 14
             ROI = Data(1).Mean.(string(bone_name)).CP(new_cp,:);
         elseif n < 26
            ROI = Data(1).Mean.(string(bone_name)).CP(new_cp,:).*[1,1,-1];
         end
     end
     NodalIndex = 0;
     %     if distcong_indx == 1
     if meanstd_indx == 1
         AllNodalData = [zeros(length(new_cp),1) zeros(length(new_cp),1) all_mean_dist all_mean_cong];
     elseif meanstd_indx == 2
         AllNodalData = [zeros(length(new_cp),1) zeros(length(new_cp),1) all_mean_dist+all_std_dist all_mean_cong+all_std_cong];
     elseif meanstd_indx == 3
         AllNodalData = [zeros(length(new_cp),1) zeros(length(new_cp),1) all_mean_dist-all_std_dist all_mean_cong-all_std_cong];
     elseif meanstd_indx == 4
         AllNodalData = [zeros(length(new_cp),1) zeros(length(new_cp),1) all_mean_dist+(2*all_std_dist) all_mean_cong+(2*all_std_cong)];
     elseif meanstd_indx == 5
         AllNodalData = [zeros(length(new_cp),1) zeros(length(new_cp),1) all_mean_dist-(2*all_std_dist) all_mean_cong-(2*all_std_cong)];
     end

     for n = 1:length(AllNodalData)
         if AllNodalData(n,3) < 0
             AllNodalData(n,3) = 0;
         elseif AllNodalData(n,4) < 0
             AllNodalData(n,4) = 0;
         end
     end

     figure()
     plot3(ROI(:,1),ROI(:,2),ROI(:,3),'.k')
     axis equal
% 
%      tesst = find(Data(1).Mean.(string(bone_name)).CP(:,3) == 20.5708);
%      temps = find(tesst == tempMeasureData(:,1));
% 
%      tesst1 = find(Data(1).Mean.(string(bone_name)).CP(:,3) == 8.272);
%      temps1 = find(tesst1 == tempMeasureData(:,1));

     Bone.Points = Data(1).Mean.(string(bone_name)).(string(bone_name)).Points;
     Bone.ConnectivityList = Data(1).Mean.(string(bone_name)).(string(bone_name)).ConnectivityList;

else
    Bone.Points = Data(1).(fns{names_indx+13}).(string(bone_name)).(string(bone_name)).Points;
    Bone.ConnectivityList = Data(1).(fns{names_indx+13}).(string(bone_name)).(string(bone_name)).ConnectivityList;
%             ROI = Data.(fns{names_indx}).(string(bone)).CP(Data.(fns{names_indx}).MeasureData(:,1),:);
    ROI = [];
    %         SPMIndex = [];
    %         MeanCP = [];
    %         perc_stance = [];
    %         part_scatter = 2;
    if distcong_indx == 1 && joint_indx >= 2
        data_values = Data.(fns{names_indx+13}).MeasureData(:,3);
    elseif distcong_indx == 2 && joint_indx >= 2
        data_values = Data.(fns{names_indx+13}).MeasureData(:,4);
    elseif distcong_indx == 1 && joint_indx == 1
        data_values = [Data(1).(fns{names_indx+13}).MeasureData(:,3); Data(2).(fns{names_indx+13}).MeasureData(:,3); Data(3).(fns{names_indx+13}).MeasureData(:,3)];
    end

    if joint_indx == 1
        NodalIndex = [Data(1).(fns{names_indx+13}).MeasureData(:,2); Data(2).(fns{names_indx+13}).MeasureData(:,2); Data(3).(fns{names_indx+13}).MeasureData(:,2)];
        AllNodalData = [Data(1).(fns{names_indx+13}).MeasureData; Data(2).(fns{names_indx+13}).MeasureData; Data(3).(fns{names_indx+13}).MeasureData];
    else
        NodalIndex = Data.(fns{names_indx+13}).MeasureData(:,2);
        AllNodalData = Data.(fns{names_indx+13}).MeasureData;
    end
end

%     RainbowFish(Bone,MeanCP,NodalIndex,AllNodalData,range,distcong_indx,SPMIndex,perc_stance,part_scatter)
RainbowFish(Bone,NodalIndex,AllNodalData,range,distcong_indx,ROI,views)
