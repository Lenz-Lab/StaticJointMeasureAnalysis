clear, clc, close all

% Load in the .mat file
[file,path] = uigetfile('*.*','Select .mat file');
addpath(path);
load(file);
%% Get all CPs
fns = fieldnames(Data);

list_names = {fns{1:length(fns)}, 'Population'};
list_distcong = {'Distance','Congruency','Coverage'};
list_joint = {'Subtalar','Talonavicular','Calcaneocuboid'};

[joint_indx,~] = listdlg('PromptString',{'Select which joint you selected'}, 'ListString', list_joint, 'SelectionMode','single');
[names_indx,~] = listdlg('PromptString',{'Select which name you want'}, 'ListString', list_names, 'SelectionMode','single');
[distcong_indx,~] = listdlg('PromptString',{'Select which calculation you want'}, 'ListString', list_distcong, 'SelectionMode','single');

% if names_indx == length(fns)+1 && distcong_indx == 3
%         Bone.Points = Data.left_05_5496366.Cuboid.Cuboid.Points;
%     Bone.ConnectivityList = Data.left_05_5496366.Cuboid.Cuboid.ConnectivityList;
%     Coverage_Surf = Data.left_05_5496366.CoverageAreaSurf.Cuboid  ;
%     if joint_indx == 1
%         Coverage_Split = Data.(fns{names_indx}).CoverageAreaSurfSplit.(string(bone));
%     else
%         Coverage_Split = [];
%     end
% 
%     RainbowFishCoverage(Bone,Coverage_Surf)%,Coverage_Split)
% end

if distcong_indx == 1
    range = [0 4];
elseif distcong_indx == 2
    range = [0 1];
end

if joint_indx == 1
    bone = 'Talus';
    views = [-185, -70];
elseif joint_indx == 2
    bone = 'Navicular';
    views = [-5, -20];
elseif joint_indx == 3
    bone = 'Cuboid';
    views = [0, 5];
end

if names_indx == length(fns)+1
    list_meanstd = {'Mean','+1_SD','-1_SD'};
    [meanstd_indx,~] = listdlg('PromptString',{'Select which you would like'}, 'ListString', list_meanstd, 'SelectionMode','single');
%     title1 = strcat(string(list_meanstd(meanstd_indx)),'_Population_',string(list_distcong(distcong_indx)),'_',string(list_joint(joint_indx)));
else
%     title1 = strcat(string(list_names(names_indx)),'_',string(list_distcong(distcong_indx)),'_',string(list_joint(joint_indx)));
end

if distcong_indx == 1 || distcong_indx == 2
    if names_indx == length(fns)+1
        %% Population Calculations
        for n = 1:length(fns)
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

        for p = 1:length(fns)
            for m = 1:length(new_cp)
                for n = 1:length(Data.(fns{p}).MeasureData)
                    if new_cp(m) == Data.(fns{p}).MeasureData(n)
                        all_dist(m,(p)) = [Data.(fns{p}).MeasureData(n,3)];
                        all_cong(m,(p)) = [Data.(fns{p}).MeasureData(n,4)];
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

        ROI = Data.Mean.(string(bone)).CP(new_cp,:);
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
        Bone.Points = Data.Mean.(string(bone)).(string(bone)).Points;
        Bone.ConnectivityList = Data.Mean.(string(bone)).(string(bone)).ConnectivityList;

        %     k = 1;
        %     tol = 6;
        %     tri_points = ROI;
        %
        %     A_A = temp_STL.(string(bone_names(bone_no_CP))).Points;
        %     C_C = temp_STL.(string(bone_names(bone_with_CP))).Points;
        %
        %     % Pair nodes with CP and calculate euclidean distance
        % %     clear temp i_surf
        %     for h = 1:length(tri_points(:,1))
        %         for j = 1:length(Data.Mean.(string(bone)).CP_Bone)
        %         % Kept the line below for legacy
        %         %     ROI = find(temp_STL.(string(bone_names(bone_no_CP))).Points(:,1) >= temp_STL.(string(bone_names(bone_with_CP))).Points(tri_points(h),1)-tol & temp_STL.(string(bone_names(bone_no_CP))).Points(:,1) <= temp_STL.(string(bone_names(bone_with_CP))).Points(tri_points(h),1)+tol & temp_STL.(string(bone_names(bone_no_CP))).Points(:,2) >= temp_STL.(string(bone_names(bone_with_CP))).Points(tri_points(h),2)-tol & temp_STL.(string(bone_names(bone_no_CP))).Points(:,2) <= temp_STL.(string(bone_names(bone_with_CP))).Points(tri_points(h),2)+tol & temp_STL.(string(bone_names(bone_no_CP))).Points(:,3) >= temp_STL.(string(bone_names(bone_with_CP))).Points(tri_points(h),3)-tol & temp_STL.(string(bone_names(bone_no_CP))).Points(:,3) <= temp_STL.(string(bone_names(bone_with_CP))).Points(tri_points(h),3)+tol);
        % %         ROI2 = find(A_A(:,1) >= C_C(tri_points(h),1)-tol & A_A(:,1) <= C_C(tri_points(h),1)+tol & A_A(:,2) >= C_C(tri_points(h),2)-tol & A_A(:,2) <= C_C(tri_points(h),2)+tol & A_A(:,3) >= C_C(tri_points(h),3)-tol & A_A(:,3) <= C_C(tri_points(h),3)+tol);
        % %         if isempty(ROI2) == 0
        %             i_CP = Data.Mean.(string(bone)).CP_Bone(find(new_cp(h,1) == Data.Mean.(string(bone)).CP_Bone(j,1)),1);
        %
        %             if isempty(i_CP) == 0
        %             i_surff(k,:) = [i_CP(1) tri_points(h,1)];
        %             end
        %             %     tempp(1) == the index of the paired node on the opposing bone surface
        %             k = k + 1;
        % %         end
        % %         clear temp tempp
        %         end
        %     end
    else
        Bone.Points = Data.(fns{names_indx}).(string(bone)).(string(bone)).Points;
        Bone.ConnectivityList = Data.(fns{names_indx}).(string(bone)).(string(bone)).ConnectivityList;
        %         ROI = Data.(fns{names_indx}).(string(bone)).CP(Data.(fns{names_indx}).MeasureData(:,1),:);
        ROI = [];
        if distcong_indx == 1
            data_values = Data.(fns{names_indx}).MeasureData(:,3);
        elseif distcong_indx == 2
            data_values = Data.(fns{names_indx}).MeasureData(:,4);
        end
        NodalIndex = Data.(fns{names_indx}).MeasureData(:,2);
        AllNodalData = Data.(fns{names_indx}).MeasureData;
    end

    RainbowFish(Bone,NodalIndex,AllNodalData,range,distcong_indx,ROI)
else
    Bone.Points = Data.(fns{names_indx}).(string(bone)).(string(bone)).Points;
    Bone.ConnectivityList = Data.(fns{names_indx}).(string(bone)).(string(bone)).ConnectivityList;
    Coverage_Surf = Data.(fns{names_indx}).CoverageAreaSurf.(string(bone));
    if joint_indx == 1
        Coverage_Split = Data.(fns{names_indx}).CoverageAreaSurfSplit.(string(bone));
    else
        Coverage_Split = [];
    end

    RainbowFishCoverage(Bone,Coverage_Surf)%,Coverage_Split)
end

%     hold on
%     figure()
%     patch('Faces',Coverage_Surf.ConnectivityList,'Vertices',Coverage_Surf.Points,...
%     'FaceColor', [0.85 0.85 0.85], ...
%     'EdgeColor','none',...
%     'FaceLighting','gouraud',...
%     'AmbientStrength', 0.15);
% camlight(views(1),views(2))
% view(views(1),views(2))
% material('dull');
% title(title1)
% grid off
% axis off
% axis equal
% hold on

% export_fig 'test.eps'

% figure()
% RainbowFish(ROI,data_values,range,100,1)
% hold on
% patch('Faces',TR.ConnectivityList,'Vertices',TR.Points,...
%     'FaceColor', [0.85 0.85 0.85], ...
%     'EdgeColor','none',...
%     'FaceLighting','gouraud',...
%     'AmbientStrength', 0.15);
% camlight(views(1),views(2))
% view(views(1),views(2))
% material('dull');
% title(title1)
% grid off
% axis off
% axis equal
% hold on

% figure()
% ColorMapPlot3_TIF(ROI,data_values,[0 1],100,1)
% hold on
% patch('Faces',TR.ConnectivityList,'Vertices',TR.Points,...
%     'FaceColor', [0.85 0.85 0.85], ...
%     'EdgeColor','none',...
%     'FaceLighting','gouraud',...
%     'AmbientStrength', 0.15);
% camlight(-185,-70)
% view(-185,-70)
% material('dull');
% title('Mean Population Congruency')
% grid off
% axis off
% axis equal

%     figure()
%     ColorMapPlot3_TIF(ROI,surf(:,3),[0 3.5],1000,2)
%     hold on
%     patch('Faces',TR.ConnectivityList,'Vertices',TR.Points,...
%         'FaceColor', [0.85 0.85 0.85], ...
%         'EdgeColor','none',...
%         'FaceLighting','gouraud',...
%         'AmbientStrength', 0.15);
%     % patch('Faces',TR.ConnectivityList,'Vertices',TR.Points)
%     camlight(-185,-70);
%     material('dull');
%     % hold on
%     axis equal
