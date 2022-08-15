%% Static Joint Measurement Analysis
% Calculates coverage area, joint space distance and congruence index between two
% different bones at correspondence particles on a particular bone surface
% throughout a dynamic activity.

% Modified By: Andrew Peterson
% Date: 7/1/2022

%% Required Files and Input
% This script requires a folder structure and files in order to process the
% data appropriately.

joint_names = {'Subtalar'};
joint = 1;
bone_names = {'Talus','Calcaneus'};
% Please update bone_names and joint_names variable with the names of the bones and joint of
% interest. Spelling is very important and must be consistent in all file
% names!

% Folder Architecture:
% Main Directory -> (folder containing each of the group folders)
%     Folders:
%     Group_A
%     Group_(...)
%     Group_(n-1) -> (contains each of the subject folders within that group)
%         Folders:
%         Subject_01
%         Subject_(...)
%         (Name)_(m-1) -> (contains files for that subject)
%             Files:
%             (Name).local.particles        (exported from ShapeWorks)
%             (Bone_Name_01)_groomed.vtk    (exported from ShapeWorks)
%             (Bone_Name_02)_groomed.vtk    (exported from ShapeWorks)
%             (Bone_Name_01).stl            (input bone model into ShapeWorks)
%             (Bone_Name_02).stl            ('opposing' bone model)

% Please read the standard operating procedure (SOP) included in the
% .github repository.
%% Clean Slate
clc, close all, clearvars -except bone_names joint_names joint
delete(gcp('nocreate'))
pool = parpool([1 100]);
clc

%% Loading Data
fprintf('Loading Data:\n')

fldr_name = uigetdir;
addpath(fldr_name)

D = dir(fullfile(sprintf('%s\\',fldr_name)));

pulled_files = [];
m = 1;
for k = 3:length(D)
    pulled_files{m} = D(k).name;
    m = m + 1;
end

temp = strsplit(fldr_name,'\');
subj_groups.(string(temp(end))).SubjectList = pulled_files;

%% Load Data for Each Subject
for m = 1:length(pulled_files)
    %%
    fprintf('   %s\n',string(pulled_files(m)))
    addpath(sprintf('%s\\%s\\',string(fldr_name),string(pulled_files(m))))

    %% Load the Bone.stl Files
    S = [];
    S = dir(fullfile(sprintf('%s\\%s\\',string(fldr_name),string(pulled_files(m))),'*.stl'));
    for b = 1:length(bone_names)
        for c = 1:length(S)
            temp_under = split(S(c).name,'_');
            temp_dot = split(temp_under(end),'.');
            temp = [temp_under(1:end-1); temp_dot];
            for d = 1:length(temp)
                % temp_check = strfind(string(bone_names(b)),string(temp(d)));
                temp_check = any(string(bone_names(b)) == string(temp(d)));
                if  temp_check == 1
                    temp_stl = [];
                    temp_stl = stlread(S(c).name);
                    Data.(string(pulled_files(m))).(string(bone_names(b))).(string(bone_names(b))) = stlread(S(c).name);

                    temp_bone = [];
                    temp_bone = Data.(string(pulled_files(m))).(string(bone_names(b))).(string(bone_names(b)));

                    % Calculate Gaussian and Mean Curvatures
                    % Meyer, M., Desbrun, M., Schröder, P., & Barr, A. H. (2003). Discrete differential-geometry operators for triangulated 2-manifolds. In Visualization and mathematics III (pp. 35-57). Springer Berlin Heidelberg.
                    [Data.(string(pulled_files(m))).(string(bone_names(b))).GaussianCurve, Data.(string(pulled_files(m))).(string(bone_names(b))).MeanCurve] = curvatures(temp_bone.Points(:,1),temp_bone.Points(:,2),temp_bone.Points(:,3),temp_bone.ConnectivityList);
                end
                % Set which side the bones are. This is important for
                % pairing the .stl points with the CP points in later
                % steps. The .ply files from ShapeWorks have a
                % different mesh than the input .stl files.
                side_check = strsplit(string(temp(d)),'.');
                if isempty(any('Right' == side_check(1))) == 0 || isempty(any('right' == side_check(1))) == 0 || isempty(any('R' == side_check(1))) == 0
                    Data.(string(pulled_files(m))).Side = 'Right';
                end
                if isempty(any('Left' == side_check(1))) == 0  || isempty(any('left' == side_check(1))) == 0 || isempty(any('L' == side_check(1))) == 0
                    Data.(string(pulled_files(m))).Side = 'Left';
                end
            end
        end
    end

    %% Load the Individual Bone Kinematics from .txt
    for b = 1:length(bone_names)
        % Assumes there is no kinematics and it is one static frame
        Data.(string(pulled_files(m))).(string(bone_names(b))).Kinematics = [1 0 0 0, 0 1 0 0, 0 0 1 0, 0 0 0 1]; % Identity Matrix
    end

    %% Load the ShapeWorks Output Groomed .vtk File
    % This file is used to align (using iterative closest point
    % algorithm) the bone with correspondence particles too in order to
    % identify CP and node pairing.
    P = [];
    P = dir(fullfile(sprintf('%s\\%s\\',string(fldr_name),string(pulled_files(m))),'*.vtk'));
    for b = 1:length(bone_names)
        for c = 1:length(P)
            temp_under = split(P(c).name,'_');
            temp_dot = split(temp_under(end),'.');
            temp = [temp_under(1:end-1); temp_dot];
            for d = 1:length(temp)
                % temp_check = strfind(string(bone_names(b)),string(temp(d)));
                temp_check = any(string(bone_names(b)) == string(temp(d)));
                if  temp_check == 1
                    temp_VTK = [];
                    temp_VTK = LoadDataFile(P(c).name);
                    Data.(string(pulled_files(m))).(string(bone_names(b))).GroomedVTK.Location   = temp_VTK;
                end
            end
        end
    end

    %% Load the Correspondence Particles (CP) from ShapeWorks
    C = [];
    C = dir(fullfile(sprintf('%s\\%s\\',string(fldr_name),string(pulled_files(m))),'*.particles'));
    for b = 1:length(bone_names)
        for c = 1:length(C)
            temp_under = split(C(c).name,'_');
            temp_dot = split(temp_under(end),'.');
            temp = [temp_under(1:end-1); temp_dot];
            for d = 1:length(temp)
                % temp_check = strfind(string(bone_names(b)),string(temp(d)));
                temp_check = any(string(bone_names(b)) == string(temp(d)));
                if  temp_check == 1
                    temp_cp = [];
                    temp_cp = importdata(C(c).name);
                    Data.(string(pulled_files(m))).(string(bone_names(b))).CP = temp_cp;
                end
            end
        end
    end
end

subjects = [];
subjects = [subjects, pulled_files];
clear pulled_files

%% Identify Indices on Bones from SSM Local Particles
g = fieldnames(Data);
for subj_count = 1:length(g)
    for bone_count = 1:length(bone_names)
        if isfield(Data.(string(subjects(subj_count))).(string(bone_names(bone_count))),'CP') == 1
            i_pair = [];
            tol = 10;
            CP = [];
            CP = Data.(string(subjects(subj_count))).(string(bone_names(bone_count))).CP;
            q = [];
            q = Data.(string(subjects(subj_count))).(string(bone_names(bone_count))).GroomedVTK.Location.Points';
            p = [];
            p = Data.(string(subjects(subj_count))).(string(bone_names(bone_count))).(string(bone_names(bone_count))).Points;

            % Will need to flip the bone if it is a right in order to align
            % properly.
            if isequal(Data.(string(subjects(subj_count))).Side,'Left')
                p = [p(:,1) p(:,2) p(:,3)]';
            end
            if isequal(Data.(string(subjects(subj_count))).Side,'Right')
                p = [-p(:,1) p(:,2) p(:,3)]';
            end

            % Shifts the .stl bone to the center of the .ply to minimize
            % errors and iterations
            p(1,:) = p(1,:) - abs(mean(q(1,:)));
            p(2,:) = p(2,:) - abs(mean(q(2,:)));
            p(3,:) = p(3,:) - abs(mean(q(3,:)));

            R = [];
            T = [];
            % Calculate the rotations and translation matrices
            [R,T] = icp(q,p,200,'Matching','kDtree');

            P = [];
            P = (R*p + repmat(T,1,length(p)))';
            q = q';

            % Troubleshooting check for proper alignment
                        figure()
                        plot3(CP(:,1),CP(:,2),CP(:,3),'.k')
                         hold on
                        plot3(p(1,:),p(2,:),p(3,:),'ob')
                        hold on
                        plot3(P(:,1),P(:,2),P(:,3),'*r')
                        hold on
                         plot3(q(:,1),q(:,2),q(:,3),'g.')
                        axis equal
                        xlabel('x')
                        ylabel('y')
                        zlabel('z')

            %% Identify Nodes and CP
            % Find the .stl nodes and their respective correspondence
            % particles and save to Data structure
            for r = 1:length(CP(:,1))
                ROI = find(P(:,1) >= CP(r,1)-tol & P(:,1) <= CP(r,1)+tol & P(:,2) >= CP(r,2)-tol & P(:,2) <= CP(r,2)+tol & P(:,3) >= CP(r,3)-tol & P(:,3) <= CP(r,3)+tol);
                found_dist = pdist2(single(CP(r,:)),single(P(ROI,:)));
                min_dist = (ROI(find(found_dist == min(found_dist))));
                i_pair(r,:) = [r min_dist];
                clear found_dist min_dist ROI
            end
            Data.(string(subjects(subj_count))).(string(bone_names(bone_count))).CP_Bone = i_pair;

            % Troubleshooting check for proper pairing
            %             c = 1000;
            %             figure()
            %             plot3(CP(i_pair(:,1),1),CP(i_pair(:,1),2),CP(i_pair(:,1),3),'.k')
            %             hold on
            %             plot3(CP(i_pair(c,1),1),CP(i_pair(c,1),2),CP(i_pair(c,1),3),'*r')
            %             hold on
            %             axis equal
        end
    end
end

%% Bone Transformations
for subj_count = 1:length(g)
    frame_count = 1;
    clear temp i_pair temp_STL
    for bone_count = 1:length(bone_names)
        if isfield(Data.(string(subjects(subj_count))).(string(bone_names(bone_count))),'CP_Bone') == 1
            bone_with_CP = bone_count;
            i_pair = Data.(string(subjects(subj_count))).(string(bone_names(bone_count))).CP_Bone;
            bone_data = Data.(string(subjects(subj_count))).(string(bone_names(bone_count))).(string(bone_names(bone_count))).Points;
        elseif isfield(Data.(string(subjects(subj_count))).(string(bone_names(bone_count))),'CP_Bone') == 0
            bone_no_CP = bone_count;
            bone_data = Data.(string(subjects(subj_count))).(string(bone_names(bone_count))).(string(bone_names(bone_count))).Points;
        end
        kine_data = Data.(string(subjects(subj_count))).(string(bone_names(bone_count))).Kinematics;
        R = [kine_data(frame_count,1:3);kine_data(frame_count,5:7);kine_data(frame_count,9:11)];
        temp{bone_count} = (R*bone_data')';
        temp{bone_count} = [temp{bone_count}(:,1)+kine_data(frame_count,4), temp{bone_count}(:,2)+kine_data(frame_count,8), temp{bone_count}(:,3)+kine_data(frame_count,12)];

        temp_STL.(string(bone_names(bone_count))) = triangulation(Data.(string(subjects(subj_count))).(string(bone_names(bone_count))).(string(bone_names(bone_count))).ConnectivityList,temp{bone_count});
        clear R bone_data kine_data
    end

    %% Find if intersecting and calculate surface area
    for bone_count = 1:length(bone_names)
        if isfield(Data.(string(subjects(subj_count))).(string(bone_names(bone_count))),'CP_Bone') == 1
            temp_center = incenter(temp_STL.(string(bone_names(bone_count))));
            temp_normal = faceNormal(temp_STL.(string(bone_names(bone_count))));

            temp_n = [];
            temp_d = [];

            % https://www.mathworks.com/help/matlab/ref/triangulation.facenormal.html
            % TriangleRayIntersection (orig, dir, vert0, vert1, vert2, varargin)

            parfor (norm_check = 1:length(temp_center),pool)
                [temp_int] = TriangleRayIntersection(temp_center(norm_check,:),temp_normal(norm_check,:),temp_STL.(string(bone_names(bone_no_CP))).Points(temp_STL.(string(bone_names(bone_no_CP))).ConnectivityList(:,1),:),temp_STL.(string(bone_names(bone_no_CP))).Points(temp_STL.(string(bone_names(bone_no_CP))).ConnectivityList(:,2),:),temp_STL.(string(bone_names(bone_no_CP))).Points(temp_STL.(string(bone_names(bone_no_CP))).ConnectivityList(:,3),:));
                if isempty(find(temp_int == 1)) == 0
                    temp_n(norm_check,:) = 1;
                end
            end

            temp_n = find(temp_n == 1);
            pool.IdleTimeout = 30; % resets pool timeout to 30 minutes

            %% Find the indices of the points and faces
            % tri_found -> faces found that intersect opposing surface
            tri_found = [];
            for tri_check = 1:length(temp_n)
                t = [];
                t = temp_STL.(string(bone_names(bone_with_CP))).ConnectivityList(temp_n(tri_check),:);
                for tri_fill = 1:length(t)
                    tri_found(end+1,:) = t(tri_fill);
                end
            end

            tri_found = unique(tri_found);

            % tri_point -> index of 'identified nodes'
            % tri_cp    -> index of correspondex particle
            tri_cp = [];
            tri_points = [];
            for n = 1:length(tri_found)
                temp = find(tri_found(n) == Data.(string(subjects(subj_count))).(string(bone_names(bone_count))).CP_Bone(:,2));
                if isempty(temp) == 0
                    tri_points(end+1,:)   = Data.(string(subjects(subj_count))).(string(bone_names(bone_count))).CP_Bone(temp(1),2);
                    tri_cp(end+1,:)       = temp(1); %Data.(string(subjects(subj_count))).(string(bone_names(bone_count))).CP_Bone(temp,1); Basically same thing since it is the index...
                end
            end

            %% Calculate Coverage Surface Area
            for n = 1:length(temp_n)
                temp_tri = temp_STL.(string(bone_names(bone_with_CP))).ConnectivityList(temp_n(n),:);
                P1 = temp_STL.(string(bone_names(bone_with_CP))).Points(temp_tri(:,1),:);
                P2 = temp_STL.(string(bone_names(bone_with_CP))).Points(temp_tri(:,2),:);
                P3 = temp_STL.(string(bone_names(bone_with_CP))).Points(temp_tri(:,3),:);
                a = P2 - P1;
                b = P3 - P1;
                c = cross(a,b,2);
                area_tri(n,:) = 1/2*sum(sqrt(sum(c.^2,2)));
                clear temp_tri
            end

            Data.(string(subjects(subj_count))).CoverageArea.(string(bone_names(bone_with_CP))) = sum(area_tri);

            % Troubleshooting
            % TR.vertices =    temp_STL.(string(bone_names(bone_with_CP))).Points;
            % TR.faces    =    temp_STL.(string(bone_names(bone_with_CP))).ConnectivityList(temp_n,:);

            % C = temp_STL.(string(bone_names(bone_with_CP))).Points;
            % figure()
            % patch(TR,'FaceColor', [0.85 0.85 0.85], ...
            % 'EdgeColor','none',...
            %   'FaceLighting','gouraud',...
            %   'AmbientStrength', 0.15);
            % camlight(0,45);
            % material('dull');
            % hold on
            % plot3(C(tri_points,1),C(tri_points,2),C(tri_points,3),'or','markersize',5)
            % hold on
            % plot3(C(iso_check,1),C(iso_check,2),C(iso_check,3),'.b','markersize',5)
            % axis equal

        else
            %% Calculate Coverage Surface Area on Opposing Bone
            temp_center = incenter(temp_STL.(string(bone_names(bone_no_CP))));
            temp_normal = faceNormal(temp_STL.(string(bone_names(bone_no_CP))));

            temp_n = [];
            temp_d = [];

            % https://www.mathworks.com/help/matlab/ref/triangulation.facenormal.html
            % TriangleRayIntersection (orig, dir, vert0, vert1, vert2, varargin)

            parfor (norm_check = 1:length(temp_center),pool)
                [temp_int] = TriangleRayIntersection(temp_center(norm_check,:),temp_normal(norm_check,:),temp_STL.(string(bone_names(bone_with_CP))).Points(temp_STL.(string(bone_names(bone_with_CP))).ConnectivityList(:,1),:),temp_STL.(string(bone_names(bone_with_CP))).Points(temp_STL.(string(bone_names(bone_with_CP))).ConnectivityList(:,2),:),temp_STL.(string(bone_names(bone_with_CP))).Points(temp_STL.(string(bone_names(bone_with_CP))).ConnectivityList(:,3),:));
                if isempty(find(temp_int == 1)) == 0
                    temp_n(norm_check,:) = 1;
                end
            end

            temp_n = find(temp_n == 1);
            pool.IdleTimeout = 30; % resets pool timeout to 30 minutes

            for n = 1:length(temp_n)
                temp_tri = temp_STL.(string(bone_names(bone_no_CP))).ConnectivityList(temp_n(n),:);
                P1 = temp_STL.(string(bone_names(bone_no_CP))).Points(temp_tri(:,1),:);
                P2 = temp_STL.(string(bone_names(bone_no_CP))).Points(temp_tri(:,2),:);
                P3 = temp_STL.(string(bone_names(bone_no_CP))).Points(temp_tri(:,3),:);
                a = P2 - P1;
                b = P3 - P1;
                c = cross(a,b,2);
                area_tri(n,:) = 1/2*sum(sqrt(sum(c.^2,2)));
                clear temp_tri
            end

            Data.(string(subjects(subj_count))).CoverageArea.(string(bone_names(bone_no_CP))) = sum(area_tri);
        end
    end

    %% Calculate Distance and Congruence Index
    k = 1;
    tol = 5;

    A_A = temp_STL.(string(bone_names(bone_no_CP))).Points;
    C_C = temp_STL.(string(bone_names(bone_with_CP))).Points;

    % Pair nodes with CP and calculate euclidean distance
    clear temp i_surf ROI
    for h = 1:length(tri_points(:,1))
        % Kept the line below for legacy
        % ROI = find(temp_STL.(string(bone_names(bone_no_CP))).Points(:,1) >= temp_STL.(string(bone_names(bone_with_CP))).Points(tri_points(h),1)-tol & temp_STL.(string(bone_names(bone_no_CP))).Points(:,1) <= temp_STL.(string(bone_names(bone_with_CP))).Points(tri_points(h),1)+tol & temp_STL.(string(bone_names(bone_no_CP))).Points(:,2) >= temp_STL.(string(bone_names(bone_with_CP))).Points(tri_points(h),2)-tol & temp_STL.(string(bone_names(bone_no_CP))).Points(:,2) <= temp_STL.(string(bone_names(bone_with_CP))).Points(tri_points(h),2)+tol & temp_STL.(string(bone_names(bone_no_CP))).Points(:,3) >= temp_STL.(string(bone_names(bone_with_CP))).Points(tri_points(h),3)-tol & temp_STL.(string(bone_names(bone_no_CP))).Points(:,3) <= temp_STL.(string(bone_names(bone_with_CP))).Points(tri_points(h),3)+tol);
        ROI = find(A_A(:,1) >= C_C(tri_points(h),1)-tol & A_A(:,1) <= C_C(tri_points(h),1)+tol & A_A(:,2) >= C_C(tri_points(h),2)-tol & A_A(:,2) <= C_C(tri_points(h),2)+tol & A_A(:,3) >= C_C(tri_points(h),3)-tol & A_A(:,3) <= C_C(tri_points(h),3)+tol);
        if isempty(ROI) == 0
            parfor n = 1:length(ROI(:,1))
                temp(n,:) = pdist2(C_C(tri_points(h),:),A_A(ROI(n),:),'euclidean');
            end
            tempp = ROI(find(temp == min(temp)));
            i_CP = Data.(string(subjects(subj_count))).(string(bone_names(bone_with_CP))).CP_Bone(find(tri_points(h,1) == Data.(string(subjects(subj_count))).(string(bone_names(bone_with_CP))).CP_Bone(:,2)),1);

            i_surf(k,:) = [i_CP(1) tri_points(h,1) min(temp) 0];
            %     tempp(1) == the index of the paired node on the opposing bone surface
            k = k + 1;
        end
        clear temp tempp
    end

    % Troubleshooting
    % Figure with 'identified nodes' in black, opposing bone nodes in blue, and
    % nodes making up the region of interest as green squares. Also will show
    % the line between the 'identified nodes' and their paired opposing surface
    % nodes with a green line.

%     A = temp_STL.(string(bone_names(bone_with_CP))).Points(tri_points,:);
%     B = temp_STL.(string(bone_names(bone_no_CP))).Points;
%     C = temp_STL.(string(bone_names(bone_with_CP))).Points;
%     
%     figure()
%     plot3(A(:,1),A(:,2),A(:,3),'.k')
%     hold on
%     plot3(B(:,1),B(:,2),B(:,3),'.b')
%     hold on
%     % plot3(C(ROI_frame,1),C(ROI_frame,2),C(ROI_frame,3),'sr') % ROI is from line 324
%     hold on
%     for q = 1:length(i_surf(:,1))
%         plot3([C(i_surf(q,2),1);B(i_surf(q,3),1)],[C(i_surf(q,2),2);B(i_surf(q,3),2)],[C(i_surf(q,2),3);B(i_surf(q,3),3)],'g')
%         hold on
%     end
%     % hold on
%     axis equal
%     axis off

    %% Pull Mean and Gaussian Curvature Data
    % These next few sections calculate the congruence index between each
    % of the paired nodes following the methods described by Ateshian et. al.
    % (CITE)
    mean1 = Data.(string(subjects(subj_count))).(string(bone_names(bone_no_CP))).MeanCurve(i_surf(:,1));
    gaus1 = Data.(string(subjects(subj_count))).(string(bone_names(bone_no_CP))).GaussianCurve(i_surf(:,1));

    mean2 = Data.(string(subjects(subj_count))).(string(bone_names(bone_with_CP))).MeanCurve(i_surf(:,2));
    gaus2 = Data.(string(subjects(subj_count))).(string(bone_names(bone_with_CP))).GaussianCurve(i_surf(:,2));

    %% Principal Curvatures
    parfor n = 1:length(mean1)
        PCMin1(n,:) = mean1(n) - sqrt(mean1(n)^2 - gaus1(n));
        PCMax1(n,:) = mean1(n) + sqrt(mean1(n)^2 - gaus1(n));
    end

    parfor n = 1:length(mean2)
        PCMin2(n,:) = mean2(n) - sqrt(mean2(n)^2 - gaus2(n));
        PCMax2(n,:) = mean2(n) + sqrt(mean2(n)^2 - gaus2(n));
    end

    %% Curvature Differences
    parfor n = 1:length(mean1)
        CD1(n,:) = PCMin1(n,:) - PCMax1(n,:);
    end

    parfor n = 1:length(mean2)
        CD2(n,:) = PCMin2(n,:) - PCMax2(n,:);
    end

    %% Relative Principal Curvatures
    parfor n = 1:length(mean1)
        u = A_A(i_surf(n,1),:);
        v = C_C(i_surf(n,2),:);
        alpha = acosd(dot(u,v)/(norm(u)*norm(v)));
        delta = sqrt(CD1(n)^2 + CD2(n)^2 + 2*CD1(n)*CD2(n)*cosd(2*alpha));
        RPCMin(n,:) = mean1(n) + mean2(n) - 0.5*delta;
        RPCMax(n,:) = mean1(n) + mean2(n) + 0.5*delta;
    end

    %% Overall Congruence Index at a Pair
    parfor n = 1:length(mean1)
        RMS(n,:) = sqrt((RPCMin(n,:)^2 + RPCMax(n,:)^2)/2);
    end

    for n = 1:length(i_surf(:,1))
        i_surf(n,4) = real(RMS(n,:));

    end

    %% Structure of the data being stored
    % i_surf(:,1) == the correspondence particle index identified to the 'identified node' .stl coordinate index
    % Plotting Example: plot3(Data.MD_Mode2_min1.Calcaneus.CP(i_surf(:,1),1),Data.MD_Mode2_min1.Calcaneus.CP(i_surf(:,1),2),Data.MD_Mode2_min1.Calcaneus.CP(i_surf(:,1),3),'*r','markersize',10)
    % i_surf(:,2) == 'identified node' .stl coordinate index
    % Plotting Example: plot3(Data.MD_Mode2_min1.Calcaneus.Calcaneus.Points(i_surf(:,2),1),Data.MD_Mode2_min1.Calcaneus.Calcaneus.Points(i_surf(:,2),2),Data.MD_Mode2_min1.Calcaneus.Calcaneus.Points(i_surf(:,2),3),'*g','markersize',5)
    % i_surf(:,3) == paired .stl coordinate index on opposing surface
    % i_surf(:,4) == euclidean distance from the .stl coordinate to the opposing bone surface
    Data.(string(subjects(subj_count))).MeasureData = i_surf;

    %%
    clearvars -except subjects bone_names Data subj_count frame_count g subj_group pool fldr_name joint_names joint bone_with_CP bone_no_CP joint

end

%% Save Data to .xls file
% An .xlsx file with distance (mm) and congruence index (mm^-1) values at 
% each correspondence particle and node, and the coverage area (mm^2) on
% each bone. Each subject is given their own sheet.
% A .mat file is also generated with similar information.
for n = 1:subj_count
    xlfilename = string(strcat(fldr_name,'\Joint_Measurements_',joint_names(joint),'.xlsx'));
    writematrix(string(joint_names(joint)),xlfilename,'Sheet',string(subjects(n)),'Range','A1');
    writematrix(["CP" "Node" "Distance" "Congruence Index" strcat(string(bone_names(bone_with_CP)),' Coverage Area') strcat(string(bone_names(bone_no_CP)),' Coverage Area')],xlfilename,'Sheet',string(subjects(n)),'Range','A2')
    writematrix(Data.(string(subjects(n))).MeasureData,xlfilename,'Sheet',string(subjects(n)),'Range','A3')
    writematrix(Data.(string(subjects(n))).CoverageArea.(string(bone_names(bone_with_CP))),xlfilename,'Sheet',string(subjects(n)),'Range','E3')
    writematrix(Data.(string(subjects(n))).CoverageArea.(string(bone_names(bone_no_CP))),xlfilename,'Sheet',string(subjects(n)),'Range','F3')

    SaveData.Data.(string(subjects(n))) = Data.(string(subjects(n)));
    save(string(strcat(fldr_name,'\',joint_names(joint),'_Joint_Measurement_SaveData.mat')),'-struct','SaveData')
end

delete(gcp('nocreate'))
