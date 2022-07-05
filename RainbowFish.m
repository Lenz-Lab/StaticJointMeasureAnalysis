function RainbowFish(Bone,NodalIndex,AllNodalData,CLimits,ColorMap_Flip,ROI)
% clearvars -except Data, clc
%
P = stlread('Particle.stl');
PP.Points = P.Points/max(max(P.Points));
T = stlread('Particle.stl');
TT.Points = T.Points/max(max(T.Points));
%
%
% Bone = Data.TAR_01.Calcaneus.Calcaneus;
% NodalIndex = Data.TAR_01.MeasureData.F_50(:,2);
% NodalData = Data.TAR_01.MeasureData.F_50(:,3);
% SPMIndex = NodalIndex; % NodalIndex([574;575;578;580],:);

% CLimits = [min(NodalData) 6]; %max(NodalData)];
% ColorMap_Flip = 1;
%
B.faces        = Bone.ConnectivityList;
B.vertices     = Bone.Points;

glyph_scale = 1.25;

Ptemp.faces     = P.ConnectivityList;
Ptemp.vertices  = PP.Points*glyph_scale;

Pspm.faces      = T.ConnectivityList;
Pspm.vertices   = TT.Points*glyph_scale;

if ColorMap_Flip == 1 % 1 = Distance
    NodalData = AllNodalData(:,3);
elseif ColorMap_Flip == 2 % 2 = Congruency
    NodalData = AllNodalData(:,4);
end

%% Create Bin Structure
% Initialize variable for while loops
k = 1;

% Set colormap, then close opened colormap
% if ColorMap_Flip == 2
    ColorMap2 = colormap(jet);
% else
%     ColorMap2 = colormap(flipud(jet));
% end
% close

% Colormap length variable
ML = length(ColorMap2(:,1));

% This section will check if the variable CLimits was an input and create
% bins for the nodal data to be binned in the future scatter3 colormap
% If CLimits does not exist will create k number of bins from 0 to length(ColorMap2)
if exist('CLimits') == 0
    while k <= ML
        if k == 1
            S.BinRange(k,:) = [0 (1/ML)];
        end
        if k > 1 && k < ML
            S.BinRange(k,:) = [S.BinRange((k-1),2) S.BinRange((k-1),2)+(1/ML)];
        end
        if k == ML
            S.BinRange(k,:) = [S.BinRange((k-1),2) inf];
        end
        k = k + 1;
    end
end

if exist('CLimits') == 1
    % If CLimits exists will create k number of bins in different ranges based
    % on different conditions.

    % Condition 1: length(CLimits) == 1
    % Creates bins from minimum of NodalData to maximum of NodalData
    if length(CLimits) == 1
        while k <= ML
            if k == 1
                S.BinRange(k) = [min(NodalData) min(NodalData)+(1/ML)*(max(NodalData)-min(NodalData))];
            end
            if k > 1 && k < ML
                S.BinRange(k) = [S.BinRange((k-1),2) S.BinRange((k-1),2)+((1/ML)*(max(NodalData)-min(NodalData)))];
            end
            if k == ML
                S.BinRange(k) = [S.BinRange((k-1),2) inf];
            end
            k = k + 1;
        end
    end

    % Condition 2: length(CLimits) > 1
    % Creates bins from CLimits(1,1) to CLimits (1,2), any above will be
    % placed in last bin
    if length(CLimits) > 1
        while k <= ML
            if k == 1
                S.BinRange(k,:) = [CLimits(1,1) CLimits(1,1)+(1/ML)*(CLimits(1,2)-CLimits(1,1))];
            end
            if k > 1 && k < ML
                S.BinRange(k,:) = [S.BinRange((k-1),2) S.BinRange((k-1),2)+((1/ML)*(CLimits(1,2)-CLimits(1,1)))];
            end
            if k == ML
                S.BinRange(k,:) = [S.BinRange((k-1),2) inf];
            end
            k = k + 1;
        end
    end

    % Condition 3: if ColorMap_Flip == 1
    % Flips the variable ColorMap2, this is what will change lower values
    % from blue to red and higher values from red to blue.
    if ColorMap_Flip == 1
        k = 0;
        for n = 1:length(ColorMap2(:,1))
            temp(n,:) = ColorMap2(end-k,:);
            k = k + 1;
        end
        clear ColorMap2
        ColorMap2 = temp;
    end
end

%% Nodal Data Index placement using Bins
% Parallel for loop takes each data point from NodalData and pairs the index
% it with a ColorMap2 value (CMap output variable) using the previously
% created bins
for n = 1:length(NodalData(:,1))
    k = 1;
    while k <= ML
        if NodalData(n,1) >= S.BinRange(k,1) && NodalData(n,1) < S.BinRange(k,2)
            CMap(n,:) = ColorMap2(k,:);
        end
        k = k + 1;
    end
end

%%
% figure()
patch(B,'FaceColor', [0.85 0.85 0.85], ...
    'EdgeColor','none',...
    'FaceLighting','gouraud',...
    'AmbientStrength', 0.15);
material('dull');
% view([180 -70]) % Subtalar
% view([0 -20]) % Talonavicular
view([10 0]) % Calcaneocuboid
camlight headlight
hold on
set(gcf,'Units','Normalized','OuterPosition',[-0.0036 0.0306 0.5073 0.9694]);

if isempty(ROI) == 0
    for n = 1:length(ROI(:,1))
        PR = Ptemp;
        PR.vertices = PR.vertices + ROI(n,:);
        if AllNodalData(n,3) >=  0 && AllNodalData(n,3) <= 6
            patch(PR,'FaceColor', CMap(n,:), ...
                'EdgeColor','none',...
                'FaceLighting','gouraud',...
                'AmbientStrength', 0.15);
            material('dull');
            hold on
            hold on
            scatter3(ROI(n,1),ROI(n,2),ROI(n,3),50,[1 1 1],'linewidth',3)
            hold on
        end
        if AllNodalData(n,3) <  0 || AllNodalData(n,3) > 6
            patch(PR,'FaceColor', [1 1 1], ...
                'EdgeColor','none',...
                'FaceLighting','gouraud',...
                'AmbientStrength', 0.15,'facealpha',1);
            material('dull');
            hold on
            %         purple [0.4940 0.1840 0.5560]
        end
    end
else
    for n = 1:length(NodalIndex(:,1))
        PR = Ptemp;
        PR.vertices = PR.vertices + B.vertices(NodalIndex(n,:),:);

        if AllNodalData(n,3) >=  0 && AllNodalData(n,3) <= 6
            %         temp = find(NodalIndex(n,:) == SPMIndex);
            %         if isempty(temp) == 1
            patch(PR,'FaceColor', CMap(n,:), ...
                'EdgeColor','none',...
                'FaceLighting','gouraud',...
                'AmbientStrength', 0.15);
            material('dull');
            hold on
            %         elseif isempty(temp) == 0
            %             TR = Pspm;
            %             TR.vertices = TR.vertices + B.vertices(NodalIndex(n,:),:);
            %             patch(TR,'FaceColor', CMap(n,:), ...
            %                 'EdgeColor','none',...
            %                 'FaceLighting','gouraud',...
            %                 'AmbientStrength', 0.15,'facealpha',1 ) %,...
            %             %'edgecolor',[0 0 0],'edgealpha',1);
            %             material('dull');
            hold on
            scatter3(Bone.Points(NodalIndex(n,:),1),Bone.Points(NodalIndex(n,:),2),Bone.Points(NodalIndex(n,:),3),50,[1 1 1],'linewidth',3)
            %             if NodalData(n,:) >= mean(CLimits)*0.33       && ColorMap_Flip == 0
            %                 scatter3(Bone.Points(NodalIndex(n,:),1),Bone.Points(NodalIndex(n,:),2),Bone.Points(NodalIndex(n,:),3),50,[0 0 0],'linewidth',3)
            %             elseif NodalData(n,:) < mean(CLimits)*0.67    && ColorMap_Flip == 0
            %                 scatter3(Bone.Points(NodalIndex(n,:),1),Bone.Points(NodalIndex(n,:),2),Bone.Points(NodalIndex(n,:),3),50,[1 1 1],'linewidth',3)
            %             elseif NodalData(n,:) >= mean(CLimits)*0.67   && ColorMap_Flip == 1
            %                 scatter3(Bone.Points(NodalIndex(n,:),1),Bone.Points(NodalIndex(n,:),2),Bone.Points(NodalIndex(n,:),3),50,[1 1 1],'linewidth',3)
            %             elseif NodalData(n,:) < mean(CLimits)*0.33    && ColorMap_Flip == 1
            %                 scatter3(Bone.Points(NodalIndex(n,:),1),Bone.Points(NodalIndex(n,:),2),Bone.Points(NodalIndex(n,:),3),50,[0 0 0],'linewidth',3)
            %             end
            hold on
        end

        if AllNodalData(n,3) <  0 || AllNodalData(n,3) > 6
            patch(PR,'FaceColor', [1 1 1], ...
                'EdgeColor','none',...
                'FaceLighting','gouraud',...
                'AmbientStrength', 0.15,'facealpha',1);
            material('dull');
            hold on
            %         purple [0.4940 0.1840 0.5560]
        end
    end
end


hold on
axis equal
grid off
set(gca,'xtick',[],'ytick',[],'ztick',[],'xcolor','none','ycolor','none','zcolor','none')
% view([-185 -70])
% camlight(0,0)
if ColorMap_Flip == 2
    colormap jet
elseif ColorMap_Flip == 1
    colormap(flipud(jet))
end
C = colorbar;
C.FontSize = 32;
caxis([CLimits(1,1),CLimits(1,2)])
set(C, 'ylim',[CLimits(1,1),CLimits(1,2)])
%             disp('Congruence')

