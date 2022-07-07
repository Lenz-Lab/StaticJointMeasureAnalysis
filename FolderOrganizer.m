% Determine the files in the folder selected
FolderPathNameParticles = uigetdir('*.*', 'Select folder with your particle files');
git FolderPathNameVTK = uigetdir('*.*', 'Select folder with your VTK files');
FolderPathNameSTL = uigetdir('*.*', 'Select folder with your STL files');

addpath(FolderPathNameParticles)
addpath(FolderPathNameSTL)
addpath(FolderPathNameVTK)

partfiles = dir(fullfile(FolderPathNameParticles, '*.*'));
partfiles = partfiles(~ismember({partfiles.name},{'.','..'}));

stlfiles = dir(fullfile(FolderPathNameSTL, '*.*'));
stlfiles = stlfiles(~ismember({stlfiles.name},{'.','..'}));

vtkfiles = dir(fullfile(FolderPathNameVTK, '*.*'));
vtkfiles = vtkfiles(~ismember({vtkfiles.name},{'.','..'}));
% 
% temp = strfind(FolderPathNameParticles,'\');
% FolderNameParticles = FolderPathNameParticles(temp(end)+1:end); % Extracts the folder name selected

subject_numbers = ["01"; "02"; "03"; "04"; "05"; "06"; "07"; "08"; "09"; "10"; "11"; "12"; "13"; "14"; "15"; "16"; "17"; "18"; "19"; "20"; "21"; "22"; "23"; "24"; "25"; "26"; "27"];
subject_folder = ["L01"; "L02"; "L03"; "L04"; "L05"; "L06"; "L07"; "L08"; "L09"; "L10"; "L11"; "L12"; "L13"; "R01"; "R02"; "R03"; "R04"; "R05"; "R06"; "R07"; "R08"; "R09"; "R10"; "R11"; "R12"; "R13"; "R14"];

temp = struct2cell(partfiles);
partfiles = temp(1,:);

temp = struct2cell(stlfiles);
stlfiles = temp(1,:);

temp = struct2cell(vtkfiles);
vtkfiles = temp(1,:);

% for m = 1:length(partfiles)
%     for n = 1:length(subject_numbers)
%         if any(strsplit(string(partfiles(m)),'_') == string("Intermediate"))
%             if any(strsplit(string(partfiles(m)),'_') == string(subject_numbers(n)))
%                 if any(strsplit(string(partfiles(m)),'_') == string('local.particles'))
%                     save([strcat('L:\Project_Data\Swiss_WBCT\Healthy\09_NC_Project\01_Joint Measurements\Navicular_Intermediate_Cuneiform_Modes\',string(subject_folder(n)),'\', string(partfiles(m)))])
%                 end
%             end
%         end
%     end
% end
% 
% 
% for m = 1:length(partfiles)
%     for n = 1:length(subject_numbers)
%         if any(strsplit(string(partfiles(m)),'_') == string("Medial"))
%             if any(strsplit(string(partfiles(m)),'_') == string(subject_numbers(n)))
%                 if any(strsplit(string(partfiles(m)),'_') == string('local.particles'))
%                     save([strcat('L:\Project_Data\Swiss_WBCT\Healthy\09_NC_Project\01_Joint Measurements\Navicular_Medial_Cuneiform_Modes\',string(subject_folder(n)),'\', string(partfiles(m)))])
%                 end
%             end
%         end
%     end
% end
% 
% 
% for m = 1:length(partfiles)
%     for n = 1:length(subject_numbers)
%         if any(strsplit(string(partfiles(m)),'_') == string("Lateral"))
%             if any(strsplit(string(partfiles(m)),'_') == string(subject_numbers(n)))
%                 if any(strsplit(string(partfiles(m)),'_') == string('local.particles'))
%                     save([strcat('L:\Project_Data\Swiss_WBCT\Healthy\09_NC_Project\01_Joint Measurements\Navicular_Lateral_Cuneiform_Modes\',string(subject_folder(n)),'\', string(partfiles(m)))])
%                 end
%             end
%         end
%     end
% end

%% Particles Files

for m = 1:length(partfiles)
    for n = 1:length(subject_numbers)
        if any(strsplit(string(partfiles(m)),'_') == string("Navicular"))
            if any(strsplit(string(partfiles(m)),'_') == string(subject_numbers(n)))
                if any(strsplit(string(partfiles(m)),'_') == string('local.particles'))
                    save([strcat('L:\Project_Data\Swiss_WBCT\Healthy\09_NC_Project\01_Joint Measurements\Navicular_Lateral_Cuneiform_Modes\',string(subject_folder(n)),'\', string(partfiles(m)))])
                    save([strcat('L:\Project_Data\Swiss_WBCT\Healthy\09_NC_Project\01_Joint Measurements\Navicular_Medial_Cuneiform_Modes\',string(subject_folder(n)),'\', string(partfiles(m)))])
                    save([strcat('L:\Project_Data\Swiss_WBCT\Healthy\09_NC_Project\01_Joint Measurements\Navicular_Intermediate_Cuneiform_Modes\',string(subject_folder(n)),'\', string(partfiles(m)))])
                end
            end
        end
    end
end

%% STL Files
for m = 1:length(stlfiles)
    for n = 1:length(subject_numbers)
        if any(strsplit(string(stlfiles(m)),'_') == string("Intermediate"))
            if any(strsplit(string(stlfiles(m)),'_') == string(subject_numbers(n)))
                if any(strsplit(string(stlfiles(m)),'_') == string('ICP.stl'))
                    save([strcat('L:\Project_Data\Swiss_WBCT\Healthy\09_NC_Project\01_Joint Measurements\Navicular_Intermediate_Cuneiform_Modes\',string(subject_folder(n)),'\', string(stlfiles(m)))])
                end
            end
        end
    end
end


for m = 1:length(stlfiles)
    for n = 1:length(subject_numbers)
        if any(strsplit(string(stlfiles(m)),'_') == string("Medial"))
            if any(strsplit(string(stlfiles(m)),'_') == string(subject_numbers(n)))
                if any(strsplit(string(stlfiles(m)),'_') == string('ICP.stl'))
                    save([strcat('L:\Project_Data\Swiss_WBCT\Healthy\09_NC_Project\01_Joint Measurements\Navicular_Medial_Cuneiform_Modes\',string(subject_folder(n)),'\', string(stlfiles(m)))])
                end
            end
        end
    end
end


for m = 1:length(stlfiles)
    for n = 1:length(subject_numbers)
        if any(strsplit(string(stlfiles(m)),'_') == string("Lateral"))
            if any(strsplit(string(stlfiles(m)),'_') == string(subject_numbers(n)))
                if any(strsplit(string(stlfiles(m)),'_') == string('ICP.stl'))
                    save([strcat('L:\Project_Data\Swiss_WBCT\Healthy\09_NC_Project\01_Joint Measurements\Navicular_Lateral_Cuneiform_Modes\',string(subject_folder(n)),'\', string(stlfiles(m)))])
                end
            end
        end
    end
end

for m = 1:length(stlfiles)
    for n = 1:length(subject_numbers)
        if any(strsplit(string(stlfiles(m)),'_') == string("Navicular"))
            if any(strsplit(string(stlfiles(m)),'_') == string(subject_numbers(n)))
                if any(strsplit(string(stlfiles(m)),'_') == string('ICP.stl'))
                    save([strcat('L:\Project_Data\Swiss_WBCT\Healthy\09_NC_Project\01_Joint Measurements\Navicular_Lateral_Cuneiform_Modes\',string(subject_folder(n)),'\', string(stlfiles(m)))])
                    save([strcat('L:\Project_Data\Swiss_WBCT\Healthy\09_NC_Project\01_Joint Measurements\Navicular_Medial_Cuneiform_Modes\',string(subject_folder(n)),'\', string(stlfiles(m)))])
                    save([strcat('L:\Project_Data\Swiss_WBCT\Healthy\09_NC_Project\01_Joint Measurements\Navicular_Intermediate_Cuneiform_Modes\',string(subject_folder(n)),'\', string(stlfiles(m)))])
                end
            end
        end
    end
end

%% VTK Files
for m = 1:length(vtkfiles)
    for n = 1:length(subject_numbers)
        if any(strsplit(string(vtkfiles(m)),'_') == string("Intermediate"))
            if any(strsplit(string(vtkfiles(m)),'_') == string(subject_numbers(n)))
                if any(strsplit(string(vtkfiles(m)),'_') == string('groomed.vtk'))
                    save([strcat('L:\Project_Data\Swiss_WBCT\Healthy\09_NC_Project\01_Joint Measurements\Navicular_Intermediate_Cuneiform_Modes\',string(subject_folder(n)),'\', string(vtkfiles(m)))])
                end
            end
        end
    end
end


for m = 1:length(vtkfiles)
    for n = 1:length(subject_numbers)
        if any(strsplit(string(vtkfiles(m)),'_') == string("Medial"))
            if any(strsplit(string(vtkfiles(m)),'_') == string(subject_numbers(n)))
                if any(strsplit(string(vtkfiles(m)),'_') == string('groomed.vtk'))
                    save([strcat('L:\Project_Data\Swiss_WBCT\Healthy\09_NC_Project\01_Joint Measurements\Navicular_Medial_Cuneiform_Modes\',string(subject_folder(n)),'\', string(vtkfiles(m)))])
                end
            end
        end
    end
end


for m = 1:length(vtkfiles)
    for n = 1:length(subject_numbers)
        if any(strsplit(string(vtkfiles(m)),'_') == string("Lateral"))
            if any(strsplit(string(vtkfiles(m)),'_') == string(subject_numbers(n)))
                if any(strsplit(string(vtkfiles(m)),'_') == string('groomed.vtk'))
                    save([strcat('L:\Project_Data\Swiss_WBCT\Healthy\09_NC_Project\01_Joint Measurements\Navicular_Lateral_Cuneiform_Modes\',string(subject_folder(n)),'\', string(vtkfiles(m)))])
                end
            end
        end
    end
end

for m = 1:length(vtkfiles)
    for n = 1:length(subject_numbers)
        if any(strsplit(string(vtkfiles(m)),'_') == string("Navicular"))
            if any(strsplit(string(vtkfiles(m)),'_') == string(subject_numbers(n)))
                if any(strsplit(string(vtkfiles(m)),'_') == string('groomed.vtk'))
                    save([strcat('L:\Project_Data\Swiss_WBCT\Healthy\09_NC_Project\01_Joint Measurements\Navicular_Lateral_Cuneiform_Modes\',string(subject_folder(n)),'\', string(vtkfiles(m)))])
                    save([strcat('L:\Project_Data\Swiss_WBCT\Healthy\09_NC_Project\01_Joint Measurements\Navicular_Medial_Cuneiform_Modes\',string(subject_folder(n)),'\', string(vtkfiles(m)))])
                    save([strcat('L:\Project_Data\Swiss_WBCT\Healthy\09_NC_Project\01_Joint Measurements\Navicular_Intermediate_Cuneiform_Modes\',string(subject_folder(n)),'\', string(vtkfiles(m)))])
                end
            end
        end
    end
end




