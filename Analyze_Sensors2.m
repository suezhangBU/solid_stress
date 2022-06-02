%%%Image segmentation & Volume Fitting program for spherical or ellipsoidal
%%%bodies.

%%%Last edits: 12/16/2021:  moved sections, adding file saving to specific folder
...change ellipsoid fitting to Douglas-rachford
%%%iterative, from Fitellipsoid paper

%%%General flow: 
%%%     1/ load tiff
%%%     2/ select Z ranges, XY bounding box for each bead
%%%     3/ For each bead:  augment data to increase z resolution, fit
%%%     sigmoidal decays to find boundary, map point cloud, and fit
%%%     ellipsoid
%%%     4/ Reports ellipsoid diameters  
%%%     5/ Plots on original image with spheroid and beads (8 bit)



%%%only inputs is cropped single channel Tiff around object of interest. X and Y must be same dimensions. 

close all; 
clear all; clc
tStart = tic;
%%Open and load tiff 
%%First open tiff:
[filename,path] = uigetfile('*.tif','Select your file');
if isequal(filename,0)
    disp('User selected cancel. Goodbye!')
    return
else
    disp(['OK, I''ll analyze ',fullfile(path,filename)])
end
%%Then load tiff:
t = Tiff(filename, 'r');
image(:,:,1) = t.read(); % Read the first image to get the array dimensions correct.
if t.lastDirectory()
    return;              %If the file only contains one page, we do not need to continue.
end
t.nextDirectory();       %Read all remaining pages (directories) in the file
while true
    image(:,:,end+1) = t.read();
    if t.lastDirectory()
        break;
    else
        t.nextDirectory();
    end
end

disp(['Great! The tiff is loaded. Let''s find our beads.'])

save_filename=erase(filename,'.tif')
% bead_number=input('What bead number is this?');
bead_number=1;
bead_name=['bead' num2str(bead_number)]

rois = meter_maid(image,path,save_filename,bead_number); 

%%
xy_res = input('Cool. Now, what is the x,y pixel resolution in um/px?\n');
z_res_OG = input('Thanks. And how many um per slice in z?\n');

%Parameterize variables:
ui='d';
%ui = input('Would you like custom settings or default? \n(default: dtheta, dphi = 5, n_r = 75\nPlease answer d/c\n','s');
if ui == 'd' || ui == 'D'
    theta_n = 5; phi_n = 5; r_n = 75; 
    fprintf('OK! default settings activated.\n')
else 
    theta_n = input('What is theta step value?\nTypical values are 1,5,10\n');
    phi_n = input('What is phi step value?\nTypical value are 1,5, or 10\n');
    r_n = input('How many r vectors would you like to try?\nCode optimization was completed with 75.\n');
end
run = 1; 


for n = 1:length(rois)
   x_range = rois{1,n}{1}(1):rois{1,n}{1}(3);
   y_range = rois{1,n}{1}(5):rois{1,n}{1}(6);
   z_range = rois{1,n}{2}:rois{1,n}{3};
   subimage = image(y_range(1):y_range(end),x_range(1):x_range(end),z_range(1):z_range(end));
   [radii,center,diameters,evecs,z_res,u,Axes_points,diameters_save] = ...
       Ellipsoid_Fitting_fxn2(n,subimage,xy_res,z_res_OG,theta_n,phi_n,r_n,run,path,save_filename);
   diameters(n,1:3) = diameters'; 
   centers_um(n,1:3) = center';
   u_all(:,n) = u; 

    centers_px(n,:) = round(center'/(xy_res.*eye(3,3))); centers_px(n,3) = round(center(3)/(4));
    absolute_points(n,:) = [...
        (y_range(1) + centers_px(n,1)),...
        (x_range(1) + centers_px(n,2)),...
        (z_range(1)-1+centers_px(n,3))];

    F_rad(n,:) = radii'; 
    F_Vec{n} = evecs;
    
    
    diameters_all(:,n)=diameters_save;
  
    
    if n == length(rois)
        fprintf('Yay! Analysis complete :D\n\n')
    end
end
% %%
% F_rad_backup = F_rad; 
% F_Vec_backup = F_Vec;
% for n = 1:length(F_Vec)
%      for m = 1:3
%          F_Vec2{n}(m,:) = F_Vec{n}(m,:)/(xy_res.*eye(3,3)); 
% %     for n = 1:2:
% %         F_Vec2{n}
% 
%     end
%     F_Vec2{n}(:,3) = F_Vec2{n}(:,3)/xy_res; 
% end
% %centers_px(n,3) = round(center(3) * xy_res);


Plotting_vecs = [];

for n = 1:length(rois);
plot_evector1 = [absolute_points(n,1) + (F_rad(n,1)/0.3551)*F_Vec{n}(1,1),absolute_points(n,2) + (F_rad(n,1)/0.3551)*F_Vec{n}(2,1),absolute_points(n,3) + F_rad(n,1)*0.3551*F_Vec{n}(3,1)];
plot_evector2 = [absolute_points(n,1) + (F_rad(n,2)/0.3551)*F_Vec{n}(1,2),absolute_points(n,2) + (F_rad(n,2)/0.3551)*F_Vec{n}(2,2),absolute_points(n,3) + F_rad(n,2)*0.3551*F_Vec{n}(3,2)];
plot_evector3 = [absolute_points(n,1) + (F_rad(n,3)/0.3551)*F_Vec{n}(1,3),absolute_points(n,2) + (F_rad(n,3)/0.3551)*F_Vec{n}(2,3),absolute_points(n,3) + F_rad(n,3)*0.3551*F_Vec{n}(3,3)];
plot_evector4 = [absolute_points(n,1) - (F_rad(n,1)/0.3551)*F_Vec{n}(1,1),absolute_points(n,2) - (F_rad(n,1)/0.3551)*F_Vec{n}(2,1),absolute_points(n,3) - F_rad(n,1)*0.3551*F_Vec{n}(3,1)];
plot_evector5 = [absolute_points(n,1) - (F_rad(n,2)/0.3551)*F_Vec{n}(1,2),absolute_points(n,2) - (F_rad(n,2)/0.3551)*F_Vec{n}(2,2),absolute_points(n,3) - F_rad(n,2)*0.3551*F_Vec{n}(3,2)];
plot_evector6 = [absolute_points(n,1) - (F_rad(n,3)/0.3551)*F_Vec{n}(1,3),absolute_points(n,2) - (F_rad(n,3)/0.3551)*F_Vec{n}(2,3),absolute_points(n,3) - F_rad(n,3)*0.3551*F_Vec{n}(3,3)];
Plotting_vecs{n} = [plot_evector1;plot_evector2;plot_evector3;plot_evector4;plot_evector5; plot_evector6];
end

% discretize vectors:
nit = 25; 
discretized_all = [];
for n = 1:length(rois)
    for m = 1:6
        x_1 = linspace(absolute_points(n,1),Plotting_vecs{n}(m,1),nit);  
        y_1 = linspace(absolute_points(n,2),Plotting_vecs{n}(m,2),nit); 
        z_1 = linspace(absolute_points(n,3),Plotting_vecs{n}(m,3),nit);
        for t = 1:nit; 
            discretized{n}{m}(1,t) = absolute_points(n,1)+(((x_1(t)-absolute_points(n,1))/nit)*t);
            discretized{n}{m}(2,t) = absolute_points(n,2)+(((y_1(t)-absolute_points(n,2))/nit)*t);
            discretized{n}{m}(3,t) = absolute_points(n,3)+(((z_1(t)-absolute_points(n,3))/nit)*t);
         end
        %discretized{n}{m}(1,t+1) = cellstr('r');
         discretized{n}{m}(3,:) = round(discretized{n}{m}(3,:));

    end
    discretized_all = vertcat(discretized_all,discretized{n}{1},discretized{n}{2},discretized{n}{3},...
        discretized{n}{4},discretized{n}{5},discretized{n}{6});
end

c_values = {[1 0 0];[1 0 0];[1 0 0];[0 1 1];[0 1 1];[0 1 1];[1 1 0];[1 1 0];[1 1 0]};
more = size(discretized_all,1); 
c_values = repmat(c_values,more/6,1);

% %%
% for n = 1:length(rois)
%     sprintf('The diameters of bead %d are: %.2f, %.2f, %.2f\n',n,diameters(n,1),diameters(n,2),diameters(n,3))
%     
%     sprintf('The centroid is at:\n')
%     disp(center)
%     
% end


%%
%save file\
date=char(datetime);
save_filename=erase(filename,'.tif')
% bead_number=input('What bead number is this?');
% bead_name=['bead' num2str(bead_number)]

save([path,save_filename, '.mat'])

% %% display image again to overlay diameters
% [filename,path] = uigetfile('*.tif','Select your file');
% if isequal(filename,0)
%     disp('User selected cancel. Goodbye!')
%     return
% else
% end
% %%then load tiff:
% t = Tiff(filename, 'r');
% fullimage(:,:,1) = t.read(); % Read the first image to get the array dimensions correct.
% if t.lastDirectory()
%     return;              %If the file only contains one page, we do not need to continue.
% end
% t.nextDirectory();       %Read all remaining pages (directories) in the file
% while true
%     fullimage(:,:,end+1) = t.read();
%     if t.lastDirectory()
%         break;
%     else
%         t.nextDirectory();
%     end
% end
% %%
% 
% 
% imagineFile_KR(fullimage,absolute_points,Plotting_vecs,discretized_all,c_values);
% 
