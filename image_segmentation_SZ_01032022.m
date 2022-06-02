%code created to manually segment OCT images. 
close all; 
clear all; clc

%%
addpath 'U:\eng_research_nialab\users\suezhang\Codes'
addpath 'U:\eng_research_nialab\users\suezhang\Codes\ForceAnalysis_V3_12-16-21'



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

%% saving information 
savepathfolder=path;
save_filename=erase(filename,'.tif')
bead_number=input('What bead number is this?');
bead_name=['bead' num2str(bead_number)]

%%

[rois,z_increment] = meter_maidSZ(image,savepathfolder,save_filename,bead_name);
%%
xy_res = input('Cool. Now, what is the x,y pixel resolution in um/px?\n');
z_res_OG = input('Thanks. And how many um per slice in z?\n');

%% i think this is just for the images
% aug_factor = z_res_OG/xy_res; 
% z_res = z_res_OG/aug_factor;

%%
xdata=[];
ydata=[];
zdata=[];

figure(2)
for i=1:length(rois)
xdata=vertcat(xdata,rois{1,i}{1,1}(:,1));
ydata=vertcat(ydata,rois{1,i}{1,1}(:,2));
zdata=vertcat(zdata,i*z_increment*ones(length(rois{1,i}{1,1}(:,1)),1));

plot3(rois{1,i}{1,1}(:,1),rois{1,i}{1,1}(:,2),i*z_increment*ones(length(rois{1,i}{1,1}(:,1)),1),'marker','.','linestyle','none')
axis equal
hold on
end
title('Unscaled point cloud')
hold off

xdata=xdata.*xy_res;
ydata=ydata.*xy_res;
zdata=zdata.*z_res_OG;
% xdata=xdata;
% ydata=ydata;
% zdata=zdata;

j=boundary(xdata,ydata,zdata); %create solid object of point cloud
figure(3)
plot3(xdata,ydata,zdata,'markersize',20,'marker','.','linestyle','none')
hold on
trisurf(j,xdata,ydata,zdata,'Facecolor','m','FaceAlpha',0.1)
title('Scaled point cloud')
axis equal
xlabel('x (um)');
ylabel('y (um)');
zlabel('z (um)');
hold off

%% ellipsoid fitting

Data = horzcat(xdata,ydata,zdata);

k = randperm(length(xdata));
Data_mat = Data(k(1:round(length(xdata)*0.6)),:);
Data_mat = rmmissing(Data_mat);

[u,CF,A,b,c]=Ellipsoid3D_Fitting_DR_SVD(Data_mat',200);


%%Find radius of ellipsoid: 
A = [u(1)   u(4)/sqrt(2)  u(5)/sqrt(2)  u(7)/2;...
    u(4)/sqrt(2)  u(2)    u(6)/sqrt(2)  u(8)/2;...
    u(5)/sqrt(2)  u(6)/sqrt(2)  u(3)    u(9)/2;...
    u(7)/2  u(8)/2  u(9)/2  u(10)];

A3 = A(1:3,1:3);
A3_inv = inv(A3); 
ofs = u(7:9)/2.0; 
center = (-A3_inv*ofs); 

Tofs = eye(4); 
Tofs(4,1:3) = center; 
R = (Tofs*(A*(Tofs)'));

[evecs,evals] = eig(R(1:3,1:3)/-R(4,4)); 

radii = sqrt(1./diag(abs(evals)));
sgns = sign(diag(evals));
radii = radii.*sgns;

[xc,yc,zc] = ellipsoid(0,0,0,radii(1),radii(2),radii(3));

a = kron(evecs(:,1),xc); 
b = kron(evecs(:,2),yc);
c = kron(evecs(:,3),zc); 

data = a+b+c; n = size(data,2); 

x_new  = data(1:n,:)+center(1); 
y_new = data(n+1:2*n,:)+center(2);
z_new = data(2*n+1:end,:)+center(3);

%%test ellipsoid plot
mind = min([ Data_mat(:,1)-10 Data_mat(:,2)-10 Data_mat(:,3)-10]); 
maxd = max([ Data_mat(:,1)+10 Data_mat(:,2)+10 Data_mat(:,3)+10]); 
nsteps = 100; 
 
[x,y,z] = meshgrid(linspace(mind(1),maxd(1),nsteps),...
    linspace(mind(2),maxd(2),nsteps),...
    linspace(mind(3),maxd(3),nsteps));

Ellipsoid = u(1)*x.^2 + u(2)*y.^2 + u(3)*z.^2 +...
    sqrt(2)*u(4)*x.*y + sqrt(2)*u(5)*x.*z + ...
    sqrt(2)*u(6)*y.*z + u(7)*x + u(8)*y + u(9)*z + u(10);

j2=boundary(Data_mat(:,1),Data_mat(:,2),Data_mat(:,3));
%%plotting ellipsoid

figure(4)
hold on
trisurf(j2,Data_mat(:,1),Data_mat(:,2),Data_mat(:,3),'Facecolor','m','FaceAlpha',0.1)
hold on
scatter3(Data_mat(:,1),Data_mat(:,2),Data_mat(:,3),300,'marker','.');
axis equal
h = surf(x_new,y_new,z_new);
alpha(0.3)
grid on
xlabel('x (um)');
ylabel('y (um)');
zlabel('z (um)');
view(90,0)
scatter3(center(1), center(2),center(3),'filled','MarkerFaceColor','m')

plot_evector1=[center(1)+radii(1)*evecs(1,1),center(2)+radii(1)*evecs(2,1),center(3) + radii(1)*evecs(3,1)];
plot_evector2=[center(1)+radii(2)*evecs(1,2),center(2)+radii(2)*evecs(2,2),center(3) + radii(2)*evecs(3,2)];
plot_evector3=[center(1)+radii(3)*evecs(1,3),center(2)+radii(3)*evecs(2,3),center(3) + radii(3)*evecs(3,3)];
plot_evector4=[center(1)-radii(1)*evecs(1,1),center(2)-radii(1)*evecs(2,1),center(3) - radii(1)*evecs(3,1)];
plot_evector5=[center(1)-radii(2)*evecs(1,2),center(2)-radii(2)*evecs(2,2),center(3) - radii(2)*evecs(3,2)];
plot_evector6=[center(1)-radii(3)*evecs(1,3),center(2)-radii(3)*evecs(2,3),center(3) - radii(3)*evecs(3,3)];

Plotting_vectors = [plot_evector1; plot_evector2; plot_evector3; plot_evector4; plot_evector5; plot_evector6];

 
plot3([center(1),plot_evector1(1)],[center(2),plot_evector1(2)],[center(3),plot_evector1(3)],'k-','LineWidth',1.5)
plot3([center(1),plot_evector2(1)],[center(2),plot_evector2(2)],[center(3),plot_evector2(3)],'r-','LineWidth',1.5)
plot3([center(1),plot_evector3(1)],[center(2),plot_evector3(2)],[center(3),plot_evector3(3)],'b-','LineWidth',1.5)
plot3([center(1),plot_evector4(1)],[center(2),plot_evector4(2)],[center(3),plot_evector4(3)],'k-','LineWidth',1.5)
plot3([center(1),plot_evector5(1)],[center(2),plot_evector5(2)],[center(3),plot_evector5(3)],'r-','LineWidth',1.5)
plot3([center(1),plot_evector6(1)],[center(2),plot_evector6(2)],[center(3),plot_evector6(3)],'b-','LineWidth',1.5)
p = patch(isosurface(x,y,z,Ellipsoid,1)); 
set(p,'FaceColor','g','EdgeColor','none'); 
alpha(p,0.3)


%fprintf('Great! We''re all done!\n\n')

%fprintf('Center of ellipsoid is at x = %.2f, y = %.2f, z = %.4f\n',center(1),center(2),center(3))
%fprintf('Final radii are:   %.2f along vector 1, %.2f along vector 2, %.2f along vector 3\n',radii(1),radii(2),radii(3))
diameters = radii.*2; 
diameters_display=diameters(:,1)
aspect_ratio_disp=diameters_display(1)/diameters_display(3)

hold off


%% saving information
clear image
saveas(figure(3),[savepathfolder,save_filename,bead_name,'pointcloud2.fig']);
saveas(figure(4),[savepathfolder,save_filename,bead_name,'ellipsoid2.fig']);
save([savepathfolder,save_filename,bead_name, 'manual.mat'])