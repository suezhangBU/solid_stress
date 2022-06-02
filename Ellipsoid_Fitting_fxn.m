function [radii,center,diameters,evecs,z_res,u,Plotting_vectors,save_coords] = Ellipsoid_Fitting_fxn(n,image,xy_res,z_res_OG,theta_n,phi_n,r_n)

%%Extend Z resolution:   increase to get similar dimensions as XY scale
image = double(image);             %convert from uint8
S = size(image);                      %get size dimensions of image stack
B = permute(image,[3 1 2]);           %permute so that dimensions go from [x y z] to [z x y]
%aug_factor = round((mean(S(1:2))/S(3))); %augment data to get it on similar scale as XY
    aug_factor = z_res_OG/xy_res; 
    
    Z = interp1(B,linspace(1,S(3),aug_factor*S(3))); %Increases to get similar dimensions as original square:  1D Interpollation along 1st dimension (which in B is our 3rd dimension)
    extended_Subimage = permute(Z, [2 3 1]);         %permute again so dimensions go from [z x y] to [x y z]
    
    z_res = z_res_OG/aug_factor;


% fidelity = input('\n\n\nDone. Would you like to check the fidelity of the permutations?\nIt might take a little while if you say yes\n','s');
% if strcmp(fidelity,'yes') == 1 || strcmp(fidelity,'Yes') == 1 || strcmp(fidelity,'y') == 1 || strcmp(fidelity,'Y') ==1
%     disp('OK! This might take a second')
%     figure()
%     subplot(1,2,1)
%     title('original image')
%     hold on
%     for n = 1:size(subimage,3)
%        n_vec = n*ones(1,length(subimage));
%         for m = 1:length(subimage)
%             scatter3(subimage(:,m,n),subimage(m,:,n),n_vec);
%         end
%     end
%     view(90,0)
%     subplot(1,2,2)
%     title('extended image')
%     hold on
%     for n = 1:5:size(extended_Subimage,3)
%         for m = 1:2:size(extended_Subimage,2)
%             scatter3(extended_Subimage(:,m,n),extended_Subimage(m,:,n),n*ones(1,100))
%         end
%     end
%     view(90,0)
% else
% end
pause(1)
%fprintf('OK cool. Now we''re working with better resolution.\nNext we''ll iterate through to map intensity values.\n')

%%Iterate through rays
%First find centroid value at center plane, in form [X,Y,Z] of total stack
centroid = [ceil(size(extended_Subimage)/2)];
centroid_um = [ceil(size(extended_Subimage)/2)]*[xy_res 0 0; 0 xy_res 0; 0 0 z_res];

dr = 1;

I = zeros(r_n/dr,360/theta_n,180/phi_n); %Intensity values as [f(r) f(theta) f(phi)]
rvec = linspace(0,r_n,r_n/dr);

for theta = theta_n:theta_n:360;
    for phi = phi_n:phi_n:180;
        for r = dr:dr:r_n
            r_x = centroid(1) + round((r * cosd(theta) * sind(phi)),0);
            r_y = centroid(2) + round((r * sind(theta) * sind(phi)),0);
            r_z = centroid(3) + round((r * cosd(phi)),0);
              if r_x <=0 || r_y <= 0 || r_z <= 0
                  I(r/dr,(theta/theta_n),(phi/phi_n)) = 0; 
              elseif r_x > size(extended_Subimage,1) ||...
                      r_y > size(extended_Subimage,2) ||...
                      r_z > size(extended_Subimage,3)
                  I(r/dr,theta/theta_n, phi/phi_n) = 0;
              else
                I(r/dr,(theta/theta_n),(phi/phi_n)) =...
                    extended_Subimage(r_x, r_y, r_z);
              end   
        end
    end
end


%fprintf('Cool. Now we''ll fit the sigmoid functions. This may take a little while...\n')    
fun = @(x,xdata)...
    (x(5) + ...
    (x(4) + x(3).*(xdata - x(2)))./...
    (1+exp(x(1).*(xdata-x(2)))));
x0 = [0.25 1 15 4 6];                     %initial guess for variables
Boundary = zeros(360/theta_n,180/phi_n);  %In format [azimuth elevation] 
r_dummy = rvec;                           %dummy r vector for plotting function
opts = optimoptions('lsqcurvefit','Display','off','Algorithm','levenberg-marquardt');         %suppress lsqcurvefit output statements

tic

for n = 1:(180/phi_n) %
    for m = 1:(360/theta_n)
        zeros_I = find(I(:,m,n) == 0);
        cropped = length(rvec)-ceil(length(zeros_I)/2);
        [solution,resnorm] =lsqcurvefit(fun,x0,rvec(1:cropped)',I(1:cropped,m,n),[],[],opts); 
    %solution = lsqcurvefit(fun,x0,rvec',I(:,m,n),[],[],opts);  %solution should minimize LS between function and r vector input
    solved = fun(solution,r_dummy);                            %solve the function with optimized coefficients 
    diff = abs(solved - (0.5*max(solved)));                   %Find ~75% value of max **USING FIT**
    edge = find(diff == min(diff));                            %find index of 75% value (rvector value)
    if length(edge) > 1
        Boundary(m,n) = NaN; 
    else
        Boundary(m,n) = r_dummy(edge);                             %store radial value
%         figure()
%         plot(rvec,I(:,m,n),'ko',r_dummy,solved,'b-'); hold on; 
%         scatter(r_dummy(edge),solved(edge),50,'k','o','filled') 
%         title(sprintf(' \\theta = %d:',m*theta_n));
    end
    end
    sprintf('phi = %d is completed',(n*phi_n))
end

toc; 

%make matrix to store the boundary values in x,y,z --> iterates through all
%the phi values and theta values for data organized as 
% [xyz coordinates, theta values, phi values] 
mapped = zeros(3,360/theta_n,180/phi_n);
for n = 1:(180/phi_n)
    for m = 1:(360/theta_n)
        mapped(1,m,n) = round(Boundary(m,n) * cosd(m*theta_n) * sind(n*phi_n))*xy_res; %x coordinates
        mapped(2,m,n) = round(Boundary(m,n) * sind(m*theta_n) * sind(n*phi_n))*xy_res; %y coordinates
        mapped(3,m,n) = round(Boundary(m,n) * cosd(n*phi_n))*z_res;                   %z coordinates
    end
end
disp(['PS:  All your data is in microns now.'])


plotYN = 'n'; %input('All done! Do you want to plot the data over an image?\nEnter Y/N\n','s');
if plotYN == 'y' || plotYN == 'Y'
    [project,path2]=uigetfile('.tif');m
    pic = read(Tiff(project))'; %beware! image is flipped!
    figure()
    imshow(pic)
    hold on
elseif plotYN == 'n' || plotYN == 'N'
    figure()
    hold on
end
for n = 1:(180/phi_n)
    %scatter3(mapped(1,:,n)*xy_res+centroid(1),mapped(2,:,n)*xy_res+centroid(2),mapped(3,:,n)*z_res+centroid(3));
    scatter3(mapped(1,:,n)+centroid_um(1),mapped(2,:,n)+centroid_um(2),mapped(3,:,n)+centroid_um(3));

end
xlabel('X axis (um)')
ylabel('Y axis (um)')
zlabel('Z axis (um)')
axis equal

Save_coords = input('Do you want me to save the raw XYZ points of the point cloud for you? 1 for yes, 0 for no')
if Save_coords == 1; 
    save_coords = zeros(1,3); 
    for n = 1:36
        save_coords = vertcat(save_coords,(mapped(:,:,n)'));
    end
        save_coords = rmmissing(save_coords);
        save_coords(:,1) = save_coords(:,1) + centroid_um(1); 
        save_coords(:,2) = save_coords(:,2) + centroid_um(2); 
        save_coords(:,3) = save_coords(:,3) + centroid_um(3);
else
end

%%Crop out some data (selective points):
selected = zeros(1,3); 

for n = 1:8                     %for phi = 5:25, reduce points by 3x
    theta_reduced = mapped(:,1:6:360/theta_n,n)';
    selected = vertcat(selected, theta_reduced);                          
end

for n = 9:10                    %for phi = 30:50, reduce points by 2x
    theta_reduced = mapped(:,1:4:360/theta_n,n)';
    selected = vertcat(selected,theta_reduced); 
end
for n = 11:25
    theta_reduced = mapped(:,1:2:360/theta_n,n)';
    selected = vertcat(selected,theta_reduced);
    %selected = vertcat(selected,mapped(:,:,n)'); 
end
for n = 26:30
    theta_reduced = mapped(:,1:4:360/theta_n,n)';
    selected = vertcat(selected,theta_reduced); 
end
for n = 31:36                   %for phi = 155:130, reduce points by 3x
    theta_reduced = mapped(:,1:6:360/theta_n,n)';
    selected = vertcat(selected,theta_reduced);
end

selected(:,1) = selected(:,1);%.*xy_res; 
selected(:,2) = selected(:,2);%.*xy_res; 
selected(:,3) = selected(:,3);%.*z_res;


% selected_fig = 'n';%linput('Would you like to see the selected data? Y/N \n','s');
% if selected_fig == 'y' || selected_fig =='Y'
% % Uncommment below if you'd like to see the selected data 
%     figure(); hold on; 
%     scatter3(selected(:,1)+centroid_um(1),selected(:,2)+centroid_um(2),selected(:,3)+centroid_um(3));
% else
%     input('That''s cool. Just press enter and I''ll continue')
% end

xdata = selected(:,1) + centroid_um(1); 
ydata = selected(:,2) + centroid_um(2); 
zdata = selected(:,3) + centroid_um(3); 

Data = horzcat(xdata,ydata,zdata);

k = randperm(length(xdata));
Data_mat = Data(k(1:round(length(xdata)*0.6)),:);
Data_mat = rmmissing(Data_mat);

%code for removing outliers
distance_data = sqrt(sum(Data_mat.^2,2));
indices = distance_data<=(mean(distance_data)+std(distance_data)); %note I use smaller than instead of bigger
Data_mat = Data_mat(indices,:);
% Uncomment below if you'd like to see the extra reduced data
% figure()
% scatter3(Data_mat(:,1), Data_mat(:,2), Data_mat(:,3))
% title(sprintf('even more reduced data:  %d points',(length(Data_mat))))

%%Ellipsoid Fitting:   uses Douglas Rachford formula from Fitellipsoid
%%paper from 2019
%%outputs from fitting function are u:  10x1 vector describing ellipsoid,
%%CF:  cost function, A/b/c = matrix, vector, scalar describing ellipsoid

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

% figure()
% scatter3(Data_mat(:,1),Data_mat(:,2),Data_mat(:,3),'.');
% hold on
% p = patch(isosurface(x,y,z,Ellipsoid,1)); 
% set(p,'FaceColor','g','EdgeColor','none'); 
% alpha(p,0.3)
% view(-70,40);
% xlabel('x axis (um)')
% ylabel('y axis (um)')
% zlabel('z axis (um)')
% % xlim([0 140])
% % ylim([0 120])
% % zlim([0 110])
% %axis equal
% title('Ellipsoid fit')

%%Test ellipsoid fit to model parameters
% figure()
% hold on
% scatter3(Data_mat(:,1),Data_mat(:,2),Data_mat(:,3),'.');
% h = surf(x_new,y_new,z_new);
% alpha(0.3)
% grid on
% xlabel('x (um)');
% ylabel('y (um)');
% zlabel('z (um)');
% view(90,0)
% scatter3(center(1), center(2),center(3),'filled','MarkerFaceColor','m')

% % figure()
% % hold on
% %scatter3(center(1), center(2),center(3),'filled','MarkerFaceColor','m')
% plot_evector1=[center(1)+radii(1)*evecs(1,1),center(2)+radii(1)*evecs(2,1),center(3) + radii(1)*evecs(3,1)];
% plot_evector2=[center(1)+radii(2)*evecs(1,2),center(2)+radii(2)*evecs(2,2),center(3) + radii(2)*evecs(3,2)];
% plot_evector3=[center(1)+radii(3)*evecs(1,3),center(2)+radii(3)*evecs(2,3),center(3) + radii(3)*evecs(3,3)];
%  
% plot3([center(1),plot_evector1(1)],[center(2),plot_evector1(2)],[center(3),plot_evector1(3)],'k-','LineWidth',1.5)
% plot3([center(1),plot_evector2(1)],[center(2),plot_evector2(2)],[center(3),plot_evector2(3)],'r-','LineWidth',1.5)
% plot3([center(1),plot_evector3(1)],[center(2),plot_evector3(2)],[center(3),plot_evector3(3)],'b-','LineWidth',1.5)
% 
% 
% 
% title(sprintf('Ellipsoid constructed from radii and eigenvectors for item %d',n))
%%
figure()
hold on
scatter3(Data_mat(:,1),Data_mat(:,2),Data_mat(:,3),'.');
h = surf(x_new,y_new,z_new);
alpha(0.3)
grid on
xlabel('x (um)');
ylabel('y (um)');
zlabel('z (um)');
view(90,0)
scatter3(center(1), center(2),center(3),'filled','MarkerFaceColor','m')

% figure()
% hold on
%scatter3(center(1), center(2),center(3),'filled','MarkerFaceColor','m')
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
%fprintf('The diameters are: %.2f, %.2f, %.2f\n',diameters(1),diameters(2), diameters(3));

%Tstop = toc(tStart);

%fprintf('The total elapsed time was %.2f minutes\n',Tstop/60)
end
