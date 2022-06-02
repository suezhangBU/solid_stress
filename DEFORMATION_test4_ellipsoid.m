% 
clear all
close all
% clc

%locEnMin function changes to allow bead to be deformed outside of the
%undeformed sphere

%To help understand the triangular intersection program downloaded:
%https://www.mathworks.com/matlabcentral/mlc-downloads/downloads/submissions/33073/versions/8/previews/html/TriangleRayIntersection_tutorial.html?access_key=
%https://www.mathworks.com/matlabcentral/fileexchange/33073-triangle-ray-intersection
%http://www.graphics.cornell.edu/pubs/1997/MT97.pdf

%% LOAD ABAQUS SPHERE POINTS BEFORE RUNNING THE PROGRAM
load('surfaceABAQUS.mat');

% %load main .mat file
% [filename_main,path_main] = uigetfile('*.mat','Select your file');
% if isequal(filename_main,0)
%     disp('User selected cancel. Goodbye!')
%     return
% else
%     disp(['OK, I''ll analyze ',fullfile(path_main,filename_main)])
% end
% 
% load([path_main,filename_main]); %load name of file with 'ABAQdeform', 'finalDeform'
%load individual bead file
[filename_ellipsoid,path] = uigetfile('*.mat','Select your file');
if isequal(filename_ellipsoid,0)
    disp('User selected cancel. Goodbye!')
    return
else
    disp(['OK, I''ll analyze ',fullfile(path,filename_ellipsoid)])
end

load([path,filename_ellipsoid]); %load name of file with 'ABAQdeform', 'finalDeform'

%% STEP 1: Hand crop and Compile Z stack images as 3D Matrix ============================================================================


% zslice=input('(Deformed Bead) What is the size of the z-slice in um/pixel?      ')
% xslice=input('(Deformed Bead) What is the size in x,y of the image in um/pixel?     ')

zslice=1;
xslice=1;
%% STEP 2: Get the Deformed Bead =========================================================================================================


%% STEP3A: Preparing Conversion to Abaqus Factor====================================================================

radiusUN=input('Radius of undeformed bead in um ')
abaqus_scale_factor=10/radiusUN

% abaqus_scale_factor=input('What is the Abaqus Scaling Factor? (e.g if Undeformed Bead radius is 30, input 1/3)  ');

%% STEP3B: Scale the Deformed Bead's vertices and centroid values by Abaqus Scaling Factor
%Scale the Deformed Beads Vertices
% isoinfo = isosurface(SM,surface_value);

isoinfo=ellipsoid_surface;
faces= isoinfo.faces;
vertices= isoinfo.vertices;
abaqus_inverse_factor=1/abaqus_scale_factor;
vertX_scaled= vertices(:,1)*xslice*abaqus_scale_factor;
vertY_scaled= vertices(:,2)*xslice*abaqus_scale_factor;
vertZ_scaled= vertices(:,3)*zslice*abaqus_scale_factor;

vertices(:,1)=vertX_scaled;
vertices(:,2)=vertY_scaled;
vertices(:,3)=vertZ_scaled;

centroid=center;
centroidXscaled= centroid(1)*xslice*abaqus_scale_factor;
centroidYscaled= centroid(2)*xslice*abaqus_scale_factor;
centroidZscaled= centroid(3)*zslice*abaqus_scale_factor;
centroid= [centroidXscaled centroidYscaled centroidZscaled];

%Ensure Abaqus Scaling is done Correctly
figure(2)
plot3(node_coord_surf_M(:,2)+centroid(1),node_coord_surf_M(:,3)+centroid(2),node_coord_surf_M(:,4)+centroid(3),'ro'); % node_coord_surf_M is the coordinates of the nodes of the sphere with radius = 10. This has been imported from ABAQUS.
axis equal
hold on
% plot3(verticesUN(:,1), verticesUN(:,2), verticesUN(:,3), 'b.');
title('If ABAQUS scaling was done correctly, the diameters of the red and blue should match up'); %this was from before when the undeformed bead was loaded in, but don't need it 

figure(3)
plot3(node_coord_surf_M(:,2)+centroid(1),node_coord_surf_M(:,3)+centroid(2),node_coord_surf_M(:,4)+centroid(3),'ro');
hold on
plot3(vertices(:,1), vertices(:,2), vertices(:,3), 'c.');
axis equal
title('Scaled Deformed Bead inside Abaqus Sphere')


% figure; 
% hold on;
% plot3(isoinfo.vertices(:,1),isoinfo.vertices(:,2),isoinfo.vertices(:,3),'b.')
% 
% disp('Press Any Key to Continue')
% pause

%% STEP4: SPHERE (Undeformed Bead) Shift and Energy of Deformation

% %Determine the boundary of the deformed bead
xmin = min(vertices(:,1));
xmax = max(vertices(:,1));
ymin = min(vertices(:,2));
ymax = max(vertices(:,2));
zmin = min(vertices(:,3));
zmax = max(vertices(:,3));

vert1 = vertices(faces(:,1),:);
vert2 = vertices(faces(:,2),:);
vert3 = vertices(faces(:,3),:);


%Calculate the energy when sphere is located minimum energy point
minCentX=centroid(1); %coordinates for where the center sphere is after energy minimization 
minCentY=centroid(2);
minCentZ=centroid(3);

xSphere= node_coord_surf_M(:,2);
ySphere= node_coord_surf_M(:,3);
zSphere= node_coord_surf_M(:,4);

xSphereMIN= xSphere+minCentX;
ySphereMIN= ySphere+minCentY;
zSphereMIN= zSphere+minCentZ;

szMin= size(xSphere(:));
xSphere_el_tM= xSphereMIN';
xSphere_elM= reshape(xSphere_el_tM, [1, szMin(1)]); %xSphere_el
            
ySphere_el_tM= ySphereMIN';
ySphere_elM= reshape(ySphere_el_tM, [1, szMin(1)]); %ySphere_el
            
zSphere_el_tM= zSphereMIN';
zSphere_elM= reshape(zSphere_el_tM, [1, szMin(1)]); %zSphere_el
Sphere1DMIN = [xSphere_elM' ySphere_elM' zSphere_elM'];

intersectstoreMIN = zeros(szMin(1),3);


orig= centroid;
for i=1:length(xSphere_elM)
    dir= [xSphere_elM(i)-centroid(1), ySphere_elM(i)-centroid(2), zSphere_elM(i)-centroid(3)];
    dir = dir ./ norm(dir);
    %% -------------------------------------------------------------------

    [INTERSECT, T, U, V, XCOOR] = TriangleRayIntersection(orig, dir, vert1, vert2, vert3, 'lineType' , 'ray');
    if max(INTERSECT)~= 1
        continue
    end
    locintersectMIN = find(INTERSECT==1); %Determines the index value of the face that has the xcoordinate of the intersection
    corintersectMIN = XCOOR(locintersectMIN(1),:);  %x,y,z coordinate value of intersection
    intersectstoreMIN(i,:) = corintersectMIN;

end

clc

ABAQint=intersectstoreMIN;

ABAQint(:,1)= ABAQint(:,1)-minCentX; % Relocate the intersection points so that the ABAQUS sphere is back to (0,0,0)
ABAQint(:,2)= ABAQint(:,2)-minCentY;
ABAQint(:,3)= ABAQint(:,3)-minCentZ;

figure(4)
plot3(xSphere, ySphere, zSphere, 'r.')
hold on
plot3(ABAQint(:,1),ABAQint(:,2), ABAQint(:,3), 'b.');
title('Abaqus Sphere at Minimum Energy');
hold on
plot3(node_coord_surf_M(19,2), node_coord_surf_M(19,3), node_coord_surf_M(19,4), 'yo') %point on sphere
hold on
plot3(ABAQint(19,1),ABAQint(19,2),ABAQint(19,3),'go')
axis equal 

ABAQdeform= [node_coord_surf_M(:,1) ABAQint]; %ABAQdeform= intersection points on deformed bead

%Get the difference between the original Abaqus and Deformed Bead points:
xDeform= node_coord_surf_M(:,2)-ABAQdeform(:,2); %xDeform= original abaqus sphere point - deformed bead position
yDeform= node_coord_surf_M(:,3)-ABAQdeform(:,3);
zDeform= node_coord_surf_M(:,4)-ABAQdeform(:,4);
%%
energytot =energymin(intersectstore, Sphere1D, xspheremin, xspheremax, yspheremin, yspheremax, zspheremin, zspheremax, vertices);
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
finalDeform= [node_coord_surf_M(:,1) xDeform yDeform zDeform]; %Distance between Abaqus and Deformed Bead
date=datestr(now, 'yyyy_mm_dd');
filename_ellipsoid=erase(filename_ellipsoid,'.mat')
save([date,filename_ellipsoid, '_ABAQUSintersection.mat'], 'ABAQdeform', 'finalDeform') % NAME OF FILE NEEDED FOR "MAKE BOUNDARY MICROGEL V1"

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
disp('ABAQdeform and finalDeform have been saved')

%ABAQdeform is the xyz coordinates of the deformed bead position 
%finalDeform is the displacement from Undeformed-->Deformed Bead needed for the Input file in ABAQUS

%% FUNCTIONS ==========================================================================================================================


function energytot =energymin(intersectstore, Sphere1D, xspheremin, xspheremax, yspheremin, yspheremax, zspheremin, zspheremax, vertices);

ucorintersect = unique(intersectstore,'rows','stable'); %Delete duplicate coordinates
usphere = unique(Sphere1D,'rows','stable'); %Delete duplicate coordinates
dist = [];
% lengthUsphere = length(usphere);
% lengthUcorinterset = length(ucorintersect);

%             inside_vertices=[];
%                 for i=1:length(vertices(:,1))
%                    if  xspheremin < vertices(i,1) && yspheremin < vertices(i,2) && zspheremin < vertices(i,3) && xspheremax >= vertices(i,1) && yspheremax >= vertices(i,2) && zspheremax >= vertices(i,3)
%                        inside_vertices(i,:)=vertices(i,:);
%                    else
%                        inside_vertices(i,:)=NaN(1,3,'single');
%                    end
%                 end
%              inside_vertices=[(1:length(inside_vertices))' inside_vertices];
%              inside_vertices(any(isnan(inside_vertices),2),:)=[];
%                
%              outside_vertices=[];
%                 for i=1:length(vertices(:,1))
%                    if  vertices(i,1)< xspheremin | vertices(i,2)< yspheremin | vertices(i,3)< zspheremin | vertices(i,1)> xspheremax |  vertices(i,2)> yspheremax | vertices(i,3)> zspheremax
%                        outside_vertices(i,:)=vertices(i,:);
%                    else
%                        outside_vertices(i,:)=NaN(1,3,'single');
%                    end
%                 end
%                 
%                 outside_vertices=[(1:length(outside_vertices))' outside_vertices];
%                  outside_vertices(any(isnan(outside_vertices),2),:)=[];
                     

    
for ii = 1:(min(length(ucorintersect),length(usphere)))
    dist(ii) = norm(usphere(ii,:) - ucorintersect(ii,:));
    k=1; %spring constant of bead To be determined by AFM
    distsqr = dist.^2;
    energydeform = distsqr * k * 0.5;
    energytot = sum(energydeform);
end

end

%Update Time Log:
%9/12/2019 6:04PM

%TASKS:
%- Overlap z slice from confocal to heat map ABAQUS z- slice