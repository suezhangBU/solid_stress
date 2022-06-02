
clear all
close all
clc

%locEnMin function changes to allow bead to be deformed outside of the
%undeformed sphere

%To help understand the triangular intersection program downloaded:
%https://www.mathworks.com/matlabcentral/mlc-downloads/downloads/submissions/33073/versions/8/previews/html/TriangleRayIntersection_tutorial.html?access_key=
%https://www.mathworks.com/matlabcentral/fileexchange/33073-triangle-ray-intersection
%http://www.graphics.cornell.edu/pubs/1997/MT97.pdf

%% LOAD ABAQUS SPHERE POINTS BEFORE RUNNING THE PROGRAM
load('surfaceABAQUS.mat');

%% STEP 1: Hand crop and Compile Z stack images as 3D Matrix ============================================================================
%Z stack images of Undeformed Bead

% [filename1, path1] = uigetfile('*.*','Multiselect','on');
% if iscell(filename1)
%     fname = filename1;
% elseif filename1 ~= 0
%     fname{1} = filename1;   
% else
%     disp('No files selected')
%     return
% end
% 
% %path1='U:\eng_research_nialab\users\suezhang\2019 Bead Deformation Codes-20210829T160056Z-001\2019 Bead Deformation Codes\SZ0004\trypsin\';
% %filename1='c4_trypsin_z108_c001.png';
% segment1_1=filename1(1:(strfind(filename1, '_z'))+1); %'c4_trypsin_z';
% segment2_1=filename1(strfind(filename1, '_c'):end); %'_c001.png';
% folderloc1= path1(1:end-1); %'U:\eng_research_nialab\users\suezhang\2019 Bead Deformation Codes-20210829T160056Z-001\2019 Bead Deformation Codes\SZ0004\trypsin';
% image_start1=25;
% image_end1=137;
% disp('Cropping Image of UNdeformed Bead');
% userInput1=input('Enter "1" To crop new image z stack of UNDEFORMED BEAD/ "0" to use previously stored z stack (0/1?) ');
% UNdeformedZstack= cropZ(userInput1, path1, filename1,folderloc1,segment1_1,segment2_1,image_start1,image_end1);
% zslice_UN=input('(UNdeformed Bead) What is the size of the z-slice in um/pixel?       ')
% xslice_UN=input('(UNdeformed Bead) What is the size in x,y of the image in um/pixel?      ')
% clc

[filename2, path2] = uigetfile('*.*','Multiselect','on');
if iscell(filename2)
    fname2 = filename2;
elseif filename2 ~= 0
    fname2{1} = filename2;   
else
    disp('No files selected')
    return
end

%Z stack images of Deformed Bead
close all
%path2='U:\eng_research_nialab\users\suezhang\2019 Bead Deformation Codes-20210829T160056Z-001\2019 Bead Deformation Codes\SZ0004\undeformed s0004\';
%filename2='c4before_z078_c001.png';
segment1_2=filename2(1:(strfind(filename2, '_z'))+1);
segment2_2=filename2(strfind(filename2, '_c'):end);;
folderloc2= path2(1:end-1);%'U:\eng_research_nialab\users\suezhang\SZ056\8-23-2021\076\ZSeries-08232021-1-000_tifs';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

image_start2=input('Image start      ');
image_end2=input('Image end      '); %changes these values 
% image_mid=round(mean([image_start2,image_end2]));

surface_value = 90; %the threshhold value of the image, determined from ImageJ analysis
sn=input('Smoothing factor (needs to be greater than 1 and odd)      ');; % nn needs to be odd

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% disp('Now cropping Image of Deformed Bead');
userInput2=1; %=input('Enter "1" To crop new image z stack of DEFORMED BEAD/ "0" to use previously stored z stack (0/1?) ');
deformedZstack= cropZ(userInput2, path2, filename2, folderloc2,segment1_2, segment2_2,image_start2,image_end2);
zslice=input('(Deformed Bead) What is the size of the z-slice in um/pixel?      ')
xslice=input('(Deformed Bead) What is the size in x,y of the image in um/pixel?     ')
clc

disp('Press any key to continue')
pause

%% STEP 2: Get the Deformed Bead =========================================================================================================
%Isosurface of Deformed Bead
SM= smooth3(deformedZstack, 'box', sn); %Smooth raw data in 3D 

%FIGURE 1: (a) Patched Isosurface (b) Plot3 of Isosurface Vertices
%Isosurface and patching
figure(10)
subplot(1,3,1)
title('DEFORMED BEAD (UNSCALED) [Units in Pixel]')
hiso= patch(isosurface(deformedZstack,surface_value), 'FaceColor', [1, .75, .65], 'EdgeColor', 'none');
view(35,30)
axis tight
daspect([1,1,.4])
lightangle(45,30);
lighting gouraud
hcap.AmbientStrength= 0.6;
hiso.SpecularColorReflectance =0;
hiso.SpecularExponent= 50;

hX= hiso.Vertices(:,1);
hY= hiso.Vertices(:,2);
hZ= hiso.Vertices(:,3);


subplot(1,3,2)
title('DEFORMED BEAD (UNSCALED) [Units in Pixel]')
hiso= patch(isosurface(SM,surface_value), 'FaceColor', [1, .75, .65], 'EdgeColor', 'none');
view(35,30)
axis tight
daspect([1,1,.4])
lightangle(45,30);
lighting gouraud
hcap.AmbientStrength= 0.6;
hiso.SpecularColorReflectance =0;
hiso.SpecularExponent= 50;

hX= hiso.Vertices(:,1);
hY= hiso.Vertices(:,2);
hZ= hiso.Vertices(:,3);

subplot(1,3,3)
plot3(hX, hY, hZ, '.');
title('DEFORMED BEAD (UNSCALED) [Units in Pixel]')
hold on

disp('Press Any Key to Continue')
pause
clc

% figure; imagesc(SM(:,:,24))

%% STEP3A: Preparing Conversion to Abaqus Factor====================================================================

%Determine the centroid of the Deformed Bead
binaryImage = SM > 0.001;
s=regionprops(binaryImage);
centroidtest = s.Centroid;

%Determine the Radius of Undeformed Bead
% binaryImage_UN = SM_Undeformed > 0.001;
% s_UN=regionprops3(binaryImage_UN, 'Centroid', 'PrincipalAxisLength');
% diameters = mean(s_UN.PrincipalAxisLength,2);
% disp('The radius of the Undeformed Bead is:')
% radiusUN=diameters/2*xslice
radiusUN=input('Radius of undeformed bead in um ')
abaqus_scale_factor=10/radiusUN

% abaqus_scale_factor=input('What is the Abaqus Scaling Factor? (e.g if Undeformed Bead radius is 30, input 1/3)  ');

%% STEP3B: Scale the Deformed Bead's vertices and centroid values by Abaqus Scaling Factor
%Scale the Deformed Beads Vertices
isoinfo = isosurface(SM,surface_value);
faces= isoinfo.faces;
vertices= isoinfo.vertices;
abaqus_inverse_factor=1/abaqus_scale_factor;
vertX_scaled= vertices(:,1)*xslice*abaqus_scale_factor;
vertY_scaled= vertices(:,2)*xslice*abaqus_scale_factor;
vertZ_scaled= vertices(:,3)*zslice*abaqus_scale_factor;

vertices(:,1)=vertX_scaled;
vertices(:,2)=vertY_scaled;
vertices(:,3)=vertZ_scaled;

binaryImage = SM > 0.001;
s=regionprops(binaryImage);
centroid = s.Centroid;
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
disp('Press Any Key to Continue')
pause

%% STEP4: SPHERE (Undeformed Bead) Shift and Energy of Deformation

%Determine the boundary of the deformed bead
xmin = min(vertices(:,1));
xmax = max(vertices(:,1));
ymin = min(vertices(:,2));
ymax = max(vertices(:,2));
zmin = min(vertices(:,3));
zmax = max(vertices(:,3));

vert1 = vertices(faces(:,1),:);
vert2 = vertices(faces(:,2),:);
vert3 = vertices(faces(:,3),:);
% 
% 
% tempv=[vert1(1,:); vert2(1,:); vert3(1,:); vert1(1,:)]
% 
% figure; plot3(tempv(:,1), tempv(:,2), tempv(:,3))
% clc
format long 
OrderCalc= 1;% input('Would you like to Calculate the Minimum Energy? [1/Y] [0/N]');
if OrderCalc
    energymat= locEnMin(1, OrderCalc, node_coord_surf_M, xmin, xmax, ymin, ymax, zmin, zmax, centroid, vert1, vert2, vert3, vertices);
    energymat_sorted=sortrows(energymat, 4)
    disp('Minimum Energy Computation: SUCCESSFUL')
    SecondOrderCalc= input('Would you like to recalculate the Minimum Energy with a finer grid size? [1/Y] [0/N]');
    if SecondOrderCalc
        energymat= locEnMin(2, SecondOrderCalc, node_coord_surf_M, xmin, xmax, ymin, ymax, zmin, zmax, centroid, vert1, vert2, vert3, vertices);
         energymat_sorted=sortrows(energymat, 4)
        disp('Minimum Energy Computation: SUCCESSFUL')
        ThirdOrderCalc= input('Would you like to recalculate the Minimum Energy with an even finer grid size? [1/Y] [0/N]');
        if ThirdOrderCalc
            energymat=locEnMin(3, ThirdOrderCalc, node_coord_surf_M, xmin, xmax, ymin, ymax, zmin, zmax, centroid, vert1, vert2, vert3, vertices);
            energymat_sorted=sortrows(energymat, 4)
            disp('Minimum Energy Computation: SUCCESSFUL')
            disp('Computation of Minimum Energy SUCCESSFUL')
        else
            disp('Proceeding to next computational steps')
        end
    else
        disp('Proceeding to next computational steps')
    end
else
    disp('Computation not executed')
end

%Calculate the energy when sphere is located minimum energy point
energyonly=energymat(:,4);
minEnergy=min(energyonly);
result=find(energyonly==minEnergy);
energymin_sphere_loc=energymat(result, :);
minCentX=energymin_sphere_loc(1); %coordinates for where the center sphere is after energy minimization 
minCentY=energymin_sphere_loc(2);
minCentZ=energymin_sphere_loc(3);

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

ABAQdeform= [node_coord_surf_M(:,1) ABAQint]; %ABAQdeform= intersection points on deformed bead

%Get the difference between the original Abaqus and Deformed Bead points:
xDeform= node_coord_surf_M(:,2)-ABAQdeform(:,2); %xDeform= original abaqus sphere point - deformed bead position
yDeform= node_coord_surf_M(:,3)-ABAQdeform(:,3);
zDeform= node_coord_surf_M(:,4)-ABAQdeform(:,4);

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
finalDeform= [node_coord_surf_M(:,1) xDeform yDeform zDeform]; %Distance between Abaqus and Deformed Bead
save('ABAQUSintersection.mat', 'ABAQdeform', 'finalDeform') % NAME OF FILE NEEDED FOR "MAKE BOUNDARY MICROGEL V1"

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
disp('ABAQdeform and finalDeform have been saved')

%ABAQdeform is the xyz coordinates of the deformed bead position 
%finalDeform is the displacement from Undeformed-->Deformed Bead needed for the Input file in ABAQUS

%% FUNCTIONS ==========================================================================================================================

function energyMinimizer= locEnMin(energystep, CalculateOrder, node_coord_surf_M, xmin, xmax, ymin, ymax, zmin, zmax, centroid, vert1, vert2, vert3, vertices)
if CalculateOrder
    counter = 0;
    energymat = []; %This is the matrix that stores the energy of deformation for every shift in the position of the sphere respect to the deformed bead
    
    if energystep==1

        stepXmax = xmin + 10.1; %radius of abaqus bead is 10, so shift by 10 in each direction 
        stepXmin = xmax - 9.9;

        stepYmax = ymin + 10.1;
        stepYmin = ymax - 9.9;

        stepZmax = zmax;
        stepZmin = zmin;
    else 
        stepXmin= input('What is the min X location for the sphere?     ');
        stepXmax= input('What is the max X location for the sphere?     ');
        stepYmin= input('What is the min Y location for the sphere?     ');
        stepYmax= input('What is the max Y location for the sphere?     ');
        stepZmin= input('What is the min Z location for the sphere?     ');
        stepZmax= input('What is the max Z location for the sphere?     ');
   
    end 
    step=input('What step size for shifting the sphere?     ');
    
 
    for xi = stepXmin:step:stepXmax
        for yi = stepYmin:step:stepYmax
            for zi = stepZmin:step:stepZmax

                xshift = xi; yshift = yi; zshift =zi;

                xSphere=node_coord_surf_M(:,2);
                ySphere=node_coord_surf_M(:,3);
                zSphere=node_coord_surf_M(:,4);

                xSphere= xSphere+xshift;
                ySphere= ySphere+yshift;
                zSphere= zSphere+zshift;

                %radi=30;
                %xSphere= xSphere*radi+xshift;
                %ySphere= ySphere*radi+yshift;
                %zSphere= zSphere*radi+zshift;

                %% STEP5: Radiate Line from CENTROID --> SPHERE (undeformed bead) --> DEFORMED BEAD and find the intersection point
                sz= size(xSphere(:));
                xSphere_el_t= xSphere';
                xSphere_el= reshape(xSphere_el_t, [1, sz(1)]); %xSphere_el

                ySphere_el_t= ySphere';
                ySphere_el= reshape(ySphere_el_t, [1, sz(1)]); %ySphere_el

                zSphere_el_t= zSphere';
                zSphere_el= reshape(zSphere_el_t, [1, sz(1)]); %zSphere_el
                Sphere1D = [xSphere_el' ySphere_el' zSphere_el'];

                %Find the boundary of the Sphere. These boundary values will help determine
                %whether the deformed bead is inside the sphere.
                xspheremin = min(xSphere_el);
                xspheremax = max(xSphere_el);
                yspheremin = min(ySphere_el);
                yspheremax = max(ySphere_el);
                zspheremin = min(zSphere_el);
                zspheremax = max(zSphere_el);
                
                

              

                    counter = counter + 1;

                    orig=centroid;


                    intersectstore = zeros(sz(1),3);
                    for i=1:length(xSphere_el)
                        dir= [xSphere_el(i)-centroid(1), ySphere_el(i)-centroid(2), zSphere_el(i)-centroid(3)];
                        dir = dir ./ norm(dir);

                        %% -------------------------------------------------------------------

                        [INTERSECT, T, U, V, XCOOR] = TriangleRayIntersection(orig, dir, vert1, vert2, vert3, 'lineType' , 'ray');
                        if max(INTERSECT)~= 1
                            continue
                        end
                        locintersect = find(INTERSECT==1); %Determines the index value of the face that has the xcoordinate of the intersection
                        corintersect = XCOOR(locintersect(1),:);  %x,y,z coordinate value of intersection
                        intersectstore(i,:) = corintersect;
                    
                    end 
            

                    %% STEP 6: Calculate the Energy of the Deformation
                    energytot = stepsix(intersectstore, Sphere1D, xspheremin, xspheremax, yspheremin, yspheremax, zspheremin, zspheremax, vertices);
                    
                        
                    disp(energytot)

                    energymat(counter,1) = xi;
                    energymat(counter,2) = yi;
                    energymat(counter,3) = zi;
                    energymat(counter,4) = energytot;                 
                    
                end
                end
    end
    
     energyMinimizer= energymat;
    end
    
   

    end 




function zStack= cropZ(userInput,path, filename_largest_image, folderloc, segment1, segment2,is,ie)
%Truncated=input('Enter "1" To crop new image z stack of DEFORMED BEAD/ "0" to use previously stored z stack (0/1?) ');
if userInput
   
    %First file that is loaded should be image where the bead has the largest area
    %path='C:\Users\Sunny\Music\ZSTACK_3dmodel\';
    %filename_largest_image='3dmodelimages_z155_c002.png';
    I=imread([path filename_largest_image]);
    
    %Location of all z stack images
    %folderloc= 'U:\eng_research_nialab\users\suezhang\2019 Bead Deformation Codes-20210829T160056Z-001\2019 Bead Deformation Codes\SZ0004\undeformed s0004';
    
    %===============================CROPPING=================================
    %Select Largest Region that contains the bead of interest
    %Use left click to draw region of interest to be cropped./
    workspace;
    imshow(I,[]);
    axis on;
    title ('Original Confocal Image', 'FontSize', 15);
    set(gcf, 'Position', get(0, 'Screensize')); %Maximize figure
    
    %Ask user to draw freehand mask.
    message = sprintf('Left click and hold to begin drawing.\nSimply lift the mouse button to finish');
    uiwait(msgbox(message));
    hFH = drawfreehand(); % Actual line of code to do the drawing.
    % Create a binary image ("mask") from the ROI object.
    binaryImage = hFH.createMask();
    xy = hFH.Position;
    Xmin=int16(round(min(xy(:,1))));
    Xmax=int16(round(max(xy(:,1))));
    Ymin=int16(round(min(xy(:,2))));
    Ymax=int16(round(max(xy(:,2))));
    
    subplot(1, 3, 1);
    imshow(I, []);
    axis on;
    drawnow;
    title('Original Confocal Image', 'FontSize', 15);
    
    % Display the freehand mask.
    subplot(1, 3, 2);
    imshow(binaryImage);
    axis on;
    title('Binary mask of the region', 'FontSize', 15);
    
    % Mask the image outside the mask, and display it.
    % Will keep only the part of the image that's inside the mask, zero outside mask.
    blackMaskedImage = I;
    blackMaskedImage(~binaryImage) = 0;
    subplot(1, 3, 3);
    imshow(blackMaskedImage);
    axis on;
    title('Masked black outside region', 'FontSize', 15);
    clc
    
    %====================== 3D MATRIX COMPILATION ===========================
%     is=38; %index of the first file you want to pick
%     ie=95; %index of the last file you want to pick
    
    sztest=I;
    sztest(~binaryImage)=0;
    sztest=sztest(Ymin:Ymax, Xmin:Xmax);
    
    rawdata_M=zeros(size(sztest,1), size(sztest,2), ie-is+1);
     
    counter= 0;
    for ii=is:ie %size(II_M,3)
        counter= counter +1;
        %filename_segment=['c4before_z' num2str(ii) '_c001.png'];
        filename_segment=[segment1 num2str(ii,'%03.f') segment2];
        imagename_segment=[folderloc '\' filename_segment];
        blackI= imread(imagename_segment);
        
        %Crop image in ROI specified by hand drawing in first file
        blackI(~binaryImage)=0;
        blackIcrop= blackI(Ymin:Ymax, Xmin:Xmax);
        
        rawdata_M(:,:,counter)= blackIcrop;
    end;
     
    save('cropresized.mat', 'rawdata_M');
    zStack= rawdata_M;
    
    disp('Cropping and Compilation of Bead Images: SUCCESSFUL')
%     disp('Press any key to continue')
% %     pause
    clc

else
    clc
end
end

function energytot =stepsix(intersectstore, Sphere1D, xspheremin, xspheremax, yspheremin, yspheremax, zspheremin, zspheremax, vertices);

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
                     

    
for ii = 1:length(usphere)
    dist(ii) = norm(usphere(ii,:) - ucorintersect(ii,:));
    k=1; %To be determined by AFM
    distsqr = dist.^2;
    energydeform = distsqr * k * 0.5;
    energytot = sum(energydeform);
end

end

%Update Time Log:
%9/12/2019 6:04PM

%TASKS:
%- Overlap z slice from confocal to heat map ABAQUS z- slice