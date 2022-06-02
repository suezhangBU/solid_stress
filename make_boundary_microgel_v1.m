clear all
close all
clc
date=char(datetime); 

[filename,path] = uigetfile('*.mat','Select your file');
if isequal(filename,0)
    disp('User selected cancel. Goodbye!')
    return
else
    disp(['OK, I''ll analyze ',fullfile(path,filename)])
end

load([path,filename]); %load name of file with 'ABAQdeform', 'finalDeform'

AA=csvread('set1_nodes.txt');
surf_node=[];
for ii=1:16
    for jj=1:size(AA,1)
        if AA(jj,ii)~=0
            surf_node=[surf_node; AA(jj,ii)];
        end;
    end;
end;

fname='spheretest.inp'; %must be in the folder 
fids=fopen(fname,'r');

for ii=1:9
    temp=fgetl(fids);
end;

temp1 = textscan(fids,'%d %f %f %f','Delimiter',',');

node_all=zeros(length(temp1{1}),4);

node_all(:,1)=temp1{1};
node_all(:,2)=temp1{2};
node_all(:,3)=temp1{3};
node_all(:,4)=temp1{4};

fclose all

fidd=fopen([filename '_deformations.txt'],'wt','n'); % the new file with the 4 first lines omitted. Can rename this file 

node_coord_surf_M=[];

for ii=1:size(finalDeform,1)
    jj=surf_node(ii);
%     zz=node_all(jj,4);
    node_coord_surf_M=[node_coord_surf_M; jj, node_all(jj,2),node_all(jj,3),node_all(jj,4)];
    
%     fwrite(fidd,['*Boundary' char(10)]);
    fprintf(fidd,['Part-1-1.' num2str( finalDeform(ii,1) ) ', 1, 1, ' num2str(-finalDeform(ii,2)) '\n']);
    fprintf(fidd,['Part-1-1.' num2str( finalDeform(ii,1) ) ', 2, 2, ' num2str(-finalDeform(ii,3)) '\n']);
    fprintf(fidd,['Part-1-1.' num2str( finalDeform(ii,1) ) ', 3, 3, ' num2str(-finalDeform(ii,4)) '\n']);
end;

fclose all;

figure(1); 
hold on;
plot3(node_coord_surf_M(:,2),node_coord_surf_M(:,3),node_coord_surf_M(:,4),'r.','MarkerSize',12)
plot3(ABAQdeform(:,2),ABAQdeform(:,3),ABAQdeform(:,4),'k.','MarkerSize',12)
axis equal


figure(1)
hold on;
for ii=1:size(ABAQdeform,1)
    xx=[node_coord_surf_M(ii,2) ABAQdeform(ii,2)];
    yy=[node_coord_surf_M(ii,3) ABAQdeform(ii,3)];
    zz=[node_coord_surf_M(ii,4) ABAQdeform(ii,4)];
    plot3(xx,yy,zz,'c');
end;


%%
save_inp=[filename,'abaqus_cellular.inp']
system(['copy spheretest_part1_cellular.inp + ' filename '_deformations.txt + spheretest_part2.inp ' save_inp]) 

filename=erase(filename,'.mat')
save([filename,'_final.mat'], 'node_coord_surf_M') %name output file here
fclose all;






