function [rois,z_increment] = meter_maidSZ(image,savepathfolder,save_filename,bead_name)
global roiPos roiStruct fROI

image=imadjustn(image);
fROI = imagineFile(image);

% --Create interupt push button
pushBtn = uicontrol('Parent',fROI,'Style','pushbutton',...
    'Units','normalized','Position',[0.025 0.05 0.175 0.05],...
    'String','Confirm Slection','Callback',@confirmROI);

z_min = input('what is the minimum z frame for this bead?\n')
z_max = input('what is the maximum z frame for this bead?\n')
totalZ=z_max-z_min;

if totalZ>=10
z_increment=round(totalZ/5);
else
z_increment=1;
end
image_slices=totalZ/z_increment;

roiPos = [];

rois = {};

z_for_roi=z_min;
n=1;
while z_for_roi< z_max
%     if isempty(roiPos)
        %uiwait(fROI)
        
        figure(n*10)
        imshow(image(:,:,z_for_roi))
        hold on
        roiStruct = drawfreehand(gca,'Color','r','Multiclick',1,'Closed',1,'Smoothing',1);
        rois{n} = [{roiStruct.Position}, z_min, z_max];
        plot(rois{1,n}{1,1}(:,1),rois{1,n}{1,1}(:,2),'marker','.','markersize',15,'linestyle','none','color','r')
        hold off
        saveas(figure(n*10),[savepathfolder,save_filename,bead_name,'_ImageSlice_',num2str(z_for_roi),'.fig']);
        %disp(get(var, 'value'));
        %uiwait(fROI)
%     else
%         roiStruct = drawrectangle(gca,'Color','y','Position',roiPos);
%         uiwait(fROI)
%     end

z_for_roi=z_for_roi+z_increment;
n=n+1;
end


end

function confirmROI(source,eventdata,fROI,n)
global roiPos roiStruct fROI
    %  --Set value of output ROI Position
    roiPos = roiStruct.Vertices;
    % --Reume execution as pushBtn callback has been completed
    uiresume(fROI)
end

