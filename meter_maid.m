function rois = meter_maid(image,savepathfolder,save_filename,bead_number)%%
global roiPos roiStruct fROI

image=imadjustn(image);
fROI = imagineFile(image);

% --Create interupt push button
pushBtn = uicontrol('Parent',fROI,'Style','pushbutton',...
    'Units','normalized','Position',[0.025 0.05 0.175 0.05],...
    'String','Confirm Slection','Callback',@confirmROI);

goon = 1; 
n = 1; 

roiPos = [];

rois = {};
while goon >= 1
%     if isempty(roiPos)
        %uiwait(fROI)
        z_min = input('what is the minimum z frame for this bead?\n')
        z_max = input('what is the maximum z frame for this bead?\n')
        z_for_roi = round((z_max + z_min)/2);
        figure(n*100)
        imshow(image(:,:,z_for_roi))
        roiStruct = drawrectangle(gca,'Color','r');
        rois{n} = [{roiStruct.Vertices}, z_min, z_max];
        
        bead_name=['_bead' num2str(bead_number)]
        title([save_filename ' ' bead_name],'Interpreter','none')
        saveas(figure(n*100),[savepathfolder,save_filename,bead_name,'_z',num2str(z_for_roi), '_ROI.fig']);
%         close(n*10)
        %disp(get(var, 'value'));
        %uiwait(fROI)
%     else
%         roiStruct = drawrectangle(gca,'Color','y','Position',roiPos);
%         uiwait(fROI)
%     end

n = n+1; 
bead_number=bead_number+1;

goon = input('would you like another bead (yes - 1, no - 0)\n');    

end


end

function confirmROI(source,eventdata,fROI,n)
global roiPos roiStruct fROI
    %  --Set value of output ROI Position
    roiPos = roiStruct.Vertices;
    % --Reume execution as pushBtn callback has been completed
    uiresume(fROI)
end

