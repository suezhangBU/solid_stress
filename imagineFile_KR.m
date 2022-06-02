% --Read MATLAB file and display scaled image
function Image_fig = imagineFile(data,absolute_points,plotting_vecs,discretized_all,c_values)%,plotting_axes)

%  if nargin == 0
%      [fname,pname] = uigetfile('*.mat');
%      data = struct2cell(load([pname,fname]));
%      data = data{1};
%      str_in = fname;
%  elseif ~isnumeric(data)
%     %  Load file
%     fname = data;
%     data = struct2cell(load(fname));
%     data = data{1};
%     str_in = strsplit(fname,'/');
%     str_in = str_in{end};
%  else
     str_in = 'Image Stack';
% end
 
 %if nargin < 2
     res = [1 1];
 %end
 
[r,c,f] = size(data);
x = linspace(0,c*res(1),c);
y = linspace(0,r*res(2),r);


    % Create a figure and axes
    fH = figure; 
    imagesc(x,y,data(:,:,1)), hold on
    axis image
    colormap(gray), colorbar
    xlabel('Position in \mu m')
    ylabel('Position in \mu m')
    title(str_in)
    set(gca,'Units','normalized');
    
    
    % Create pop-up menu
%     popup_contounr = uicontrol(fH,'Style', 'popup','Units','normalized',...
%         'String', title_str,...
%         'Position', [0.1 0.90 0.1 0.05],...
%         'Callback', @ContLine); 
%     popup_bckgrnd = uicontrol(fH,'Style', 'popup','Units','normalized',...
%         'String', {'Original','Thresholded','Dual D-E'},...
%         'Position', [0.25 0.90 0.1 0.05],...
%         'Callback', @Bckgrnd); 
    
   % Create push button
%     btn = uicontrol('Style', 'pushbutton', 'String', 'Clear',...
%         'Position', [0.1 0.15 0.1 0.2],...
%         'Callback', 'cla');       

   % Create slider
    sld = uicontrol(fH,'Style', 'slider','Units','normalized',...
        'Min',1,'Max',f,'Value',1,...
        'Position', [0.35 0.15 0.3 0.05],...
        'Callback', @StackSlide); 
					
    % Add a text uicontrol to label the slider.
    txt1 = uicontrol(fH,'Style','text','Units','normalized',...
        'Position',[0.50 0.1 0.075 0.05],...
        'String','Slice #: ');
    txt2 = uicontrol('Style','text','Units','normalized',...
        'Position',[0.575 0.1 0.035 0.05],...
        'String',num2str(1));
    
    % Make figure visble after adding all components
%     f.Visible = 'on';
    % This code uses dot notation to set properties. 
    % Dot notation runs in R2014b and later.
    % For R2014a and earlier: set(f,'Visible','on');
    
%     function ContLine(source,event)
%         val = source.Value;
%         str = title_str{floor(val)};
%         % For R2014a and earlier: 
%         % val = get(source,'Value');
%         % maps = get(source,'String'); 
%         lower_cont = y_l_c{floor(val)};
%         imagineFile(stack,lower_cont,str)
%     end
%     function Bckgrnd(source,event)
%         val = floor(source.Value);
%         if val == 1
%             imagineEndPlate(c,slice_array,lo_cont,str_in)
%         elseif val == 2
%             imagineEndPlate(c,slice_array_de,lo_cont,str_in)
%         elseif val == 3
%             imagineEndPlate(c,slice_array_inv_de,lo_cont,str_in)
%         end
%     end

    Image_fig = fH; 
    
    function StackSlide(source,~)
        val = floor(source.Value);
        % For R2014a and earlier:
        % val = get(source,'Value');
        txt2.String = num2str(val);
        fH;imagesc(x,y,data(:,:,val)); 
        any_points = find(val == absolute_points(:,3));
%         axes_points = find(val == round(plotting_vecs(:,3)));
        disc_points = find(val == discretized_all);
        disc_points_col = ceil((disc_points/size(discretized_all,1)));
        disc_points_row = round((rem(disc_points/size(discretized_all,1),1))/(1/size(discretized_all,1)),1);
        hold on
        if length(any_points) == 0;
        else
            for n = 1:length(any_points)
                scatter(absolute_points(any_points(n),2),absolute_points(any_points(n),1),'co','filled')
            end
        end
        if length(disc_points) == 0;
        else
            for t = 1:length(disc_points)
                if disc_points_row(t) == 0
                    disc_points_row(t) = 18; 
                else
                end
                scatter(discretized_all(disc_points_row(t)-1,disc_points_col(t)),discretized_all(disc_points_row(t)-2,disc_points_col(t)),7,'MarkerEdgeColor','none','MarkerFaceColor',c_values{disc_points_row(t)});
            end
        end
       
    end
 end