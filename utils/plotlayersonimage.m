function plotlayersonimage(imLayer,labelLayer,im,ysrf,ybtm, lineColor,holdon,flip)

    if isempty(labelLayer)
        labelLayer = 1:max(imLayer(:));
    end

    if nargin >=8
        if strcmp(flip,'flip')
        imLayer = fliplr(imLayer);
        im = fliplr(im);
        else
        fprintf('%s is not accepted',flip);
        end
    end
     
    if nargin < 7
        holdon = 0;
    else
        if strcmp(holdon,'hold on')
            holdon = 1;
        end
    end 
    
    
    if nargin < 6 || strcmp(lineColor,'rand')
        mclrflag =1;
    else
        if isnumeric(lineColor) && (sum(size(lineColor)==[1 3])==2 || ...
      sum(size(lineColor)==[3 1])==2) && sum((lineColor<=[1 1 1] & ...
      lineColor>=[0 0 0]))==3 || sum(strcmpi({'y','m','c','r','g','b','w','k',...
      'yellow','magenta','cyan','red','green','blue','white','black'},lineColor))>0
            mclrflag = 0;
        else
            return;
        end
    end
            
    if holdon == 0
        % Create figure
        figure;
%         figure1 = figure('InvertHardcopy','off','Color',[1 1 1]);
%         set(gcf, 'PaperUnits', 'centimeters');
%         set(gcf, 'PaperSize', [10 6]);
%         % % Create axes
%         axes1 = axes('Parent',figure1,...
%             'Position',[0.13 0.215605749486653 0.643504273504274 0.709394250513347]);
% %         hold(axes1,'on'); axis ij;
        if ~isempty(im), imagesc(im); end
        hold on; axis on;axis ij;
%       h = colorbar; 
%         ylabel(h,'Radar Amplitude (dB)');
%         set(get(h,'YLabel'),'Rotation',-90);
%       set(gca,'XTick',[],'YTick',[]);
    else
        fig1 = gcf;
        hold on;
    end

    if ~isempty(ysrf) && ~isempty(ybtm)

    plot(1:length(ysrf),ysrf,'-k','linewidth',3);
    plot(1:length(ybtm),ybtm,'-k','linewidth',3);

    end

    if mclrflag == 0
        [y,x] = find(imLayer > 0);        
        plot(x,y,'.','color',lineColor);
        
    else
        if mclrflag == 1
            cmap = colormap('jet');
            step = (max(labelLayer) - min(labelLayer))/size(cmap,1);
            for i = 1:length(labelLayer)
                ncolor = floor((labelLayer(i) - min(labelLayer))...
                /step) ;
            ncolor = size(cmap,1) - ncolor ;
            if ncolor > size(cmap,1) 
                ncolor = size(cmap,1);
            else
                if ncolor < 1
                    ncolor = 1;
                end
            end
                
                [y,x] = find(imLayer == labelLayer(i));
                xline = nan * size(imLayer,2);
                yline = nan * size(imLayer,2);
                xline(x) = x;
                yline(x) = y;
                plot(xline,yline,'.','color',cmap(ncolor,:),'markersize',3);
            end
        end
    end
    ylim([0 1100]); xlim([0 size(imLayer,2)]);
    plot([100 100],[30 100],'k-','linewidth',3); % 200 m
    plot([100 478],[30 30],'k-','linewidth',3); % 500 m
%     axis on; set(gca,'box','on','XTick',[],'YTick',[]);
    set(findall(gcf,'-property','FontSize'),'FontSize',14);
end
