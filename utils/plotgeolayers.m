function plotgeolayers(geolayers,lineColor,holdon,flip)

    if nargin < 4 || strcmp(flip,'flip')
        geolayers = flipgeolayers(geolayers);
    else
        fprintf('%s is not accepted\n',flip);
    end
     
    if nargin < 3
        holdon = 0;
    else
        if strcmp(holdon,'hold on')
            holdon = 1;
        end
    end 
    
    
    if nargin < 2 || strcmp(lineColor,'rand')
        lineColor = rand(geolayers.num_layer,3);
        mclrflag = 1;
    else
        if isnumeric(lineColor) && (sum(size(lineColor)==[1 3])==2 || ...
      sum(size(lineColor)==[3 1])==2) && sum((lineColor<=[1 1 1] & ...
      lineColor>=[0 0 0]))==3 || sum(strcmpi({'y','m','c','r','g','b','w','k',...
      'yellow','magenta','cyan','red','green','blue','white','black'},lineColor))>0
           mclrflag = 0;
        else
            disp('No input color!\n');
            return;
        end
    end
            
    if holdon == 0
        % Create figure
        figure1 = figure('InvertHardcopy','off','Color',[1 1 1]);
        set(gcf, 'PaperUnits', 'centimeters');
        set(gcf, 'PaperSize', [18 6]);
        % % Create axes
        axes1 = axes('Parent',figure1,...
            'Position',[0.13 0.215605749486653 0.643504273504274 0.709394250513347]);
        hold(axes1,'on'); hold on;
        
    else
        hold on;
    end
    
    %n = 0; n2 = 0;
    % Create multiple lines using matrix input to plot
    x = geolayers.x; y = geolayers.y;
    x = x - x(1); y = y - y(1);
    distance = sqrt(x.^2 + y.^2);
    if distance(end) > 1000
        distance = distance/1000;
    end
    maxdist = max(distance);
    plot(distance,geolayers.elevation_surface/1000,'-','color',[0.5,0,0],...
        'markersize',1.5,'linewidth',3,'displayname','Surface');
    plot(distance,geolayers.elevation_bed/1000,'-','color',[0,0,0.56],...
        'markersize',1.5,'linewidth',3,'displayname','Bed');
    if geolayers.elevation_surface(1) < geolayers.elevation_surface(end)
        legend('Location','northwest');
        legend('Surface','Bed');
        legend('show');
    else
        legend('show');
    end
    
    %% Initial age
    if isfield(geolayers.layer,'age')
        age = extractfield(geolayers.layer,'age');
%           ageuncert = geolayers.layer(:).age_uncert;
        logage = log10(age);
        colorbar;
%         h = colorbar; set(h,'YDir','reverse');
        cmap = colormap('jet');
        caxis([min(logage) max(logage)]);
        stepage = (max(logage) - min(logage))/(size(cmap,1)-1);
        a = unique(age); nconfirm = length(a); disp(nconfirm);
    else
        age = [];
    end
    
    nage = 0; 
    edgeHei = zeros(geolayers.num_layer,1);
    for i = 1:geolayers.num_layer
        ielevation = geolayers.layer(i).elevation/1000;
        
        if isempty(find(~isnan(ielevation),1))
            continue;
        end
        
        
        if ~isempty(age) 

            idx = find(~isnan(fliplr(ielevation)),1); 
            
            iage = geolayers.layer(i).age;
            
            if isempty(iage) || isnan(iage)
                
                wide = length(find(ielevation > 0));
                if wide > 0
                    plot(distance,ielevation,'k-',...
                        'linewidth',1.5);
                end
                continue;
            else
                edgeHei(i) = ielevation(geolayers.num_trace - idx + 1);
                nage = nage + 1;
            end
            
            
            ncolor = ceil((log10(iage) - min(logage))...
                /stepage) ;
            ncolor = size(cmap,1) - ncolor ;
            if ncolor > size(cmap,1) 
                ncolor = size(cmap,1);
            else
                if ncolor < 1
                    ncolor = 1;
                end
            end
            try
            plot(distance,ielevation,'-',...
                    'color',cmap(ncolor,:),'linewidth',3);
%                 
            plot([maxdist+2 maxdist + 7],[edgeHei(i) edgeHei(i)],'-',...
                'color',cmap(ncolor,:),'linewidth',5);
            catch
                ncolor
            end
        else
            if mclrflag
                plot(distance,ielevation,'-',...
                    'color',lineColor(i,:),'linewidth',1.5);
                %n = n + 1;
            else
                plot(distance,ielevation,'-',...
                    'color',lineColor,'linewidth',1.5);
            end
        end
%         n2 = n2 + 1;
%     end
    end
    disp(nage);
    xlim([0 maxdist+ 7]);
    ylim([-0.5 3]);
    set(gca,'box','on');
    % Create xlabel
    xlabel({'Distance (km)'});

    % Create ylabel
    ylabel({'Elevation (km)'});
     

        
    xtickformat('%.1f');ytickformat('%.1f');
    
    yyaxis right;
    
    if ~isempty(edgeHei)
        return
    end

    edgeHei(edgeHei == 0) = [];
    age(isnan(age) | isempty(age)) = [];
    edgeHeiage = [edgeHei, age'];
    edgeHeiage = sortrows(edgeHeiage);
    edgeHei = edgeHeiage(:,1);
    age = edgeHeiage(:,2);
    
    edgeHei1 = round(edgeHei * 10);
    
    idx = zeros(5,1);
    elearray = 4:25;
    for i = 1:length(elearray)
        if ~isempty(find(edgeHei1 == elearray(i),1))
            idx(i) = find(edgeHei1 == elearray(i),1);
        end
    end
    idx(idx ==0)=[];
    
    set(gca,'YTick',edgeHei(idx),'YTickLabel',cellstr(num2str(round(age(idx)/1000))));
    ylim([-0.5 3]); ylabel('Ice age (kyr)');
    hold off;
    set(findall(gcf,'-property','FontSize'),'FontSize',18);
end

