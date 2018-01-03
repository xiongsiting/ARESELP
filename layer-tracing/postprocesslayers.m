function [newpeaks, label] = postprocesslayers(ipeaks,DIST,HDIST)

if nargin < 3
    HDIST = inf;
end

if nargin < 2
    DIST = 7;
end
    
    
label = 1:max(ipeaks(:));
points = zeros(length(label),5);
% layer length, y1, x1, y2, x2

newpeaks = ipeaks;
fullwidth = zeros(size(newpeaks));

% Initial setting
for i = 1:length(label)
    [y,x] = find(newpeaks == label(i));
    points(i,1) = length(x);
    
    % keep the full length layers as the original labels
    if length(x) == size(newpeaks,2)
        indices = sub2ind(size(newpeaks),y,x);
        fullwidth(indices) = label(i);
    end
    [xmin,xind] = min(x); [xmax,xind2] = max(x);
    points(i,2:5) = [y(xind) xmin y(xind2) xmax];

end

for i = 1:length(label) % is the label array same as i???

    % fullwidth or no layer has assigned ith label
    
    if points(i,1) == size(newpeaks,2) || isempty(find(newpeaks==i,1))
        continue;
    end

    % to the layer shorter than image width
    % connect left ones
    while points(i,3) > 1
        % connect label(i) layer to label(n)?
        n = connectleft(i,newpeaks,points,DIST,HDIST);
        if n == -1
            break;
        end

        % check backwards
        checki = connectright(n,newpeaks,points,DIST,HDIST);


        if checki > 0 && label(checki) == label(i)
            % connect the label(n) layer to label (i)
            newpeaks(newpeaks == label(n)) = label(i);
            
            % % they have the same label now
            label(label == n) = label(i);
            
            % % change start point to the extended ones
            points(label == i,2) = points(n,2);
            points(label == i,3) = points(n,3);
            points(n,4) = points(i,4);
            points(n,5) = points(i,5);
        else
            break;
        end
        
    end
    
    % connect right ones
    while points(i,5) < size(newpeaks,2)
        n = connectright(i,newpeaks,points,DIST,HDIST);
        if n == -1
            break;
        end
        checki = connectleft(n,newpeaks,points,DIST,HDIST);
        if checki > 0 && label(checki) == label(i)
            newpeaks(newpeaks == label(n)) = label(i);
            label(label == n) = label(i);
            %% change end point to the extended ones
            points(label == i,4) = points(n,4);
            points(label == i,5) = points(n,5);
            points(n,2) = points(i,2);
            points(n,3) = points(i,3);
        else
            break;
        end
    end
    
end

label = unique(label);
for i = 1:length(label)
    if isempty(find(newpeaks == label(i),1))
        newpeaks(newpeaks == label(i)) = 0;
        label(i) = -1;
    end
end
label(label == -1) = [];


% figure; imagesc(newpeaks); hold on;
% for i = 1:length(label)
% [y,x] = find(newpeaks == label(i));
% plot(x,y,'.','color',rand(1,3));
% end

end

function n = connectleft(i,newpeaks,points,DIST,HDIST)

n = -1;
starty = points(i,2); startx = points(i,3);

% in the point ends array ,find ones left to current
ind = find(points(:,5) <= startx);
if isempty(ind)
    return;
else
    if length(ind) > 1
        % sort from nearest to furtherst segments
        tmp = sortrows([points(ind,5) ind]);
        tmp = flipud(tmp);
        ind = tmp(:,2); % get their location in points array
    end
end

% cutting the layers at the start x, and find ys of other layers
pos1 = find(newpeaks(:,startx) > 0);
% get other labels of other layers
pos1_lbl = newpeaks(pos1,startx);

% search every layer segment which has end left to the current start
for j = 1:length(ind)
    searchedy = points(ind(j),4); searchedx = points(ind(j),5);
    if abs(searchedx - startx) > HDIST
        n = -1;
        return;
    end
    try
    % find the common layer between current start and searched end
    pos2 = find(newpeaks(:,searchedx) > 0);
    pos2_lbl = newpeaks(pos2,searchedx);
    [~,ind1,~] = intersect(pos1_lbl,pos2_lbl); % ind can be an array
    catch
        searchedx
    end
    if isempty(ind1)
        return;
    end
    % find the nearest common layer to the current start point
    [~,pos_ind] = min(abs(starty - pos1(ind1)));
    % pos_ind is for the common list

    pos_lbl = pos1_lbl(ind1(pos_ind)); % % label of the common layer
    dist1 = starty - pos1(pos1_lbl == pos_lbl);
    dist2 = searchedy - pos2(pos2_lbl == pos_lbl);
    try
    if abs(dist1 - dist2) < DIST && dist1*dist2 > 0
        n = ind(j); % n is row no. in points array
        return;
    end
    catch 
        dist1 = dist1;
    end
end

end

function n = connectright(i,newpeaks,points,DIST,HDIST)
n = -1;

endy = points(i,4); endx = points(i,5);
ind = find(points(:,3) >= endx);
if isempty(ind)
    return;
else
    if length(ind) > 1
        tmp = sortrows([points(ind,3) ind]);
        ind = tmp(:,2);
    end
end

pos1 = find(newpeaks(:,endx) > 0);
pos1_lbl = newpeaks(pos1,endx);

for j = 1:length(ind)
    searchedy = points(ind(j),2); searchedx = points(ind(j),3);
    if abs(searchedx - endx) > HDIST
        n = -1;
        return;
    end
    
    pos2 = find(newpeaks(:,searchedx) > 0);
    pos2_lbl = newpeaks(pos2,searchedx);
    [~,ind1,~] = intersect(pos1_lbl,pos2_lbl);
    if isempty(ind1)
        return;
    end
    [~,pos_ind] = min(abs(endy - pos1(ind1)));
    pos_lbl = pos1_lbl(ind1(pos_ind));
    try
    dist1 = endy - pos1(pos1_lbl == pos_lbl);
    dist2 = searchedy - pos2(pos2_lbl == pos_lbl);
    if abs(dist1 - dist2) < DIST && dist1*dist2 > 0
        n = ind(j);
        return;
    end
    catch
        dist1 = dist1;
    end
end
end

% 
% function n = findleft(i,newpeaks,points)
% DIST = 16;
% n = -1;
% starty = points(i,2); startx = points(i,3);
% ind = find(points(:,5) <= startx);
% if isempty(ind)
%     return;
% end
% 
% tmp = sort([points(ind,5) ind],'descend');
% ind = tmp(:,2);
% 
% pos1 = find(newpeaks(:,startx) > 0);
% pos1_lbl = newpeaks(pos1,startx);
% 
% for j = 1:length(ind)
%     searchedy = points(ind(j),4); searchedx = points(ind(j),5);
%     pos2 = find(newpeaks(:,searchedx) > 0);
%     pos2_lbl = newpeaks(pos2,searchedx);
%     [~,ind1,~] = intersect(pos1_lbl,pos2_lbl);
%     [~,pos_ind] = min(abs(starty - pos1(ind1)));
%     pos_lbl = pos1_lbl(ind1(pos_ind));
%     dist1 = starty - pos1(pos1_lbl == pos_lbl);
%     dist2 = searchedy - pos2(pos2_lbl == pos_lbl);
%     if abs(dist1 - dist2) < DIST
%         n = ind(j);
%         return;
%     end
% end
% % 
% % 
% % pos = find(fullwidth(:,startx) > 0);
% % [~,pos_ind] = min(abs(starty-pos));
% % dist1 = starty - pos(pos_ind);
% % pos_lbl = fullwidth(pos(pos_ind),startx);
% % [y,~] = find(fullwidth == pos_lbl);
% % 
% % 
% % 
% % if  isempty(ind)
% %     return;
% % end
% % searched_y = points(ind,4); searched_x = points(ind,5);
% % pos_y = y(searched_x);
% % searched_y((searched_y - pos_y).*dist1 < 0) = -9999;
% % searched_y(abs(searched_y - pos_y - dist1) > DIST) = -9999;
% % searched_x(searched_y == -9999) = -inf;
% % [~,mind] = max(searched_x);
% % 
% % % [dist,mind] = min(abs(searched_y - pos_y - dist1));
% % % if dist > DIST
% % %     return;
% % % end
% % 
% % n = ind(mind);
% 
% end
% 
% function n = findright(i,fullwidth,points)
% DIST = 16;
% n = -1;
% 
% endy = points(i,4); endx = points(i,5);
% pos = find(fullwidth(:,endx) > 0);
% [~,pos_ind] = min(abs(endy-pos));
% dist1 = endy - pos(pos_ind);
% pos_lbl = fullwidth(pos(pos_ind),endx);
% [y,~] = find(fullwidth == pos_lbl);
% 
% ind = find(points(:,3) >= endx);
% if isempty(ind)
%     return;
% end
% searched_y = points(ind,2); searched_x = points(ind,3);
% pos_y = y(searched_x);
% searched_y((searched_y - pos_y).*dist1 < 0) = -9999;
% searched_y(abs(searched_y - pos_y - dist1) > DIST) = -9999;
% searched_x(searched_y == -9999) = inf;
% [~,mind] = min(searched_x);
% 
% % [dist,mind] = min(abs(searched_y - pos_y - dist1));
% % if dist > DIST
% %     return;
% % end
% 
% n = ind(mind);
% end
