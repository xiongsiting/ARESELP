function imLayer = traceonelayer(peakim,seedpoints,imLayer,labelOffset,params)

DIST = 7;
BLOCKSIZE = 51;
DELTA_ORI_THRESH = inf;
gapfill = ceil(0.25 * BLOCKSIZE);
buffdist = 2 * DIST;
HDIST = 3 * BLOCKSIZE;

if ~isempty(params)
    DIST = params{1};
    if DIST < 3
        DIST = 3;
    end

    BLOCKSIZE = params{2};

    if max(size(params)) > 2
        DELTA_ORI_THRESH = params{3};
    end

    if max(size(params)) > 3
        gapfill = params{4};
    end  

    if max(size(params)) > 4
        buffdist = params{5};
    end
    
    if max(size(params)) > 5
        HDIST = params{6};
    end
end
        
    if ~mod(BLOCKSIZE,2)
        hfblk = BLOCKSIZE /2;
    else
        hfblk = (BLOCKSIZE - 1)/2;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    jlocs = seedpoints(1,2); % col no.
    ilocs = seedpoints(1,1); % row no.
   
    first_angle = 0;
    pre_angle = 0; % every check set initial angle = 0
    % this is used for smoothing
    imLayer(ilocs,jlocs) = labelOffset + 1;

    j = jlocs;
    while j <= size(peakim,2)
        directright = 1; % right search
        % extract block image
        [blkim,xrgt] = extractblockim(peakim,ilocs,j,hfblk,directright);
        if isempty(blkim) || xrgt < 1
            break; end
        
        % calculate HT slope angle
%         [hough_angle,~] = htSlopeAngle(blkim,gapfill,buffdist); %,HT_N,gapfill);
        hough_angle = rdSlopeAngle(blkim);
        if isnan(hough_angle)
            break; end  
%         if abs(pre_angle - hough_angle) > DELTA_ORI_THRESH, break;end

        % right propogation
        xx = 0 : 1 : xrgt;           
        yy = hough_angle * xx;
        
        % bound in block
        idx = find(yy > hfblk | yy < -hfblk);
        xx(idx) = []; yy(idx) = [];

        % location in image
        xx = j + xx; 
        yy = ilocs + round(yy);

        condition = struct('mindist',DIST,'crossing','',...
            'deltaangle',[pre_angle,DELTA_ORI_THRESH]);
        stopflag = stoptracing(xx,yy,imLayer,labelOffset,...
                directright,condition);
        if stopflag == 1, 
            break; end
       
       % inverse the line coordinates of the right and left lines
%         indices = sub2ind(size(imLayer),yy,xx);
        
        indices1 = sub2ind(size(imLayer),yy(2:end),xx(2:end));
        imLayer(indices1) = 1 + labelOffset;
        
        pre_angle = hough_angle;
        if j == jlocs
            first_angle = pre_angle;
        end

        j = xx(end);
        ilocs = yy(end); % this is not yy(end)???
        
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % go back to seed point
    jlocs = seedpoints(1,2); % col no.
    ilocs = seedpoints(1,1); % row no.
    
    pre_angle = first_angle; % every check set initial pre_angle = 0
    
    j = jlocs;
    while j >= 1
        
        % search left
        directright = 0;
        
        [blkim,xlft] = extractblockim(peakim,ilocs,j,hfblk,directright);
        if isempty(blkim) || xlft < 1
            break;end
        
%        [hough_angle,~] = htSlopeAngle(blkim,gapfill,buffdist); %,HT_N,gapfill);
       hough_angle = rdSlopeAngle(blkim);
       if isnan(hough_angle)
           break; end 
%        if abs(pre_angle - hough_angle) > DELTA_ORI_THRESH,break;end

        xx = -xlft: 1 : 0;
        yy = hough_angle *xx; % if angle larger than +- 45 degrees,
                               % yy would be out the bound of block
                               
        % bound in block
        idx = find(yy > hfblk | yy < -hfblk);
        xx(idx) = []; yy(idx) = [];
        
        xx = j + xx; %yy = yy + floor(0.3* shift);
        yy = ilocs + round(yy);         
        
        condition = struct('mindist',DIST,'crossing','',...
            'deltaangle',[pre_angle,DELTA_ORI_THRESH]);
        stopflag = stoptracing(xx,yy,imLayer,labelOffset,...
            directright,condition);
        if stopflag, 
            break;end    
        % inverse the line coordinates of the right and left lines
%         try
%         indices = sub2ind(size(imLayer),yy,xx);
        indices1 = sub2ind(size(imLayer),yy(1:end-1),xx(1:end-1));
        imLayer(indices1) = 1 + labelOffset;
        
        j = xx(1);
        ilocs = yy(1);
        pre_angle = hough_angle;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % clean up
    
    ilength = sum(sum(imLayer == labelOffset + 1));
    if ilength < HDIST
        imLayer(imLayer == labelOffset + 1) = 0;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
% % % %         
% % % %         
% % % %     for j = jlocs : hfblk : size(peakim,2)
% % % %         directright = 1;
% % % %         
% % % %         
% % % % 
% % % %         xx = 0 : 1 : xprop;           
% % % %         yy = hough_angle * xx;
% % % %         if sum(yy > hfblk)~=0 || sum(yy < -hfblk)~=0
% % % %             break;
% % % %         end
% % % %         
% % % %         % uncomment this for adjusting the layers         
% % % % %         yy(end) = yy(end) + shift;
% % % % %         hough_angle = (yy(end) - yy(1)) / (xprop);
% % % % %         yy = hough_angle * xx;
% % % % 
% % % %         xx = j + xx; %yy = yy + floor(0.3* shift);
% % % %         yy = ilocs + round(yy);
% % % %   
% % % %         condition = struct('mindist',DIST,'crossing','','maxangle',maxangle);
% % % %         stopflag = stoptracing(xx,yy,imLayer,labelOffset,...
% % % %             directright,condition);
% % % %         if stopflag == 1, break; end
% % % % 
% % % %             
% % % % 
% % % %     end
% % % %     
% % % %     ilocs = seedpoints(1,1);
% % % %     angle = 0;
% % % %     for j = jlocs : -hfblk : 1
% % % % 
% % % %         directright = 0;
% % % %         [blkim,xprop] = extractblockim(peakim,ilocs,j,hfblk,directright);
% % % %         if isnan(blkim), break;end
% % % %         
% % % %        [hough_angle,~] = htSlopeAngle(blkim,gapfill,buffdist); %,HT_N,gapfill);
% % % %        if isnan(hough_angle), break; end 
% % % %        if abs(angle - hough_angle) > DELTA_ORI_THRESH,break;end
% % % % 
% % % %         xx = -xprop: 1 : 0;
% % % %         yy = hough_angle *xx;
% % % %         if sum(yy > hfblk)~=0 || sum(yy < -hfblk)~=0
% % % %             break;
% % % %         end
% % % %         
% % % %         % shift value used to adjust the layers during tracing
% % % % %         yy(end) = yy(end) + shift;
% % % % %         hough_angle = (yy(end)-yy(1))/ (xx(end) - xx(1));
% % % % %         yy = hough_angle * xx;
% % % %         
% % % %         xx = j + xx; %yy = yy + floor(0.3* shift);
% % % %         yy = ilocs + round(yy);         
% % % %         
% % % %         condition = struct('mindist',DIST,'crossing','','maxangle',maxangle); 
% % % %         stopflag = stoptracing(xx,yy,imLayer,labelOffset,...
% % % %             directright,condition);
% % % %         if stopflag, break;end    
% % % %         % inverse the line coordinates of the right and left lines
% % % % %         try
% % % %             indices = sub2ind(size(imLayer),yy,xx);
% % % %             imLayer(indices) = 1 + labelOffset;
% % % %             ilocs = ilocs - round(hough_angle * hfblk);
% % % %             angle = hough_angle;
% % % %     end
% % % %     

   
end
   
% end

function [blkim,headingwidth] = extractblockim(peakim,ilocs,j,hfblk,directright)

blkim = [];headingwidth = 0;

if j == size(peakim,2) || j == 1
    blkim = [];
    headingwidth = 0;
    return;
end

if directright
    
    ymin = ilocs - hfblk; ymax = ilocs + hfblk;
%     ycen = hfblk + 1;xcen = hfblk + 1; 
        
        if j > size(peakim,2) - hfblk 
            xmin = j - hfblk; xmax = size(peakim,2);
            xprop = size(peakim,2) - j; 
            % if j == size(peakim,2) - hfblk + 1
            % xprop = hfblk -1;
        else
            if j <= size(peakim,2) - hfblk && j > hfblk
                xmin = j - hfblk; xmax = j + hfblk; 
                xprop = hfblk;
            else
                xmin = 1; xmax = j + hfblk;
%                 xcen = j; 
                xprop = hfblk;
            end
        end
        if ymin < 1 || ymax > size(peakim,1) || xmin < 1 || xmax > size(peakim,2)
            return;
        else
            headingwidth = xprop;
            blkim = peakim(ymin:ymax, xmin:xmax);
        end
        
else
    
    ymin = ilocs - hfblk; ymax = ilocs + hfblk;
%     ycen = hfblk + 1; xcen = hfblk + 1; 
    if j < hfblk + 1
        xmin = 1; xmax = j + hfblk;
%         xcen = j; 
        xprop = j-1;
    else
        if j <= size(peakim,2)- hfblk && j > hfblk
            xmin = j - hfblk; xmax = j + hfblk;
            xprop = hfblk;
        else
            xmin = j - hfblk; xmax = size(peakim,2);
            xprop = hfblk;
        end
    end
    
    if ymin < 1 || ymax > size(peakim,1) || xmin < 1 || xmax > size(peakim,2)
        return;
    else
    
        blkim = peakim(ymin:ymax,xmin:xmax);
        headingwidth = xprop;
    end
    
end

end

function [stopflag, imLayer] = stoptracing(xx,yy,imLayer,labelOffset...
    ,directright,condition)

stopflag = 0;

if length(xx) == 1 || length(yy) == 1
    stopflag = 1;
    return;
end

if directright
    idxhead = length(yy);
    idxend = 1;
    idxend2 = idxend + 1;
else
    idxhead = 1;
    idxend = length(yy);
    idxend2 = idxend -1;
end

%% Do not change the order of the following codes
% not used 
% if isfield(condition,'maxangle')
%     agtan = (yy(idxhead) - yy(idxend))/(xx(idxhead) - xx(idxend));
%     if abs(agtan) > tan(condition.maxangle)
%         stopflag = 1;
%         return;
%     end
% end

if isfield(condition,'deltaangle')
    pre_angle = condition.deltaangle(1);
    thresh = condition.deltaangle(2);
    agtan = (yy(idxhead) - yy(idxend))/(xx(idxhead) - xx(idxend));
    if abs(atan(agtan) - atan(pre_angle)) > thresh /180 *pi
        stopflag = 1;
        return;
    end
end


if isfield(condition,'crossing')
    % inverse the line coordinates of the right and left lines
    % prevent from cross-lines, add on 03-07-2014
    pt1x = xx(min(idxend2,idxhead));
    pt2x = xx(max(idxend2,idxhead));
    pt1y = yy(min(idxend2,idxhead));
    pt2y = yy(max(idxend2,idxhead));
    
    yrange = min(pt1y,pt2y): max(pt1y,pt2y);
    xrange = min(pt1x,pt2x): max(pt1x,pt2x);

    subim = imLayer(yrange,xrange);
    if sum(subim(:)) ~= 0
        imLayer(imLayer == 1 + labelOffset) = 0; % remove this layer
        stopflag = 1;
        return;
    end
end

    
if isfield(condition,'mindist')

    ind = find(imLayer(:,xx(idxhead))~=0);
    if ~isempty(ind)
        [minval,minidx] = min(abs(ind-yy(idxhead)));
        
        elabel = imLayer(ind(minidx),xx(idxhead));
        [~, xexist] = find(imLayer == elabel);

        % when two layers are about to intersect
        if minval < condition.mindist
        % all points in the closest layer locate the opposite side
        % to the current layer, then connect
            if directright && isempty(xexist < xx(idxend)) || ...
                    ~directright && isempty(xexist > xx(idxend))
                  % leave to post-processing, just stop tracing here
%                 if  minval < floor(condition.mindist/4)
%                     % connect two of these layers
%                     imLayer(imLayer == 1 + labelOffset) = elabel;
%                     indices = sub2ind(size(imLayer),yy,xx);
%                     imLayer(indices) = elabel;
%                 end
                stopflag = 1;
                return;
            end

            % previous one has priority, so delete this one
            imLayer(imLayer == 1 + labelOffset) = 0; 
            stopflag = 1;
            return;
        end
    end
%     % It seems this can be deleted, not useful???
    if length(xx) > 1
        ind = find(imLayer(:,xx(idxend2))~=0);
        if ~isempty(ind)
            [minval,~] = min(abs(ind-yy(idxend2)));
            if minval <= condition.mindist
                imLayer(imLayer == 1 + labelOffset) = 0;
                stopflag = 1;
                return;
            end
        end
    end
end
end



