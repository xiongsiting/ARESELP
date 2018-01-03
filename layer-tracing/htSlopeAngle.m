function [ht_agtan,shift] = htSlopeAngle(im,gapfill,buffdist)%, N, K)

    if nargin < 3
        buffdist = ceil(size(im,1)/3);
    end
    if nargin < 2 || gapfill == 0
        gapfill = ceil(size(im,2) * 0.25);
    end
    ht_agtan = nan;shift = 0;
    
    % keep minlen not change
    minlen = ceil(size(im,2) * 0.25);
    
    ht_agtan_initial = slopestimate(im,gapfill,minlen);
    if isnan(ht_agtan_initial),return; end

    [newim, shift] = adjustslopeim(im,ht_agtan_initial,buffdist);
    % newim might be invalide due to indices cutting
    if isnan(newim), ht_agtan = ht_agtan_initial; return; end
%     imagesc(newim); hold on;
%     plot(51:101,round((1:51)*ht_agtan_initial) + 51,'r-','linewidth',3);
%     gapfill = ceil(sqrt(2) * gapfill); % change gapfill for the second time
    ht_agtan = slopestimate(newim,gapfill,minlen);
    
% % % % % % % % % %     ht_agtan_new = houghreal(im);
% % %     ht_agtan = ht_agtan_initial;
%     plot(51:101,round((1:51)*ht_agtan) + 51,'w-','linewidth',3);

end

function ht_agtan = slopestimate(im, gapfill,minlen)
    if nargin < 2
        gapfill = ceil(size(im,2) * 0.25);
        minlen = ceil(size(im,2) * 0.25);
    end

    ht_agtan = nan;
%     gapfill = 10;
    
    % im = edge(im,'canny');
    % thetares = 2 / size(im,2);  % theta resolution for increasing one y axis
    % when resolution increases, the computing time increase a lot
    % for a random image of 51 * 51 pixels
    % [H, theta,rho] = hough(im) uses 0.005 seconds,
    % while [H, theta,rho] = hough(a,'Theta',-90:thetares:89) uses 0.17 seconds
    % maybe use thetares = [-90:thetares:-45 45:thetares:89];
    [H, theta,rho] = hough(im);
    % [H, theta,rho] = hough(im,'Theta',[-90:thetares:-45 45:thetares:89]); 
    % less than 45 degrees
    P = houghpeaks(H); % default is pick the largest peak in HT domain

    if isempty(P)
        return;
    else
    % % %     %% usng Matlab houghlines
    lines = houghlines(im,theta,rho, P, 'FillGap',gapfill,'MinLength',minlen);
    % lines = houghlines(im,theta,rho, P, 'FillGap',gapfill,'MinLength',minlen);

    % figure; imagesc(im); hold on;
    % for i = 1:max(size(lines))
    %     xy = [lines(i).point1; lines(i).point2];
    %     plot(xy(:,1),xy(:,2),'w-','linewidth',3);
    % end

    if isempty(lines)
        return;
    else
        ht_angle = 0;
        nline = 0;
        for k = 1:length(lines)
            xy = [lines(k).point1; lines(k).point2]; % column 1 is x, 2 is y
            if xy(2,1) - xy(1,1) == 0
                nline = nline - 1;
                continue;
            else
                % angle is 0 at east
                % positive clockwise
                angle = (xy(2,2)-xy(1,2))/(xy(2,1)-xy(1,1));
                ht_angle = ht_angle + angle;
                nline = nline + 1;
            end                                        
        end
        if nline < 1
            return;
        else
            ht_angle = ht_angle / nline;
            ht_agtan = ht_angle;
        end

    end
    end
end



function [newim, shift] = adjustslopeim(im,ht_agtan,buffdist)

% % % buffdist = 7;
% % % if nargin < 2 || K == 0
% % %     Thresh = max(im(:)) * 0.1 * min(size(im));
% % % %     Thresh = size(im,2) * 0.5;
% % % %     gapfill = floor(size(im,2) * 0.25);
% % % else
% % %     Thresh = max(im(:)) * K * min(size(im));
% % % %     Thresh = size(im,2) * K;
% % % %     gapfill = K;
% % % end

newim = nan; shift = 0; 
% ratio = nan; % not used
%% check buffer pixels, very slow
    ycen = ceil(size(im,1)/2);
    hfwin = floor(size(im,2)/2);

    x = -hfwin:hfwin;
    y = round(x * ht_agtan); 

    x = x + ceil(size(im,2)/2);
    y = y + ycen;

    if ~mod(size(im,2),2)
        x(1) = [];
        y(1) = [];
    end
    
    % x, y should be the same size as im width as this moment

    idx = find(y <= 0 | y > size(im,1));
    x(idx) = [];
    y(idx) = []; 
    
    nbuff = buffdist*2 + 1;
    if abs(ht_agtan) <= 1
    
        xx = repmat(x,nbuff,1);
        yy = zeros(size(xx));

        for i = -buffdist: buffdist
            yy(buffdist + 1 + i,:) = y + i;
        end
    else
        yy = repmat(y,nbuff,1);
        xx = zeros(size(x));
        for i = -buffdist: buffdist
            xx(buffdist + 1 + i,:) = x + i;
        end
    end

        xx = reshape(xx,[size(xx,1)*size(xx,2),1]);
        yy = reshape(yy,[size(yy,1)*size(yy,2),1]);
        
        idx = find(xx <=0 | xx > size(im,2) | yy <=0 | yy > size(im,1));
        if ~isempty(idx), xx(idx) = []; ...
        yy(idx) = []; end
        
        indices = sub2ind(size(im),yy,xx);
    newim = zeros(size(im));
    newim(indices) = im(indices);
    
%     bim = reshape(im(indices),[nbuff,size(xx,2)/nbuff]);
%     [~,idx] = max(sum(bim,2));
%     shift = idx - buffdist -1;
    
    % uncomment to calculate the ratio of buffer area over overall pixels
    % in block, sum(buffer)/sum(blockim)
    % not used because this function is also needed to calculate shift
%     bsum = sum(bim,2); 
%     ratio = sum(bsum(buffdist:buffdist + 2))/sum(bsum);
%     ratio = sum(bsum)/sum(im(:));
    
    
    
end
% % %     %% working slower than houghlines
% % %     angles = theta(P(:,2));
% % %     ind1 = find(angles <= 0);
% % %     ind2 = find(angles > 0);
% % %     if ~isempty(ind1)
% % %         angles(ind1) = 90 + angles(ind1);
% % %     else
% % %         if ~isempty(ind2)
% % %             angles(ind2) = angles(ind2) - 90;
% % %         
% % %         else
% % %             return;
% % %         end
% % %     end
% % %     
% % %     if angles >= 45 || angles <= -45
% % %         return ;
% % %     else
% % %         ht_angle = mean(angles)/180*pi;
% % %         ht_agtan = tan(ht_angle);
% % %         shift = 0;
% % %         return;
% % %     end




% % %     %     val = sum(im(indices) > 0);
% % %     bim = reshape(im(indices),[nbuff,size(xx,2)/nbuff]);
% % %     bsum = sum(bim,2); 
% % %     ratio = sum(bsum(buffdist:buffdist + 2))/sum(bsum);
% % % 
% % %     %     val = sum(im(indices));
% % %     %     
% % %     %       if abs(ht_agtan) >= 1
% % %     %         nline;
% % %     %        
% % %     %         imagesc(im); hold on;
% % %     %         plot(x,y,'w-','linewidth',3);
% % %     %         hold on; plot(xx,yy,'r.');
% % %     %       end
% % % 
% % %     if ratio > 0
% % % %         bim = reshape(im(indices),[nbuff,size(xx,2)/nbuff]);
% % %         [~,idx] = max(sum(bim,2));
% % %         shift = idx - buffdist -1;
% % % %         y = y - ycen;
% % % %         [~,idx2] = max(bim(:,end));
% % % %         yend = y(end) + idx2 - buffdist - 1;
% % % %         [~,idx1] = max(bim(:,1));
% % % %         yfrst = y(1) + idx1 - buffdist - 1;
% % % %         
% % % %         
% % % %         ht_agtan = (yend - yfrst) / size(im,2);
% % % %         ht_agtan = (y(end) + shift) / x(end);
% % %         ht_agtan = ht_angle;
% % % %         shift = 0;
% % % 
% % %         else
% % %         ht_agtan = nan;
% % %     end
% % % end
