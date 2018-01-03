function [imLayer,debuginfo] = tracelayers(peakim,seedpoints,params)

    if nargin < 2
        N = 10000;
    else 
        if size(seedpoints,2) == 1
            % input seedpoints are not an array
            N = seedpoints;
        end
    end
    
    if nargin < 3
        params = {7,51,15};    
    end
    
    if size(seedpoints,2) == 1 && ~isempty(N)
        % input seedpoints are not an array
        [rows,cols] = find(peakim > 0);
        val = peakim(peakim > 0);
        seedpoints = [val rows cols];
        seedpoints = flipud(sortrows(seedpoints,1));
        seedpoints(N + 1:end,:) = [];
        seedpoints(:,1) = [];
    else
        if size(seedpoints,2) ~= 2
            return;
        end
    end
    
    debuginfo = [];
    %% FIND PEAKS FROM PEAK IMAGE
    labelOffset = 0;
    imLayer = zeros(size(peakim));
    mindist = params{1};
    fprintf('Searching for %d seedpoints\n',size(seedpoints,1));
    while ~isempty(seedpoints)
%         tic;
%         if ~mod(i,5000)
%             fprintf('%d......',i);toc;
%         end
%
        locs = seedpoints(1,1);
        cols = seedpoints(1,2);
        seed_point = [locs,cols];
        pos = find(imLayer(:,cols) ~=0);
        [minval,~] = min(abs(pos-locs));
        if isempty(pos) || minval >=mindist
            imLayer = traceonelayer(peakim,seed_point,imLayer,...
                labelOffset,params);
            
            if max(imLayer(:)) == labelOffset + 1
                [y,x] = find(imLayer == labelOffset + 1);
                xx = repmat(x',mindist*2 + 1,1);
                yy = zeros(size(xx));
                for i = 1:2*mindist + 1
                    yy(i,:) = y + i - mindist - 1;
                end
                N = length(x) * (2 * mindist + 1);
                buffer = [reshape(yy,[N,1]),reshape(xx,[N,1])];
                [~,~,idx] = intersect(buffer,seedpoints,'rows');
                seedpoints(idx,:) = [];
                
                labelOffset = labelOffset + 1;

                debuginfo = vertcat(debuginfo,[labelOffset,length(x),...
                    size(seedpoints,1)]);
                fprintf('%d th layer\tlength: %d\tpoints left: %d ...\n',...
                    labelOffset, length(x),size(seedpoints,1));
            else
                seedpoints(1,:) = [];
            end
        else
            seedpoints(1,:) = [];
        end
        
    end
    % else
    %     display('no seed points have been chosen');
    % end



end