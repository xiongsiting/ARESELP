function [icoefs,inums,ikurto, iratio ] = findpeakscwt(idata,scales,wavelet,...
    idxStart,idxEnd,bgSkip,fplot)

    if nargin < 7
        % uncomment this to plot scalogram
        fplot = [];
    end

    if nargin < 6
        bgSkip = 50;
    end
     
    if nargin < 5
        idxEnd = length(idata);
    end
    
    if nargin < 4
        idxStart = 1;
    end

%     ipeaks = zeros(size(idata));
    icoefs = zeros(size(idata)); % sum of coeffcients
    inums = zeros(size(idata)); % sum of peak numbers on scalogram
%     istd = zeros(size(idata));
%     iskew = zeros(size(idata));
    if nargout > 2
        ikurto = zeros(size(idata));
        iratio = zeros(size(idata));
    end

    cwtData = cwt(idata,scales,wavelet);
%     if strcmp(fplot,'plot')
%         imagesc(cwtData); hold on;
%         colormap('gray'); colorbar; 
%         axis on; 
%     end
    
    for r = 1:size(cwtData,1)
        cwtRow = cwtData(r,:);
          bgmax = 0;
%         trim end
%         cwtRow(length(idata) - bgSkip * r:end) = 0;
%         background maximum, trim down the last 150 pixels
        if idxEnd > length(idata) - bgSkip 
            bgmax = 0;
        else
            bgSig = cwtRow(idxEnd + bgSkip:size(cwtRow,2) - bgSkip - scales(r));
            bgmax = max(bgSig);
            try
            if isempty(find(cwtRow(idxStart:idxEnd) > bgmax,1))
                % edge affect the wavelet transform
                continue;
            end
            catch
                bgmax;
            end
            
% %             bgSkipStep = 5;
% %             while isempty(find(cwtRow(idxStart:idxEnd) > bgmax,1))
% %                 if idxEnd + bgSkip <= size(cwtRow,2)- bgSkip
% %                     bgSkip = bgSkip + bgSkipStep;
% %                     fprintf('No peak found in %d of CWT scalogram,bgmax:%f\n',r,bgmax);
% %                     % if the bgSkip has been narrow to be no signal
% %                     % then bgmax is infinity which makes no locs to be found
% %                     bgSig = cwtRow(idxEnd + bgSkip:size(cwtRow,2)- bgSkip);
% %                     bgmax = max(bgSig);
% %                 else
% %                     bgmax = inf;
% %                     break;
% %                 end
% %             end
        end
%                 
      
        % findpeaks along the cwt row
        [~,locs] = findpeaks(cwtRow,'minpeakheight',bgmax);
        % if there is no peaks found in this cwt row, continue in next row
        if isempty(locs)
            continue;
        end
        % remove peaks above surface and below bottom
        % remove also surface and bottom layers
%         locs(locs <= idxStart + 5 | locs >= idxEnd - 2 * bgSkip) = [];

        if ~strcmp(fplot,'plot')
            locs(locs < idxStart| locs > idxEnd) = [];

            halfScale = floor(scales(r)/2);
            for j = 1:length(locs)
                sigSegment = idata(locs(j) - halfScale: locs(j) + halfScale);
                [~,ind] = max(sigSegment);
                if locs(j) - halfScale - 1  + ind < idxEnd
                    locs(j) = locs(j) - halfScale - 1  + ind;
                end
            end
        end

        icoefs(locs) = icoefs(locs) + cwtRow(locs)';
        inums(locs) = inums(locs) + 1; 
        
        if nargout > 2
            % Higher order statistics
            if scales(1) == 1
                shwin = floor(scales(2)/2);
            else
                shwin = floor(scales(1)/2);
            end
            for i = 1:length(locs)
                lhwin = floor(scales(r)/2);
    %             rstd = std(idata(locs(i) - hfwin:locs(i) + hfwin));
    %             skew = skewness(idata(locs(i) - hfwin:locs(i) + hfwin));
                kurto = kurtosis(idata(locs(i) - lhwin:locs(i) + lhwin));
                ratio = mean(idata(locs(i) - shwin:locs(i) + shwin))/...
                    mean(idata(locs(i) - lhwin:locs(i) + lhwin));

    %             if isnan(rstd) || isnan(skew) || isnan(kurto)
                if isnan(ratio) || isnan(kurto)
                    continue;
                end

    %             istd(locs(i)) = istd(locs(i)) + rstd;
    %                 
    %             iskew(locs(i)) = iskew(locs(i)) + skew;

                ikurto(locs(i)) = ikurto(locs(i)) + kurto;
                iratio(locs(i)) = iratio(locs(i)) + ratio;

            end
        
        end
        % uncomment this to plot the local maxima on scalogram
        if strcmp(fplot,'plot')
%             yyaxis left;
%             bar(inums,'edgecolor','c','facecolor','c');
%             ylabel('Nos. of CWT Maxima');
%             set(gca,'YColor','k','Ydir','normal');
            
%             yyaxis right;
%             ylim([1 23]);
            if strcmp(wavelet,'mexh')
                colorline = 'b.';
            else
                colorline = 'b.';
            end
            plot(locs - 1,ones(size(locs)) * r,colorline);
%             scaletitle = ['scales (1-' num2str(scales(end)) ')']; 
%             ylabel(scaletitle);
%             set(gca,'YColor','k');
        end
        
    end
    
    
    % uncomment code below to calculate the prominence of peaks
    % this codes are not working very well
%     % calculate prominence
%     locs = find(inums > 0);
%     mmdiff = zeros(size(locs));
%     for j = 1:length(locs)
% 
%         if j < 2
%             between = idata(idxStart:locs(3));
%         else
%             if j > length(locs) - 1
%                 between = idata(locs(j-1):idxEnd);
%             else
%                 between = idata(locs(j-1):locs(j+1));
%             end
%         end
%         try
%         [minval,~] = min(between);
%         mmdiff(j) = idata(locs(j)) - minval;  
%         catch
%             j=j;
%         end
%     end
%     ipeaks(locs) = mmdiff;
    
end