function convert2mat(filename)
    mat = netcdf_to_mat(filename);
    Bottom = mat.Bottom';
    Surface = mat.Surface';
    Time = mat.fasttime'*10e-7;
    Elevation = mat.altitude';
    Latitude = mat.lat';
    Longitude = mat.lon';
    GPS_time = mat.time';
    Data = mat.amplitude;
    %% find surface and Bottom
    [ysrf,ybtm] = picksurfacebottom(Data);
    for i = 1:length(ysrf)
        if isnan(Surface(i))
            Surface(i) = Time(ysrf(i));
        end
        if isnan(Bottom(i))
           Bottom(i) = Time(ybtm(i));
        end
    end
    clear mat;
    filename_prefix = split(filename,'.');
    outfilename = filename_prefix(1)+'.mat';
    clear filename;
    save(char(outfilename));
end

function [ysrf,ybtm] = picksurfacebottom(data)
SKIP1 = 50; SKIP2 = 100;
ysrf = zeros(length(size(data,2)),1);
ybtm = ysrf;
data = mat2gray(data);
for i = 1:size(data,2)
    profile = data(:,i);
    a = zeros(length(profile),1);
    for j = 31:length(profile)
        a(j) = profile(j).^2/mean(profile(j-30:j-1).^2);
    end
    [~,locs] = findpeaks(a(1:end-200),'sortstr','descend');
    ind2 = find(locs > 700);
    locs2 = locs(ind2(1));
    ind1 = find(locs < 200);
    locs1 = locs(ind1(1));
    
    if locs1 > locs2 
        ysrf(i) = locs2;
        ybtm(i) = locs1;
    else
        ysrf(i) = locs1;
        ybtm(i) = locs2;
    end
end
ysrf = medfilt1(ysrf,SKIP1);
ybtm = medfilt1(ybtm,SKIP2);
msrf = floor(median(ysrf));
mbtm = floor(median(ybtm));

for i = 1:size(data,2)
    profile1 = data(msrf-SKIP1:msrf+SKIP1,i);
    profile2 = data(mbtm-SKIP2:mbtm+SKIP2,i);
    [~,maxind] = max(profile1);
    ysrf(i) = maxind + msrf - SKIP1 + 1;
    [~,maxind] = max(profile2);
    ybtm(i) = maxind + mbtm - SKIP2 + 1;   
end

for i = 1:size(data,2)
    if i <= 50
        ysrf(i) = mean(ysrf(i:i+50));
        ybtm(i) = mean(ybtm(i:i+50));
    else
        if i > size(data,2) - 50
            ysrf(i) = mean(ysrf(i-50:i));
            ybtm(i) = mean(ybtm(i-50:i));
        else
            ysrf(i) = mean(ysrf(i-25:i+25));
            ybtm(i) = mean(ybtm(i-25:i+25));
        end
    end
end
ysrf = round(ysrf); ybtm = round(ybtm);
end