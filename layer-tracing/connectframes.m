function geoinfo = connectframes(geoinfo1,geoinfo2)

DIST = 10;C = 3e8; epsilon = 3.15;
DIST = DIST*(geoinfo1.time_range(2) - geoinfo1.time_range(1))*C/2/sqrt(epsilon);

geoinfo = geoinfo1;
%% Calculate the common part
lat1 = geoinfo1.latitude;
lat2 = geoinfo2.latitude;

[~,ind1] = min(abs(lat1 - lat2(1))); % cut point in the 1st radar echogram
[~,ind2] = min(abs(lat2 - lat1(end))); % cut point in the 2th radar echogram

midlat = 0.5*lat1(ind1) + 0.5*lat2(ind2);

[~,ind1] = min(abs(lat1 - midlat));
[~,ind2] = min(abs(lat2 - midlat));
ind2 = ind2 + 1;

geoinfo.elevation_bed = [geoinfo1.elevation_bed(1:ind1) geoinfo2.elevation_bed(ind2:end)];
geoinfo.elevation_surface = [geoinfo1.elevation_surface(1:ind1) ...
    geoinfo2.elevation_surface(ind2:end)];
geoinfo.latitude = [geoinfo1.latitude(1:ind1) geoinfo2.latitude(ind2:end)];
geoinfo.longitude = [geoinfo1.longitude(1:ind1) geoinfo2.longitude(ind2:end)];
geoinfo.num_trace = ind1 + geoinfo2.num_trace - ind2 + 1;
geoinfo.thickness = [geoinfo1.thickness(1:ind1) geoinfo2.thickness(ind2:end)];
geoinfo.time_gps = [geoinfo1.time_gps(1:ind1) geoinfo2.time_gps(ind2:end)];
geoinfo.traveltime_surface = [geoinfo1.traveltime_surface(1:ind1) ...
    geoinfo2.traveltime_surface(ind2:end)];
geoinfo.x = [geoinfo1.x(1:ind1) geoinfo2.x(ind2:end)];
geoinfo.y = [geoinfo1.y(1:ind1) geoinfo2.y(ind2:end)];
geoinfo.time_range = geoinfo1.time_range;

ele2 = zeros(geoinfo2.num_layer,2);
for j = 1:geoinfo2.num_layer
   ele2(j,1) = geoinfo2.layer(j).elevation(ind2);
end

num_layer = 1;
for i = 1:geoinfo1.num_layer
   ele1 = geoinfo1.layer(i).elevation(ind1);
   [minele,ind] = min(abs(ele1 - ele2(:,1)));
   if minele < DIST && ~isnan(ele1) && ele2(ind,2) == 0 && ~isnan(ele2(ind,1))
       dist = ele2(ind,1) - ele1;
       try
       depth = geoinfo2.layer(ind).depth(ind2) - geoinfo1.layer(i).depth(ind1);
       catch
           depth;
       end
       geoinfo.layer(num_layer).elevation = [geoinfo1.layer(i).elevation(1:ind1) ...
           geoinfo2.layer(ind).elevation(ind2:end)-dist];
       geoinfo.layer(num_layer).echo_intensity = [geoinfo1.layer(i).echo_intensity(1:ind1) ...
           geoinfo2.layer(ind).echo_intensity(ind2:end)];
       geoinfo.layer(num_layer).depth = [geoinfo1.layer(i).depth(1:ind1) ...
           geoinfo2.layer(ind).depth(ind2:end)-depth];
       ele2(ind,2) = i;
       num_layer = num_layer + 1;
   else
       geoinfo.layer(num_layer).elevation = [geoinfo1.layer(i).elevation(1:ind1) ...
           nan*zeros(1,geoinfo2.num_trace-ind2+1)];
       geoinfo.layer(num_layer).echo_intensity = [geoinfo1.layer(i).echo_intensity(1:ind1) ...
           nan*zeros(1,geoinfo2.num_trace-ind2+1)];
       geoinfo.layer(num_layer).depth = [geoinfo1.layer(i).depth(1:ind1) ...
           nan*zeros(1,geoinfo2.num_trace-ind2+1)];
       num_layer = num_layer + 1;
   end
end
for j = 1:geoinfo2.num_layer
   if ele2(j,2) == 0      
       geoinfo.layer(num_layer).elevation = [nan*zeros(1,ind1) ...
           geoinfo2.layer(j).elevation(ind2:end)];
       geoinfo.layer(num_layer).echo_intensity = [nan*zeros(1,ind1) ...
           geoinfo2.layer(j).echo_intensity(ind2:end)];
       geoinfo.layer(num_layer).depth = [nan*zeros(1,ind1) ...
           geoinfo2.layer(j).depth(ind2:end)];
       num_layer = num_layer + 1;
   end
end   

geoinfo.num_layer = num_layer - 1;

nempty = 0;
% remove those layer with all nan data
i = 1;
while i <= geoinfo.num_layer - nempty 
elevation = geoinfo.layer(i).elevation;
elevation(isnan(elevation)) = [];
if isempty(elevation)
    geoinfo.layer(i) = [];
    nempty = nempty + 1;
else
    i = i + 1;
end
end
geoinfo.num_layer = geoinfo.num_layer - nempty;

% figure; imagesc(peaks); hold on;

ipeaks = zeros(length(geoinfo.time_range),length(geoinfo.elevation_bed));
for i = 1:geoinfo.num_layer
time_y = geoinfo.layer(i).depth*sqrt(epsilon)*2/C+ geoinfo.traveltime_surface;
y = nan*zeros(length(geoinfo.elevation_bed),1);
for j = 1:length(time_y)
[~,y(j)] = min(abs(geoinfo.time_range - time_y(j)));
if y(j) == 1
    y(j) = nan;
end
end
x = 1:length(geoinfo.elevation_surface);
x(isnan(y)) = [];
y(isnan(y)) = [];
indices = sub2ind(size(ipeaks),y,x');
ipeaks(indices) = i;
end
[newpeaks, label] = postprocesslayers(ipeaks);

geoinfo.layer(length(label)+1:end) = [];
geoinfo.num_layer = 0;
for i = 1:length(label)
    [y,x] = find(newpeaks == label(i));
    if isempty(x)
        continue;
    end
    geoinfo.num_layer = geoinfo.num_layer + 1;
%     indices = sub2ind(size(newpeaks),y,x);
    
    depth = ones(1,geoinfo.num_trace)*nan;
    elevation = depth;
    depth(x) = (geoinfo.time_range(y) - geoinfo.traveltime_surface(x))*C/2/sqrt(epsilon);
    elevation(x) = geoinfo.elevation_surface(x) - depth(x);
    geoinfo.layer(geoinfo.num_layer).depth = depth;
    geoinfo.layer(geoinfo.num_layer).elevation = elevation;
%     geoinfo1.layer(i)=struct(field1,depth,field2,echo_intensity,field3,elevation);
end



end