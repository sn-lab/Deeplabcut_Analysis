function analyze_objectinteraction(h5_filename)
% FUNCTION analyze_objectinteraction(h5_filename)
%
% analyzes h5 files output from deeplabcut tracking of object recognition
% videos with multiple objects placed in the red box, using the 
% "Redbox_2" training dataset. Analysis includes plotting the summary of 
% mouse and box tracking data, and calculating the time the mouse spent 
% "interacting" with each detected object. In order to detect objects, this 
% script opens the original video file in addition to using the deeplabcut 
% tracking results h5 file. 
%
% Note: this script assumes that the original ideo file is in AVI format,
% that it is in the same folder as the h5 file, and that it shares the same
% base filename with the h5 file. For example: 
% C:/video1.avi
% C:/video1DLC_mobnet_100_Redbox_2Oct13shuffle1_10000.h5
%
%INPUTS
%h5_filename: full path/filename of a h5 file to analyze

%set parameters
fps = 30;
nearobject_threshold = 0.04; %near object positions
verynearobject_threshold = 0.02; %very near object positions
box_length_in_meters = 0.4445;
figure_size = [10 10 25 10];

%load tracking data
[save_dir, filename, ~] = fileparts(h5_filename);
ind = strfind(filename,'DLC_mobnet');
filename = filename(1:ind-1);
output_filename = fullfile(save_dir,filename);

data = h5read(h5_filename,'/df_with_missing/table');
ary = data.values_block_0';
%columns: top-right, top-left, bottom-left, bottome-right, right-ear, left-ear, nose, tail-base; each have 3 rows, for [x, y, likelihood]
colnames = {'TR_x','TR_y','TR_l','TL_x','TL_y','TL_l','BL_x','BL_y','BL_l','BR_x','BR_y','BR_l','RE_x','RE_y','RE_l','LE_x','LE_y','LE_l','NO_x','NO_y','NO_l','TB_x','TB_y','TB_l'};
[num_frames, num_cols] = size(ary);

%vertically flip all y values
ystrings = repmat({'_y'},1,num_cols);
col_is_y = cellfun(@strfind,colnames,ystrings,'UniformOutput',false);
col_is_y = find(~cellfun(@isempty,col_is_y));
ary(:,col_is_y) = 1080-ary(:,col_is_y);

%load video subset (~100 frames) for object detection
avi_filename = [output_filename '.avi'];
assert(exist(avi_filename,'file')==2,['Cannot find avi file associated with the h5 file. ',...
    'Make sure the original video file is in the same folder as the h5 file, and shares',...
    ' the same filename before the "DLC_mobnet_..." (e.g. video1.avi and video1DLC_mobnet_100....h5']);
v = VideoReader(avi_filename);
time_incr = v.Duration/100;
sample_frame_times = 0:time_incr:v.Duration;
sample_frames = nan([v.height, v.width, 3, length(sample_frame_times)]);
switch v.VideoFormat
    case 'Mono8'
        scale_factor = 1/255;
    case 'RGB24'
        scale_factor = 1/255;
    case 'RGB48'
        scale_factor = 1/65535;
    otherwise
        scale_factor = 1/255;
        warning('unexpected video format, object detection may not work');
end
for t = 1:length(sample_frame_times)
    v.CurrentTime = sample_frame_times(t);
    if hasFrame(v)
        sample_frames(:,:,:,t) = readFrame(v);
    end
end
sample_frames = flipud(sample_frames);
median_frame = median(sample_frames,4,'omitnan')*scale_factor;


%% calculate box corner locations
topright_x_ind = find(strcmp(colnames,'TR_x'));
topright_y_ind = find(strcmp(colnames,'TR_y'));
topleft_x_ind = find(strcmp(colnames,'TL_x'));
topleft_y_ind = find(strcmp(colnames,'TL_y'));
botleft_x_ind = find(strcmp(colnames,'BL_x'));
botleft_y_ind = find(strcmp(colnames,'BL_y'));
botright_x_ind = find(strcmp(colnames,'BR_x'));
botright_y_ind = find(strcmp(colnames,'BR_y'));

topright_x = ary(:,topright_x_ind);
topright_y = ary(:,topright_y_ind);
topright_inds = topright_x>800 & topright_x<1440 & topright_y>540 & topright_y<1080;
topright_x = topright_x(topright_inds);
topright_x = median(topright_x,'omitnan');
topright_y = topright_y(topright_inds);
topright_y = median(topright_y,'omitnan');

topleft_x = ary(:,topleft_x_ind);
topleft_y = ary(:,topleft_y_ind);
topleft_inds = topleft_x<800 & topleft_x>0 & topleft_y>540 & topleft_y<1080;
topleft_x = topleft_x(topleft_inds);
topleft_x = median(topleft_x,'omitnan');
topleft_y = topleft_y(topleft_inds);
topleft_y = median(topleft_y,'omitnan');

botleft_x = ary(:,botleft_x_ind);
botleft_y = ary(:,botleft_y_ind);
botleft_inds = botleft_x<800 & botleft_x>0 & botleft_y<540 & botleft_y>0;
botleft_x = botleft_x(botleft_inds);
botleft_x = median(botleft_x,'omitnan');
botleft_y = botleft_y(botleft_inds);
botleft_y = median(botleft_y,'omitnan');

botright_x = ary(:,botright_x_ind);
botright_y = ary(:,botright_y_ind);
botright_inds = botright_x>800 & botright_x<1440 & botright_y<540 & botright_y>0;
botright_x = botright_x(botright_inds);
botright_x = median(botright_x,'omitnan');
botright_y = botright_y(botright_inds);
botright_y = median(botright_y,'omitnan');

corners_x = [topright_x topleft_x botleft_x botright_x topright_x];
corners_y = [topright_y topleft_y botleft_y botright_y topright_y];

box_lengths(1) = get_dist(topright_x,topright_y,topleft_x,topleft_y);
box_lengths(2) = get_dist(topleft_x,topleft_y,botleft_x,botleft_y);
box_lengths(3) = get_dist(botleft_x,botleft_y,botright_x,botright_y);
box_lengths(4) = get_dist(botright_x,botright_y,topright_x,topright_y);
box_length_in_pixels = mean(box_lengths,'omitnan');
assert(all((box_lengths./box_length_in_pixels)>0.95),'detected box is out of shape - tracking is probably too poor to continue')
pixels_per_meter = box_length_in_pixels/box_length_in_meters;

%draw box for validation
figure('Units','centimeters','Position',figure_size);
ax = subplot(2,4,1);
image(median_frame);
set(ax,'YDir','normal');
xlim([-100 1550]);
ylim([-280 1370])
title('median frame')

%calculate some points along the box edge (for later, to track mice distance from edges)
num_edgepts = round(box_length_in_pixels/10);

topedge_xstep = (topleft_x-topright_x)/(num_edgepts-1);
leftedge_xstep = (botleft_x-topleft_x)/(num_edgepts-1);
botedge_xstep = (botright_x-botleft_x)/(num_edgepts-1);
rightedge_xstep = (topright_x-botright_x)/(num_edgepts-1);
edgepts_x = [topright_x:topedge_xstep:topleft_x, topleft_x:leftedge_xstep:botleft_x, botleft_x:botedge_xstep:botright_x, botright_x:rightedge_xstep:topright_x];

topedge_ystep = (topleft_y-topright_y)/(num_edgepts-1);
leftedge_ystep = (botleft_y-topleft_y)/(num_edgepts-1);
botedge_ystep = (botright_y-botleft_y)/(num_edgepts-1);
rightedge_ystep = (topright_y-botright_y)/(num_edgepts-1);
edgepts_y = [topright_y:topedge_ystep:topleft_y, topleft_y:leftedge_ystep:botleft_y, botleft_y:botedge_ystep:botright_y, botright_y:rightedge_ystep:topright_y];

%side 1 = top, side 2 = left, side 3 = bottom, side 4 = right
side_ind = [1*ones(1, num_edgepts) 2*ones(1, num_edgepts) 3*ones(1, num_edgepts) 4*ones(1, num_edgepts)];
num_edgepts = length(edgepts_x);


%% detect object locations (hue and dark-value method)
min_object_size = 500; %in pixels
object_filt = 20;

%convert median frame to hsv (h and v may be most useful to find objects, since the background is hue-stable and not dark)
within_box_mask = poly2mask(corners_x,corners_y,1080,1440);
median_frame_hsv = rgb2hsv(median_frame);
median_frame_hue = 360*median_frame_hsv(:,:,1);
median_frame_hue(median_frame_hue>180) = median_frame_hue(median_frame_hue>180)-360; %scaled as a difference from the red hue
median_frame_v = median_frame_hsv(:,:,3);

%calculate the median hue of the background, and calculate the every pixel's hue difference from that median
full_bkgd_mask = poly2mask(corners_x,corners_y,size(median_frame_hue,1),size(median_frame_hue,2));
full_bkgd_hue = median_frame_hue(full_bkgd_mask);
median_bkgd_hue = median(full_bkgd_hue,'omitnan');
median_frame_hue = median_frame_hue-median_bkgd_hue; %scaled as a difference from the background hue (probably close to red)
median_frame_hue(median_frame_hue>180) = median_frame_hue(median_frame_hue>180)-360;
median_frame_hue(median_frame_hue<-180) = median_frame_hue(median_frame_hue<-180)+360;
median_frame_hue_abs = abs(median_frame_hue);

%using the smallest 90% difference hues (confident these are background pixels, since objects are small), 
%calculate the S.D. of the background hue and convert hue differences to hue z-scores
hue_80th = prctile(median_frame_hue_abs(full_bkgd_mask),80);
bkgd_mask = median_frame_hue_abs<=hue_80th & full_bkgd_mask;
bkgd_hue_mean = mean(median_frame_hue_abs(bkgd_mask),'omitnan');
bkgd_hue_std = std(median_frame_hue_abs(bkgd_mask),'omitnan');
median_frame_hue_zscr = full_bkgd_mask.*abs(median_frame_hue_abs-bkgd_hue_mean)/bkgd_hue_std;

%using the pixels that are confidently background by hue, calculate the
%background v range median and S.D., to convert v values to z-scores
bkgd_v_mean = mean(median_frame_v(bkgd_mask),'omitnan');
bkgd_v_std = std(median_frame_v(bkgd_mask),'omitnan');
median_frame_vdark_zscr = full_bkgd_mask.*(bkgd_v_mean-median_frame_v)/bkgd_v_std;
median_frame_vdark_zscr(median_frame_vdark_zscr<0) = 0;

%add z-scores together, median filter to smooth pixels, then threshold to
%find objects mask
median_frame_comb_zscr = (2*median_frame_vdark_zscr)+median_frame_hue_zscr;
median_frame_cz_filt = medfilt2(median_frame_comb_zscr,[object_filt object_filt]);
objects_mask = median_frame_cz_filt>10;

%find connected components, remove all that are too small, find connected
%components again
objects_label = bwlabel(objects_mask); 
num_labels = max(max(objects_label));
for i=1:num_labels
    label_size = sum(sum(objects_label==i));
    if label_size<min_object_size
        objects_label(objects_label==i) = 0;
    end
end
objects_label = bwlabel(logical(objects_label));
num_objects = max(max(objects_label));
assert(num_objects>0,'No objects detected');
assert(num_objects<5,'Too many objects detected, video quality may be too poor');

%calculate object center location
objects_x = nan(1,num_objects);
objects_y = nan(1,num_objects);
for o = 1:num_objects
    object_mask = objects_label==o; %current object mask
    %calculate object center
    pixels_x = repmat(1:v.width,[v.height 1]);
    pixels_y = repmat((1:v.height)',[1 v.width]);
    objects_x(o) = mean(pixels_x(object_mask),'omitnan');
    objects_y(o) = mean(pixels_y(object_mask),'omitnan');
end

%make sure objects are sorted roughly left-to-right, top-to-bottom
order_val = (objects_x-topleft_x) + 2*(topleft_y-objects_y); %give extra weight to top-to-bottom, so left-right ordering takes priority
[~, order] = sort(order_val);
objects_x = objects_x(order);
objects_y = objects_y(order);
for o = 1:num_objects
    objects_label(objects_label==order(o)) = -o; %re-order labels with unused negative numbers
end
objects_label = -objects_label; %make labels positive again

%calculate each object's size, outline, and color
objects_radius = nan(1,num_objects);
objects_edgepts_x = nan(num_objects,100); %pre-allocating plenty of space
objects_edgepts_y = nan(num_objects,100);
colors = linspecer(4);
objects_by_color = repmat(objects_label,[1 1 3]);
for o = 1:num_objects
    object_mask = objects_label==o; %current object mask
    objects_radius(o) = sqrt(sum(sum(object_mask,'omitnan'),2,'omitnan')/pi);
    
    %calculate object center
    pixels_x = repmat(1:v.width,[v.height 1]);
    pixels_y = repmat((1:v.height)',[1 v.width]);
    objects_x(o) = mean(pixels_x(object_mask),'omitnan');
    objects_y(o) = mean(pixels_y(object_mask),'omitnan');
    
    %draw object outline
    vertical = any(object_mask,2); %rows that have object pixels
    bottom = find(vertical,1,'first'); %first (bottom) row that has object pixels
    top = find(vertical,1,'last'); %last (top) row that ""
    cur_object_y = [bottom:5:top-3 top]; %sampling points from bottom to top
    num_sidepts = length(cur_object_y);
    cur_object_y = [cur_object_y fliplr(cur_object_y) cur_object_y(1)]; %sampling points from bottom to top, to bottom
    num_object_edgepts = 1+(num_sidepts*2);
    cur_object_x = nan([1 num_object_edgepts]);
    for i = 1:num_sidepts
        cur_object_x(i) = find(object_mask(cur_object_y(i),:),1,'first'); %left side x
        cur_object_x(i+num_sidepts) = find(object_mask(cur_object_y(i+num_sidepts),:),1,'last'); %right side x
    end
    cur_object_x(end) = cur_object_x(1);
    objects_edgepts_x(o,1:num_object_edgepts) = cur_object_x;
    objects_edgepts_y(o,1:num_object_edgepts) = cur_object_y;
    
    %add the objects color to an RGB image
    for c = 1:3
        tmp = objects_by_color(:,:,c);
        tmp(tmp==o) = colors(o,c);
        objects_by_color(:,:,c) = tmp;
    end
end

figure(); %new "object detection details" figure
ax = subplot(3,3,1);
imshow(median_frame)
set(ax,'YDir','normal');
title('median frame')
ax = subplot(3,3,2);
imshow(full_bkgd_mask)
set(ax,'YDir','normal');
title('within-box mask')
ax = subplot(3,3,3);
imshow(bkgd_mask)
set(ax,'YDir','normal');
title('80th percentile hue mask')
ax = subplot(3,3,4);
imshow(median_frame_hue_zscr/20)
set(ax,'YDir','normal');
title('hue z-score')
ax = subplot(3,3,5);
imshow(median_frame_vdark_zscr/10)
set(ax,'YDir','normal');
title('(darker) value z-score')
ax = subplot(3,3,6);
imshow(median_frame_cz_filt/30)
set(ax,'YDir','normal');
title('combined z-score, median filtered')
ax = subplot(3,3,7);
imshow(objects_mask)
set(ax,'YDir','normal');
title('objects mask')
ax = subplot(3,3,8);
imshow(objects_by_color)
set(ax,'YDir','normal');
title('detected objects')
ax = subplot(3,3,9);
imshow(median_frame/2)
set(ax,'YDir','normal');
hold on
title('object centers and outlines')
for o = 1:num_objects
    plot(objects_edgepts_x(o,:),objects_edgepts_y(o,:),'Color',colors(o,:),'LineWidth',1)
    scatter(objects_x(o),objects_y(o),'kx','MarkerEdgeColor',colors(o,:))
    text(objects_x(o)+(1.5*objects_radius(o)),objects_y(o)+(1.5*objects_radius(o)),num2str(o),'Color',colors(o,:))
end
saveas(gcf,[output_filename '_DetectionDetails.fig'])
close(gcf)

%% incorporate object detection into summary plot
figure(gcf); %current figure is the summary figure again
subplot(2,4,2)
plot(corners_x,corners_y);
hold on
for o = 1:num_objects
    plot(objects_edgepts_x(o,:),objects_edgepts_y(o,:),'Color',colors(o,:),'LineWidth',1)
    text(objects_x(o)+(1.5*objects_radius(o)),objects_y(o)+(1.5*objects_radius(o)),num2str(o),'Color',colors(o,:))
    xlim([-100 1550]);
    ylim([-280 1370])
    title('detected walls and object locations');
end

%% delete all mouse tracking that's outside of the box (change to NaNs)
leftear_x_ind = find(strcmp(colnames,'LE_x'));
leftear_y_ind = find(strcmp(colnames,'LE_y'));
rightear_x_ind = find(strcmp(colnames,'RE_x'));
rightear_y_ind = find(strcmp(colnames,'RE_y'));
nose_x_ind = find(strcmp(colnames,'NO_x'));
nose_y_ind = find(strcmp(colnames,'NO_y'));
tail_x_ind = find(strcmp(colnames,'TB_x'));
tail_y_ind = find(strcmp(colnames,'TB_y'));

full_edgepts_x = repmat(edgepts_x,[num_frames 1]);
full_edgepts_y = repmat(edgepts_y,[num_frames 1]);

nose_x = ary(:,nose_x_ind);
nose_y = ary(:,nose_y_ind);
edge_dists = get_dist(repmat(nose_x,[1 num_edgepts]),repmat(nose_y,[1 num_edgepts]),full_edgepts_x,full_edgepts_y);
[~, closest_edgept_inds] = min(edge_dists,[],2,'omitnan');
closest_edgepts_x = full_edgepts_x(sub2ind(repmat([num_frames num_edgepts],[num_frames 1]),(1:num_frames)', closest_edgept_inds));
closest_edgepts_y = full_edgepts_y(sub2ind(repmat([num_frames num_edgepts],[num_frames 1]),(1:num_frames)', closest_edgept_inds));
closest_sides = side_ind(closest_edgept_inds)';
out_of_box_idx = false(size(nose_x));
out_of_box_idx(closest_sides==1) = nose_y(closest_sides==1)>closest_edgepts_y(closest_sides==1);
out_of_box_idx(closest_sides==2) = nose_x(closest_sides==2)<closest_edgepts_x(closest_sides==2);
out_of_box_idx(closest_sides==3) = nose_y(closest_sides==3)<closest_edgepts_y(closest_sides==3);
out_of_box_idx(closest_sides==4) = nose_x(closest_sides==4)>closest_edgepts_x(closest_sides==4);
nose_x(out_of_box_idx) = nan;
nose_y(out_of_box_idx) = nan;

leftear_x = ary(:,leftear_x_ind);
leftear_y = ary(:,leftear_y_ind);
edge_dists = get_dist(repmat(leftear_x,[1 num_edgepts]),repmat(leftear_y,[1 num_edgepts]),full_edgepts_x,full_edgepts_y);
[~, closest_edgept_inds] = min(edge_dists,[],2,'omitnan');
closest_edgepts_x = full_edgepts_x(sub2ind(repmat([num_frames num_edgepts],[num_frames 1]),(1:num_frames)', closest_edgept_inds));
closest_edgepts_y = full_edgepts_y(sub2ind(repmat([num_frames num_edgepts],[num_frames 1]),(1:num_frames)', closest_edgept_inds));
closest_sides = side_ind(closest_edgept_inds)';
out_of_box_idx = false(size(nose_x));
out_of_box_idx(closest_sides==1) = leftear_y(closest_sides==1)>closest_edgepts_y(closest_sides==1);
out_of_box_idx(closest_sides==2) = leftear_x(closest_sides==2)<closest_edgepts_x(closest_sides==2);
out_of_box_idx(closest_sides==3) = leftear_y(closest_sides==3)<closest_edgepts_y(closest_sides==3);
out_of_box_idx(closest_sides==4) = leftear_x(closest_sides==4)>closest_edgepts_x(closest_sides==4);
leftear_x(out_of_box_idx) = nan;
leftear_y(out_of_box_idx) = nan;

rightear_x = ary(:,rightear_x_ind);
rightear_y = ary(:,rightear_y_ind);
edge_dists = get_dist(repmat(rightear_x,[1 num_edgepts]),repmat(rightear_y,[1 num_edgepts]),full_edgepts_x,full_edgepts_y);
[~, closest_edgept_inds] = min(edge_dists,[],2,'omitnan');
closest_edgepts_x = full_edgepts_x(sub2ind(repmat([num_frames num_edgepts],[num_frames 1]),(1:num_frames)', closest_edgept_inds));
closest_edgepts_y = full_edgepts_y(sub2ind(repmat([num_frames num_edgepts],[num_frames 1]),(1:num_frames)', closest_edgept_inds));
closest_sides = side_ind(closest_edgept_inds)';
out_of_box_idx = false(size(nose_x));
out_of_box_idx(closest_sides==1) = rightear_y(closest_sides==1)>closest_edgepts_y(closest_sides==1);
out_of_box_idx(closest_sides==2) = rightear_x(closest_sides==2)<closest_edgepts_x(closest_sides==2);
out_of_box_idx(closest_sides==3) = rightear_y(closest_sides==3)<closest_edgepts_y(closest_sides==3);
out_of_box_idx(closest_sides==4) = rightear_x(closest_sides==4)>closest_edgepts_x(closest_sides==4);
rightear_x(out_of_box_idx) = nan;
rightear_y(out_of_box_idx) = nan;

tail_x = ary(:,tail_x_ind);
tail_y = ary(:,tail_y_ind);
edge_dists = get_dist(repmat(tail_x,[1 num_edgepts]),repmat(tail_y,[1 num_edgepts]),full_edgepts_x,full_edgepts_y);
[~, closest_edgept_inds] = min(edge_dists,[],2,'omitnan');
closest_edgepts_x = full_edgepts_x(sub2ind(repmat([num_frames num_edgepts],[num_frames 1]),(1:num_frames)', closest_edgept_inds));
closest_edgepts_y = full_edgepts_y(sub2ind(repmat([num_frames num_edgepts],[num_frames 1]),(1:num_frames)', closest_edgept_inds));
closest_sides = side_ind(closest_edgept_inds)';
out_of_box_idx = false(size(nose_x));
out_of_box_idx(closest_sides==1) = tail_y(closest_sides==1)>closest_edgepts_y(closest_sides==1);
out_of_box_idx(closest_sides==2) = tail_x(closest_sides==2)<closest_edgepts_x(closest_sides==2);
out_of_box_idx(closest_sides==3) = tail_y(closest_sides==3)<closest_edgepts_y(closest_sides==3);
out_of_box_idx(closest_sides==4) = tail_x(closest_sides==4)>closest_edgepts_x(closest_sides==4);
tail_x(out_of_box_idx) = nan;
tail_y(out_of_box_idx) = nan;

%delete all mouse tracking that's too far from the mouse's estimated position
%estimate the mouse's size
mouse_size_in_pixels = median(get_dist(nose_x,nose_y,tail_x,tail_y),'omitnan');
% mouse_size_in_meters = mouse_size_in_pixels/pixels_per_meter;

%create a time-smoothed, body-center estimate
head_x = median([nose_x leftear_x rightear_x],2,'omitnan');
head_y = median([nose_y leftear_y rightear_y],2,'omitnan');
body_x =  mean([head_x tail_x],2,'omitnan');
body_y =  mean([head_y tail_y],2,'omitnan');
smooth_body_x = rolling_average(body_x,1,round(fps/2),'median');
smooth_body_y = rolling_average(body_y,1,round(fps/2),'median');

%find all mouse labels that are very far from the body center estimate, and nan them
dist_tmp = get_dist(nose_x,nose_y,smooth_body_x,smooth_body_y);
too_far_idx = dist_tmp>mouse_size_in_pixels*3;
nose_x(too_far_idx) = nan;
nose_y(too_far_idx) = nan;

dist_tmp = get_dist(leftear_x,leftear_y,smooth_body_x,smooth_body_y);
too_far_idx = dist_tmp>mouse_size_in_pixels*3;
leftear_x(too_far_idx) = nan;
leftear_y(too_far_idx) = nan;

dist_tmp = get_dist(rightear_x,rightear_y,smooth_body_x,smooth_body_y);
too_far_idx = dist_tmp>mouse_size_in_pixels*3;
rightear_x(too_far_idx) = nan;
rightear_y(too_far_idx) = nan;

dist_tmp = get_dist(tail_x,tail_y,smooth_body_x,smooth_body_y);
too_far_idx = dist_tmp>mouse_size_in_pixels*3;
tail_x(too_far_idx) = nan;
tail_y(too_far_idx) = nan;


%% estimate the mouse's body position again
smooth_tail_x = rolling_average(tail_x,1,round(fps/2),'median');
smooth_tail_y = rolling_average(tail_y,1,round(fps/2),'median');
smooth_nose_x = rolling_average(nose_x,1,round(fps/4),'median'); %head moves faster, so smooth less
smooth_nose_y = rolling_average(nose_y,1,round(fps/4),'median');
smooth_rightear_x = rolling_average(rightear_x,1,round(fps/4),'median');
smooth_rightear_y = rolling_average(rightear_y,1,round(fps/4),'median');
smooth_leftear_x = rolling_average(leftear_x,1,round(fps/4),'median');
smooth_leftear_y = rolling_average(leftear_y,1,round(fps/4),'median');
smooth_head_x = median([smooth_nose_x smooth_leftear_x smooth_rightear_x],2,'omitnan');
smooth_head_y = median([smooth_nose_y smooth_leftear_y smooth_rightear_y],2,'omitnan');
smooth_body_x =  mean([smooth_head_x smooth_tail_x],2,'omitnan');
smooth_body_y =  mean([smooth_head_y smooth_tail_y],2,'omitnan');

subplot(2,4,5)
plot(corners_x,corners_y);
hold on
scatter(smooth_body_x,smooth_body_y,1,[0 0 1]);%colors)
xlim([-100 1550]);
ylim([-280 1370])
title('all mouse body positions')

subplot(2,4,6)
plot(corners_x,corners_y);
hold on
xlim([-100 1550]);
ylim([-280 1370])
title('first 30 s trajectory')
plot(smooth_body_x(1:fps*30),smooth_body_y(1:fps*30),'b','LineWidth',0.5);

%calculate head direction
headang_ears = atan2((smooth_rightear_y-smooth_leftear_y),(smooth_rightear_x-smooth_leftear_x))+pi/2;
headang_lnose = atan2((smooth_nose_y-smooth_leftear_y),(smooth_nose_x-smooth_leftear_x))+pi/4;
headang_rnose = atan2((smooth_nose_y-smooth_rightear_y),(smooth_nose_x-smooth_rightear_x))-pi/4;

%get vector sum head direction (giving more weight to ears)
vs_x = sum([2*cos(headang_ears) cos(headang_lnose) cos(headang_rnose)],2,'omitnan');
vs_y = sum([2*sin(headang_ears) sin(headang_lnose) sin(headang_rnose)],2,'omitnan');

%get smoother vector sum head angle
smooth_vs_x = rolling_average(vs_x,1,round(fps/4),'median');
smooth_vs_y = rolling_average(vs_y,1,round(fps/4),'median');
smooth_headang = atan2(smooth_vs_y,smooth_vs_x);


%% analyze results for detected objects
%calculate direction of all objects to mouse head
objects_ang = atan2((repmat(objects_y,[num_frames 1])-repmat(smooth_head_y,[1 num_objects])),(repmat(objects_x,[num_frames 1])-repmat(smooth_head_x,[1 num_objects])));
objects_head_ang = get_angular_dist(objects_ang,repmat(smooth_headang,[1 num_objects]),'radians');
angled_toward_objects = objects_head_ang<(pi/2);

%calculate distances of mouse body to object edges
full_obj_edge_x = repmat(permute(objects_edgepts_x,[3 1 2]),[num_frames 1 1]); %new format is: [num_frames num_objects num_object_edgepts]
full_obj_edge_y = repmat(permute(objects_edgepts_y,[3 1 2]),[num_frames 1 1]);
full_head_x = repmat(smooth_head_x,[1 num_objects 100]); %there are 100 possible object edgepoints (from pre-allocation)
full_head_y = repmat(smooth_head_y,[1 num_objects 100]);
objects_dist = get_dist(full_head_x,full_head_y, full_obj_edge_x, full_obj_edge_y);
objects_dist = min(objects_dist,[],3,'omitnan'); %find closest edgepoint
nearobjects = (objects_dist/pixels_per_meter)<nearobject_threshold; %format is [num_frames num_objects]
verynearobjects = (objects_dist/pixels_per_meter)<verynearobject_threshold;

%calculate near and verynear interactions (orientation towards + near or verynear)
objects_nearinteract = nearobjects & angled_toward_objects;
objects_verynearinteract = verynearobjects & angled_toward_objects;
objects_nearinteraction_time = sum(objects_nearinteract)/fps;
objects_verynearinteraction_time = sum(objects_verynearinteract)/fps; 


%% calculate summary statistics
steps = get_dist(smooth_body_x(1:end-1),smooth_body_y(1:end-1),smooth_body_x(2:end),smooth_body_y(2:end));
distance_travelled = sum(abs(steps))/pixels_per_meter;

results.number_of_frames = num_frames;
results.assumed_camera_fps = fps;
results.total_time_in_seconds = num_frames/fps;
results.distance_travelled_in_meters = distance_travelled;
results.number_of_objects = num_objects;
results.object_interaction_time_4cm = objects_nearinteraction_time;
results.object_interaction_time_2cm = objects_verynearinteraction_time;

%print results to figure
subplot(2,4,3)
set(gca,'Units','Normalized','Position',[0.6 0.05 0.4 0.95])
axis off
text(0,0.9,['Number of frames: ' num2str(results.number_of_frames)])
text(0,0.82,['Assumed camera fps: ' num2str(results.assumed_camera_fps)])
text(0,0.74,['Total time: ' num2str(results.total_time_in_seconds, '%.1f') ' s'])
text(0,0.66,['Distance travelled: ' num2str(results.distance_travelled_in_meters, '%.1f') ' m'])
text(0,0.58,'Object locations: ')
text(0,0.53,['[' num2str(1:num_objects) ']'],'FontSize',8)
text(0,0.45,'Time interacting with objects - within 4 cm: ')
text(0,0.40,['[' num2str(results.object_interaction_time_4cm, '%.1f  ') '] s,   [' num2str(round(1000*results.object_interaction_time_4cm/results.total_time_in_seconds)/10,'%.1f  ') '] %']);
text(0,0.32,'Time interacting with objects - within 2 cm: ')
text(0,0.27,['[' num2str(results.object_interaction_time_2cm, '%.1f  ') '] s,   [' num2str(round(1000*results.object_interaction_time_2cm/results.total_time_in_seconds)/10, '%.1f  ') '] %']);

saveas(F1,[output_filename '.png'])
save([output_filename '.mat'],'results');

%write custom xls file of results
C = cell(4,2);
C(:,1) = {'number of frames', 'assumed camera fps', 'total time (s)', 'distance travelled (m)'};
C(:,2) = {num_frames, fps, results.total_time_in_seconds, distance_travelled};
writecell(C,[output_filename '.xls'],'Range','A1');
C = cell(5,1+num_objects);
C(:,1) = {'object location', 'Time interacting with objects - within 4 cm (s)', ...
    'Time interacting with objects - within 4 cm (%)', 'Time interacting with objects - within 2 cm (s)', ...
    'Time interacting with objects - within 2 cm (%)'};
C(1,2:1+num_objects) = num2cell(1:num_objects);
C(2,2:1+num_objects) = num2cell(results.object_interaction_time_4cm);
C(3,2:1+num_objects) = num2cell(100*results.object_interaction_time_4cm/results.total_time_in_seconds);
C(4,2:1+num_objects) = num2cell(results.object_interaction_time_2cm);
C(5,2:1+num_objects) = num2cell(100*results.object_interaction_time_2cm/results.total_time_in_seconds);
writecell(C,[output_filename '.xls'],'Range','A6');
