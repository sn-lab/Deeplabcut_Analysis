function analyze_objectinteraction(h5_filename,param)
% FUNCTION analyze_objectinteraction(h5_filename,param)
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
%param (optional):
%     fps = framerate of the original video recordings
%     down_fps = effective fps after rolling average performed to remove some tracking errors
%     bodypart = bodypart to use for distance measurements (e.g. mouse-object distance)
%     view_cone = viewing cone angle (i.e. horizontal FOV, in degrees) to consider an object inside the mouse's visual field
%     close_threshold = minimum distance (cm) from mouse to object to count as interacting
%     show_figure = whether to show (1) or close (0) the results figure when finished
%     hue_weight = %the weight of hue (vs dark-luminance) in object detection
%     bkgd_fraction = fraction of box background to use as "non-object pixels"
%     zscore_threshold = z-score threshold to define object pixels
%     min_obj_size = minimum size (in cm^2) to identify an object
%     med_filt = size of median filter, in cm

if nargin<2
    %set default parameters
    fps = 30; %framerate of the original video recordings
    down_fps = 10; %effective fps after rolling average performed to remove some tracking errors
    bodypart = 'head'; %bodypart to use for distance measurements (e.g. mouse-object distance)
    view_cone = 90; %viewing cone angle (i.e. horizontal FOV, in degrees) to consider an object inside the mouse's visual field
    close_threshold = 4; %minimum distance (cm) from mouse to object to count as interacting
    show_figure = 1; %whether to show (1) or close (0) the results figure when finished
    hue_weight = 0.3; %the weight of hue (vs dark-luminance) in object detection
    bkgd_fraction = 0.85; %fraction of box background to use as "non-object pixels"
    zscore_threshold = 3; %z-score threshold to define object pixels
    min_obj_size = 5; %minimum size (in cm^2) to identify an object
    med_filt = 1; %size of median filter, in cm
else
    fps = param.fps;
    down_fps = param.down_fps;
    bodypart = param.bodypart;
    view_cone = param.view_cone;
    close_threshold = param.close_threshold;
    show_figure = param.show_figure;
    hue_weight = param.hue_weight;
    bkgd_fraction = param.bkgd_fraction;
    zscore_threshold = param.zscore_threshold;
    min_obj_size = param.min_obj_size;
    med_filt = param.med_filt;
end

box_length_in_meters = 0.4445;
figure_size = [10 10 25 10];

%load tracking data
[save_dir, filename, ~] = fileparts(h5_filename);
ind = strfind(filename,'DeepCut');
if isempty(ind)
    ind = strfind(filename,'DLC');
end
if isempty(ind)
    ind = strfind(filename,'.h5');
end
filename = filename(1:ind-1);
output_filename = fullfile(save_dir,filename);
data = h5read(h5_filename,'/df_with_missing/table');
ary = data.values_block_0';
%columns: top-right, top-left, bottom-left, bottom-right, right-ear, left-ear, nose, tail-base; each have 3 columns, for [x, y, likelihood]
colnames = {'TR_x','TR_y','TR_l','TL_x','TL_y','TL_l','BL_x','BL_y','BL_l','BR_x','BR_y','BR_l','RE_x','RE_y','RE_l','LE_x','LE_y','LE_l','NO_x','NO_y','NO_l','TB_x','TB_y','TB_l'};
[num_frames, num_cols] = size(ary);

%load video information
avi_filename = [output_filename '.avi'];
assert(exist(avi_filename,'file')==2,['Cannot find avi file associated with the h5 file. ',...
    'Make sure the original video file is in the same folder as the h5 file, and shares',...
    ' the same filename before "DeepCut" (e.g. video1.avi and video1DeepCutAp....h5']);
v = VideoReader(avi_filename);

%vertically flip all y values
ystrings = repmat({'_y'},1,num_cols);
col_is_y = cellfun(@strfind,colnames,ystrings,'UniformOutput',false);
col_is_y = find(~cellfun(@isempty,col_is_y));
ary(:,col_is_y) = v.height-ary(:,col_is_y);

%get image/plot limits
largest_range = max([v.width v.height]);
xlimits = [(v.width-largest_range) (v.width+largest_range)]/2;
ylimits = [(v.height-largest_range) (v.height+largest_range)]/2;

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
topright_inds = topright_x>(v.width/2) & topright_x<v.width & topright_y>(v.height/2) & topright_y<v.height;
topright_x = topright_x(topright_inds);
topright_x = median(topright_x,'omitnan');
topright_y = topright_y(topright_inds);
topright_y = median(topright_y,'omitnan');

topleft_x = ary(:,topleft_x_ind);
topleft_y = ary(:,topleft_y_ind);
topleft_inds = topleft_x<(v.width/2) & topleft_x>0 & topleft_y>(v.height/2) & topleft_y<v.height;
topleft_x = topleft_x(topleft_inds);
topleft_x = median(topleft_x,'omitnan');
topleft_y = topleft_y(topleft_inds);
topleft_y = median(topleft_y,'omitnan');

botleft_x = ary(:,botleft_x_ind);
botleft_y = ary(:,botleft_y_ind);
botleft_inds = botleft_x<(v.width/2) & botleft_x>0 & botleft_y<(v.height/2) & botleft_y>0;
botleft_x = botleft_x(botleft_inds);
botleft_x = median(botleft_x,'omitnan');
botleft_y = botleft_y(botleft_inds);
botleft_y = median(botleft_y,'omitnan');

botright_x = ary(:,botright_x_ind);
botright_y = ary(:,botright_y_ind);
botright_inds = botright_x>(v.width/2) & botright_x<v.width & botright_y<(v.height/2) & botright_y>0;
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
while any((box_lengths./box_length_in_pixels)<0.95 | (box_lengths./box_length_in_pixels)>1.05)
    %find the worst angle and replace that corner
    warning('box appears to be out of shape -- tracking quality may be poor')
    coords = [corners_x' corners_y'];
    corner_angles = nan(1,4);
    for c = 1:4
        c1 = coords(1+mod(c-2,4),:) - coords(c,:);
        angle1 = atan2(c1(2),c1(1));
        c2 = coords(1+mod(c,4),:) - coords(c,:);
        angle2 = atan2(c2(2),c2(1));
        corner_angles(c) = get_angular_dist(angle1,angle2,'radians');
    end
    [~,bad_corner] = max(abs(corner_angles-(pi/2)));
    good_corners = 1+mod([bad_corner-1, bad_corner+1]-1,4);
    diagonal_len = get_dist(corners_x(good_corners(1)),corners_y(good_corners(1)),corners_x(good_corners(2)),corners_y(good_corners(2)));
    side_len = diagonal_len/sqrt(2);
    diagonal_ang = atan2(diff(corners_y(good_corners)),diff(corners_x(good_corners)));
    bad_ang = atan2(diff(corners_y([good_corners(1) bad_corner])),diff(corners_x([good_corners(1) bad_corner])));
    ang_d = bad_ang-diagonal_ang;
    if (ang_d>0 && ang_d<pi) || ang_d<-pi
        new_ang = pi/4 + diagonal_ang;
    else
        new_ang = -pi/4 + diagonal_ang;
    end
    %calculate new corner from new angle and length
    corners_x(bad_corner) = corners_x(good_corners(1)) + side_len*cos(new_ang);
    corners_y(bad_corner) = corners_y(good_corners(1)) + side_len*sin(new_ang);
    corners_x(5) = corners_x(1);
    corners_y(5) = corners_y(1);
    topright_x = corners_x(1);
    topright_y = corners_y(1);
    topleft_x = corners_x(2);
    topleft_y = corners_y(2);
    botleft_x = corners_x(3);
    botleft_y = corners_y(3);
    botright_x = corners_x(4);
    botright_y = corners_y(4);
    box_lengths(1) = get_dist(topright_x,topright_y,topleft_x,topleft_y);
    box_lengths(2) = get_dist(topleft_x,topleft_y,botleft_x,botleft_y);
    box_lengths(3) = get_dist(botleft_x,botleft_y,botright_x,botright_y);
    box_lengths(4) = get_dist(botright_x,botright_y,topright_x,topright_y);
    box_length_in_pixels = mean(box_lengths,'omitnan');
end
pixels_per_meter = box_length_in_pixels/box_length_in_meters;

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

%% delete all mouse tracking that's too far from the mouse's estimated position
%estimate the mouse's size
mouse_size_in_pixels = median(get_dist(nose_x,nose_y,tail_x,tail_y),'omitnan');

%create time-smoothed bodypart estimates (downsampled to 5 Hz)
verysmooth_fps = 5;
head_x = median([nose_x leftear_x rightear_x],2,'omitnan');
head_y = median([nose_y leftear_y rightear_y],2,'omitnan');
smooth_tail_x = rolling_average(tail_x,1,round(fps/verysmooth_fps),'median');
smooth_tail_y = rolling_average(tail_y,1,round(fps/verysmooth_fps),'median');
smooth_head_x = rolling_average(head_x,1,round(fps/verysmooth_fps),'median');
smooth_head_y = rolling_average(head_y,1,round(fps/verysmooth_fps),'median');
smooth_body_x = mean([smooth_head_x smooth_tail_x],2,'omitnan');
smooth_body_y = mean([smooth_head_y smooth_tail_y],2,'omitnan');

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


%% estimate the mouse's body position again (downsampled)
smooth_tail_x = rolling_average(tail_x,1,round(fps/down_fps),'median');
smooth_tail_y = rolling_average(tail_y,1,round(fps/down_fps),'median');
smooth_nose_x = rolling_average(nose_x,1,round(fps/down_fps),'median');
smooth_nose_y = rolling_average(nose_y,1,round(fps/down_fps),'median');
smooth_rightear_x = rolling_average(rightear_x,1,round(fps/down_fps),'median');
smooth_rightear_y = rolling_average(rightear_y,1,round(fps/down_fps),'median');
smooth_leftear_x = rolling_average(leftear_x,1,round(fps/down_fps),'median');
smooth_leftear_y = rolling_average(leftear_y,1,round(fps/down_fps),'median');
smooth_head_x = median([smooth_nose_x smooth_leftear_x smooth_rightear_x],2,'omitnan');
smooth_head_y = median([smooth_nose_y smooth_leftear_y smooth_rightear_y],2,'omitnan');
smooth_body_x =  mean([smooth_head_x smooth_tail_x],2,'omitnan');
smooth_body_y =  mean([smooth_head_y smooth_tail_y],2,'omitnan');

%% calculate median frame (mouse excluded) for object detection
%sample frames where the mouse was in different positions
[~,xI] = sort(smooth_body_x);
[~,yI] = sort(smooth_body_y);
sample_frames_inds = [xI(ceil(0.02*v.NumFrames)), xI(floor(0.98*v.NumFrames)), yI(ceil(0.02*v.NumFrames)), yI(floor(0.98*v.NumFrames))];
sample_frames_inds = unique(sample_frames_inds);
sample_frames = nan([v.height, v.width, 3, length(sample_frames_inds)]);

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
for i = 1:length(sample_frames_inds)
    fidx = sample_frames_inds(i);
    try
        frame = double(read(v,fidx));
        frame = flipud(frame);
        %mask mouse location out of current frame
        mouse_mask = false(size(frame));
        mouse_mask(round(smooth_body_y(fidx)-(mouse_size_in_pixels/1.5)):round(smooth_body_y(fidx)+(mouse_size_in_pixels/1.5)), round(smooth_body_x(fidx)-(mouse_size_in_pixels/1.5)):round(smooth_body_x(fidx)+(mouse_size_in_pixels/1.5)), :) = true;
        frame(mouse_mask) = nan;
        sample_frames(:,:,:,i) = frame;
    end
end
median_frame = median(sample_frames,4,'omitnan')*scale_factor;

%% if the mouse mask covered any pixels for all frames, replace by background color
if any(isnan(median_frame(:)))
    %calculate the median hsv of the background
    full_bkgd_mask = repmat(poly2mask(corners_x,corners_y,size(median_frame,1),size(median_frame,2)),[1 1 3]);
    masked_median_frame = median_frame;
    masked_median_frame(~full_bkgd_mask) = nan;
    median_bkgd_rgb = median(median(masked_median_frame,1,'omitnan'),2,'omitnan');

    nanidx = isnan(median_frame);
    full_bkgd_rgb = repmat(median_bkgd_rgb,[v.Height,v.Width,1]);
    median_frame(nanidx) = full_bkgd_rgb(nanidx);
end

%draw median frame and box for visual validation
figure('Units','centimeters','Position',figure_size);
ax = subplot(2,4,1);
image(median_frame);
set(ax,'YDir','normal');
xlim(xlimits);
ylim(ylimits)
title('median frame')

%% detect object locations (hue and dark-value method)
%convert median frame to hsv (h and v may be most useful to find objects, since the background is hue-stable and not dark)
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

%using the smallest % difference hues (confident these are background pixels, since objects are small), 
%calculate the S.D. of the background hue and convert hue differences to hue z-scores
hue_nonobj = prctile(median_frame_hue_abs(full_bkgd_mask),100*bkgd_fraction);
bkgd_mask = median_frame_hue_abs<=hue_nonobj & full_bkgd_mask;
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
med_filt_pix = max([1 round(100*med_filt/pixels_per_meter)]); %size of median filter, in pixels (from cm)
median_frame_comb_zscr = (hue_weight*median_frame_hue_zscr) + ((1-hue_weight)*median_frame_vdark_zscr);
median_frame_cz_filt = medfilt2(median_frame_comb_zscr,[med_filt_pix med_filt_pix]);
objects_mask = median_frame_cz_filt>zscore_threshold;

%find connected components, remove all that are too small, find connected
%components again
objects_label = bwlabel(objects_mask); 
num_labels = max(objects_label(:));
for i=1:num_labels
    label_size = sum(objects_label(:)==i); %size in pixels
    label_size = label_size/(pixels_per_meter^2); %size in m^2
    label_size = label_size*(100^2); %size in cm^2
    if label_size<min_obj_size
        objects_label(objects_label==i) = 0;
    end
end
objects_label = bwlabel(logical(objects_label));

figure(); %new "object detection details" figure
ax = subplot(3,3,1);
imshow(median_frame)
set(ax,'YDir','normal');
title('panel 1: median frame')
ax = subplot(3,3,2);
imshow(full_bkgd_mask)
set(ax,'YDir','normal');
title('panel 2: box mask')
ax = subplot(3,3,3);
imshow(bkgd_mask)
set(ax,'YDir','normal');
title('panel 3: bkgd pixels')
ax = subplot(3,3,4);
imshow(median_frame_hue_zscr/20)
set(ax,'YDir','normal');
title('panel 4: hue z-score')
ax = subplot(3,3,5);
imshow(median_frame_vdark_zscr/10)
set(ax,'YDir','normal');
title('panel 5: dark z-score')
ax = subplot(3,3,6);
imshow(median_frame_cz_filt/30)
set(ax,'YDir','normal');
title('panel 6: combined z-score')
ax = subplot(3,3,7);
imshow(objects_mask)
set(ax,'YDir','normal');
title('panel 7: objects mask')

num_objects = max(max(objects_label));
if num_objects<1 || num_objects>4 %save figure and exit
    saveas(gcf,[output_filename '_DetectionDetails.png'])
    close(gcf) %close detectiondetails
    close(gcf) %close results figure
    assert(num_objects>0,'No objects detected');
    assert(num_objects<5,'Too many objects detected (4 max is assumed), video quality may be too poor');
end

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

%finish detection details figure
ax = subplot(3,3,8);
imshow(objects_by_color)
set(ax,'YDir','normal');
title('panel 8: detected objects')
ax = subplot(3,3,9);
imshow(median_frame/2)
set(ax,'YDir','normal');
hold on
title('panel 9: object outlines')
for o = 1:num_objects
    plot(objects_edgepts_x(o,:),objects_edgepts_y(o,:),'Color',colors(o,:),'LineWidth',1)
    scatter(objects_x(o),objects_y(o),'kx','MarkerEdgeColor',colors(o,:))
    text(objects_x(o)+(1.5*objects_radius(o)),objects_y(o)+(1.5*objects_radius(o)),num2str(o),'Color',colors(o,:))
end
saveas(gcf,[output_filename '_DetectionDetails.png'])
close(gcf)

%% incorporate object detection into summary plot
figure(gcf); %current figure is the summary figure again
subplot(2,4,2)
plot(corners_x,corners_y);
hold on
for o = 1:num_objects
    plot(objects_edgepts_x(o,:),objects_edgepts_y(o,:),'Color',colors(o,:),'LineWidth',1)
    text(objects_x(o)+(1.5*objects_radius(o)),objects_y(o)+(1.5*objects_radius(o)),num2str(o),'Color',colors(o,:))
    xlim(xlimits);
    ylim(ylimits)
    title('detected walls and object locations');
end

subplot(2,4,5)
plot(corners_x,corners_y);
hold on
scatter(smooth_body_x,smooth_body_y,1,[0 0 1]);%colors)
xlim(xlimits);
ylim(ylimits)
title('all mouse body positions')

subplot(2,4,6)
plot(corners_x,corners_y);
hold on
xlim(xlimits);
ylim(ylimits)
plotframes = min([num_frames 30*fps]);
title(['first ' num2str(ceil(plotframes/fps)) ' s trajectory'])
plot(smooth_body_x(1:plotframes),smooth_body_y(1:plotframes),'b','LineWidth',0.5);

%calculate distance travelled (first remove nans; i.e. bad tracking)
nonan_sbx = smooth_body_x;
nonan_sbx(isnan(nonan_sbx)) = [];
nonan_sby = smooth_body_y;
nonan_sby(isnan(nonan_sby)) = [];
steps = get_dist(nonan_sbx(1:end-1),nonan_sby(1:end-1),nonan_sbx(2:end),nonan_sby(2:end));
distance_travelled = sum(abs(steps))/pixels_per_meter;
total_time = num_frames/fps;


%% calculate mouse distance to objects
%use primary bodypart
full_obj_edge_x = repmat(permute(objects_edgepts_x,[3 1 2]),[num_frames 1 1]);
full_obj_edge_y = repmat(permute(objects_edgepts_y,[3 1 2]),[num_frames 1 1]);
switch bodypart
    case 'head'
        bodypart_x = repmat(smooth_head_x,[1 num_objects 100]); %there are 100 possible object edgepoints (from pre-allocation)
        bodypart_y = repmat(smooth_head_y,[1 num_objects 100]);
        
    case 'body'
        bodypart_x = repmat(smooth_body_x,[1 num_objects 100]); 
        bodypart_y = repmat(smooth_body_y,[1 num_objects 100]);
        
    case 'nose'
        bodypart_x = repmat(smooth_nose_x,[1 num_objects 100]); 
        bodypart_y = repmat(smooth_nose_y,[1 num_objects 100]);
        
    case 'tail'
        bodypart_x = repmat(smooth_tail_x,[1 num_objects 100]);
        bodypart_y = repmat(smooth_tail_y,[1 num_objects 100]);
        
    case 'head and tail' %first do head, then repeat for tail
        bodypart_x = repmat(smooth_head_x,[1 num_objects 100]); 
        bodypart_y = repmat(smooth_head_y,[1 num_objects 100]);
end

%calculate distance between mouse and nearest object edge
objects_edge_dist = get_dist(bodypart_x, bodypart_y, full_obj_edge_x, full_obj_edge_y);
objects_edge_dist = min(objects_edge_dist,[],3,'omitnan'); %find closest edgepoint

switch bodypart
    case 'head and tail' %repeat process for tail
        bodypart_x = repmat(smooth_tail_x,[1 num_objects 100]); 
        bodypart_y = repmat(smooth_tail_y,[1 num_objects 100]);
        objects_edge_dist2 = get_dist(bodypart_x, bodypart_y, full_obj_edge_x, full_obj_edge_y);
        objects_edge_dist2 = min(objects_edge_dist2,[],3,'omitnan'); %find closest edgepoint
        objects_edge_dist(:,:,2) = objects_edge_dist2;
        objects_edge_dist = max(objects_edge_dist,[],3,'omitnan'); %use bodypart that is farther away
end

%find timepoints where mouse is close to object
close_to_objects = (objects_edge_dist/pixels_per_meter)<(close_threshold/100); %format is [num_frames num_objects]


%% find timepoints where the mouse body is on top of each object
on_top_of_objects = nan([num_frames num_objects]);
location_inds = sub2ind(size(objects_label),round(smooth_body_y),round(smooth_body_x));
nan_inds = isnan(location_inds);
location_inds(nan_inds) = 1;
for o = 1:num_objects
    objects_mask = objects_label==o;
    on_top_of_objects(:,o) = objects_mask(location_inds);
end
on_top_of_objects(repmat(nan_inds,[1 num_objects])) = 0;
objects_ontopof_time = sum(on_top_of_objects,1,'omitnan')/fps;
objects_ontopof_pct = 100*sum(on_top_of_objects,1,'omitnan')/sum(~nan_inds);


%% calculate when mouse is oriented toward objects
head_ang_ears = mod(atan2((smooth_rightear_y-smooth_leftear_y),(smooth_rightear_x-smooth_leftear_x))+pi/2,2*pi);
head_ang_lnose = mod(atan2((smooth_nose_y-smooth_leftear_y),(smooth_nose_x-smooth_leftear_x))+pi/4,2*pi);
head_ang_rnose = mod(atan2((smooth_nose_y-smooth_rightear_y),(smooth_nose_x-smooth_rightear_x))-pi/4,2*pi);

%get vector sum head direction (giving more weight to ears)
vs_x = sum([2*cos(head_ang_ears) cos(head_ang_lnose) cos(head_ang_rnose)],2,'omitnan');
vs_y = sum([2*sin(head_ang_ears) sin(head_ang_lnose) sin(head_ang_rnose)],2,'omitnan');

%get smoother vector sum head angle
smooth_vs_x = rolling_average(vs_x,1,round(fps/down_fps),'median');
smooth_vs_y = rolling_average(vs_y,1,round(fps/down_fps),'median');
smooth_head_ang = mod(atan2(smooth_vs_y,smooth_vs_x),2*pi);

%calculate direction of all objects to mouse head
full_obj_center_x = repmat(objects_x,[num_frames 1]);
full_obj_center_y = repmat(objects_y,[num_frames 1]);
full_head_x = repmat(smooth_head_x,[1 num_objects]);
full_head_y = repmat(smooth_head_y,[1 num_objects]);
objects_ang = mod(atan2((full_obj_center_y-full_head_y),((full_obj_center_x-full_head_x))),2*pi);
objects_head_ang = get_angular_dist(objects_ang,repmat(smooth_head_ang,[1 num_objects]),'radians'); %angle from straight ahead (from the mouse's perspective)
angled_toward_objects = objects_head_ang<deg2rad(view_cone/2);


%% calculate interactions
objects_interaction = close_to_objects & angled_toward_objects & ~on_top_of_objects; %nice
objects_interaction_time = sum(objects_interaction,1,'omitnan')/fps;
objects_interaction_pct = 100*objects_interaction_time/sum(objects_interaction_time);


%% summarize results
results.number_of_frames = num_frames;
results.camera_fps = fps;
results.downsampled_fps = down_fps;
results.primary_bodypart = bodypart;
results.view_cone = view_cone;
results.close_threshold = close_threshold;
results.total_time_in_seconds = total_time;
results.distance_travelled_in_meters = distance_travelled;
results.number_of_objects = num_objects;
results.object_interaction_time = objects_interaction_time;
results.object_interaction_percent = objects_interaction_pct;
results.object_ontopof_time = objects_ontopof_time;
results.object_ontopof_percent = objects_ontopof_pct;
results.object_interaction_vectors = objects_interaction;
results.object_detection.hl_ratio = hue_weight;
results.object_detection.bkgd_fraction = bkgd_fraction;
results.object_detection.zscore_threshold = zscore_threshold;
results.object_detection.min_obj_size = min_obj_size;
results.object_detection.med_filt = med_filt;

%print results to figure
subplot(2,4,3)
set(gca,'Units','Normalized','Position',[0.6 0.05 0.4 0.95])
axis off
text(0,0.95,['Number of frames: ' num2str(num_frames)])
text(0,0.88,['Camera fps: ' num2str(fps)])
text(0,0.81,['Downsampled fps: ' num2str(down_fps)])
text(0,0.74,['Primary bodypart: ' bodypart])
text(0,0.67,['View cone angle: ' num2str(view_cone) ' deg'])
text(0,0.60,['Close threshold: ' num2str(close_threshold) ' cm'])
text(0,0.53,['Total time: ' num2str(total_time, '%.1f') ' s'])
text(0,0.46,['Distance travelled: ' num2str(distance_travelled, '%.1f') ' m'])
text(0,0.39,'Object labels: ')
text(0,0.34,['[' num2str(1:num_objects) ']'],'FontSize',8)
text(0,0.27,'Time interacting with objects: ')
text(0,0.22,['[' num2str(objects_interaction_time, '%.1f  ') '] s,   [' num2str(objects_interaction_pct,'%.1f  ') '] %']);
text(0,0.15,'Time on top of objects: ')
text(0,0.10,['[' num2str(objects_ontopof_time, '%.1f  ') '] s,   [' num2str(objects_ontopof_pct,'%.1f  ') '] %']);

saveas(gcf,[output_filename '.png'])
save([output_filename '.mat'],'results');
if ~show_figure
    close(gcf)
end

%write custom xls file of results
C = cell(8,2);
C(:,1) = {'number of frames', 'camera fps', 'downsampled fps', 'primary bodypart', 'view cone angle (deg)', 'close threshold (cm)', 'total time (s)', 'distance travelled (m)'};
C(:,2) = {num_frames, fps, down_fps, bodypart, view_cone, close_threshold, total_time, distance_travelled};
writecell(C,[output_filename '.xls'],'Range','A1');

C = cell(5,2);
C(:,1) = {'h-l ratio', 'background fraction', 'z-score threshold', 'minimum object size (cm)', 'median filter (cm)'};
C(:,2) = {hue_weight, bkgd_fraction, zscore_threshold, min_obj_size, med_filt};
writecell(C,[output_filename '.xls'],'Range','D1');

C = cell(5,1+num_objects);
C(:,1) = {'object labels', 'Time interacting with objects (s)', 'Time interacting with objects (%)', 'Time on top of objects (s)', 'Time on top of objects (%)'};
C(1,2:1+num_objects) = num2cell(1:num_objects);
C(2,2:1+num_objects) = num2cell(objects_interaction_time);
C(3,2:1+num_objects) = num2cell(objects_interaction_pct);
C(4,2:1+num_objects) = num2cell(objects_ontopof_time);
C(5,2:1+num_objects) = num2cell(objects_ontopof_pct);
writecell(C,[output_filename '.xls'],'Range','A10');
