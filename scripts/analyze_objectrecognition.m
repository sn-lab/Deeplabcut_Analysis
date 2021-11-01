function analyze_objectrecognition(csv_filename, avi_filename)
% FUNCTION analyze_objectrecognition(csv_filename, avi_filename)
%
% analyzes csv files output from deeplabcut tracking of object recognition
% videos with multiple objects placed in the red box, using the 
% "Redbox_2" training dataset. Analysis includes plotting the summary of 
% mouse and box tracking data, and calculating the time the mouse spent 
% "interacting" with each detected object. In order to detect objects, this 
% script opens the original video file in addition to using the deeplabcut 
% tracking results xls file.
%
%INPUTS
%csv_filename: full path/filename of a csv file to analyze
%avi_filename: full path/filename of the avi file associated with the csv

[save_dir, filename, ~] = fileparts(csv_filename);
output_filename = fullfile(save_dir,filename);
fps = 30;
nearobject_threshold = 0.04; %near object positions
verynearobject_threshold = 0.02; %very near object positions
box_length_in_meters = 0.4445;
figure_size = [10 10 25 10];

%load tracking data
opts = detectImportOptions(csv_filename);
opts.VariableNamesLine = 2;
tbl = readtable(csv_filename, opts);
colnames = tbl.Properties.VariableNames;
ary = table2array(tbl);
[num_frames, num_cols] = size(ary);

%vertically flip all y values
ystrings = repmat({'_1'},1,num_cols);
col_is_y = cellfun(@strfind,colnames,ystrings,'UniformOutput',false);
col_is_y = find(~cellfun(@isempty,col_is_y));
ary(:,col_is_y) = 1080-ary(:,col_is_y);

%load video subset (~100 frames) for object detection
assert(exist(avi_filename,'file')==2,'Cannot find avi file');
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
topright_x_ind = find(strcmp(colnames,'toprightcorner'));
topright_y_ind = find(strcmp(colnames,'toprightcorner_1'));
topleft_x_ind = find(strcmp(colnames,'topleftcorner'));
topleft_y_ind = find(strcmp(colnames,'topleftcorner_1'));
botleft_x_ind = find(strcmp(colnames,'botleftcorner'));
botleft_y_ind = find(strcmp(colnames,'botleftcorner_1'));
botright_x_ind = find(strcmp(colnames,'botrightcorner'));
botright_y_ind = find(strcmp(colnames,'botrightcorner_1'));

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
F1 = figure('Units','centimeters','Position',figure_size);
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


%% detect object locations
object_search_width = 200;
object_gaussfilt = 20;

%estimate possible object locations
object_location_names = {'top-right','top-left','bottom-left','bottom-right'};
topright_object_x = botleft_x + (2*(topright_x - botleft_x)/3);
topleft_object_x = topleft_x + (botright_x - topleft_x)/3;
botleft_object_x = botleft_x + (topright_x - botleft_x)/3;
botright_object_x = topleft_x + 2*(botright_x - topleft_x)/3;

topright_object_y = botleft_y + 2*(topright_y - botleft_y)/3;
topleft_object_y = botright_y + 2*(topleft_y - botright_y)/3;
botleft_object_y = botleft_y + (topright_y - botleft_y)/3;
botright_object_y = botright_y + (topleft_y - botright_y)/3;

possible_objects_x = round([topright_object_x topleft_object_x botleft_object_x botright_object_x]);
possible_objects_y = round([topright_object_y topleft_object_y botleft_object_y botright_object_y]);
num_objects = 4;

%estimate similar locations where no object should be (for background subtraction step of object detection)
nonright_object_x = (topright_object_x + botright_object_x)/2;
nonright_object_y = (topright_object_y + botright_object_y)/2;
nonleft_object_x = (topleft_object_x + botleft_object_x)/2;
nonleft_object_y = (topleft_object_y + botleft_object_y)/2;
non_objects_x = [nonright_object_x nonleft_object_x nonleft_object_x nonright_object_x ];
non_objects_y = [nonright_object_y nonleft_object_y nonleft_object_y nonright_object_y];

% get background color value
xrange = round(non_objects_x(1)-(object_search_width/4)):round(non_objects_x(1)+(object_search_width/4));
yrange = round(non_objects_y(1)-(object_search_width/4)):round(non_objects_y(1)+(object_search_width/4));
bkgd_right = median(median(median_frame(yrange,xrange,:),'omitnan'),2,'omitnan'); % get backgroun substraction val

xrange = round(non_objects_x(2)-(object_search_width/4)):round(non_objects_x(2)+(object_search_width/4));
yrange = round(non_objects_y(2)-(object_search_width/4)):round(non_objects_y(2)+(object_search_width/4));
bkgd_left = median(median(median_frame(yrange,xrange,:),'omitnan'),2,'omitnan'); % get backgroun substraction val

bkgd_colors = [bkgd_right bkgd_left bkgd_left bkgd_right];

%time to figure out the actual object locations
num_objects_edgepts = round(object_search_width/2);
objects_edgepts_x = nan(num_objects_edgepts,4); %pre-allocating plenty of space
objects_edgepts_y = nan(num_objects_edgepts,4);
objects_x = nan(1,4);
objects_y = nan(1,4);
F2 = figure();
for o = 1:num_objects
    %get median object location pic
    xrange = round(-object_search_width/2):round(object_search_width/2);
    yrange = round(-object_search_width/2):round(object_search_width/2);
    object_pic = median_frame(possible_objects_y(o)+yrange,possible_objects_x(o)+xrange,:);
    pic_h = size(object_pic,1);
    pic_w = size(object_pic,2);
    ax = subplot(4,4,1+4*(o-1));
    imshow(object_pic)
    set(ax,'YDir','normal');
    ylabel(object_location_names{o},'FontSize',10,'Rotation',0,'Position',[-130,80])
    if o==1
        title('estimated location','FontSize',9)
    end
    
    %subtract background color, threshold color difference to detect object pixels
    object_pic = object_pic-repmat(bkgd_colors(1,o,:),[pic_h pic_w 1]);
    object_pic_diffs = sqrt(sum(object_pic.^2,3,'omitnan'));
    object_pic_diffs = object_pic_diffs>0.5;
    ax = subplot(4,4,2+4*(o-1));
    imshow(object_pic_diffs)
    set(ax,'YDir','normal');
    if o==1
        title('bkgd subtracted','FontSize',9);
    end
    
    %gaussian filter object pixels, threshold to get object
    object_pic_filt = imgaussfilt(object_pic_diffs*255,object_gaussfilt);
    object_mask = object_pic_filt>60;
    
    %only use the single largest blob
    object_label = bwlabel(object_mask); 
    num_labels = max(max(object_label));
    label_sizes = zeros(1,num_labels);
    for i=1:num_labels
        label_sizes(i) = sum(sum(object_label==i));
    end
    if ~isempty(label_sizes)
        object_mask = object_label==find(max(label_sizes));
    end
    ax = subplot(4,4,3+4*(o-1));
    imshow(object_mask)
    set(ax,'YDir','normal');
    if o==1
        title('object mask','FontSize',9);
    end
    
    %get approximate center of object
    [r, c] = find(object_mask == 1);
    objects_x(o) = mean(c) + possible_objects_x(o) - round(object_search_width/2);
    objects_y(o) = mean(r) + possible_objects_y(o) - round(object_search_width/2);
    
    %get edge coordinate of any detected objects
    if sum(sum(object_mask))>500
        vertical = any(object_mask,2);
        bottom = find(diff(vertical)==1,1,'first') + 1; %index from "find(diff())" is the pixel before the change
        top = find(diff(vertical)==-1,1,'first'); 
        cur_object_y = [bottom:5:top-3 top];
        num_sidepts = length(cur_object_y);
        cur_object_y = [cur_object_y fliplr(cur_object_y) cur_object_y(1)];
        num_object_edgepts = 1+(num_sidepts*2);
        cur_object_x = nan([1 num_object_edgepts]);
        for i = 1:num_sidepts
            cur_object_x(i) = possible_objects_x(o) + find(object_mask(cur_object_y(i),:)==1,1,'first') - round(object_search_width/2);
            cur_object_x(i+num_sidepts) = possible_objects_x(o) + find(object_mask(cur_object_y(i+num_sidepts),:)==1,1,'last') - round(object_search_width/2);
        end
        cur_object_x(end) = cur_object_x(1);
        cur_object_y = possible_objects_y(o) + cur_object_y - round(object_search_width/2);
        objects_edgepts_x(1:num_object_edgepts,o) = cur_object_x;
        objects_edgepts_y(1:num_object_edgepts,o) = cur_object_y;
        
        ax = subplot(4,4,4*o);
        imshow(median_frame)
        hold on
        set(ax,'YDir','normal');
        plot(objects_edgepts_x(:,o),objects_edgepts_y(:,o),'g','LineWidth',3)
        scatter(objects_x(o),objects_y(o),'ob','SizeData',9)
    else
        subplot(4,4,4*o)
        text(0.5, 0.6, ['no ' object_location_names{o}],'FontSize',9,'HorizontalAlignment','center')
        text(0.5, 0.4, 'object detected','FontSize',9,'HorizontalAlignment','center')
        axis off
        box off
    end
    if o==1
        title('object outline','FontSize',9);
    end
    
end  

figure(F1);
subplot(2,4,2)
plot(corners_x,corners_y);
hold on
for o = 1:num_objects
    plot(objects_edgepts_x,objects_edgepts_y,'g','LineWidth',1)
    xlim([-100 1550]);
    ylim([-280 1370])
    title('detected walls and object locations');
end


%% delete all mouse tracking that's outside of the box (change to NaNs)
nose_x_ind = find(strcmp(colnames,'nose'));
nose_y_ind = find(strcmp(colnames,'nose_1'));
leftear_x_ind = find(strcmp(colnames,'leftear'));
leftear_y_ind = find(strcmp(colnames,'leftear_1'));
rightear_x_ind = find(strcmp(colnames,'rightear'));
rightear_y_ind = find(strcmp(colnames,'rightear_1'));
tail_x_ind = find(strcmp(colnames,'tailbase'));
tail_y_ind = find(strcmp(colnames,'tailbase_1'));

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
objects_ang = atan2((repmat(objects_y,[num_frames 1])-smooth_head_y),(repmat(objects_x,[num_frames 1])-smooth_head_x));

%calculate distances of mouse body to object edges
full_obj_edge_x = repmat(permute(objects_edgepts_x,[3 2 1]),[num_frames 1 1]);
full_obj_edge_y = repmat(permute(objects_edgepts_y,[3 2 1]),[num_frames 1 1]);
objects_dist = get_dist(repmat(smooth_body_x,[1 num_objects num_objects_edgepts]),repmat(smooth_body_y,[1 num_objects num_objects_edgepts]), full_obj_edge_x, full_obj_edge_y);
objects_dist = min(objects_dist,[],3,'omitnan'); %find closest edgepoint
nearobjects = (objects_dist/pixels_per_meter)<nearobject_threshold;
verynearobjects = (objects_dist/pixels_per_meter)<verynearobject_threshold;

neartopright = nearobjects(:,1);
neartopleft = nearobjects(:,2);
nearbotleft = nearobjects(:,3);
nearbotright = nearobjects(:,4);
veryneartopright = verynearobjects(:,1);
veryneartopleft = verynearobjects(:,2);
verynearbotleft = verynearobjects(:,3);
verynearbotright = verynearobjects(:,4);

%calculate orientation to each object position
topright_object_head_ang = get_angular_dist(objects_ang(:,1),smooth_headang,'radians');
topleft_object_head_ang = get_angular_dist(objects_ang(:,2),smooth_headang,'radians');
botleft_object_head_ang = get_angular_dist(objects_ang(:,3),smooth_headang,'radians');
botright_object_head_ang = get_angular_dist(objects_ang(:,4),smooth_headang,'radians');

%calculate near and verynear interactions (orientation towards + near or verynear)
nearinteract_topright_object = neartopright & topright_object_head_ang<(pi/2);
verynearinteract_topright_object = veryneartopright & topright_object_head_ang<(pi/2);
nearinteract_topleft_object = neartopleft & topleft_object_head_ang<(pi/2);
verynearinteract_topleft_object = veryneartopleft& topleft_object_head_ang<(pi/2);
nearinteract_botleft_object = nearbotleft & botleft_object_head_ang<(pi/2);
verynearinteract_botleft_object = verynearbotleft & botleft_object_head_ang<(pi/2);
nearinteract_botright_object = nearbotright & botright_object_head_ang<(pi/2);
verynearinteract_botright_object = verynearbotright & botright_object_head_ang<(pi/2);

near_object_interaction_time = [sum(nearinteract_topright_object) sum(nearinteract_topleft_object) sum(nearinteract_botleft_object) sum(nearinteract_botright_object)]/fps;
verynear_object_interaction_time = [sum(verynearinteract_topright_object) sum(verynearinteract_topleft_object) sum(verynearinteract_botleft_object) sum(verynearinteract_botright_object)]/fps;


%% calculate summary statistics
steps = get_dist(smooth_body_x(1:end-1),smooth_body_y(1:end-1),smooth_body_x(2:end),smooth_body_y(2:end));
distance_travelled = sum(abs(steps))/pixels_per_meter;

results.number_of_frames = num_frames;
results.assumed_camera_fps = fps;
results.total_time_in_seconds = num_frames/fps;
results.distance_travelled_in_meters = distance_travelled;
results.object_location_names = object_location_names;
results.object_interaction_time_4cm = near_object_interaction_time;
results.object_interaction_time_2cm = verynear_object_interaction_time;

%print results to figure
subplot(2,4,3)
set(gca,'Units','Normalized','Position',[0.6 0.05 0.4 0.95])
axis off
text(0,0.9,['Number of frames: ' num2str(results.number_of_frames)])
text(0,0.82,['Assumed camera fps: ' num2str(results.assumed_camera_fps)])
text(0,0.74,['Total time: ' num2str(results.total_time_in_seconds, '%.1f') ' s'])
text(0,0.66,['Distance travelled: ' num2str(results.distance_travelled_in_meters, '%.1f') ' m'])
text(0,0.58,'Object locations: ')
text(0,0.53,['[' cell2str(results.object_location_names) ']'],'FontSize',8)
text(0,0.45,'Time interacting with objects - within 4 cm: ')
text(0,0.40,['[' num2str(results.object_interaction_time_4cm, '%.1f  ') '] s,   [' num2str(round(1000*results.object_interaction_time_4cm/results.total_time_in_seconds)/10,'%.1f  ') '] %']);
text(0,0.32,'Time interacting with objects - within 2 cm: ')
text(0,0.27,['[' num2str(results.object_interaction_time_2cm, '%.1f  ') '] s,   [' num2str(round(1000*results.object_interaction_time_2cm/results.total_time_in_seconds)/10, '%.1f  ') '] %']);

saveas(gcf,[output_filename '.png'])
save([output_filename '.mat'],'results');

%write custom xls file of results
C = cell(4,2);
C(:,1) = {'number of frames', 'assumed camera fps', 'total time (s)', 'distance travelled (m)'};
C(:,2) = {num_frames, fps, results.total_time_in_seconds, distance_travelled};
writecell(C,[output_filename '.xls'],'Range','A1');
C = cell(5,5);
C(:,1) = {'object location', 'Time interacting with objects - within 4 cm (s)', 'Time interacting with objects - within 4 cm (%)', 'Time interacting with objects - within 2 cm (s)', 'Time interacting with objects - within 2 cm (%)'};
C(1,2:5) = {results.object_location_names{1}, results.object_location_names{2}, results.object_location_names{3}, results.object_location_names{4}};
C(2,2:5) = {results.object_interaction_time_4cm(1), results.object_interaction_time_4cm(2), results.object_interaction_time_4cm(3), results.object_interaction_time_4cm(4)};
C(3,2:5) = {100*results.object_interaction_time_4cm(1)/results.total_time_in_seconds, 100*results.object_interaction_time_4cm(2)/results.total_time_in_seconds, 100*results.object_interaction_time_4cm(3)/results.total_time_in_seconds, 100*results.object_interaction_time_4cm(4)/results.total_time_in_seconds};
C(4,2:5) = {results.object_interaction_time_2cm(1), results.object_interaction_time_2cm(2), results.object_interaction_time_2cm(3), results.object_interaction_time_2cm(4)};
C(5,2:5) = {100*results.object_interaction_time_2cm(1)/results.total_time_in_seconds, 100*results.object_interaction_time_2cm(2)/results.total_time_in_seconds, 100*results.object_interaction_time_2cm(3)/results.total_time_in_seconds, 100*results.object_interaction_time_2cm(4)/results.total_time_in_seconds};
writecell(C,[output_filename '.xls'],'Range','A6');