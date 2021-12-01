function analyze_openfield(h5_filename, param)
% FUNCTION analyze_openfield(h5_filename, param)
%
% analyzes h5 files output from deeplabcut tracking of open field (aka
% "red box") videos using the "Redbox_2" training dataset. Analysis includes
% plotting the summary of mouse and box tracking data, and calculating the 
% mouse's time and travel distance near walls and in the box center
%
%INPUTS
%h5_filename: full path/filename of a h5 file to analyze
%param (optional):
%     fps = framerate of the original video recordings
%     down_fps = effective fps after rolling average performed to remove some tracking errors
%     bodypart = bodypart to use for distance measurements (e.g. mouse-wall distance)
%     close_threshold = minimum distance (cm) from mouse bodypart to the wall to count as "near wall"

if nargin<2
    %set default parameters
    fps = 30; %framerate of the original video recordings
    down_fps = 10; %effective fps after rolling average performed to remove some tracking errors
    bodypart = 'head'; %bodypart to use for distance measurements (e.g. mouse-wall distance)
    close_threshold = 10; %minimum distance (cm) from mouse bodypart to the wall to count as "near wall"
else
    fps = param.fps;
    down_fps = param.down_fps;
    bodypart = param.bodypart;
    close_threshold = param.close_threshold;
end

box_length_in_meters = 0.4445;
figure_size = [10 10 25 10];

%load tracking data
[save_dir, filename, ~] = fileparts(h5_filename);
ind = strfind(filename,'DeepCut');
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

%estimate original video size
ystrings = repmat({'_y'},1,num_cols);
col_is_y = cellfun(@strfind,colnames,ystrings,'UniformOutput',false);
col_is_y = find(~cellfun(@isempty,col_is_y));
xstrings = repmat({'_x'},1,num_cols);
col_is_x = cellfun(@strfind,colnames,xstrings,'UniformOutput',false);
col_is_x = find(~cellfun(@isempty,col_is_x));
arymax = max(ary,[],'omitnan');
v.width = max(arymax(col_is_x),[],'omitnan');
v.height = max(arymax(col_is_y),[],'omitnan');

%vertically flip all y values
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
assert(all((box_lengths./box_length_in_pixels)>0.95),'detected box is out of shape - tracking is probably too poor to continue')
pixels_per_meter = box_length_in_pixels/box_length_in_meters;

%draw box for validation
figure('Units','centimeters','Position',figure_size);
subplot(2,4,1)
plot(corners_x,corners_y,'r');
xlim(xlimits);
ylim(ylimits)
title('detected box edges')

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

pc = pixels_per_meter*(close_threshold/100); %how many pixels "close" is
subplot(2,4,2)
hold on
patch([corners_x corners_x(1) corners_x(1)-pc corners_x(4)-pc corners_x(3)+pc corners_x(2)+pc corners_x(1)-pc corners_x(1)],...
    [corners_y corners_y(1) corners_y(1)-pc corners_y(4)+pc corners_y(3)+pc corners_y(2)-pc corners_y(1)-pc corners_y(1)],'b');
xlim(xlimits);
ylim(ylimits)
title('edge location');


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


%% estimate the mouse's body position again
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

%% calculate results
%calculate distances of mouse bodypart to edge
switch bodypart
    case 'head'
        bodypart_x = repmat(smooth_head_x,[1 num_edgepts]);
        bodypart_y = repmat(smooth_head_y,[1 num_edgepts]);
        edge_dists = get_dist(bodypart_x,bodypart_y, full_edgepts_x, full_edgepts_y);
        edge_dist = min(edge_dists,[],2,'omitnan'); %find closest edgepoint

    case 'body'
        bodypart_x = repmat(smooth_body_x,[1 num_edgepts]);
        bodypart_y = repmat(smooth_body_y,[1 num_edgepts]);
        edge_dists = get_dist(bodypart_x,bodypart_y, full_edgepts_x, full_edgepts_y);
        edge_dist = min(edge_dists,[],2,'omitnan'); %find closest edgepoint

    case 'nose'
        bodypart_x = repmat(smooth_nose_x,[1 num_edgepts]);
        bodypart_y = repmat(smooth_nose_y,[1 num_edgepts]);
        edge_dists = get_dist(bodypart_x,bodypart_y, full_edgepts_x, full_edgepts_y);
        edge_dist = min(edge_dists,[],2,'omitnan'); %find closest edgepoint

    case 'tail'
        bodypart_x = repmat(smooth_tail_x,[1 num_edgepts]);
        bodypart_y = repmat(smooth_tail_y,[1 num_edgepts]);
        edge_dists = get_dist(bodypart_x,bodypart_y, full_edgepts_x, full_edgepts_y);
        edge_dist = min(edge_dists,[],2,'omitnan'); %find closest edgepoint

    case 'head and tail'
        bodypart_x = repmat(smooth_head_x,[1 num_edgepts]);
        bodypart_y = repmat(smooth_head_y,[1 num_edgepts]);
        edge_dists = get_dist(bodypart_x,bodypart_y, full_edgepts_x, full_edgepts_y);
        edge_dist1 = min(edge_dists,[],2,'omitnan'); %find closest edgepoint
        bodypart_x = repmat(smooth_tail_x,[1 num_edgepts]);
        bodypart_y = repmat(smooth_tail_y,[1 num_edgepts]);
        edge_dists = get_dist(bodypart_x,bodypart_y, full_edgepts_x, full_edgepts_y);
        edge_dist = min(edge_dists,[],2,'omitnan'); %find closest edgepoint
        edge_dist(:,1) = edge_dist1;
        edge_dist(:,2) = edge_dist2;
        edge_dist = max(edge_dist,[],2,'omitnan'); %which bodypart is farther away
end

closewall = (edge_dist/pixels_per_meter)<=(close_threshold/100);
incenter = (edge_dist/pixels_per_meter)>(close_threshold/100);

time_closewall = sum(closewall,'omitnan')/fps;
time_incenter = sum(incenter,'omitnan')/fps;
steps = [0; get_dist(smooth_body_x(1:end-1),smooth_body_y(1:end-1),smooth_body_x(2:end),smooth_body_y(2:end))];
distance_travelled = sum(abs(steps),'omitnan')/pixels_per_meter;
distance_closewall = sum(abs(steps(closewall)),'omitnan')/pixels_per_meter;
distance_incenter = sum(abs(steps(incenter)),'omitnan')/pixels_per_meter;
total_time = num_frames/fps;

results.number_of_frames = num_frames;
results.camera_fps = fps;
results.downsampled_fps = down_fps;
results.primary_bodypart = bodypart;
results.close_threshold = close_threshold;
results.total_time_in_seconds = total_time;
results.distance_travelled_in_meters = distance_travelled;
results.time_close_to_wall = time_closewall;
results.distance_travelled_close_to_wall = distance_closewall;
results.time_in_center = time_incenter;
results.distance_travelled_in_center = distance_incenter;

%print results to figure
subplot(2,4,3)
set(gca,'Units','Normalized','Position',[0.6 0.05 0.4 0.95])
axis off
text(0,0.95,['Number of frames: ' num2str(num_frames)])
text(0,0.88,['Camera fps: ' num2str(fps)])
text(0,0.81,['Downsampled fps: ' num2str(down_fps)])
text(0,0.74,['Primary bodypart: ' bodypart])
text(0,0.67,['Close threshold: ' num2str(close_threshold) ' cm'])
text(0,0.60,['Total time: ' num2str(total_time, '%.1f') ' s'])
text(0,0.53,['Distance travelled: ' num2str(distance_travelled, '%.1f') ' m'])
text(0,0.46,['Time close to walls: ' num2str(time_closewall, '%.1f'), ' s'])
text(0,0.39,['Distance travelled close to walls: ' num2str(distance_closewall, '%.1f'), ' m'])
text(0,0.32,['Time in center: ' num2str(time_incenter, '%.1f'), ' s'])
text(0,0.25,['Distance travelled in center: ' num2str(distance_incenter, '%.1f'), ' m'])

saveas(gcf,[output_filename '.png'])
save([output_filename '.mat'],'results');

%write custom xls file of results
C = cell(11,2);
C(:,1) = {'number of frames', 'camera fps', 'downsampled fps', 'primary bodypart', 'close threshold', 'total time (s)', 'distance travelled (m)', 'time near wall (s)', 'distance travelled near wall (m)', 'time in center (s)', 'distance travelled in center (m)'};
C(:,2) = {num_frames, fps, down_fps, bodypart, close_threshold, total_time, distance_travelled, time_closewall, distance_closewall, time_incenter, distance_incenter};
writecell(C,[output_filename '.xls'],'Range','A1');
C = cell(5,1);
C(:) = {'% of total time/distance', 100*time_closewall/total_time, 100*distance_closewall/distance_travelled, 100*time_incenter/total_time, 100*distance_incenter/distance_travelled};
writecell(C,[output_filename '.xls'],'Range','C7');