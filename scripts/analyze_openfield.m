function analyze_openfield(h5_filename)
% FUNCTION analyze_openfield(h5_filename)
%
% analyzes h5 files output from deeplabcut tracking of open field (aka
% "red box") videos using the "Redbox_2" training dataset. Analysis includes
% plotting the summary of mouse and box tracking data, and calculating the 
% mouse's time and travel distance near walls and in the box center
%
%INPUTS
%h5_filename: full path/filename of a h5 file to analyze

%set parameters
fps = 30;
nearedge_threshold = .10; %how close to the edge (in m) a mouse has to be to be call "near"
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
subplot(2,4,1)
plot(corners_x,corners_y,'r');
xlim([-100 1550]);
ylim([-280 1370])
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

subplot(2,4,2)
scatter(edgepts_x,edgepts_y,200,'MarkerFaceColor','flat');
hold on
patch([corners_x 1550 1550 -100 -100 1550 corners_x(1)],[corners_y 1370 -280 -280 1370 1370 corners_y(1)],'w','EdgeColor','None');
xlim([-100 1550]);
ylim([-280 1370])
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

%% calculate results
%calculate distances of mouse body to edge
edge_dists = get_dist(repmat(smooth_head_x,[1 num_edgepts]), repmat(smooth_head_y,[1 num_edgepts]), full_edgepts_x, full_edgepts_y);
edge_dist = min(edge_dists,[],2,'omitnan')/pixels_per_meter;
nearedge = edge_dist<=nearedge_threshold;
incenter = edge_dist>nearedge_threshold;

time_near_edge = sum(nearedge,'omitnan')/fps;
time_in_center = sum(incenter,'omitnan')/fps;
steps = [0; get_dist(smooth_body_x(1:end-1),smooth_body_y(1:end-1),smooth_body_x(2:end),smooth_body_y(2:end))];
distance_travelled = sum(abs(steps),'omitnan')/pixels_per_meter;
distance_near_edge = sum(abs(steps(nearedge)),'omitnan')/pixels_per_meter;
distance_in_center = sum(abs(steps(incenter)),'omitnan')/pixels_per_meter;
total_time = num_frames/fps;

results.number_of_frames = num_frames;
results.assumed_camera_fps = fps;
results.total_time_in_seconds = total_time;
results.distance_travelled_in_meters = distance_travelled;
results.near_wall_time = time_near_edge;
results.near_wall_distance = distance_near_edge;
results.in_center_time = time_in_center;
results.in_center_distance = distance_in_center;

%print results to figure
subplot(2,4,3)
set(gca,'Units','Normalized','Position',[0.6 0.05 0.4 0.95])
axis off
text(0,0.9,['Number of frames: ' num2str(results.number_of_frames)])
text(0,0.82,['Assumed camera fps: ' num2str(results.assumed_camera_fps)])
text(0,0.74,['Total time: ' num2str(results.total_time_in_seconds, '%.1f') ' s'])
text(0,0.66,['Distance travelled: ' num2str(results.distance_travelled_in_meters, '%.1f') ' m'])
text(0,0.58,['Time near walls: ' num2str(results.near_wall_time, '%.1f'), ' s'])
text(0,0.50,['Distance near walls: ' num2str(results.near_wall_distance, '%.1f'), ' s'])
text(0,0.42,['Time in center: ' num2str(results.in_center_time, '%.1f'), ' s'])
text(0,0.34,['Distance in center: ' num2str(results.in_center_time, '%.1f'), ' s'])

saveas(gcf,[output_filename '.png'])
save([output_filename '.mat'],'results');

%write custom xls file of results
C = cell(8,2);
C(:,1) = {'number of frames', 'assumed camera fps', 'total time (s)', 'distance travelled (m)', 'time near wall (s)', 'distance near wall (m)', 'time in center (s)', 'distance in center (m)'};
C(:,2) = {num_frames, fps, results.total_time_in_seconds, distance_travelled, time_near_edge, distance_near_edge, time_in_center, distance_in_center};
writecell(C,[output_filename '.xls'],'Range','A1');
C = cell(5,1);
C(:) = {'% of total time/distance', 100*time_near_edge/total_time, 100*distance_near_edge/distance_travelled, 100*time_in_center/total_time, 100*distance_near_edge/distance_travelled};
writecell(C,[output_filename '.xls'],'Range','A3');