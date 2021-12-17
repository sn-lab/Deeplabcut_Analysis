function analyze_ymaze(h5_filename, param)
% FUNCTION analyze_ymaze(h5_filename, param)
%
% analyzes h5 files output from deeplabcut tracking of Y-maze videos using 
% the "Ymaze_2" training dataset. Analysis includes plotting the summary of 
% mouse and ymaze tracking data, calculating time the mouse spent within 
% each arm, the sequence of arm entries, and the spontaneous alternation %
%
%INPUTS
%h5_filename: full path/filename of a h5 file to analyze
%param (optional):
%     fps = framerate of the original video recordings
%     down_fps = effective fps after rolling average performed to remove some tracking errors
%     bodypart = bodypart to use for distance measurements (e.g. mouse-arm distance)
%     arm_threshold = minimum distance (cm) down arm a mouse bodypart has to be to count as in that arm
%     show_figure = whether to show (1) or close (0) the results figure when finished

if nargin<2
    %set default parameters
    fps = 30; %framerate of the original video recordings
    down_fps = 10; %effective fps after rolling average performed to remove some tracking errors
    bodypart = 'head'; %bodypart to use for distance measurements (e.g. mouse-arm distance)
    arm_threshold = 8; %minimum distance (cm) down arm a mouse bodypart has to be to count as in that arm
    show_figure = 1; %whether to show (1) or close (0) the results figure when finished
else
    fps = param.fps;
    down_fps = param.down_fps;
    bodypart = param.bodypart;
    arm_threshold = param.arm_threshold;
    show_figure = param.show_figure;
end

arm_length_in_meters = 0.365;
arm_width_in_meters = 0.07; 
figure_size = [10 10 25 10];
insidearm_threshold_fraction = (arm_length_in_meters-(arm_threshold/100))/arm_length_in_meters; %how close to the center a mouse can be and still be identified as in an arm

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
%columns: right-arm, left-arm, middle-arm, right-ear, left-ear, nose, tail-base; each have 3 columns, for [x, y, likelihood]
colnames = {'RA_x','RA_y','RA_l','LA_x','LA_y','LA_l','MA_x','MA_y','MA_l','RE_x','RE_y','RE_l','LE_x','LE_y','LE_l','NO_x','NO_y','NO_l','TB_x','TB_y','TB_l'};
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

%% calculate arm end locations
arm1_x_ind = find(strcmp(colnames,'LA_x'));
arm1_y_ind = find(strcmp(colnames,'LA_y'));
arm2_x_ind = find(strcmp(colnames,'RA_x'));
arm2_y_ind = find(strcmp(colnames,'RA_y'));
arm3_x_ind = find(strcmp(colnames,'MA_x'));
arm3_y_ind = find(strcmp(colnames,'MA_y'));
arms_x = [ary(:,arm1_x_ind);ary(:,arm2_x_ind);ary(:,arm3_x_ind)];
arms_y = [ary(:,arm1_y_ind);ary(:,arm2_y_ind);ary(:,arm3_y_ind)];

left_inds = arms_x>0 & arms_x<(v.width/3); %start with a wide search for the left arm
leftarm_x = median(arms_x(left_inds),'omitnan'); %estimate the left arm x with median
left_inds = arms_x>(leftarm_x-50) & arms_x<(leftarm_x+50); %narrow the search range
leftarm_x = median(arms_x(left_inds),'omitnan'); %re-estimate

middle_inds = arms_x>(v.width/3) & arms_x<(2*v.width/3); %repeat the above for the center and right arms
middlearm_x = median(arms_x(middle_inds),'omitnan');
middle_inds = arms_x>(middlearm_x-50) & arms_x<(middlearm_x+50);
middlearm_x = median(arms_x(middle_inds),'omitnan');

right_inds = arms_x>(2*v.width/3) & arms_x<v.width;
rightarm_x = median(arms_x(right_inds),'omitnan');
right_inds = arms_x>(rightarm_x-50) & arms_x<(rightarm_x+50);
rightarm_x = median(arms_x(right_inds),'omitnan');

leftarm_y = median(arms_y(left_inds),'omitnan'); %using the narrow search ranges, find the y's
middlearm_y = median(arms_y(middle_inds),'omitnan');
rightarm_y = median(arms_y(right_inds),'omitnan');

center_x = mean([leftarm_x middlearm_x rightarm_x]);
center_y = mean([leftarm_y middlearm_y rightarm_y]);

leftarm_length = get_dist(center_x,center_y,leftarm_x,leftarm_y);
middlearm_length = get_dist(center_x,center_y,middlearm_x,middlearm_y);
rightarm_length = get_dist(center_x,center_y,rightarm_x,rightarm_y);
arm_length_in_pixels = mean([leftarm_length middlearm_length rightarm_length]);
pixels_per_meter = arm_length_in_pixels/arm_length_in_meters;
arm_width = arm_width_in_meters*pixels_per_meter;
assert(all(([leftarm_length middlearm_length rightarm_length]./arm_length_in_pixels)>0.95),...
    'detected y-maze arm lengths are out of shape - tracking is probably too poor to continue')


%draw y-maze for validation
figure('Units','centimeters','Position',figure_size);
subplot(2,4,1)
mazeends_x = [leftarm_x center_x middlearm_x center_x rightarm_x];
mazeends_y =  [leftarm_y center_y middlearm_y center_y rightarm_y];
plot(mazeends_x,mazeends_y,'r');
xlim(xlimits);
ylim(ylimits)
title('detected ymaze skeleton')

%calculate some points along ymaze arms (for later, to only use mouse tracking near arms)
middlearm_ysign = (middlearm_y-center_y)/abs(middlearm_y-center_y);
num_armpts = length(center_y:middlearm_ysign*arm_width/4:middlearm_y);
leftarm_xstep = (center_x-leftarm_x)/(num_armpts-1);
middlearm_xstep = (center_x-middlearm_x)/(num_armpts-1);
rightarm_xstep = (center_x-rightarm_x)/(num_armpts-1);
leftarm_ystep = (center_y-leftarm_y)/(num_armpts-1);
middlearm_ystep = (center_y-middlearm_y)/(num_armpts-1);
rightarm_ystep = (center_y-rightarm_y)/(num_armpts-1);
mazepts_x = [leftarm_x:leftarm_xstep:center_x, middlearm_x:middlearm_xstep:center_x, rightarm_x:rightarm_xstep:center_x];
mazepts_y = [leftarm_y:leftarm_ystep:center_y, middlearm_y:middlearm_ystep:center_y, rightarm_y:rightarm_ystep:center_y];
num_mazepts = length(mazepts_x);

subplot(2,4,2)
plot(mazeends_x,mazeends_y);
hold on
scatter(mazepts_x,mazepts_y,70,'MarkerFaceColor','flat')
xlim(xlimits);
ylim(ylimits)
title('ymaze area');

%% delete all mouse tracking that's too far from the ymaze (change to NaNs)
leftear_x_ind = find(strcmp(colnames,'LE_x'));
leftear_y_ind = find(strcmp(colnames,'LE_y'));
rightear_x_ind = find(strcmp(colnames,'RE_x'));
rightear_y_ind = find(strcmp(colnames,'RE_y'));
nose_x_ind = find(strcmp(colnames,'NO_x'));
nose_y_ind = find(strcmp(colnames,'NO_y'));
tail_x_ind = find(strcmp(colnames,'TB_x'));
tail_y_ind = find(strcmp(colnames,'TB_y'));

nose_x = ary(:,nose_x_ind);
nose_y = ary(:,nose_y_ind);
nose_maze_dists = get_dist(repmat(nose_x,[1 num_mazepts]),repmat(nose_y,[1 num_mazepts]),repmat(mazepts_x,[num_frames 1]),repmat(mazepts_y,[num_frames 1]));
out_of_arm_idx = all(nose_maze_dists>arm_width/1.5,2);
nose_x(out_of_arm_idx) = nan;
nose_y(out_of_arm_idx) = nan;

leftear_x = ary(:,leftear_x_ind);
leftear_y = ary(:,leftear_y_ind);
leftear_maze_dists = get_dist(repmat(leftear_x,[1 num_mazepts]),repmat(leftear_y,[1 num_mazepts]),repmat(mazepts_x,[num_frames 1]),repmat(mazepts_y,[num_frames 1]));
out_of_arm_idx = all(leftear_maze_dists>arm_width/1.5,2);
leftear_x(out_of_arm_idx) = nan;
leftear_y(out_of_arm_idx) = nan;

rightear_x = ary(:,rightear_x_ind);
rightear_y = ary(:,rightear_y_ind);
rightear_maze_dists = get_dist(repmat(rightear_x,[1 num_mazepts]),repmat(rightear_y,[1 num_mazepts]),repmat(mazepts_x,[num_frames 1]),repmat(mazepts_y,[num_frames 1]));
out_of_arm_idx = all(rightear_maze_dists>arm_width/1.5,2);
rightear_x(out_of_arm_idx) = nan;
rightear_y(out_of_arm_idx) = nan;

tail_x = ary(:,tail_x_ind);
tail_y = ary(:,tail_y_ind);
tail_maze_dists = get_dist(repmat(tail_x,[1 num_mazepts]),repmat(tail_y,[1 num_mazepts]),repmat(mazepts_x,[num_frames 1]),repmat(mazepts_y,[num_frames 1]));
out_of_arm_idx = all(tail_maze_dists>arm_width/1.5,2);
tail_x(out_of_arm_idx) = nan;
tail_y(out_of_arm_idx) = nan;


%% estimate the mouse's size
mouse_size_in_pixels = median(get_dist(nose_x,nose_y,tail_x,tail_y),'omitnan');
% mouse_size_in_meters = mouse_size_in_pixels/pixels_per_meter;


%% delete all mouse tracking that's too far from the mouse's estimated position
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
plot(mazeends_x,mazeends_y);
hold on
scatter(smooth_body_x,smooth_body_y,1,[0 0 1]);%colors)
xlim(xlimits);
ylim(ylimits)
title('all mouse body positions')

subplot(2,4,6)
% plot(mazeends_x,mazeends_y);
hold on
xlim(xlimits);
ylim(ylimits)
plotframes = min([num_frames 30*fps]);
title(['first ' num2str(ceil(plotframes/fps)) ' s trajectory'])
plot(smooth_body_x(1:plotframes),smooth_body_y(1:plotframes),'b','LineWidth',0.5);

%% calculate results
%calculate mouse distance from 3 arms
num_arms = 3;
full_arms_x = repmat([leftarm_x middlearm_x rightarm_x],[num_frames 1]);
full_arms_y = repmat([leftarm_y middlearm_y rightarm_y],[num_frames 1]);
switch bodypart
    case 'head'
        bodypart_x = repmat(smooth_head_x,[1 num_arms]);
        bodypart_y = repmat(smooth_head_y,[1 num_arms]);
        arm_dists = get_dist(bodypart_x, bodypart_y, full_arms_x, full_arms_y);
        
    case 'body'
        bodypart_x = repmat(smooth_body_x,[1 num_arms]);
        bodypart_y = repmat(smooth_body_y,[1 num_arms]);
        arm_dists = get_dist(bodypart_x, bodypart_y, full_arms_x, full_arms_y);
        
    case 'nose'
        bodypart_x = repmat(smooth_nose_x,[1 num_arms]);
        bodypart_y = repmat(smooth_nose_y,[1 num_arms]);
        arm_dists = get_dist(bodypart_x, bodypart_y, full_arms_x, full_arms_y);
        
    case 'tail'
        bodypart_x = repmat(smooth_tail_x,[1 num_arms]);
        bodypart_y = repmat(smooth_tail_y,[1 num_arms]);
        arm_dists = get_dist(bodypart_x, bodypart_y, full_arms_x, full_arms_y);
        
    case 'head and tail'
        bodypart_x = repmat(smooth_head_x,[1 num_arms]);
        bodypart_y = repmat(smooth_head_y,[1 num_arms]);
        arm_dists1 = get_dist(bodypart_x, bodypart_y, full_arms_x, full_arms_y);
        bodypart_x = repmat(smooth_tail_x,[1 num_arms]);
        bodypart_y = repmat(smooth_tail_y,[1 num_arms]);
        arm_dists2 = get_dist(bodypart_x, bodypart_y, full_arms_x, full_arms_y);
        arm_dists(:,:,1) = arm_dists1;
        arm_dists(:,:,2) = arm_dists2;
        arm_dists = max(arm_dists,[],3,'omitnan'); %which bodypart is farther away
end

inleft = (arm_dists(:,1)/leftarm_length)<insidearm_threshold_fraction;
inmiddle = (arm_dists(:,2)/middlearm_length)<insidearm_threshold_fraction;
inright = (arm_dists(:,3)/rightarm_length)<insidearm_threshold_fraction;
inmultiple = inleft&inright | inleft&inmiddle | inright&inmiddle;

inleft(inmultiple) = 0;
inright(inmultiple) = 0;
inmiddle(inmultiple) = 0;
innone = ~inleft&~inright&~inmiddle;
assert((sum(inmultiple)/num_frames)<0.05, 'mouse tracking is too poor, need to make analysis even more robust');

time_in_locations = [sum(inleft) sum(inright) sum(inmiddle) sum(innone)]/fps;
sequence_of_locations = inleft + 2*inright + 3*inmiddle;
sequence_of_arms = sequence_of_locations;
sequence_of_arms(sequence_of_arms==0) = [];
sequence_of_arm_entries = sequence_of_arms;
sequence_of_arm_entries(logical([0; diff(sequence_of_arm_entries)==0])) = [];
num_arm_entries = length(sequence_of_arm_entries);

num_alternations = 0;
for i = 3:num_arm_entries
    if length(unique(sequence_of_arm_entries(i-2:i)))==3 %if the arms in this 3-arm sequence are all unique, it's an alternation
        num_alternations = num_alternations + 1;
    end
end
spontaneous_alternation_percent = 100*num_alternations/(num_arm_entries-2);

%calculate distance travelled (first remove nans; i.e. bad tracking)
nonan_sbx = smooth_body_x;
nonan_sbx(isnan(nonan_sbx)) = [];
nonan_sby = smooth_body_y;
nonan_sby(isnan(nonan_sby)) = [];
steps = get_dist(nonan_sbx(1:end-1),nonan_sby(1:end-1),nonan_sbx(2:end),nonan_sby(2:end));
distance_travelled = sum(abs(steps),'omitnan')/pixels_per_meter;
total_time = num_frames/fps;

results.number_of_frames = num_frames;
results.camera_fps = fps;
results.downsampled_fps = down_fps;
results.primary_bodypart = bodypart;
results.arm_threshold = arm_threshold;
results.total_time_in_seconds = total_time;
results.distance_travelled_in_meters = distance_travelled;
results.location_names = {'left arm (1)','right arm (2)','middle arm (3)','center (0)'};
results.seconds_in_locations = time_in_locations;
results.sequence_of_locations = sequence_of_locations';
results.sequence_of_arm_locations = sequence_of_arms';
results.sequence_of_arm_entries = sequence_of_arm_entries';
results.number_of_arm_entries = num_arm_entries;
results.number_of_alternations = num_alternations;
results.spontaneous_alternation_percent = spontaneous_alternation_percent;

%print results to figure
subplot(2,4,3)
set(gca,'Units','Normalized','Position',[0.6 0.05 0.4 0.95])
axis off
text(0,0.95,['Number of frames: ' num2str(num_frames)])
text(0,0.89,['Camera fps: ' num2str(fps)])
text(0,0.83,['Downsampled fps: ' num2str(down_fps)])
text(0,0.77,['Primary bodypart: ' bodypart])
text(0,0.71,['Arm threshold: ' num2str(arm_threshold) ' cm'])
text(0,0.65,['Total time: ' num2str(total_time, '%.1f') ' s'])
text(0,0.59,['Distance travelled: ' num2str(distance_travelled, '%.1f') ' m'])
text(0,0.53,'Names of locations: ')
text(0,0.48,['[' cell2str(results.location_names) ']'],'FontSize',8)
text(0,0.42,'Time in locations: ')
text(0,0.37,['[' num2str(time_in_locations, '%.1f  ') '] s'])
text(0,0.31,['Number of arm entries: ' num2str(num_arm_entries)])
text(0,0.25,'Sequence of arm entries: ')
text(0,0.20,num2str(sequence_of_arm_entries'),'FontSize',7)
text(0,0.14,['Number of alternations: ' num2str(num_alternations)])
text(0,0.08,['Spontaneous alternation: ' num2str(spontaneous_alternation_percent, '%.1f') ' %'])

saveas(gcf,[output_filename '.png'])
save([output_filename '.mat'],'results');
if ~show_figure
    close(gcf)
end

%write custom xls file of results
C = cell(10,2);
C(:,1) = {'number of frames', 'camera fps', 'downsampled fps', 'primary bodypart', 'arm threshold', 'total time (s)', 'distance travelled (m)', 'number of arm entries', 'number of alternations', 'spontaneous alternation (%)'};
C(:,2) = {num_frames, fps, down_fps, bodypart, arm_threshold, total_time, distance_travelled, num_arm_entries, num_alternations, spontaneous_alternation_percent};
writecell(C,[output_filename '.xls'],'Range','A1');
C = cell(2,5);
C(:,1) = {'location name', 'time in location (s)'};
C(1,2:5) = {results.location_names{1}, results.location_names{2}, results.location_names{3}, results.location_names{4}};
C(2,2:5) = {time_in_locations(1), time_in_locations(2), time_in_locations(3), time_in_locations(4)};
writecell(C,[output_filename '.xls'],'Range','A12');
writecell({'sequence of arm entries'},[output_filename '.xls'],'Range','A15');
writecell({sequence_of_arm_entries},[output_filename '.xls'],'Range','B15');