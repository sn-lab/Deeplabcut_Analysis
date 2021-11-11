function analyze_ymaze(h5_filename)
% FUNCTION analyze_ymaze(h5_filename)
%
% analyzes h5 files output from deeplabcut tracking of Y-maze videos using 
% the "Ymaze_2" training dataset. Analysis includes plotting the summary of 
% mouse and ymaze tracking data, calculating time the mouse spent within 
% each arm, the sequence of arm entries, and the spontaneous alternation %
%
%INPUTS
%h5_filename: full path/filename of a h5 file to analyze

fps = 30;
arm_length_in_meters = 0.365;
armnotch_length = 0.0762;
arm_width_in_meters = 0.07; 
figure_size = [10 10 25 10];
insidearm_threshold_fraction = (arm_length_in_meters-armnotch_length)/arm_length_in_meters; %how close to the center a mouse can be and still be identified as in an arm

%load tracking data
[save_dir, filename, ~] = fileparts(h5_filename);
ind = strfind(filename,'DLC_mobnet');
filename = filename(1:ind-1);
output_filename = fullfile(save_dir,filename);

data = h5read(h5_filename,'/df_with_missing/table');
ary = data.values_block_0';
%columns: top-right, top-left, bottom-left, bottome-right, right-ear, left-ear, nose, tail-base; each have 3 rows, for [x, y, likelihood]
colnames = {'RA_x','RA_y','RA_l','LA_x','LA_y','LA_l','BA_x','BA_y','BA_l','RE_x','RE_y','RE_l','LE_x','LE_y','LE_l','NO_x','NO_y','NO_l','TB_x','TB_y','TB_l'};
[num_frames, num_cols] = size(ary);

%vertically flip all y values
ystrings = repmat({'_y'},1,num_cols);
col_is_y = cellfun(@strfind,colnames,ystrings,'UniformOutput',false);
col_is_y = find(~cellfun(@isempty,col_is_y));
ary(:,col_is_y) = 1080-ary(:,col_is_y);


%% calculate arm end locations
arm1_x_ind = find(strcmp(colnames,'LA_x'));
arm1_y_ind = find(strcmp(colnames,'LA_y'));
arm2_x_ind = find(strcmp(colnames,'RA_x'));
arm2_y_ind = find(strcmp(colnames,'RA_y'));
arm3_x_ind = find(strcmp(colnames,'BA_x'));
arm3_y_ind = find(strcmp(colnames,'BA_y'));
arms_x = [ary(:,arm1_x_ind);ary(:,arm2_x_ind);ary(:,arm3_x_ind)];
arms_y = [ary(:,arm1_y_ind);ary(:,arm2_y_ind);ary(:,arm3_y_ind)];

left_inds = arms_x>25 & arms_x<275; %start with a wide search for the left arm
leftarm_x = median(arms_x(left_inds),'omitnan'); %estimate the left arm x with median
left_inds = arms_x>(leftarm_x-50) & arms_x<(leftarm_x+50); %narrow the search range
leftarm_x = median(arms_x(left_inds),'omitnan'); %re-estimate

center_inds = arms_x>575 & arms_x<825; %repeat the above for the center and right arms
centerarm_x = median(arms_x(center_inds),'omitnan');
center_inds = arms_x>(centerarm_x-50) & arms_x<(centerarm_x+50);
centerarm_x = median(arms_x(center_inds),'omitnan');

right_inds = arms_x>1175 & arms_x<1425;
rightarm_x = median(arms_x(right_inds),'omitnan');
right_inds = arms_x>(rightarm_x-50) & arms_x<(rightarm_x+50);
rightarm_x = median(arms_x(right_inds),'omitnan');

leftarm_y = median(arms_y(left_inds),'omitnan'); %using the narrow search ranges, find the y's
centerarm_y = median(arms_y(center_inds),'omitnan');
rightarm_y = median(arms_y(right_inds),'omitnan');

center_x = mean([leftarm_x centerarm_x rightarm_x]);
center_y = mean([leftarm_y centerarm_y rightarm_y]);

leftarm_length = get_dist(center_x,center_y,leftarm_x,leftarm_y);
centerarm_length = get_dist(center_x,center_y,centerarm_x,centerarm_y);
rightarm_length = get_dist(center_x,center_y,rightarm_x,rightarm_y);
arm_length_in_pixels = mean([leftarm_length centerarm_length rightarm_length]);
pixels_per_meter = arm_length_in_pixels/arm_length_in_meters;
arm_width = arm_width_in_meters*pixels_per_meter;
assert(all(([leftarm_length centerarm_length rightarm_length]./arm_length_in_pixels)>0.95),...
    'detected y-maze arm lengths are out of shape - tracking is probably too poor to continue')


%draw y-maze for validation
figure('Units','centimeters','Position',figure_size);
subplot(2,4,1)
mazeends_x = [leftarm_x center_x centerarm_x center_x rightarm_x];
mazeends_y =  [leftarm_y center_y centerarm_y center_y rightarm_y];
plot(mazeends_x,mazeends_y,'r');
xlim([-100 1550]);
ylim([-280 1370])
title('detected ymaze skeleton')

%calculate some points along ymaze arms (for later, to only use mouse tracking near arms)
centerarm_ysign = (centerarm_y-center_y)/abs(centerarm_y-center_y);
num_armpts = length(center_y:centerarm_ysign*arm_width/4:centerarm_y);
leftarm_xstep = (center_x-leftarm_x)/(num_armpts-1);
centerarm_xstep = (center_x-centerarm_x)/(num_armpts-1);
rightarm_xstep = (center_x-rightarm_x)/(num_armpts-1);
leftarm_ystep = (center_y-leftarm_y)/(num_armpts-1);
centerarm_ystep = (center_y-centerarm_y)/(num_armpts-1);
rightarm_ystep = (center_y-rightarm_y)/(num_armpts-1);
mazepts_x = [leftarm_x:leftarm_xstep:center_x, centerarm_x:centerarm_xstep:center_x, rightarm_x:rightarm_xstep:center_x];
mazepts_y = [leftarm_y:leftarm_ystep:center_y, centerarm_y:centerarm_ystep:center_y, rightarm_y:rightarm_ystep:center_y];
num_mazepts = length(mazepts_x);

subplot(2,4,2)
plot(mazeends_x,mazeends_y);
hold on
scatter(mazepts_x,mazepts_y,70,'MarkerFaceColor','flat')
xlim([-100 1550]);
ylim([-280 1370])
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
head_x = median([nose_x leftear_x rightear_x],2,'omitnan');
head_y = median([nose_y leftear_y rightear_y],2,'omitnan');
body_x =  mean([head_x tail_x],2,'omitnan');
body_y =  mean([head_y tail_y],2,'omitnan');
smooth_body_x = rolling_average(body_x,1,round(fps/2),'median');
smooth_body_y = rolling_average(body_y,1,round(fps/2),'median');

subplot(2,4,5)
plot(mazeends_x,mazeends_y);
hold on
scatter(body_x,body_y,1,[0 0 1]);%colors)
xlim([-100 1550]);
ylim([-280 1370])
title('all mouse body positions')

subplot(2,4,6)
% plot(mazeends_x,mazeends_y);
hold on
xlim([-100 1550]);
ylim([-280 1370])
title('first 30 s trajectory')
plot(body_x(1:fps*30),body_y(1:fps*30),'b','LineWidth',0.5);

% p1 = scatter(body_x(1),body_y(1),3,[1 0 0]);
% for i = 1:10:min([1000 num_frames])
%     p1.XData = median(smooth_body_x(i:i+4),'omitnan');
%     p1.YData = median(smooth_body_y(i:i+4),'omitnan');
%     
%     % Capture the plot as an image 
%     frame = getframe(gcf); 
%     im = frame2im(frame); 
%     [imind,cm] = rgb2ind(im,256); 
%     % Write to the GIF File 
%     if i == 1 
%       imwrite(imind,cm,[output_filename '.gif'],'gif', 'Loopcount',inf,'DelayTime',0.033); 
%     else 
%       imwrite(imind,cm,[output_filename '.gif'],'gif','WriteMode','append'); 
%     end 
% end

%% calculate results
%calculate mouse distance from 3 arms
leftarm_dist = get_dist(smooth_body_x, smooth_body_y, repmat(leftarm_x,[num_frames 1]), repmat(leftarm_y,[num_frames 1]));
rightarm_dist = get_dist(smooth_body_x, smooth_body_y, repmat(rightarm_x,[num_frames 1]), repmat(rightarm_y,[num_frames 1]));
centerarm_dist = get_dist(smooth_body_x, smooth_body_y, repmat(centerarm_x,[num_frames 1]), repmat(centerarm_y,[num_frames 1]));

inleft = (leftarm_dist/leftarm_length)<insidearm_threshold_fraction;
inright = (rightarm_dist/rightarm_length)<insidearm_threshold_fraction;
incenter = (centerarm_dist/centerarm_length)<insidearm_threshold_fraction;

inmultiple = inleft&inright | inleft&incenter | inright&incenter;

inleft(inmultiple) = 0;
inright(inmultiple) = 0;
incenter(inmultiple) = 0;
innone = ~inleft&~inright&~incenter;

assert((sum(inmultiple)/num_frames)<0.05, 'mouse tracking is too poor, need to make analysis even more robust');

time_in_locations = [sum(inleft) sum(inright) sum(incenter) sum(innone)]/fps;

sequence_of_locations = inleft + 2*inright + 3*incenter;

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

steps = get_dist(smooth_body_x(1:end-1),smooth_body_y(1:end-1),smooth_body_x(2:end),smooth_body_y(2:end));
distance_travelled = sum(abs(steps),'omitnan')/pixels_per_meter;

results.number_of_frames = num_frames;
results.assumed_camera_fps = fps;
results.total_time_in_seconds = num_frames/fps;
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
text(0,0.9,['Number of frames: ' num2str(results.number_of_frames)])
text(0,0.82,['Assumed camera fps: ' num2str(results.assumed_camera_fps)])
text(0,0.74,['Total time: ' num2str(results.total_time_in_seconds, '%.1f') ' s'])
text(0,0.66,['Distance travelled: ' num2str(results.distance_travelled_in_meters, '%.1f') ' m'])
text(0,0.58,'Names of locations: ')
text(0,0.53,['[' cell2str(results.location_names) ']'],'FontSize',8)
text(0,0.45,'Time in locations: ')
text(0,0.4,['[' num2str(results.seconds_in_locations, '%.1f  ') '] s'])
text(0,0.32,['Number of arm entries: ' num2str(results.number_of_arm_entries)])
text(0,0.24,'Sequence of arm entries: ')
text(0,0.19,num2str(results.sequence_of_arm_entries),'FontSize',7)
text(0,0.11,['Number of alternations: ' num2str(results.number_of_alternations)])
text(0,0.03,['Spontaneous alternation: ' num2str(results.spontaneous_alternation_percent, '%.1f') ' %'])


saveas(gcf,[output_filename '.png'])
save([output_filename '.mat'],'results');


%write custom xls file of results
C = cell(7,2);
C(:,1) = {'number of frames', 'assumed camera fps', 'total time (s)', 'distance travelled (m)', 'number of arm entries', 'number of alternations', 'spontaneous alternation (%)'};
C(:,2) = {num_frames, fps, num_frames/fps, distance_travelled, num_arm_entries, num_alternations, spontaneous_alternation_percent};
writecell(C,[output_filename '.xls'],'Range','A1');
C = cell(2,5);
C(:,1) = {'location name', 'time in location (s)'};
C(1,2:5) = {results.location_names{1}, results.location_names{2}, results.location_names{3}, results.location_names{4}};
C(2,2:5) = {results.seconds_in_locations(1), results.seconds_in_locations(2), results.seconds_in_locations(3), results.seconds_in_locations(4)};
writecell(C,[output_filename '.xls'],'Range','A9');
writecell({'sequence of arm entries'},[output_filename '.xls'],'Range','A12');
writecell({sequence_of_arm_entries},[output_filename '.xls'],'Range','B12');