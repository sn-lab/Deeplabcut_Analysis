clear all
% function analyze_fixedwalking(h5_filename, param)
% FUNCTION analyze_fixedwalking(h5_filename, param)
%
% analyzes h5 files output from deeplabcut tracking of spine-fixed (or 
% tail-held) videos using the "CAM3" training dataset. `Analysis includes
% ________________
%
%INPUTS
%h5_filename: full path/filename of a h5 file to analyze
%param (optional):
%     
%     fps = framerate of the original video recordings
%     down_fps = effective fps after rolling average performed to remove some tracking errors

% if nargin<2
    %set default parameters
    h5_filename = 'D:\SN Lab\Spinal Cord\walking wireframe\TF\';
%     nnn = '062321_mouse419008_TF_CAM3DeepCut_resnet50_062121_roundTM_CAM3_TRAINJun21shuffle1_709000.h5';
    nnn = '062321_mouse419008_TF_CAM3DeepCut_resnet50_062121_roundTM_CAM3_TRAINJun21shuffle1_709000.h5';
    h5_filename = [h5_filename nnn];
    bodypart = 'toe'; %bodypart used to find gait cycle start/stop times
    fps = 122.5; %framerate of the original video recordings
    down_fps = 30; %effective fps after rolling average performed to remove some tracking errors
    show_figure = 1; %whether to show (1) or close (0) the results figure when finished
    manually_inspect = false;
% else
%     bodypart = param.bodypart;
%     fps = param.fps;
%     down_fps = param.down_fps;
%     show_figure = param.show_figure;
% end

%more parameters:
max_stride_speed = 30; %max strides per second (anything faster is assumed to be a tracking error)
min_stride_range = 0.5; %approx. fraction of max stride range to be still consisdered a stride
min_stance_speed = 0.5; %min strides lengths per second of the stance phase (anything slower is ignored)
min_swing_speed = 4; %min stride lengths per second of swing phase (which is always quick)
max_standing_speed = 0.25; %max stride lengths per second to count as standing still
max_error_duration = 0.25; %how long a tracking error can be and still be detected/omitted
quality_threshold = 0.75; %minimum quality of tracking (0-1) to include
num_gait_bins = 11; %number of data bins to align the all gaits to

%secondary parameters
max_stance_frames = round(fps/min_stance_speed);
max_swing_frames = round(fps/min_swing_speed);
label_names = {'UpperHip','Hip','Knee','Ankle','Midway','Toe'}';
colnames = {'UpperHipX','UpperHipY','UpperHipQ','HipX','HipY','HipQ','KneeX','KneeY','KneeQ','AnkleX','AnkleY','AnkleQ','MidX','MidY','MidQ','ToeX','ToeY','ToeQ'};
gait_var_name = [bodypart 'X'];
gait_var_name(1) = upper(gait_var_name(1));
gait_var_col = find(strcmp(colnames,gait_var_name));
assert(any(contains({'AnkleX','ToeX'},gait_var_name)),'Only "toe" and "ankle" bodyparts should be used for fixedwalking analysis')

figure_pos = [10 10];
figure_size = [25 10];

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

%load tracking data
data = h5read(h5_filename,'/df_with_missing/table');
ary = data.values_block_0';
[num_frames, num_datatypes] = size(ary);
timevec = (1:num_frames)/fps;

%get columns subscripts of certain labels
x_col = find(contains(colnames,'X'));
y_col = find(contains(colnames,'Y'));

%get upper limits of label coordinates
arymax = max(ary,[],'omitnan');
v.width = max(arymax(x_col),[],'omitnan')+1;
v.height = max(arymax(y_col),[],'omitnan')+1;

%vertically flip all y values (video coordinates are upside-down from plot coordinates)
ary(:,y_col) = v.height-ary(:,y_col);

%if mouse is facing left, flips all x values (so that x local minima are the starts of swing phases) 
high_quality_frames = ary(:,strcmp(colnames,'ToeQ'))>quality_threshold;
toeX = ary(:,strcmp(colnames,'ToeX'));
toeX = mean(toeX(high_quality_frames));
high_quality_frames = ary(:,strcmp(colnames,'AnkleQ'))>quality_threshold;
ankleX = ary(:,strcmp(colnames,'AnkleX'));
ankleX = mean(ankleX(high_quality_frames));
if toeX<ankleX
    ary(:,x_col) = v.width-ary(:,x_col);
end

%get image/plot limits
largest_range = max([v.width v.height]);
xlimits = [(v.width-largest_range) (v.width+largest_range)]/2;
ylimits = [(v.height-largest_range) (v.height+largest_range)]/2;

%find poorly tracked data
ary_errors = false(size(ary));
poor_inds = ary(:,contains(colnames,'Q'))<quality_threshold;
ary_errors(:,x_col) = poor_inds;
ary_errors(:,y_col) = poor_inds;
percent_error = 100*(sum(poor_inds)./size(poor_inds,1));
fprintf('Percent errors: %2.1f upperhip, %2.1f hip, %2.1f knee, %2.1f ankle, %2.1f midway, %2.1f toe\n',percent_error);
% assert(all(percent_error<50),'tracking quality is too poor to continue (error>50%)');

%check whether ankle or toe should be used
ankle_ind = find(strcmpi(label_names,'Ankle'));
toe_ind = find(strcmpi(label_names,'Toe'));
gait_var_ind = find(strcmpi(label_names,bodypart));
gait_error = percent_error(gait_var_ind);
best_error = min(percent_error([ankle_ind toe_ind]));
if gait_error>best_error
    %%%%%ask for user input to decide to switch?
    gait_var_ind = -gait_var_ind + ankle_ind + toe_ind;
    bodypart = label_names{gait_var_ind};
    gait_var_name = [bodypart 'X'];
    gait_var_col = find(strcmp(colnames,gait_var_name));
end

%create short raw data wireframe movie to verify that tracking is OK
if manually_inspect
    vid = VideoWriter(fullfile(save_dir,'rawwireframe.avi'));
    open(vid);
    figure();
    title('raw wireframe')
    p1 = plot(ary(1,x_col),ary(1,y_col),'-ok');
    for i = 1:round(2.5*fps)%1:round(fps/30):round(5*fps)
        p1.XData = ary(i,x_col);
        p1.YData = ary(i,y_col);
        xlim(xlimits)
        ylim(ylimits)
        pause(1/fps);
        frame = getframe(gcf);
        writeVideo(vid,frame);
    end
    close(vid)
    answer = input('Is tracking OK? (type yes to continue with analysis): ','s');
    assert(strncmpi('y',answer,1),'Analysis cancelled.')
    close(gcf)
end


%plot raw gait variable data
tmpX = ary(:,gait_var_col)';
figure('Units','centimeters','Position',[figure_pos figure_size]);
subplot(3,2,1)
plot(timevec,tmpX,'b');
xlabel('time (s)')
ylabel('coordinate (pixels)')
title([gait_var_name ' raw data'])

%remove poorly tracked data
tmpX = ary(:,x_col);
tmpX(poor_inds) = nan;
ary(:,x_col) = tmpX;
tmpY = ary(:,y_col);
tmpY(poor_inds) = nan;
ary(:,y_col) = tmpY;

%use gait var ind to estimate stride range (in pixels)
tmpX = ary(:,gait_var_col)'; %recalculate gait variable data
%%%%%%%make this relative to upperhip? (or something that doesn't move)
% stable_name = 'UpperHipX';
% stable_col = find(strcmp(colnames,stable_name));
% tmpX = ary(:,gait_var_col)' - ary(:,stable_col)'; %recalculate gait variable data
pixels_per_stride = prctile(tmpX,98) - prctile(tmpX,2); 

%look for additional tracking errors using max speed constraint (doesn't remove too many more) 
err_window = round(max_error_duration*fps);
for k = 1:length(x_col) %loop for every label type
    x_ind = x_col(k);
    tmpX = ary(:,x_ind)';

    %calculate instantaneous velocity of label to find when label velocity is faster than possible
    vel = diff(tmpX)*fps/pixels_per_stride; %strides per second in x dimension
    errors = [false abs(vel)>max_stride_speed];
    directions = vel./abs(vel);
    inds = 1:length(vel);
    
    %find longer errors (tracking stays lost for many frames in a row)
    error_dirs = directions(errors);
    error_inds = inds(errors);
    returns = find(abs(diff(error_dirs))==2);
    return_inds = [error_inds(returns); error_inds(returns+1)];
    span_lengths = diff(return_inds);
    too_long_spans = find(span_lengths>err_window);
    return_inds(:,too_long_spans) = [];
    for e = 1:size(return_inds,2)
        errors(return_inds(1,e):return_inds(2,e)) = true;
    end
    
    %record tracking errors
    ary_errors(:,x_ind) = ary_errors(:,x_ind) | errors';
    ary_errors(:,x_ind+1) = ary_errors(:,x_ind+1) | errors';
end

%remove tracking errors from data array
ary(ary_errors) = nan;
percent_errors_per_label = 100*sum(ary_errors(:,x_col))/num_frames;

%plot gait variable data with errors removed
subplot(3,2,3)
plot(timevec,ary(:,gait_var_col)','b');
xlabel('time (s)')
ylabel('coordinate (pixels)')
title('errors removed')
    
%replace nans with interpolated values
ary = fillmissing(ary,'makima');
tmpX = ary(:,gait_var_col)'; %re-calculate
pixels_per_stride = prctile(tmpX,98) - prctile(tmpX,2); %re-calculate

%plot gait variable data with errors interpolated
subplot(3,2,5)    
plot(timevec,ary(:,gait_var_col)','b');
xlabel('time (s)')
ylabel('coordinate (pixels)')
title('errors interpolated')

%calculate joint angles
jointnames = {'hip','knee','ankle'};
up = {'UpperHip', 'Hip', 'Knee'}; %upper point
mp = {'Hip', 'Knee', 'Ankle'}; %middle point
lp = {'Knee', 'Ankle', 'Mid'}; %lower point

num_joints = length(jointnames);
jointangles = nan(num_frames,num_joints);
for a = 1:num_joints 
    upx = ary(:,strcmp(colnames,[up{a} 'X'])); %x-coord of upper point
    upy = ary(:,strcmp(colnames,[up{a} 'Y']));
    mpx = ary(:,strcmp(colnames,[mp{a} 'X']));
    mpy = ary(:,strcmp(colnames,[mp{a} 'Y']));
    lpx = ary(:,strcmp(colnames,[lp{a} 'X']));
    lpy = ary(:,strcmp(colnames,[lp{a} 'Y']));

    rupx = upx-mpx; %upper x-coord relative to mp
    rupy = upy-mpy; 
    rlpx = lpx-mpx; %lower x-coord relative to mp
    rlpy = lpy-mpy;

    upa = (atan2(rupy,rupx))*180/pi; %angle of mp-rp from horizontal
    lpa = (atan2(rlpy,rlpx))*180/pi;

    jointangles(:,a) = upa-lpa;
end
jointangles(:,1) = 90-jointangles(:,1); %hip: positive = flex from 90deg
jointangles(:,3) = 90-jointangles(:,3); %ankle: positive = dorsiflex from 90deg
tmpj = jointangles(:,2); %knee: positive = flex from 180deg
tmpj(tmpj<-180) = tmpj(tmpj<-180)+360;
tmpj(tmpj>180) = tmpj(tmpj>180)-360;
jointangles(:,2) = 180+tmpj; 

%find gait cycles using local minimums and maximums with constraints
min_gait_separation = round(fps/(0.5*max_stride_speed));
max_gait_separation = round(fps/(0.5*min_stance_speed)); %the slowest allowed gait cycle (aka minimum separation between bouts of gaits)
min_stride_size = min_stride_range*pixels_per_stride; %the smallest movement (in pixels) that can be classified as a good stride
mins = find(islocalmin(tmpX,'MinSeparation',min_gait_separation));
maxs = find(islocalmax(tmpX,'MinSeparation',min_gait_separation));

%plot all local minima/maxima
subplot(3,2,2)
plot(timevec,tmpX,'b');
hold on
scatter(timevec(mins),tmpX(mins),'og');
scatter(timevec(maxs),tmpX(maxs),'or');
xlabel('time (s)')
ylabel('coordinate (pixels)')
title('local minima/maxima')

%%%%%%%%%%%possible change of strategy: only use mins to find swings and go
%%%%%%%%%%%from there, don't use maxs at all

%use local maxima/minima to find good gaits
%sort mins and maxs together, in chronological order
minmax_inds = [mins maxs]; %combined mins and maxs
minmax_dir_log = [false(size(mins)) true(size(maxs))]; %order of combined (0=mins, 1=maxs)
[minmax_inds, sorts] = sort(minmax_inds); %sort mins and maxs together
minmax_dir_log = minmax_dir_log(sorts); %keep track of which is a max and which is a min
minmax_tmpX = tmpX(minmax_inds); %value of all minima/maxima

%if there are 2 mins in a row, close together, delete the higher one
mintomin_first = [diff(minmax_dir_log)==0 false] & minmax_dir_log==0;
bigger_than_next = [diff(minmax_tmpX)<0 false];
close_to_next = [diff(minmax_inds)<(max_swing_frames-1) false];
mm_inds_to_erase = mintomin_first & bigger_than_next & close_to_next;
mintomin_second = [false diff(minmax_dir_log)==0] & minmax_dir_log==0;
bigger_than_prev = [false diff(minmax_tmpX)>0];
close_to_previous = [false diff(minmax_inds)<(max_swing_frames-1)];
mm_inds_to_erase = mm_inds_to_erase | (mintomin_second & bigger_than_prev & close_to_previous);
minmax_inds(mm_inds_to_erase) = [];
minmax_dir_log(mm_inds_to_erase) = [];
minmax_tmpX = tmpX(minmax_inds); %recalculate

%%%%%%%%%%%%%%possible to-do: if there are 2 swings in a row, close together, and the 2nd has both points higher than the 1st, delete the middle 2 points 

%look for good swings (min to max)
minmax_diff_pixels = abs(diff(minmax_tmpX));
minmax_diff_frames = diff(minmax_inds);
swing_inds = find(diff(minmax_dir_log)==1);
swing_range = minmax_diff_pixels(swing_inds);
swing_frames = minmax_diff_frames(swing_inds);
good_swings = swing_range>=min_stride_size & swing_frames<=max_swing_frames;
good_swing_inds = swing_inds(good_swings);
minmax_inds_good_swings = minmax_inds(good_swing_inds); %indices of start of swings

%look for good stances preceding all swings
%in window before all good swings, look for highest max?
bad_swings = false(size(good_swing_inds));
gait_cycle_inds = nan(3,length(good_swing_inds));
for s = 1:length(good_swing_inds)
    swing_start_ind = good_swing_inds(s); %ind of minmax corredponding to the current swing start

    if swing_start_ind>1 %if it =1, then there is no preceding stance phase start and this isn't a good gait
        %get all inds of minmax that are within the possible preceding stance phase 
        swing_start_frame = minmax_inds(swing_start_ind);
        earliest_stance_frame = swing_start_frame-max_stance_frames;
        if s>1
            previous_swing_start_frame = minmax_inds(good_swing_inds(s-1));
            earliest_stance_frame = max([earliest_stance_frame previous_swing_start_frame]);
        end
        stance_start_ind = find(minmax_inds>earliest_stance_frame & minmax_inds<swing_start_frame);
        
        if ~isempty(stance_start_ind)
            %pick the single biggest stance start ind
            stance_start_tmpX = minmax_tmpX(stance_start_ind);
            [~, biggest_max] = max(stance_start_tmpX);
            stance_start_ind = stance_start_ind(biggest_max);
            
            %get the stance after this swing (use differential if the next max isn't close enough)
            latest_next_swing_frame = swing_start_frame+max_stance_frames+max_swing_frames;
            if s<length(good_swing_inds)
                next_swing_start_frame = minmax_inds(good_swing_inds(s+1));
            	latest_next_swing_frame = min([latest_next_swing_frame next_swing_start_frame]);
            end
            swing_stop_ind =  find(minmax_inds<latest_next_swing_frame & minmax_inds>swing_start_frame);
            if ~isempty(swing_stop_ind)
                %pick the [next] single biggest stance start ind
                stance_start_tmpX = minmax_tmpX(swing_stop_ind);
                [~, biggest_max] = max(stance_start_tmpX);
                swing_stop_ind = swing_stop_ind(biggest_max);
                
                gait_cycle_inds(1,s) = minmax_inds(stance_start_ind);
                gait_cycle_inds(2,s) = minmax_inds(swing_start_ind);
                gait_cycle_inds(3,s) = minmax_inds(swing_stop_ind);
            end
        end
    end
end
gait_cycle_inds(isnan(gait_cycle_inds)) = [];
num_gaits = size(gait_cycle_inds,2);

%plot all good gait cycle timepoints
subplot(3,2,4)
plot(timevec,tmpX,'b');
hold on
scatter(timevec(gait_cycle_inds(1,:)),tmpX(gait_cycle_inds(1,:)),'xg');
scatter(timevec(gait_cycle_inds(2,:)),tmpX(gait_cycle_inds(2,:)),'og');
scatter(timevec(gait_cycle_inds(3,:)),tmpX(gait_cycle_inds(3,:)),'xr');
xlabel('time (s)')
ylabel('coordinate (pixels)')
title([num2str(num_gaits) ' gait cycles'])

%create vector of gait cycle bins (0-nbins during gaits, nan otherwise)
gait_phase_vector = nan([1 num_frames]);
stance_log = false([1 num_frames]);
swing_log = false([1 num_frames]);
for g = 1:num_gaits
    relative_inds = (1:num_frames)-gait_cycle_inds(1,g); %[-3 -2 -1 0 1 2 3 4 5]
    phase = relative_inds/diff(gait_cycle_inds([1 3],g)); %[-1 -0.7 -0.3 0 0.3 0.7 1 1.3 1.7]
    current_bin = round(1+(phase*(num_gait_bins-1)));
    current_bin(current_bin<1 | current_bin>num_gait_bins) = nan;
    gait_phase_vector = min([gait_phase_vector; current_bin],[],'omitnan');
    stance_log(gait_cycle_inds(1,g):gait_cycle_inds(2,g)-1) = true;
    swing_log(gait_cycle_inds(2,g):gait_cycle_inds(3,g)-1) = true;
end

%get logical indices during gaits
gait_log = ~isnan(gait_phase_vector);

%get logical indices of standing still
vel = ([0 diff(tmpX)])*fps; %instantaneous velocity of the gait variable (pixels/sec)
vel = rolling_average(vel,2,round(fps/down_fps),'mean'); 
standing_log = vel<=(max_standing_speed*pixels_per_stride);
standing_log = standing_log & ~gait_log; %omit any slow moving times during gaits

%combine logical indices of standing, stance, and swing into a single state vector
state_vector = standing_log + 2*stance_log + 3*swing_log;

subplot(3,2,6)
standing_color = [0.7 0.7 0.7];
gait_color = [1 0 1];
patch_starts = find(diff([0 standing_log])==1);
patch_ends = find(diff([standing_log 0])==-1);
ylims = [0 v.height];
for p = 1:length(patch_starts)
    patch(timevec([patch_starts(p) patch_ends(p) patch_ends(p) patch_starts(p)]),[ylims(1) ylims ylims(2)],'k','FaceColor',standing_color,'EdgeColor','none')
    if p==1; hold on; end
end
for p = 1:num_gaits
    patch(timevec([gait_cycle_inds(1,p) gait_cycle_inds(3,p) gait_cycle_inds(3,p) gait_cycle_inds(1,p)]),[ylims(1) ylims ylims(2)],'k','FaceColor',gait_color,'EdgeColor','none')
end
ylim(ylims);
plot(timevec,tmpX,'b');
xlabel('time (s)')
ylabel('coordinate (pixels)')
title('standing vs walking')
% saveas(gcf,[output_filename '_AnalysisDetails.png']);
% close(gcf)

%align raw data and joint angle arrays for all gaits
jointangles_gaits = nan([num_gaits, num_gait_bins, num_joints]);
rawdata_gaits = nan([num_gaits, num_gait_bins, num_datatypes]);
alignedTime = 0:(1/(num_gait_bins-1)):1;
for g = 1:num_gaits
    start_ind = max([gait_cycle_inds(1,g)-10, 1]);
    stop_ind = min([gait_cycle_inds(3,g)+10 num_frames]);
    num_inds = diff(gait_cycle_inds([1 3],g));
    inds = start_ind:stop_ind;
    unalignedTime = (inds-gait_cycle_inds(1,g))/num_inds;
    for d = 1:num_datatypes
        unalignedData = ary(inds,d)';
        rawdata_gaits(g,:,d) = align_data(unalignedData,unalignedTime,alignedTime);
    end
    for j = 1:num_joints
        unalignedData = jointangles(inds,j)';
        jointangles_gaits(g,:,j) = align_data(unalignedData,unalignedTime,alignedTime);
    end
end

%calculate frequency/duration of each gait
gait_duration = diff(gait_cycle_inds([1 3],:))/fps;
gait_frequency = 1./gait_duration;
stance_range = -diff(tmpX(gait_cycle_inds([1 2],:))); %pixels
stance_duration = diff(gait_cycle_inds([1 2],:))/fps; %seconds
stance_speed = stance_range./stance_duration; %pixels per second
swing_range = diff(tmpX(gait_cycle_inds([2 3],:)));
swing_duration = diff(gait_cycle_inds([2 3],:))/fps; 
swing_speed = swing_range./swing_duration; 

%calculate walking velocity vector (excluding the swing phase)
walking_speed_vector = [0 -diff(tmpX)]*fps; %pixel/sec
walking_speed_vector(walking_speed_vector<0) = 0;
for g = 1:num_gaits
    walking_speed_vector(gait_cycle_inds(1,g):gait_cycle_inds(3,g)) = stance_speed(g);
end
walking_speed_vector = rolling_average(walking_speed_vector,2,round(fps/down_fps),'mean');
walking_speed_vector = fillmissing(walking_speed_vector,'makima');

%save/plot final results
total_time = num_frames/fps;
total_standing_time = sum(standing_log)/fps;
total_gait_time = sum(gait_log)/fps;

results.number_of_frames = num_frames;
results.camera_fps = fps;
results.downsampled_fps = down_fps;
results.primary_bodypart = bodypart;
results.total_time_in_seconds = total_time;
results.total_gait_time = total_gait_time;
results.total_standing_time = total_standing_time;
results.label_names = label_names';
results.percent_errors_per_label = percent_errors_per_label;
results.data_names = colnames;
results.data_array = ary;
results.data_array_errors = ary_errors;
results.jointnames = jointnames;
results.jointangle_array = jointangles;
results.state_vector = state_vector;
results.walking_speed_vector = walking_speed_vector;
results.gait_phase_vector = gait_phase_vector;
results.num_gaits = num_gaits;
results.gait_cycle_inds = gait_cycle_inds;
results.gait_freq = gait_frequency;
results.stance_speed = stance_speed;
results.data_in_gait_bins = rawdata_gaits;
results.jointangles_in_gait_bins = jointangles_gaits;

parameters.max_stride_speed = max_stride_speed;
parameters.min_stride_range = min_stride_range;
parameters.min_stance_speed = min_stance_speed;
parameters.min_swing_speed = min_swing_speed;
parameters.max_standing_speed = max_standing_speed;
parameters.max_error_duration = max_error_duration; 
parameters.quality_threshold = quality_threshold; 
parameters.num_gait_bins = num_gait_bins;
results.parameters = parameters;

%print results to figure
figure('Units','centimeters','Position',[figure_pos 1.5*figure_size]);
subplot(2,4,1)
% plot of wireframe of gait cycle; average of all gaits, wireframe per bin
mean_rawdata_gaits = squeeze(mean(rawdata_gaits));
for b = 1:num_gait_bins
    x = mean_rawdata_gaits(b,x_col);
    y = mean_rawdata_gaits(b,y_col);
    colors = [0.7-0.7*(b/num_gait_bins) 0.7-0.7*(b/num_gait_bins) 1];
    plot(x,y,'-ok','Color',colors)
    if b==1; hold on; end
end
xlabel('x coord (pixels)')
xlabel('y coord (pixels)')
title('mean gait wireframe')

subplot(2,4,2)
% plot of joint angle by gait cycle; average/sd of all
mean_jointangles_gaits = squeeze(mean(jointangles_gaits));
std_jointangles_gaits = squeeze(std(jointangles_gaits));
x = 0:100/(num_gait_bins-1):100;
colors = linspecer(3);
for j = 1:3
    y = mean_jointangles_gaits(:,j)';
    s = std_jointangles_gaits(:,j)';
    patch([x fliplr(x)],[y-s fliplr(y+s)],'k','FaceColor',colors(j,:),'EdgeColor','none','FaceAlpha',0.3)
    hold on
    p(j) = plot(x,y,'-k','Color',colors(j,:));
end
legend(p,{'hip','knee','ankle'})
xlabel('gait cycle (%)')
ylabel('joint angle (deg)')
title('avg gait joint angles')

subplot(2,4,5)
% plot histograms of gait cycle frequency/duration
histogram(gait_frequency,0.5:15)
xlabel('frequency (Hz)')
ylabel('number of gaits')
title('gait cycle frequencies')

%plot walking speed over time
subplot(2,4,6)
plot(timevec,walking_speed_vector,'b');
xlabel('time (s)')
ylabel('speed (pixels/sec)')
title('walking/running speed')
% %plot histogram of stance phase speed
% subplot(2,4,6)
% histogram(stance_speed)
% xlabel('speed (pixels/sec)')
% ylabel('number of gaits')
% title('stance phase speeds')

subplot(2,4,3)
set(gca,'Units','Normalized','Position',[0.6 0.05 0.4 0.95]);
axis off
text(0,0.95,['Number of frames: ' num2str(num_frames)]);
text(0,0.89,['Camera fps: ' num2str(fps)]);
text(0,0.83,['Downsampled fps: ' num2str(down_fps)]);
text(0,0.77,['Primary bodypart: ' bodypart]);
text(0,0.71,['Total time: ' num2str(total_time, '%.1f') ' s']);
text(0,0.65,['Number of gaits: ' num2str(num_gaits)]);
text(0,0.59,['Median gait frequency: ' num2str(median(gait_frequency))]);
text(0,0.53,['Median stance speed: ' num2str(median(stance_speed)) ' pixels/sec']);

% saveas(gcf,[output_filename '.png'])
% save([output_filename '.mat'],'results');
% if ~show_figure
%     close(gcf)
% end