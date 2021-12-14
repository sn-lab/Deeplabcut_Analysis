function data = rolling_average(data,dim,window_size,average_type)
%FUNCTION rolling_average(data,dim,window_size,average_type)
%
%Performs a rolling-window average on an N-dimensional data array (nans
%omitted)
%
%INPUTS:
%data: data array
%dim: dimensions to perform rolling average
%window_size: number of data points to average, centered on each datapoint
%               (window_size<2 returns unchaged data)
%average_type: 'mean', 'median', or 'mode'

dsize = size(data);
ndims = length(dsize);
dims = 1:ndims;

if nargin<2
    dim = find(dsize>1);
end
if nargin<3
    window_size = min([5 round(0.25*dsize(dim))]);
end
if nargin<4
    average_type = 'mean';
end

assert(ndims<11,'number of dimensions in data array cannot exceed 10')
assert(dim<=ndims & dim>0,'chosen dimension is outside data array')
assert(dsize(dim)>2,'data array is not large enough along the chosen dimension')
assert(rem(window_size,1)==0,'window_size must be an integer')
assert(window_size>=0,'window_size must be a positive integer')

if window_size>1 %window_sizes of 0 or 1 will be ignored
    %change the data array dimensions so that rolling average is performed on
    %the 1st dimension, and the 2nd dimensions is unused (singleton)
    order = dims;
    if ndims>1
        order(dim) = 1;
        order(1) = dim;
        order = [order(1) ndims+1 order(2:end)];
    end
    data = permute(data,order);

    %perform rolling average by shifting data forward/backward half of the
    %window size
    empty_rows_size = size(data);
    repmat_input = ones(size(1,ndims+1));
    repmat_input(2) = window_size;
    data = repmat(data,repmat_input);

    forward_shifts = floor(window_size/2);
    backward_shifts = ceil(window_size/2)-1;
    for i = 1:forward_shifts
        empty_rows_size(1) = i;
        data(:,i+1,:,:,:,:,:,:,:,:,:) = [nan(empty_rows_size); data(1:end-i,1,:,:,:,:,:,:,:,:,:)];
    end
    for i = 1:backward_shifts
        empty_rows_size(1) = i; 
        data(:,i+1+forward_shifts,:,:,:,:,:,:,:,:,:) = [data(1+i:end,1,:,:,:,:,:,:,:,:,:); nan(empty_rows_size)];
    end

    switch average_type
        case 'mean'
            data = mean(data,2,'omitnan');
        case 'median'
            data = median(data,2,'omitnan');
        case 'mode'
            data = mode(data,2,'omitnan');
        otherwise
            error('unrecognized "average type" (muse be "mean", "median", or "mode"')
    end

    %change data array dimensions back to original
    reverse_order = order;
    reverse_order(2) = []; %take out unused dimension
    reverse_order(1) = dim; %switch dim back to original position
    reverse_order(dim) = 1;
    data = permute(data,reverse_order);
end
    