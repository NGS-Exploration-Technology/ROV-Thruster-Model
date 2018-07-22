function [filtered] = median_filter(samples, median_size)
%MEDIAN_FILTER samples to remove outliers
%function [filtered] = median_filter(samples, median_size)
%inputs
%   samples = vector of samples to be filtered
%   size = size of median filter
%outputs
%   filtered = output vector

index = 1;

while(index<(length(samples)-median_size))
    median_vector = samples(index:index+median_size-1);
    filtered(index) = median(median_vector);
    index = index+1;
end

%pad end of filtered with median
while(index<(length(samples)+1))
    filtered(index)=median(median_vector);
    index = index+1;
end

filtered = filtered';
