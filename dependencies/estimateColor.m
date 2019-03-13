function estimate = estimateColor(color)
% a weird function that attempts to give a vague human description to a 3
% value color code...
load('colorcodes.mat','cols','names')

distances = Inf(1,length(cols));
color = color * 255;
for n = 1:length(cols)
    distances(n) = sum(abs(cols(n,:) - color));
end
[~,ind] = min(distances);
estimate = names{ind};