cpalette = imread('recolour/results/PROP+CENTER/Magnets_1__Decoded/0707.png');
palette = imread('recolour/results/PROP+CENTER/Magnets_1__Decoded/1307.png');
target = imread('data/Magnets_1__Decoded/SAI_14_07.png');

cmatches = dlmread('sparse_matches/mixed/Magnets_1__Decoded/07071407.txt');
matches = dlmread('sparse_matches/mixed/Magnets_1__Decoded/13071407.txt');

cmatch_size = size(cmatches, 1);
match_size = size(matches, 1);

col_palette(1:match_size+cmatch_size, 3) = double(0);
col_target(1:match_size+cmatch_size, 3) = double(0);

for k = 1:match_size
    col_palette(k,:) = palette(matches(k,2) + 1, matches(k,1) + 1, :);
    col_target(k,:) = target(matches(k,4) + 1, matches(k,3) + 1, :);
end

for k = match_size+1:match_size+cmatch_size
    col_palette(k,:) = cpalette(cmatches(k-match_size,2) + 1, cmatches(k-match_size) + 1, :);
    col_target(k,:) = target(cmatches(k-match_size,4) + 1, cmatches(k-match_size) + 1, :);
end

result = fullfile('recolour/results/', '13071407.png');

ctfunction_corr(palette, target, 8, 2, col_palette, col_target, result, 'CIELab');

%plotting colours in CT
%scatter3(config.model(:,1), config.model(:,2), config.model(:,3), 10, [config.model(:,1), config.model(:,2), config.model(:,3)]./255, 'filled')
%scatter3(config.scene(:,1), config.scene(:,2), config.scene(:,3), 10, [config.scene(:,1), config.scene(:,2), config.scene(:,3)]./255, 'filled')
%