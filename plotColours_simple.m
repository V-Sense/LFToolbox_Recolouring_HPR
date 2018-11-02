palette = imread('data/Magnets_1__Decoded/SAI_00_10.png');
target = imread('data/Magnets_1__Decoded/SAI_00_11.png');

matches = dlmread('sparse_matches/(prop_ccol+rows)/Magnets_1__Decoded/00100011.txt');

match_size = size(matches, 1);

col_palette(1:match_size+cmatch_size, 3) = double(0);
col_target(1:match_size+cmatch_size, 3) = double(0);

for k = 1:match_size
    col_palette(k,:) = palette(matches(k,2) + 1, matches(k,1) + 1, :);
    col_target(k,:) = target(matches(k,4) + 1, matches(k,3) + 1, :);
end

result = fullfile('recolour/results/', '00100011.png');

ctfunction_corr(palette, target, 8, 2, col_palette, col_target, result, 'RGB');

%plotting colours in CT
%scatter3(config.model(:,1), config.model(:,2), config.model(:,3), 10, [config.model(:,1), config.model(:,2), config.model(:,3)]./255, 'filled')
%scatter3(config.scene(:,1), config.scene(:,2), config.scene(:,3), 10, [config.scene(:,1), config.scene(:,2), config.scene(:,3)]./255, 'filled')
%