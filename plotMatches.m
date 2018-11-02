im01 = imread('data/deco+HPR/EPFL/Color_Chart_1__Decoded/SAI_07_07.png');
im02 = imread('data/deco+HPR/EPFL/Color_Chart_1__Decoded/SAI_00_05.png');

matches = dlmread('sparse_matches/EPFL/Color_Chart_1__Decoded/07070005.txt');

match_size = size(matches, 1);

matches_pixels01(1:match_size, 2) = uint16(0);
matches_pixels02(1:match_size, 2) = uint16(0);

for k = 1:match_size
    matches_pixels01(k,1) = matches(k,1);
    matches_pixels01(k,2) = matches(k,2);

    matches_pixels02(k,1) = matches(k,3);
    matches_pixels02(k,2) = matches(k,4);
end

figure;
ax=axes;

title(ax,'Candidate Point Matches');

showMatchedFeatures(im01, im02, matches_pixels01(1300:1330,:), matches_pixels02(1300:1330,:), 'montage', 'Parent', ax)