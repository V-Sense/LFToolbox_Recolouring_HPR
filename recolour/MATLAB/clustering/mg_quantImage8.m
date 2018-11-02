function [X] = mg_quantImage8( RGB, n)
%apply matlab's rgb2ind to the RGB values given in 'RGB' to find n colors.
[im,quant] = rgb2ind(RGB, n);
X = unique(255*quant, 'rows');

end

