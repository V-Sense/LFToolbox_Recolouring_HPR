function [X] = mg_quantImage16( RGB, n)
%apply matlab's rgb2ind to the RGB values given in 'RGB' to find n colors.
[im,quant] = rgb2ind(RGB, n);
X = unique(65535*quant, 'rows');

end

