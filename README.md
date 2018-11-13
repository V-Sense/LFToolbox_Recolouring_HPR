# Light Field toolbox - Recolouring step

This code is mostly meant to correct colour inconsistencies in Light Field datasets taken using plenoptic cameras [1].
This will allow you to propagate the colours from the centre of a Light Field to all the images, which will ensure you get consistent sub-aperture images accross the whole set.

If you use or adapt any part of this code in your work, please remember to cite the appropriate paper [1].
* Authors   : Pierre Matysiak (matysiap@scss.tcd.ie), Mairéad Grogan (mgrogan@tcd.ie)
* Copyright : (c) 2018 [V-Sense](https://v-sense.scss.tcd.ie)
* License   : GPL v3+, see GPLv3.txt

For additional details and visual results, see [our webpage](https://v-sense.scss.tcd.ie/research/light-fields/a-pipeline-for-lenslet-light-field-quality-enhancement/).

The additional post-processing steps presented in [1] (RAW data demultiplexing and denoising) are not present in this repository.
You can find them respectively [here](https://github.com/V-Sense/LFToolbox-CLIM_VSENSE) and [here](https://github.com/V-Sense/LFBM5D).

The proper order in which to use these tools is : RAW data demultiplexing --> Recolouring --> Denoising.

# Usage

For simplicity's sake, we made it so that you only have to use the one file, everything else is done automatically :
- put any number of datasets in the folder called "data"
- at the beginning of the file "recolour_centre_neighbour.m", add/replace the paths to your datasets, assuming that the root folder is ".../data/"
- launch the code by calling "recolour_centre_neighbour.m" ; no arguments are necessary.

Please note that the recolouring code uses OpenMP, please make sure you have the proper dependancies installed and that your compiler supports it.

## References

[1] P. Matysiak, M. Grogan, M. Le Pendu, M. Alain and A. Smolic, ["A Pipeline for Lenslet Light Field Quality Enhancement"](https://v-sense.scss.tcd.ie/research/light-fields/a-pipeline-for-lenslet-light-field-quality-enhancement/), International Conference on Image Processing (ICIP) 2018.

[2] M. Grogan and R. Dahyot, [“Robust registration of gaussian mixtures for colour transfer”](https://arxiv.org/abs/1705.06091), ArXiv e-prints (May 2017). arXiv:cs.CV/1705.06091.

[3] Y. Hu, R. Song, and Y. Li, [“Efficient coarse-to-fine patchmatch for large displacement optical flow”](https://www.cv-foundation.org/openaccess/content_cvpr_2016/papers/Hu_Efficient_Coarse-To-Fine_PatchMatch_CVPR_2016_paper.pdf), in Proc. CVPR, 2016.
