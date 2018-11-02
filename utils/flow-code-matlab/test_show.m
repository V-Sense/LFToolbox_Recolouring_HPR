clear;

addpath('.');

flow_path = 'dragon1.flo';

flow_img_path = 'dragon1.png';

flow = readFlowFile(flow_path);

imwrite(flowToColor(flow),flow_img_path);
