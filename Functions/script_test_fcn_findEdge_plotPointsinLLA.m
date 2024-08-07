%script_fcn_findEdge_plotPointsinLLA
%script for fcn_findEdge_plotPointsinLLA
%
% Written by:
% 7_11_2024 -- Aleksandr Goncharov
% Revision: 
% 7_17_2024 -- Aleksandr Goncharov
% Added more cases 
% 7_29_2024 -- Aleksandr Goncharov
% Updated cases and made sure no legend case worked
% 2024_08_07 - Jiabao Zhao
% -- pull this code from geometry class and rename it

%% Basic Example
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%   ____            _        ______                           _      
%  |  _ \          (_)      |  ____|                         | |     
%  | |_) | __ _ ___ _  ___  | |__  __  ____ _ _ __ ___  _ __ | | ___ 
%  |  _ < / _` / __| |/ __| |  __| \ \/ / _` | '_ ` _ \| '_ \| |/ _ \
%  | |_) | (_| \__ \ | (__  | |____ >  < (_| | | | | | | |_) | |  __/
%  |____/ \__,_|___/_|\___| |______/_/\_\__,_|_| |_| |_| .__/|_|\___|
%                                                      | |           
%                                                      |_|          
% See: https://patorjk.com/software/taag/#p=display&f=Big&t=Basic%20Example
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%§

%% Example 1 - All inputs using sample ENU_data;

ENU_data =[418.5194  230.3028 -344.2069;
    419.1155  229.8859 -344.2069;
    419.8043  229.5207 -344.2069;
    420.3828  229.0694 -344.2070;
    421.0570  228.7612 -344.2070;
    421.5912  228.2227 -344.2070;
    422.2313  227.8078 -344.2070;
    422.9040  227.4922 -344.2071;
    423.4936  227.0386 -344.2071;
    424.0961  226.5475 -344.2071;
    424.7295  226.1757 -344.2071];

ENU_data2=ENU_data+1;
ENU_data3=ENU_data-1;

marker_size=30;
marker_type="*";
RGB_triplet=[0.2 0.8 0.2];
RGB_triplet2=[.8 .1 .1];
RGB_triplet3=[0 .2 .2];
fig_num=1134;
legend_option=1;
reference_longitude=[];
reference_latitude=[];
reference_altitude=[];

legend_position=[.75 .8 0 0];

legend_name='Legend Test';
legend_name2='Legend Test 2';
legend_name3='Legend Test 3';

%fcn_geometry_plotGridCentersBoundaryPoints(ENU_data,marker_size,RGB_triplet,(legend_option),(legend_name),(fig_num))
fcn_findEdge_plotPointsinLLA(ENU_data,marker_size,RGB_triplet,marker_type,legend_option,legend_name,legend_position,reference_longitude,reference_latitude,reference_altitude,fig_num);
fcn_findEdge_plotPointsinLLA(ENU_data2,marker_size,RGB_triplet2,marker_type,legend_option,legend_name2,legend_position,reference_longitude,reference_latitude,reference_altitude,fig_num);
fcn_findEdge_plotPointsinLLA(ENU_data3,marker_size,RGB_triplet3,marker_type,legend_option,legend_name3,legend_position,reference_longitude,reference_latitude,reference_altitude,fig_num);

%% Example 2 - All inputs using sample ENU_data, different colors, different figure number

ENU_data =[418.5194  230.3028 -344.2069;
    419.1155  229.8859 -344.2069;
    419.8043  229.5207 -344.2069;
    420.3828  229.0694 -344.2070;
    421.0570  228.7612 -344.2070;
    421.5912  228.2227 -344.2070;
    422.2313  227.8078 -344.2070;
    422.9040  227.4922 -344.2071;
    423.4936  227.0386 -344.2071;
    424.0961  226.5475 -344.2071;
    424.7295  226.1757 -344.2071];

ENU_data2=ENU_data+1;
ENU_data3=ENU_data-1;
marker_size=30;
marker_type=".";
RGB_triplet=[0.2 0.8 0.2];
RGB_triplet2=[.8 .1 .1];
RGB_triplet3=[0 .2 .2];
fig_num=1134;
legend_option=0;

legend_position=[.75 .8 0 0];

legend_name='Legend Test';
legend_name2='Legend Test 2';
legend_name3='Legend Test 3';

%fcn_geometry_plotGridCentersBoundaryPoints(ENU_data,marker_size,RGB_triplet,(legend_option),(legend_name),(fig_num))
fcn_findEdge_plotPointsinLLA(ENU_data,marker_size,RGB_triplet,marker_type,legend_option,legend_name,legend_position,reference_longitude,reference_latitude,reference_altitude,fig_num);
fcn_findEdge_plotPointsinLLA(ENU_data2,marker_size,RGB_triplet2,marker_type,legend_option,legend_name2,legend_position,reference_longitude,reference_latitude,reference_altitude,fig_num);
fcn_findEdge_plotPointsinLLA(ENU_data3,marker_size,RGB_triplet3,marker_type,legend_option,legend_name3,legend_position,reference_longitude,reference_latitude,reference_altitude,fig_num);

%% Example 3 - All inputs using sample ENU_data, no legend

ENU_data =[418.5194  230.3028 -344.2069;
    419.1155  229.8859 -344.2069;
    419.8043  229.5207 -344.2069;
    420.3828  229.0694 -344.2070;
    421.0570  228.7612 -344.2070;
    421.5912  228.2227 -344.2070;
    422.2313  227.8078 -344.2070;
    422.9040  227.4922 -344.2071;
    423.4936  227.0386 -344.2071;
    424.0961  226.5475 -344.2071;
    424.7295  226.1757 -344.2071];

ENU_data2=ENU_data+1;
ENU_data3=ENU_data-1;

marker_size=30;
marker_type="o";
RGB_triplet=[0.2 0.8 0.2];
RGB_triplet2=[.8 .1 .1];
RGB_triplet3=[0 .2 .2];
fig_num=1134;
legend_option=[];

legend_position=[.75 .8 0 0];

legend_name='Legend Test';
legend_name2='Legend Test 2';
legend_name3='Legend Test 3';

%fcn_geometry_plotGridCentersBoundaryPoints(ENU_data,marker_size,RGB_triplet,(legend_option),(legend_name),(fig_num))
fcn_findEdge_plotPointsinLLA(ENU_data,marker_size,RGB_triplet,marker_type,legend_option,legend_name,legend_position,reference_longitude,reference_latitude,reference_altitude,fig_num);
fcn_findEdge_plotPointsinLLA(ENU_data2,marker_size,RGB_triplet2,marker_type,legend_option,legend_name2,legend_position,reference_longitude,reference_latitude,reference_altitude,fig_num);
fcn_findEdge_plotPointsinLLA(ENU_data3,marker_size,RGB_triplet3,marker_type,legend_option,legend_name3,legend_position,reference_longitude,reference_latitude,reference_altitude,fig_num);

%% Example 4 - All inputs using sample ENU_data, legend moved
ENU_data =[418.5194  230.3028 -344.2069;
    419.1155  229.8859 -344.2069;
    419.8043  229.5207 -344.2069;
    420.3828  229.0694 -344.2070;
    421.0570  228.7612 -344.2070;
    421.5912  228.2227 -344.2070;
    422.2313  227.8078 -344.2070;
    422.9040  227.4922 -344.2071;
    423.4936  227.0386 -344.2071;
    424.0961  226.5475 -344.2071;
    424.7295  226.1757 -344.2071];

ENU_data2=ENU_data+1;
ENU_data3=ENU_data-1;
marker_type=".";
marker_size=30;
RGB_triplet=[0.2 0.8 0.2];
RGB_triplet2=[.8 .1 .1];
RGB_triplet3=[0 .2 .2];
fig_num=1134;
legend_option=1;

legend_position=[.2 .8 0 0];

legend_name='Legend Test';
legend_name2='Legend Test 2';
legend_name3='Legend Test 3';

%fcn_geometry_plotGridCentersBoundaryPoints(ENU_data,marker_size,RGB_triplet,(legend_option),(legend_name),(fig_num))
fcn_findEdge_plotPointsinLLA(ENU_data,marker_size,RGB_triplet,marker_type,legend_option,legend_name,legend_position,reference_longitude,reference_latitude,reference_altitude,fig_num);
fcn_findEdge_plotPointsinLLA(ENU_data2,marker_size,RGB_triplet2,marker_type,legend_option,legend_name2,legend_position,reference_longitude,reference_latitude,reference_altitude,fig_num);
fcn_findEdge_plotPointsinLLA(ENU_data3,marker_size,RGB_triplet3,marker_type,legend_option,legend_name3,legend_position,reference_longitude,reference_latitude,reference_altitude,fig_num);

%% Example 5 - Max Speed
ENU_data =[418.5194  230.3028 -344.2069;
    419.1155  229.8859 -344.2069;
    419.8043  229.5207 -344.2069;
    420.3828  229.0694 -344.2070;
    421.0570  228.7612 -344.2070;
    421.5912  228.2227 -344.2070;
    422.2313  227.8078 -344.2070;
    422.9040  227.4922 -344.2071;
    423.4936  227.0386 -344.2071;
    424.0961  226.5475 -344.2071;
    424.7295  226.1757 -344.2071];

ENU_data2=ENU_data+1;
ENU_data3=ENU_data-1;
marker_type=".";
marker_size=30;
RGB_triplet=[0.2 0.8 0.2];
RGB_triplet2=[.8 .1 .1];
RGB_triplet3=[0 .2 .2];
fig_num=-1;
legend_option=1;

legend_position=[.2 .8 0 0];

legend_name='Legend Test';
legend_name2='Legend Test 2';
legend_name3='Legend Test 3';

%fcn_geometry_plotGridCentersBoundaryPoints(ENU_data,marker_size,RGB_triplet,(legend_option),(legend_name),(fig_num))
fcn_findEdge_plotPointsinLLA(ENU_data,marker_size,RGB_triplet,marker_type,legend_option,legend_name,legend_position,reference_longitude,reference_latitude,reference_altitude,fig_num);
fcn_findEdge_plotPointsinLLA(ENU_data2,marker_size,RGB_triplet2,marker_type,legend_option,legend_name2,legend_position,reference_longitude,reference_latitude,reference_altitude,fig_num);
fcn_findEdge_plotPointsinLLA(ENU_data3,marker_size,RGB_triplet3,marker_type,legend_option,legend_name3,legend_position,reference_longitude,reference_latitude,reference_altitude,fig_num);

