clear all
close all
%**************************************************************************
% Author of program: Karol Vegso
% Affiliation: Institute of Physics, Slovak Academy of Sciences
% Department: Dep. of multilayers and nanostructures
%**************************************************************************
% This program performs time binning of images, pixel binning of
% images, rotation of images and 
% then reshaping of images to final GISAXS patterns
%**************************************************************************
% path to folder with images
path_to_input_folder='d:\Pilatus_images\';
% path to folder with result images
path_to_output_folder='d:\output_folder\';
% path to calibration file
path_to_calib_file='d:\calibration\calibration_text_file.txt';
% specify image name - root
image_name_root='image';
% image name extension
image_name_ext='.tif';
% image number width
image_no_width=5;
% image number - down limit, start point
no_image_start=0;
% image number - upper limit, stop point
no_image_stop=1200;
% specify time binning
time_binning=1;
%**************************************************************************
% rotation images
%**************************************************************************
% rotate images in degrees by angle
rot_angle=0.0; % [deg]
%**************************************************************************
% pixel binning - parameters
%**************************************************************************
% ROW DIRECTION - y direction (y axis of image, vertical axis of image)
%**************************************************************************
% specify pixel binning in row direction
pixel_binning_row=1;
% start pixel in row direction
start_pixel_row=1; % [pixels]
% stop pixel in row direction
stop_pixel_row=407; % [pixels]
% number of rows in new image
no_rows_new_image=length(start_pixel_row:pixel_binning_row:stop_pixel_row);
%**************************************************************************
% COLUMN DIERCTION - x direction (x axis of image, horizontal axis of image)
%**************************************************************************
% specify pixel binning in column direction
pixel_binning_col=1;
% start pixel in column direction
start_pixel_col=1; % [pixels]
% stop pixel in col direction
stop_pixel_col=487; % [pixels]
% number of columns in new image
no_cols_new_image=length(start_pixel_col:pixel_binning_col:stop_pixel_col);
%**************************************************************************
% read calibration text file
%**************************************************************************
% read calibration text file
calib_info=readtable(path_to_calib_file, 'Delimiter', ' ', 'Format', '%s%u%u%f%f%f%f%u%u%u%f%f%f%f%f%f');
% type of detector
type_of_detector=char(table2array(calib_info(1,1)));
% image dimension
image_dimension=[table2array(calib_info(1,2)), table2array(calib_info(1,3))];
% pixel size
pixel_size=[table2array(calib_info(1,4)), table2array(calib_info(1,5))];
% X-ray energy
E=table2array(calib_info(1,6)); % [keV]
% X-ray wavelength
lambda=1.24/E; % [nm]
% Sample to detector distance
SDD=table2array(calib_info(1,7)); % [in pixels]
% geometry, 1 for transmission, 2 for reflection
geometry=table2array(calib_info(1,8));
% PhiMode, 1 for [-180°, 180°], 2 for [0°, 360°], 3 for [-270°, 90°], 4 for [-90°, 270°]
phi_mode=table2array(calib_info(1,9));
% PolarizationMode, 1, 2, 3, 4 for none, horiozntal, vertical, unpolarized
polarization_mode=table2array(calib_info(1,10));
% horizontal polarization fraction, 1 for fully horizontally polarized beam,
% 0 for fully vertically polarized beam
horizontal_polarization_fraction=table2array(calib_info(1,11));
% incident angle
incident_angle=table2array(calib_info(1,12)); % [deg]
% Beam0
beam0=[table2array(calib_info(1,13)), table2array(calib_info(1,14))];
% Specular
specular=[table2array(calib_info(1,15)), table2array(calib_info(1,16))];
%**************************************************************************
% reshape parameters
% parameters of reshaping image to final GISAXS pattern
%**************************************************************************
% reshape in qy direction to number of points in x direction
no_points_x=round(250/pixel_binning_col);
% reshape in qz direction to number of points in y direction
no_points_y=round(120/pixel_binning_row);
% reshape along qy axis from qy minimal value
qy_min=-0.2; % [A^(-1)]
% reshape along qy axis to qy maximal value
qy_max=0.2; % [A^(-1)]
% reshape along qz axis from qz minimal value
qz_min=0.0; % [A^(-1)]
% reshape along qz axis to qz maximal value
qz_max=0.2; % [A^(-1)]
%**************************************************************************
% main loop
%**************************************************************************
% initialize counter
counter=0.0;
for index_0=no_image_start:time_binning:no_image_stop
    %**********************************************************************
    % perform time binning
    %**********************************************************************
    % initialzie integrate_image corresponding to time binning
    integrate_image=zeros(image_dimension(1,2), image_dimension(1,1));
    for index_1=index_0:(index_0+time_binning-1)
        %******************************************************************
        % calculate image number
        no_digits=0;
        current_image_no=index_1;
        if (index_1 == 0)
            no_digits=1;
        else
            while (current_image_no > 0)
                no_digits=no_digits+1;
                current_image_no=floor(current_image_no/10);
            end
        end
        %******************************************************************
        image_no_str=num2str(index_1);
        for index_2=1:(image_no_width-no_digits)
            image_no_str=strcat('0', image_no_str);
        end
        full_image_name=strcat(image_name_root, '_', image_no_str, image_name_ext);
        full_path_to_image=strcat(path_to_input_folder, full_image_name)
        %******************************************************************
        % create existing gixsdata obejct and calibrate it
        obj=gixsdata(full_path_to_image);
        obj.Camera=type_of_detector;
        obj.PixelSize=pixel_size;
        %obj.ImDim=image_dimension;
        obj.SDD=SDD;
        obj.XEnergy=E;
        obj.Geometry=geometry;
        obj.PhiMode=phi_mode;
        obj.PolarizationMode=polarization_mode;
        obj.HorizontalPolarizationFraction=horizontal_polarization_fraction;
        obj.Beam0=beam0;
        obj.Specular=specular;
        obj.IncidentAngle=incident_angle;
        %******************************************************************
        obj.ImFile=char(full_path_to_image);
        %**************************************************************************
        if (rot_angle == 0)
        else
            obj.RawData=imrotate(obj.RawData,rot_angle,'bilinear','crop');
        end
        %**************************************************************************
        integrate_image(:,:)=integrate_image(:,:)+obj.RawData(:,:);
        %**************************************************************************
        % destroy/delete existing gixsdata object
        delete(obj);
        %**************************************************************************
    end
    %**********************************************************************
    % perform pixel binning
    %**********************************************************************
    % initialize of pixel binning matrix/array
    new_image=zeros(no_rows_new_image, no_cols_new_image);
    index_5=0;
    for index_3=start_pixel_row:pixel_binning_row:stop_pixel_row
        index_5=index_5+1;
        index_6=0;
        for index_4=start_pixel_col:pixel_binning_col:stop_pixel_col
            index_6=index_6+1;
            new_image(index_5,index_6)=...
                sum(sum((integrate_image(index_3:(index_3+pixel_binning_row-1),index_4:(index_4+pixel_binning_col-1)))));
        end
    end
    %**************************************************************************
    % calculate horizontal (x) and vertical dimensions (y) of new image
    size_new_image=size(new_image(:,:));
    size_new_image_x=size_new_image(1,2);
    size_new_image_y=size_new_image(1,1);
    image_dimension_new=[size_new_image_x, size_new_image_y];
    %**************************************************************************
    % create new path to new image
    %**************************************************************************
    % calculate image number
    no_digits=0;
    current_image_no=counter;
    if (counter == 0)
        no_digits=1;
    else
        while (current_image_no > 0)
            no_digits=no_digits+1;
            current_image_no=floor(current_image_no/10);
        end
    end
    %******************************************************************
    image_no_str=num2str(counter);
    for index_2=1:(image_no_width-no_digits)
        image_no_str=strcat('0', image_no_str);
    end
    full_image_name_create=strcat(image_name_root, '_', image_no_str, image_name_ext);
    full_path_to_image_create=strcat(path_to_output_folder, full_image_name_create)
    %**********************************************************************
    % create new gixsdata obejct and calibrate it
    obj=gixsdata(full_path_to_image_create);
    obj.Camera='Other';
    obj.PixelSize=[pixel_size(1,1)*pixel_binning_col, pixel_size(1,2)*pixel_binning_row];
    %obj.ImDim=image_dimension_new;
    obj.SDD=SDD;
    obj.XEnergy=E;
    obj.Geometry=geometry;
    obj.PhiMode=phi_mode;
    obj.PolarizationMode=polarization_mode;
    obj.HorizontalPolarizationFraction=horizontal_polarization_fraction;
    obj.Beam0=[beam0(1,1)/pixel_binning_col, beam0(1,2)/pixel_binning_row];
    obj.Specular=[specular(1,1)/pixel_binning_col, specular(1,2)/pixel_binning_row];
    obj.IncidentAngle=incident_angle;
    obj.RawData=new_image(:,:);
    %**************************************************************************
    % 2 for q representation
    obj.PlotAxisLabel = 2;
    %**************************************************************************
    % Define reshaping parameter for (qy vs qz)
    param_reshape.X = 5;  % qy
    param_reshape.Y = 3;  % qz
    param_reshape.XNOfPts = no_points_x; % number of points for X axis
    param_reshape.YNOfPts = no_points_y; % number of points for y axis
    param_reshape.XRange = [qy_min, qy_max];  % range for x
    param_reshape.YRange = [qz_min, qz_max];  % range for y
    %**************************************************************************
    % reshape
    dataflag=2;       % 2 for corrected data; 1 for masked rawdata
    [x,y,img_reshaped,countdata] = reshape_image(obj,param_reshape,dataflag);
    %**************************************************************************
    full_image_name_save=strcat(image_name_root, '_', image_no_str, '.txt');
    full_path_to_image_save=strcat(path_to_output_folder, full_image_name_save)
    save(full_path_to_image_save, 'img_reshaped', '-ascii');
    %**************************************************************************
    % destroy / delete gixsdata object
    delete(obj);
    %**********************************************************************
    % increase value in counter by one
    counter=counter+1;
    %**************************************************************************
end
%**************************************************************************
% end of program
%**************************************************************************