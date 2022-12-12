clear all
close all
%**************************************************************************
% Author: Karol Vegso
% Affiliation: Institute of Physics, Slovak Aacdemy of Sciences
% Department: Dep. of multilayers and nanostructures
%**************************************************************************
% This program performs time binning of images, then pixel binning of
% images, rotation of images and 
% then takes q-cuts integrated as a function of azimuthal angle phi
% This algorithm is designed for WAXS measurement in transmission mode.
%**************************************************************************
% path to folder with images
path_to_input_folder='d:\Pilatus_images\';
% path to folder with result q-cuts matrix
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
% ROW DIRECTION - y direction (y image axis, vertical axis)
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
% COLUMN DIERCTION -  x direction (x image axis, horizontal axis)
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
% linecut parameters
%**************************************************************************
% number of q points to cut
no_of_q_points=300;
% q - minimum limit
q_min=0; % [A^(-1)]
% q - maximum limit
q_max=3.0; % [A^(-1)]
% phi - minimum limit
phi_min=-180; % [deg]
% phi - maximum limit
phi_max=180; % [deg]
% number of temporal linecuts
no_linecuts=length(no_image_start:time_binning:no_image_stop);
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
        % destroy/delete existing object
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
    %**********************************************************************
    obj.PlotAxisLabel = 2;
    %**********************************************************************
    % precised using constraints
    % define constraint
    constr = [1 2 phi_min phi_max; 1 1 q_min q_max];
    % perform linecut
    xflag = 1; % q linecut
    nofpts = no_of_q_points; % number of points in the linecut
    dataflag = 2; % 2 for corrected data; 1 for masked rawdata
    [x,y] = linecut(obj,xflag,constr,nofpts,dataflag);
    %**************************************************************************
    delete(obj);
    %**********************************************************************
    % increase value in counter by one
    counter=counter+1;
    %**********************************************************************
    % save x vector
    x=transpose(x(:));
    % insert new linecut to buffer matrix
    buffer_linecut(:,counter)=y(:,1);
    %**************************************************************************
end
% save x vector
save(strcat(path_to_output_folder, strcat(image_name_root, '_results_x.txt')), 'x', '-ascii');
% save q cuts matrix
save(strcat(path_to_output_folder, strcat(image_name_root, '_results.txt')), 'buffer_linecut', '-ascii');
%**************************************************************************
% end of program
%**************************************************************************