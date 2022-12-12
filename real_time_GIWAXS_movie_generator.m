%**************************************************************************
% GIWAXS evaluation
%**************************************************************************
clear all
%**************************************************************************
% read calibration text file
%**************************************************************************
path_to_calib_file='d:\laboratory\calibration.txt';
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
% root name of image
image_name_root='image';
% number of first image
no_start_image=0;
% number of last image
no_stop_image=400;
% time binning
time_binning=1;
% number of digits
no_digits=5;
% image name extension
image_name_extension='tif';
% path to folder
path_to_folder='d:\folder_with_images\';
% direction of evaluation, forward/backward
direction='forward';
%******************************************************************************
% remove bad pixels, 1 is True, 0 is False
remove_bad_pixels=1;
bad_pixels=[[119, 168]]; % [52 445]
size_bad_pixels=size(bad_pixels);
%**************************************************************************
myVideo = VideoWriter('D:\myVideoAlica.avi', 'Motion JPEG AVI');
open(myVideo);
%**************************************************************************
min_intensity=0; % [Counts]
max_intensity=5; % [Counts]
%**************************************************************************
% reshape parameters
no_of_points_x=400;
no_of_points_y=300;
% q_r minimal value
qr_min=-1.8; % [A-1]
% q_r maximal value
qr_max=1.8; % [A-1]
% q_z minimal value
qz_min=0.0; % [A-1]
% q_z maximal value
qz_max=2.5; % [A-1]
%**************************************************************************
% time constant between successive frames
time_step=100; % [ms]
time_const=time_binning*time_step;
if strcmp(direction, 'forward')
    time_actual=no_start_image*time_step;
    for index_0=no_start_image:time_binning:no_stop_image
        sum_matrix=zeros(no_of_points_y, no_of_points_x);
        for index_1=index_0:(index_0+time_binning-1)
           %**********************************************************************
           % calculate number of digits
           %**********************************************************************
           number=index_1;
           count_digits=0;
           if (number == 0)
               count_digits=1;
           else
               while (number > 0)
                   count_digits=count_digits+1;
                   number=floor(number/10);
               end
           end
           %**********************************************************************
           if (count_digits == no_digits)
               image_name_number=int2str(index_1);
           else
               image_name_number=int2str(index_1);
               for index_2=1:(no_digits-count_digits)
                   image_name_number=strcat('0', image_name_number);
               end
           end
           %**********************************************************************
           image_name=strcat(image_name_root, '_', image_name_number, '.', image_name_extension);
           %**********************************************************************
           path_to_image=strcat(path_to_folder, '\',image_name)
           %**************************************************************************
           % evaluation script
           %**************************************************************************
           obj=gixsdata(path_to_image);
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
           %**************************************************************************
           obj.ImFile=char(path_to_image);
           %**************************************************************************
           % remove bad pixels
           %**************************************************************************
           if (remove_bad_pixels == 1)
               image_matrix=obj.RawData;
               for index_2=1:size_bad_pixels(1)
                   image_matrix(bad_pixels(index_2,1), bad_pixels(index_2,2))=0.0;
               end
               obj.RawData=image_matrix(:,:);
           end
           %**************************************************************************
           obj.PlotAxisLabel = 2;
           %**************************************************************************
           % Define reshaping parameter for (qz vs qr)
           param_reshape.X = 6;  % qr
           param_reshape.Y = 3;  % qz
           param_reshape.XNOfPts = no_of_points_x; % number of points for X axis
           param_reshape.YNOfPts = no_of_points_y; % number of points for y axis
           param_reshape.XRange = [qr_min, qr_max];  % range for x
           param_reshape.YRange = [qz_min, qz_max];  % range for y
           %**************************************************************************
           % reshape
           dataflag=2;       % 2 for corrected data; 1 for masked rawdata
           [x,y,img_reshaped,countdata] = reshape_image(obj,param_reshape,dataflag);
           %**************************************************************************
           sum_matrix(:,:)=sum_matrix(:,:)+img_reshaped(:,:);
           %**************************************************************************
           delete(obj);
           %**************************************************************************
        end
        %**************************************************************************
        time_actual=time_actual+time_const;
        time_to_string=num2str(time_actual);
        time_string=strcat(time_to_string, ' ms');
        % plot reshaped image
        figure
        imagesc(x,y,sum_matrix, [min_intensity max_intensity]);
        axis ij;
        set(gca,'ydir','norm');
        xlabel('q_r (A^{-1})');
        ylabel('q_z (A^{-1})');
        title(time_string);
        c=colorbar;
        colormap parula;
        c.Label.String='Intensity (counts)';
        c.Label.FontSize=12;
        daspect(gca, [1 1 1]);
        %xticks(table2array(current_info(1,7)):qr_step:table2array(current_info(1,8)));
        %yticks(table2array(current_info(1,9)):qz_step:table2array(current_info(1,10)));
        %**************************************************************************
        frame=getframe(gcf);
        writeVideo(myVideo, frame);
        %**************************************************************************
    end
    % close video
    close(myVideo);
    %**************************************************************************
else
    if strcmp(direction, 'backward')
        time_actual=no_stop_image*time_step;
        for index_0=no_stop_image:((-1)*time_binning):no_start_image
            sum_matrix=zeros(no_of_points_y, no_of_points_x);
            for index_1=index_0:(-1):(index_0-time_binning+1)
                %**********************************************************************
                % calculate number of digits
                %**********************************************************************
                number=index_1;
                count_digits=0;
                if (number == 0)
                    count_digits=1;
                else
                    while (number > 0)
                        count_digits=count_digits+1;
                        number=floor(number/10);
                    end
                end
                %**********************************************************************
                if (count_digits == no_digits)
                    image_name_number=int2str(index_1);
                else
                    image_name_number=int2str(index_1);
                    for index_2=1:(no_digits-count_digits)
                        image_name_number=strcat('0', image_name_number);
                    end
                end
                %**********************************************************************
                image_name=strcat(image_name_root, '_', image_name_number, '.', image_name_extension);
                %**********************************************************************
                path_to_image=strcat(path_to_folder, '\',image_name)
                %**************************************************************************
                % evaluation script
                %**************************************************************************
                obj=gixsdata(path_to_image);
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
                %**************************************************************************
                obj.ImFile=char(path_to_image);
                %**************************************************************************
                % remove bad pixels
                %**************************************************************************
                if (remove_bad_pixels == 1)
                    image_matrix=obj.RawData;
                    for index_2=1:size_bad_pixels(1)
                        image_matrix(bad_pixels(index_2,1), bad_pixels(index_2,2))=0.0;
                    end
                    obj.RawData=image_matrix(:,:);
                end
                %**************************************************************************
                obj.PlotAxisLabel = 2;
                %**************************************************************************
                % Define reshaping parameter for (qz vs qr)
                param_reshape.X = 6;  % qr
                param_reshape.Y = 3;  % qz
                param_reshape.XNOfPts = no_of_points_x; % number of points for X axis
                param_reshape.YNOfPts = no_of_points_y; % number of points for y axis
                param_reshape.XRange = [qr_min, qr_max];  % range for x
                param_reshape.YRange = [qz_min, qz_max];  % range for y
                %**************************************************************************
                % reshape
                dataflag=2;       % 2 for corrected data; 1 for masked rawdata
                [x,y,img_reshaped,countdata] = reshape_image(obj,param_reshape,dataflag);
                %**************************************************************************
                sum_matrix(:,:)=sum_matrix(:,:)+img_reshaped(:,:);
                %**************************************************************************
                delete(obj);
                %**************************************************************************
            end
            %**************************************************************************
            time_actual=time_actual-time_const;
            time_to_string=num2str(time_actual);
            time_string=strcat(time_to_string, ' ms');
            % plot reshaped image
            figure
            imagesc(x,y,sum_matrix, [min_intensity max_intensity]);
            axis ij;
            set(gca,'ydir','norm');
            xlabel('q_r (A^{-1})');
            ylabel('q_z (A^{-1})');
            title(time_string);
            c=colorbar;
            colormap parula;
            c.Label.String='Intensity (counts)';
            c.Label.FontSize=12;
            daspect(gca, [1 1 1]);
            %xticks(table2array(current_info(1,7)):qr_step:table2array(current_info(1,8)));
            %yticks(table2array(current_info(1,9)):qz_step:table2array(current_info(1,10)));
            %**************************************************************************
            frame=getframe(gcf);
            writeVideo(myVideo, frame);
            %**************************************************************************
        end
    end
    %**************************************************************************
    % close video
    myVideo.FrameRate=2;
    close(myVideo);
    %******************************************************************************
end
%**************************************************************************
% End
%**************************************************************************