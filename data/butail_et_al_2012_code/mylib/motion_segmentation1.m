function retval = motion_segmentation1(method, img_info, varargin)
%function retval = motion_segmentation1(method, img_info, varargin)
%
%
% methods are: 'fwback'


switch method
    case {'fwback'}
        k=varargin{1};
        ball_radius=varargin{2};
        filter_noise_dev=varargin{3};
        max_min=varargin{4};
        file_list=varargin{5};
        file_loc=varargin{6};
        ud=varargin{7};
        camstr=varargin{8};
        retval=fwback_diff(img_info, k, ball_radius, filter_noise_dev, max_min, file_list, file_loc, ud, camstr);
    otherwise
        disp('Method does not exist');
end