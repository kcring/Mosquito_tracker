function checkstatus

run config
Xi=strXi;
if ~strcmp(floc(end), '/')
    floc=[floc, '/'];
end

fprintf('\n');
% menu
menu(1).str='Configure (setupdir.m)';
% menu(2).str='Get bg params (get_fwback_diff_params.m)';
menu(2).str='Auto Track (trackemall.m)';
menu(3).str='Verify/Join tracks (trackone.m)';
% menu(4).str='Filtering tracks (fwd_smoothing_kf.m)';
% menu(5).str='Convert to world frame (cam2world.m)';


% status
nm=size(menu,2);
for ii=1:nm
	menu(ii).status=' ';
end


% info
[expname image_id]=scan_expfile(floc);
if isempty(expname)
    return
end
for ii=1:nm
	menu(ii).info=[];
    menu(ii).info.floc=floc;
    menu(ii).info.frmloc=frmloc;
    menu(ii).info.imageid=image_id;  
    menu(ii).info.expname=expname;
end

fprintf('\t [I] Information, [?] Warning, [!] Error\n\n');
nitems=size(menu,2);
for ii=1:nitems
	menu(ii)=checkmenu(menu(ii), ii, Xi);
end
dispmenu(menu)



function menuitem=checkmenu(menuitem, ii, Xi)

switch ii
	case 1
        if isempty(menuitem.info.frmloc)
            menuitem.info.frmloc=[menuitem.info.floc, 'frames/'];
        end	
        fprintf('[I] floc=%s\n', menuitem.info.floc);
        fprintf('[I] frames loc=%s\n', menuitem.info.frmloc);
        if isempty(menuitem.info.imageid)
            menuitem.info.imageid='???';
            stat=0;
        end
        fprintf('[I] imageid=%s\n', menuitem.info.imageid);   
        
        
		% check folders and set status as 'x' if done
		stat=1;
		if isempty(menuitem.info)
			stat=0;
            fprintf('[!] Could not find info\n');
		end

        if ~exist(menuitem.info.floc, 'dir')
            stat=0;
            fprintf('[!] Could not find %s\n', menuitem.info.floc);
        end

    	if ~exist([menuitem.info.floc, '/calib'], 'dir')
            stat=0;
            fprintf('[!] Could not find %s/calib\n', menuitem.info.floc);                
        end
	
	    metafile=dir([menuitem.info.floc, '/calib/metadata*.txt']);
        if ~size(metafile, 1)
            stat=0;
            fprintf('[!] Could not find %s/calib/metadata*.txt\n', menuitem.info.floc);                
        else
            fprintf('[I] Found metadata (%s)\n', metafile(1).name);                
            md=csvread([menuitem.info.floc, '/calib/', metafile(1).name]);
            menuitem.info.azimuth=md(1);
            menuitme.info.incl=md(2);
        end


        calibfun=dir([menuitem.info.floc 'calib/' '*calib*.m']);
        if ~size(calibfun,1)
            stat=0;
            fprintf('\n[?] Could not find calibration for the cameras ...\n');
        end

        climate_data=dir([menuitem.info.floc 'calib/' 'Kestrel*.*']);
        if ~size(climate_data,1)
            stat=0;
            fprintf('[?] Did not find climate data (Kestrel*.csv file) in calib folder ...\n');
        else
            fprintf('[I] Found climate data (%s) ...\n', climate_data(1).name);
        end
        
        
        % backward compatibility (going to add expfile myself...)
%   		if stat
% 			if ~exist([menuitem.info.floc, '/calib/expfile.txt'], 'file')
% 				stat=0;
%                 fprintf('[!] Could not find %s/calib/expfile.txt\n', menuitem.info.floc);                
%                 fprintf('Adding a file myself... *** Run checkstatus again ***..\n');
%                 expfile=[menuitem.info.floc, '/calib/', 'expfile.txt'];
%                 fid=fopen(expfile, 'w');
%                 fprintf(fid, '%s', menuitem.info.imageid(1:end-1));
%                 fclose(fid);
% 			end
%         end

        
        if ~exist([menuitem.info.frmloc], 'dir')
            stat=0;
            fprintf('[?] Could not find %s\n', menuitem.info.frmloc);                
        end

        if ~exist([menuitem.info.floc, '/output'], 'dir')
            stat=0;
            fprintf('[!] Could not find %s/output\n', menuitem.info.floc);                
        end

        if ~exist([menuitem.info.floc, '/output/data'], 'dir')
            stat=0;
            fprintf('[!] Could not find %s/output/data\n', menuitem.info.floc);                
        end


        if ~exist([menuitem.info.floc, '/output/movies'], 'dir')
            stat=0;
            fprintf('[!] Could not find %s/output/movies\n', menuitem.info.floc);                
        end



        % display information


        fprintf('[I] Found calibration (%s) ...\n', calibfun(1).name);
        nimgs=0;
        if exist(menuitem.info.frmloc, 'dir')
            flist=dir([menuitem.info.frmloc, menuitem.info.imageid, '*']);
            nimgs=size(flist,1);
        end

        if nimgs ==0
            fprintf('[?] %d images found...\n', nimgs);
        else
            fprintf('[I] %d images found...\n', nimgs);
        end
            

        if stat
            menuitem.status='x';
        end
% 	case 2
%         
%         stat=0;
% 
%         bglist=dir([menuitem.info.floc, '/output/data/', 'bg_params_', ...
%             menuitem.info.imageid, '*.mat']);
% 
%         for ii =1:size(bglist,1)
%             fprintf('[I] Background params (%s) found...\n', bglist(ii).name);
%         end
% 
%         if size(bglist,1) == 2
%             stat=1;
%         else
%             fprintf('[!] No background params found... check image_id in config.m also.. \n');
%         end
%    
%         
%         if stat
%             menuitem.status='x';
%         end
		% check fg files for each camera
        

	case 2
        stat=1;

        
        %*********** Backward compatibility -- backing up files
%         if stat
%             old_datfile=dir([menuitem.info.floc, '/output/data/', expname, '__datfile1*.mat']);
% 
%             if size(old_datfile,1)
%                 if ~exist([menuitem.info.floc, '/output/data/bkp'], 'dir')
%                     mkdir([menuitem.info.floc, '/output/data/bkp']);
%                 end
%                 for ii=1:size(old_datfile,1)
%                     if strcmp(old_datfile(ii).name, [expname, '__datfile1.mat'])
%                         copyfile([menuitem.info.floc, '/output/data/', old_datfile(ii).name], ...
%                             [menuitem.info.floc, '/output/data/', 'data_mq_', expname, '.mat']);
%                     end
%                         movefile([menuitem.info.floc, '/output/data/', old_datfile(ii).name], ...
%                             [menuitem.info.floc, '/output/data/bkp/', old_datfile(ii).name]);
%                 end
%             end
%         end
        %****************Backwardn compatibility.......
        
        if stat
            datfile=dir([menuitem.info.floc, '/output/data/', 'data_mq_auto_', ...
                menuitem.info.expname, '.mat']);
            if size(datfile,1)
%                 fprintf('[I] Data (%s) found ...\n', datfile(ii).name);
                stat=1;
            else
                stat=0;
            end
            
%             if size(datfile,1)
%                load([menuitem.info.floc, '/output/data/', datfile(1).name], 'Xh');
% 
%                 m_flags=Xh(Xi.fi(1):nx:end,:);
%                 m_ind=find(m_flags==1);
%                 [m_mqid k]=ind2sub(size(m_flags), m_ind);
% 
%                 if isempty(m_mqid)
%                     fprintf('[?] Could not find swarm data....\n');
%                     stat=0;
%                 else
%                     fprintf('[I] Found swarm data ....\n');
%                 end
%             else 
%                 fprintf('[!] Could not find any data ...\n');
%                 stat=0;
%             end
        end
        if stat
            fprintf('[I] Found automatically tracked mosquito data ....\n');
            menuitem.status='x';
        end
        
    case 3
        stat=1;
    
        if stat
            datfile=dir([menuitem.info.floc, '/output/data/', 'data_mq_', menuitem.info.expname, '.mat']);
%             for ii=1:size(datfile,1)
%                 fprintf('[I] %s found ...\n', datfile(ii).name);
%             end
            if size(datfile,1)
                fprintf('[I] Found manually joined tracks ....\n');
                stat=1;
            else
                stat=0;
            end
        end
        
    
        if stat
            menuitem.status='x';
        end
%     case 5
%         stat=0;
% 
%         
%         if stat
%             datfile=dir([menuitem.info.floc, '/output/data/', 'data_mq_', menuitem.info.expname, '.mat']);        
%             datfile_F=dir([menuitem.info.floc, '/output/data/', 'data_mq_', menuitem.info.expname, '_F.mat']);
%             for ii=1:size(datfile_F,1)
%                 fprintf('[I] Filtered data (%s) found ...\n', datfile_F(ii).name);
%             end
% 
%             if size(datfile_F,1)
%     %            load([menuitem.info.floc, '/output/data/', datfile_F(1).name], 'Xh');
%     %            Xh_F=Xh;
%     %            Xh_F=Xh_F(any(Xh_F,2), any(Xh_F,1));
%     %            clear('Xh');
%     %            load([menuitem.info.floc, '/output/data/', datfile(1).name], 'Xh'); 
%     %            Xh=Xh(any(Xh,2), any(Xh,1));
%     %            if numel(Xh_F) ~= numel(Xh)
%     %                fprintf('[!] Filtered data is old..regenerate filtered data \n');
%     %            else
%     %                stat=1;
%     %            end
%                 if (datenum(datfile_F.date)-datenum(datfile.date) < 0)
%                     fprintf('[!] Filtered data is old..regenerate filtered data \n');
%                     stat=0;
%                 end
%             else
%                 fprintf('[!] Could not find filtered data...\n');
%                 stat=0;
%             end
%         end
%         if stat
%             menuitem.status='x';
%         end
%         
%     case 6
%         stat=0;
%        
%         if stat
%         
%             datfile_F=dir([menuitem.info.floc, '/output/data/', 'data_mq_', menuitem.info.expname, '_F.mat']);
%             datfile_F_W=dir([menuitem.info.floc, '/output/data/', 'data_mq_', menuitem.info.expname, '_F_W.mat']);        
%             for ii=1:size(datfile_F_W,1)
%                 fprintf('[I] World ref. data (%s) found ...\n', datfile_F_W(ii).name);
%             end
% 
%             if size(datfile_F_W,1)
%     %            load([menuitem.info.floc, '/output/data/', datfile_F_W(1).name], 'Xh');
%     %            Xh_F_W=Xh;
%     %            Xh_F_W=Xh_F_W(any(Xh_F_W,2), any(Xh_F_W,1));
%     %            clear('Xh');
%     %            load([menuitem.info.floc, '/output/data/', datfile_F(1).name], 'Xh'); 
%     %            Xh_F=Xh(any(Xh,2), any(Xh,1));
%     %            if numel(Xh_F_W) ~= numel(Xh_F)
%     %                fprintf('[!] World reference data is old... regenerate world reference data \n');
%     %            else
%     %                stat=1;
%     %            end
%                 if (datenum(datfile_F_W.date)-datenum(datfile_F.date) < 0)
%                     fprintf('[!] World reference data is old... regenerate world reference data  \n');
%                     stat=0;
%                 end
% 
%             else
%                 fprintf('[!] Could not find World reference data...\n');
%                 stat=0;
%             end
%         end        
%         if stat
%             menuitem.status='x';
%         end
	otherwise
end


function dispmenu(menu)
fprintf('_______________________________________________\n');
fprintf('      ========== Tracker checklist ========\n');
fprintf('\n');
nitems=size(menu,2);
for ii=1:nitems
	fprintf('%d.[%s] %s.\n', ii, menu(ii).status, menu(ii).str);
end
fprintf('\n');