function disp_climate_data
% climate data

addpath ../tracker
run initproc


kestrel_file_list=dir([floc, 'calib/Kestrel_*.csv']);
if(size(kestrel_file_list,1))
    climate_data=csvread([floc, 'calib/', kestrel_file_list(1).name], 2, 1);
else
    fprintf('[!] Climate data not found in %s/calib\n', floc);
end


% input the time from the first line in the kestrel file
fid=fopen([floc, 'calib/', kestrel_file_list(1).name], 'r');
for ii=1:3
    ll=fgetl(fid);
end
fclose(fid);
start_climate=datenum(ll(11:19));

% input the time of the experiment for the first file of the data
start_exp=datenum('19:04:53');
end_exp=start_exp+datenum('0:20:00');


time=climate_data(:,1)-climate_data(1,1);
start_exp=dv2sec(datevec(start_exp-start_climate));
end_exp=dv2sec(datevec(end_exp-start_climate));

% wind speed
plot_data(expname, floc, time, start_exp, end_exp, climate_data(:,4), 'windspeed', '(m/s)', 1);
% wind direction
plot_data(expname, floc, time, start_exp, end_exp, climate_data(:,2), 'winddir', '(deg.)', 2);

% temperature
plot_data(expname, floc, time, start_exp, end_exp, climate_data(:,7), 'temperature', '(F)', 3);

% relative humidity
plot_data(expname, floc, time, start_exp, end_exp, climate_data(:,9), 'relativhumidity', '(%)', 4);

% barometric pressure
plot_data(expname, floc, time, start_exp, end_exp, climate_data(:,13), 'barometricpressure', '(inHg)', 5);

function plot_data(expname, floc, time, start_exp, end_exp, data, desc, units, id)

figure(id); gcf;
subplot(2,1,1);
gca; cla;

plot(time, data);
hold on;
xv=[start_exp*ones(2,1);
    end_exp*ones(2,1)];
yv=[0; max(data); max(data); 0];
fill(xv,yv, 'r', 'EdgeColor', 'none', 'FaceAlpha', .25);
ylabel([ desc, units]);
xlabel('time (s)');
% plot(dv2sec(datevec(start_exp-start_climate))*ones(1,2), [0 max(ws)], , ');
% plot(dv2sec(datevec(end_exp-start_climate))*ones(1,2), [0 max(ws)], 'k');


subplot(2,1,2); 
gca; cla;
[exp_data_start_time exp_start_idx]=min(abs(time-start_exp));
[exp_data_end_time exp_end_idx]=min(abs(time-end_exp));

plot(time(exp_start_idx:exp_end_idx), data(exp_start_idx:exp_end_idx));
ylabel([desc, units]);
xlabel(['time (s) ', expname], 'Interpreter', 'none');
print('-dpng', sprintf('%s/output/movies/climate_%s_%s.png', floc, desc, expname));
expdata=[time(exp_start_idx:exp_end_idx), data(exp_start_idx:exp_end_idx)];
csvwrite(sprintf('%s/output/data/%s.csv', floc, desc), expdata);

function secs=dv2sec(dv)

secs=dv(6)+dv(5)*60+dv(4)*3600;