function cdata = get_climate_data(floc)
cdata.windspeed=[];
cdata.winddir=[];
if exist([floc, 'output/data/windspeed.csv'], 'file')
    fprintf('[I] No climate data found ...\n');
    cdata.windspeed=csvread([floc, 'output/data/windspeed.csv'],0,1);
end

if exist([floc, 'output/data/winddir.csv'], 'file')
    cdata.winddir=csvread([floc, 'output/data/winddir.csv'],0,1);
end
    