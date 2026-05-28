clear variables

addpath ../tracker/

% olddata
oloc='/Volumes/DUMP01/2016/nimr/20160326_192011/';

t1file=dir([oloc, '/output/data/data_mq_EXP_*.mat']);

load([oloc, 'output/data/', t1file(1).name]);

Xh_old=Xh(any(Xh,2),:);
cams_old=cams;

% newdata
nloc='/Users/sachit/Dropbox/EASeL/marker_exp/output_data/20160326_192011/';
t1file=dir([nloc, '/output/data/data_mq_EXP_*.mat']);

load([nloc, 'output/data/', t1file(1).name]);

Xh_new=Xh;
cams_new=cams;
rw=Xh_old*0;

for ii=1:size(Xh_old,1)/6
    pix1=w2cam_nd(Xh_old(6*(ii-1)+1:6*(ii-1)+3,:), cams_old(1));
    pix2=w2cam_nd(Xh_old(6*(ii-1)+1:6*(ii-1)+3,:), cams_old(2));
     
    for jj=1:size(pix1,2)
        if ~isnan(pix1(1,jj))
            rw(6*(ii-1)+1:6*(ii-1)+3,jj)=lsTriangulate([pix1(:,jj),pix2(:,jj)], cams_new);
        end
    end
end
rw(rw==0)=nan;
Xh_old(Xh_old==0)=nan;
figure(1); gcf;clf;
for ii=1:size(Xh_old,1)/6
    plot3(Xh_old(6*(ii-1)+1,:), Xh_old(6*(ii-1)+2,:), Xh_old(6*(ii-1)+3,:));
    hold on;
    fnan=find(~isnan(Xh_old(6*(ii-1)+1,:)));
    text(Xh_old(6*(ii-1)+1,fnan(1)), Xh_old(6*(ii-1)+2,fnan(1)), Xh_old(6*(ii-1)+3,fnan(1)), ...
            sprintf('%d', ii), 'fontsize', 24);

end
mkrid=input('marker id is: ');
valmq=input('mqs to migrate e.g. [2, 3, 5]: ');

% update the marker as 99
Xh_new(6*(99-1)+1:6*(99-1)+3,:)=rw(6*(mkrid-1)+1:6*(mkrid-1)+3,:);


% update all valid mosquitoes from 98 onwards
nxt=99-1;
for jj=valmq
    Xh_new(6*(nxt-1)+1:6*(nxt-1)+3,:)=rw(6*(jj-1)+1:6*(jj-1)+3,:);
end

figure(1); gcf; clf;
Xhplot=Xh_new;
Xhplot(Xhplot==0)=nan;
plot3(Xhplot(1:6:end,:)', Xhplot(2:6:end,:)', Xhplot(3:6:end,:)');

input('ok to save: ');

Xh_new(isnan(Xh_new))=0;
Xh=Xh_new;
save([nloc, 'output/data/', t1file(1).name], 'Xh', '-append');



