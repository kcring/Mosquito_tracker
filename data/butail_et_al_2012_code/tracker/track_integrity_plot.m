function track_integrity_plot
close all

run initproc

datafile=[dataloc, sprintf('data_mq_auto_%s.mat', expname)];
fs=10;

if exist(datafile, 'file')
    load(datafile, 'Xh')
else
    fprintf('[!] auto datafile not found...');
end

Xh=Xh(:,any(Xh,1));

nzvals=Xh(Xi.ri(3):Xi.nX:end,:);


nzvals(nzvals~=0)=1;
ntracked=sum(nzvals); % at each time step
tracklengths=sum(nzvals,2);
tracklengths=tracklengths(tracklengths~=0);
avg_track_length=mean(tracklengths);

nzvals(nzvals==0)=nan;
nzvals=nzvals.*([1:size(nzvals,1)]'*ones(1,size(nzvals,2)));
figure(1); gcf; clf;

plot((ones(size(nzvals,1),1)*[1:size(nzvals,2)])', nzvals', 'LineWidth',2);
hold on;

axis tight
ylim=get(gca,'Ylim');
xlim=get(gca,'Xlim');
xtl=get(gca,'Xtick');
set(gca, 'XtickLabel', xtl/25);
nmq_plot_range=[50:50:floor(size(Xh,2)/50)*50];
plot(ones(2,1)*nmq_plot_range, ylim'*ones(1,numel(nmq_plot_range)), '--', 'Color', [.5 .5 .5]);
text(nmq_plot_range-2, (ylim(2)+8)*ones(1,numel(nmq_plot_range)), num2str(ntracked(nmq_plot_range)', '%d\n'), 'FontSize', fs);
text((xlim(1)+xlim(2))/3, ylim(2)+20, 'Number of tracks', 'FontSize', fs);

set(gca, 'FontSize', fs);
xlabel('Time(s)');
ylabel(sprintf('Track ID (Avg. track length=%.2f)', avg_track_length));
box off;

print('-dpng', sprintf('%s/output/movies/track_integrity_%s.png', floc, expname));
