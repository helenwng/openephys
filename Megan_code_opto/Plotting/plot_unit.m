function plot_unit(animal,subfold,exp_type,unit)

% animal = animal name (e.g., 'M21')
% subfold = possible subfolder within animal name (e.g., 'LP' or 'pen1')
% exp_type (e.g., 'intensities' or 'trains')
% unit = name of unit to plot (eg., 6457)

% load appropriate results file
main_dir = 'H:\LPproject\LPresults';
cd(main_dir)
all_figs = 'H:\LPproject\LPresults\UnitFigures';

s=dir;
names = {s.name};
which_fold = find(cellfun(@(x) contains(x,animal,'ignorecase',1),names));
if length(which_fold)>1
    which_fold = find(cellfun(@(x) contains(x,animal,'ignorecase',1)&&contains(x,subfold,'ignorecase',1),names));
end
animal_dir = fullfile(main_dir,names{which_fold});
ss=dir(animal_dir);
ss_names = {ss.name};
which_subfold = find(cellfun(@(x) contains(x,exp_type,'ignorecase',1),ss_names));
exp_dir = fullfile(animal_dir,ss_names{which_subfold});
sss = dir(exp_dir);
sss_names = {sss.name};
results_file = sss_names{cellfun(@(x) contains(x,'results','ignorecase',1),sss_names)};
load(fullfile(exp_dir,results_file))

% find desired unit
all_units = [unitinfo.name];
unit_ind = find(all_units==unit);
unit_rast = unitinfo(unit_ind).rast;
unit_psthV = FRs(unit_ind).psthVisual;
unit_psthB = FRs(unit_ind).psthBlank;

% example psth - creON ai32
fig_name = sprintf('%s\\%s_%s_clust%d_fig',all_figs,animal,exp_type,unit);
figure;
subplot(1,3,1)
light_start = params.av_light_start;
pulse_dur = params.pulse_dur;
totaltime = params.prestim+params.stimtime+params.poststim;
make_raster_plot_v2(unit_rast,params.prestim,totaltime,params.trial_type(:,contains(params.IVs,'light_bit')),light_start,pulse_dur)
xlabel('Time from stimulus onset (s)','fontsize',14)
ylabel('Trial (by condition)','fontsize',14)
set(gca,'FontSize',14,'fontname','arial','fontweight','bold','linewidth',2);

subplot(1,3,2)
binsize = .025;
make_psth_plot_v2(unit_psthV,binsize,params.prestim,params.stimtime,totaltime,params.trial_type(:,contains(params.IVs,'light_bit')),params.av_light_start,params.light_dur);
xlabel('Time from stimulus onset (s)','fontsize',14)
ylabel('Average FR (spks/s)','fontsize',14)
set(gca,'FontSize',14,'fontname','arial','fontweight','bold','linewidth',2);

subplot(1,3,3)
binsize = .025;
make_psth_plot_v2(unit_psthB,binsize,params.prestim,params.stimtime,totaltime,params.trial_type(:,contains(params.IVs,'light_bit')),params.av_light_start,params.light_dur);
xlabel('Time from stimulus onset (s)','fontsize',14)
ylabel('Average FR (spks/s)','fontsize',14)
set(gca,'FontSize',14,'fontname','arial','fontweight','bold','linewidth',2);

set(gcf, 'Position', [100, 100, 1500, 400])
print(gcf,'-dpng',fig_name)
% print2eps(save_fig_name,clust304_fig)
print(gcf,'-painters','-depsc',fig_name)

if strcmpi(exp_type,'trains') || strcmpi(subfold,'inLP')
    figure;
    make_raster_plot_v2(unit_rast(:,round(min(light_start))*1000-100:round(min(light_start))*1000+299),-((round(min(light_start))-.1)-params.prestim),length(round(min(light_start))*1000-100:round(min(light_start))*1000+299)/1000,params.all_light,.1*ones(length(light_start)),pulse_dur)
    xlabel('Time from stimulus onset (s)','fontsize',14)
    ylabel('Trial (by condition)','fontsize',14)
    set(gca,'FontSize',14,'fontname','arial','fontweight','bold','linewidth',2);
    print(gcf,'-dpng',strcat(fig_name,'_zoom'))
% print2eps(save_fig_name,clust304_fig)
    print(gcf,'-painters','-depsc',strcat(fig_name,'_zoom'))

end

return


