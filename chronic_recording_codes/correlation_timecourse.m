
%%% Load analyzed correlation data in the workspace %%%

% in this order: light, dark (every correlation value corresponds to a 5
% min unit (for the x-axis time indices)

%% light-dark (12-12h)
bin_num = 1:(size(spike_corr_light{5}{1},1)*3);
bin_ct = size(spike_corr_light{5}{1},1);
sw = 5; %sliding window size in minutes

%first col stores light session, second col stores dark session, thrid col
%stores light+cpp session
all_corr_by_ani = cell(1,numel(spike_corr_dark));
corr_cat_by_ani = cell(1,numel(spike_corr_dark));

for ai = 1:numel(spike_corr_dark)

    if isempty(spike_corr_dark{ai})
        continue
    else
        curr_time_dark(1:size(spike_corr_dark{ai},2)) = NaN;
        curr_time_light(1:size(spike_corr_light{ai},2)) = NaN;
        curr_time_light_cpp(1:size(spike_corr_light_cpp{ai},2)) = NaN;
        
        all_corr_by_ani{ai}(size(spike_corr_dark{ai}{1},1),1) = NaN;
    
        for ti = 1:size(spike_corr_dark{ai}{1},1)
            for pi = 1:size(spike_corr_dark{ai},2)
           
                 curr_time_dark(pi) = spike_corr_dark{ai}{pi}(ti,1);
                 curr_time_light(pi) = spike_corr_light{ai}{pi}(ti,1);
                 curr_time_light_cpp(pi) = spike_corr_light_cpp{ai}{pi}(ti,1);
    
            end

%                  if ti > 110 
%                      curr_time_light_cpp
%                  end
    
            all_corr_by_ani{ai}(ti,1) = mean(curr_time_light,'omitnan');
            all_corr_by_ani{ai}(ti,2) = mean(curr_time_dark,'omitnan'); 
            all_corr_by_ani{ai}(ti,3) = mean(curr_time_light_cpp,'omitnan');  

        end
    end

    %concatenante light and dark, cpp at the end
    corr_cat_by_ani{ai}(bin_num,1) = ...,
        cat(1,all_corr_by_ani{ai}(1:end,1),all_corr_by_ani{ai}(1:end,2),all_corr_by_ani{ai}(1:end,3));
end
    
meta_corr = NaN(numel(bin_num),numel(find(animal_hit))); % move per animals values from cells to columns
meta_corr_ave = NaN(numel(bin_num),1);

%normalize all correlation values to the mean of the dark session for each
%animal (and then for get the meta average)
mean_corr_dark = NaN(1,numel(find(animal_hit)));
meta_corr_norm = NaN(numel(bin_num),numel(find(animal_hit)));
meta_corr_ave_norm = NaN(numel(bin_num),1);
meta_corr_sem_norm = NaN(numel(bin_num),1);

ani_ct = 0;
for ani = 1:numel(animal_hit)
    if animal_hit(ani) == 0
        continue
    else
        ani_ct = ani_ct + 1;
        for bi = 1:size(corr_cat_by_ani{ani},1)
            meta_corr(bi,ani_ct) = corr_cat_by_ani{ani}(bi,1);
            
        end

         mean_corr_dark(ani_ct) = mean(meta_corr(bin_ct+1:bin_ct*2,ani_ct),'omitnan');
         meta_corr_norm(:,ani_ct) = meta_corr(:,ani_ct)./mean_corr_dark(ani_ct);
    end
end

for bii = 1:size(meta_corr_ave,1)
    meta_corr_ave(bii,1) = mean(meta_corr(bii,1:size(meta_corr,2)),'omitnan');
    meta_corr_ave_norm(bii,1) = mean(meta_corr_norm(bii,1:size(meta_corr_norm,2)),'omitnan');
    meta_corr_sem_norm(bii,1) = nansem(meta_corr_norm(bii,1:size(meta_corr,2)));
end

% %% calculating slopes
% slope_all = NaN(1,3)
% corr_norm_all = NaN(bin_ct,3);
% fit_all = cell(1,3);
% r_square = NaN(1,3);
% rmse = NaN(1,3);
% 
% bin_session = (0.5:0.5:24*0.5)';
% for si = 1:3
%     corr_norm_all(:,si) = meta_corr_ave_norm(bin_ct*(si-1)+1:bin_ct*si);
%     [f,gof] = fit(bin_session,corr_norm_all(1:24,si),'poly1');
%     slope_all(si) = f.p1;
%     r_square(si) = gof.rsquare;
%     rmse(si) = gof.rmse;
%     
% %     figure()
% % 
% %     plot(f,bin_session,corr_norm_all(1:24,si))
% end

%% percentage of change- for light/light+CPP
%baseline- first hour after light on
% light_bl = meta_corr_ave_norm(1:2);
% light_cpp_bl = meta_corr_ave_norm(bin_ct*2+1:bin_ct*2+2);
% 
% %values after CPP inj- 4 hours after light on (starting from bin 9)
% light_val = meta_corr_ave_norm(9:bin_ct);
% light_cpp_val = meta_corr_ave_norm(bin_ct*2+9:bin_ct*3);
% 
% light_val_prc = light_val./mean(light_bl,'omitnan')*100;
% light_cpp_prc = light_cpp_val./mean(light_cpp_bl,'omitnan')*100;
% 
% Light.prc = light_val_prc(1:10);
% Light_CPP.prc = light_cpp_prc(1:10);


%% plotting light/dark raw
time_indx = cat(2,sw:sw:bin_ct*sw,(bin_ct+sw+1)*sw:sw:(bin_ct*2+sw)*sw); %in minutes
colors = parula(5);

figure()
hold on
y_lim = [0 0.13];
%set(gca,'yscale','log');


rectangle('Position',[0 y_lim(1) 720 y_lim(2)-y_lim(1)],...
    'FaceColor','w','EdgeColor','none')

rectangle('Position',[720 y_lim(1) 720 y_lim(2)-y_lim(1)],...
    'FaceColor',[.85 .85 .85] ,'EdgeColor','none')

for afi = 1:size(meta_corr,2)

    plot(time_indx(1:bin_ct),meta_corr(1:bin_ct,afi),'Color',colors(afi,1:3))
    plot(time_indx(bin_ct+1:bin_ct*2),meta_corr(bin_ct+1:bin_ct*2,afi),'Color',colors(afi,1:3))

    plot(time_indx(1:bin_ct),meta_corr_ave(1:bin_ct,1),'k','LineWidth',3)
    plot(time_indx(bin_ct+1:bin_ct*2),meta_corr_ave(bin_ct+1:bin_ct*2,1),'k','LineWidth',3)
end
 
ax = gca;
ax.XLim = [0 1440];
ax.YLim = y_lim;
ax.XLabel.String = 'Experiment Time in ZT (minute)';
ax.YLabel.String = 'Correlation';
ax.FontSize = 14;
ax.LineWidth = 2;
ax.YLabel.FontSize = 16;
ax.XLabel.FontSize = 16;
title('Light/Dark')

hold off

%% light/dark normalized
time_indx = cat(2,sw:sw:bin_ct*sw,(bin_ct+sw+1)*sw:sw:(bin_ct*2+sw)*sw); %in minutes
sem_CI = generate_sem_CI(meta_corr_ave_norm,meta_corr_sem_norm);

%colors = parula(5);

fig1 = figure('Position',[500 100 800 500]);
hold on
y_lim = [0 3];
%set(gca,'yscale','log');


rectangle('Position',[0 y_lim(1) 720 y_lim(2)-y_lim(1)],...
    'FaceColor','w','EdgeColor','none')

rectangle('Position',[720 y_lim(1) 720 y_lim(2)-y_lim(1)],...
    'FaceColor',[.85 .85 .85] ,'EdgeColor','none')


for di = 1:(bin_ct-1)
    x_co = [time_indx(di) time_indx(di+1) time_indx(di+1) time_indx(di)];
    y_co = [sem_CI(di,1) sem_CI(di+1,1) sem_CI(di+1,2) sem_CI(di,2)];
    v = [x_co(1) y_co(1); x_co(2) y_co(2);x_co(3) y_co(3);x_co(4) y_co(4)];
    f = [1 2 3 4];


    patch('Faces',f,'Vertices',v,...
        'FaceColor',[0.5 0.5 0.5],'EdgeColor','none','FaceAlpha',0.3)
end

for di = bin_ct+1:(bin_ct*2-1)
    x_co = [time_indx(di) time_indx(di+1) time_indx(di+1) time_indx(di)];
    y_co = [sem_CI(di,1) sem_CI(di+1,1) sem_CI(di+1,2) sem_CI(di,2)];
    v = [x_co(1) y_co(1); x_co(2) y_co(2);x_co(3) y_co(3);x_co(4) y_co(4)];
    f = [1 2 3 4];


    patch('Faces',f,'Vertices',v,...
        'FaceColor',[0.5 0.5 0.5],'EdgeColor','none','FaceAlpha',0.4)
end

for afi = 1:size(meta_corr_norm,2)

%     plot(time_indx(1:bin_ct),meta_corr_norm(1:bin_ct,afi),'Color',colors(afi,1:3))
%     plot(time_indx(bin_ct+1:bin_ct*2),meta_corr_norm(bin_ct+1:bin_ct*2,afi),'Color',colors(afi,1:3))

    plot(time_indx(1:bin_ct),meta_corr_ave_norm(1:bin_ct,1),'k','LineWidth',2)
    plot(time_indx(bin_ct+1:bin_ct*2),meta_corr_ave_norm(bin_ct+1:bin_ct*2,1),'k','LineWidth',2)

end

 
ax = gca;
ax.Layer = 'top';
ax.XLim = [0 1440];
ax.XTick = 0:180:1440;
ax.YLim = y_lim;
ax.XLabel.String = 'Experiment Time in ZT (minute)';
ax.YLabel.String = 'Normalized Correlation';
ax.FontSize = 14;
ax.LineWidth = 2;
ax.YLabel.FontSize = 16;
ax.XLabel.FontSize = 16;
title('Light/Dark')

hold off

exportgraphics(fig1,'light_dark_normalized_correlation.pdf',...,
    'ContentType','vector','BackgroundColor','none')


%% plotting Light+CPP (normalized to dark)
time_indx_cpp = (5:5:139*5);
cpp_region = (279:417);
% meta_corr_cpp_norm = NaN(139,5);
% meta_corr_cpp_norm_ave = NaN(139,1);

% for api = 1:5
%     for cpi = 1:size(meta_corr_cpp_norm,1)
%         meta_corr_cpp_norm(cpi,api) = meta_corr(cpi+278,api)/meta_corr(cpi+139,api);
%     end
% end
% 
% for cpii = 1:139
%     meta_corr_cpp_norm_ave(cpii,1) = mean(meta_corr_cpp_norm(cpii,:),'omitnan');
% end

% % normalized
% figure()
% hold on

%  for afi = 1:size(meta_corr,2)
% 
%     plot(time_indx_cpp,meta_corr_cpp_norm(1:139,afi),'Color',colors(afi,1:3))
%     plot(time_indx_cpp,meta_corr_cpp_norm_ave(:,1),'k','LineWidth',3)
%  end
% 
% ax = gca;
% ax.XLim = [0 720];
% %ax.YLim = y_lim;
% ax.XLabel.String = 'Experiment Time in ZT (minute)';
% ax.YLabel.String = 'Normalzied Correlation';
% ax.FontSize = 14;
% ax.LineWidth = 2;
% ax.YLabel.FontSize = 16;
% ax.XLabel.FontSize = 16;
% 
% title('Light+CPP')
% 
% hold off

%raw
figure()
hold on
 
 for afi = 1:size(meta_corr,2)

    plot(time_indx_cpp,meta_corr(cpp_region,afi),'Color',colors(afi,1:3))
    plot(time_indx_cpp,meta_corr_ave(cpp_region,1),'k','LineWidth',3)
 end

ax = gca;
ax.XLim = [0 720];
%ax.YLim = y_lim;
ax.XLabel.String = 'Experiment Time in ZT (minute)';
ax.YLabel.String = 'Correlation';
ax.FontSize = 14;
ax.LineWidth = 2;
ax.YLabel.FontSize = 16;
ax.XLabel.FontSize = 16;

title('Light+CPP')

hold off

%normalized
figure()
hold on
 
 for afi = 1:size(meta_corr,2)

    plot(time_indx_cpp,meta_corr_norm(cpp_region,afi),'Color',colors(afi,1:3))
    plot(time_indx_cpp,meta_corr_ave_norm(cpp_region,1),'k','LineWidth',3)
 end

ax = gca;
ax.XLim = [0 720];
%ax.YLim = y_lim;
ax.XLabel.String = 'Experiment Time in ZT (minute)';
ax.YLabel.String = 'Normalized Correlation';
ax.FontSize = 14;
ax.LineWidth = 2;
ax.YLabel.FontSize = 16;
ax.XLabel.FontSize = 16;

title('Light+CPP')

hold off

%% all three session together, normalized
time_indx_all = cat(2,sw:sw:bin_ct*sw,(bin_ct+1)*sw:sw:(bin_ct*2+1)*sw,...
    (bin_ct*2+1)*sw:sw:bin_ct*3*sw); %in minutes
colors = parula(5);

%normalized
figure();
hold on
y_lim = [0 2];
%set(gca,'yscale','log');


rectangle('Position',[0 y_lim(1) 720 y_lim(2)-y_lim(1)],...
    'FaceColor','w','EdgeColor','none');

rectangle('Position',[720 y_lim(1) 720 y_lim(2)-y_lim(1)],...
    'FaceColor',[.85 .85 .85] ,'EdgeColor','none')

rectangle('Position',[1440 y_lim(1) 720 y_lim(2)-y_lim(1)],...
    'FaceColor','w' ,'EdgeColor','none')


for afi = 1:size(meta_corr_norm,2)

    plot(time_indx_all(1:bin_ct),meta_corr_norm(1:bin_ct,afi),'Color',colors(afi,1:3))
    plot(time_indx_all(bin_ct+1:bin_ct*2),meta_corr_norm(bin_ct+1:bin_ct*2,afi),'Color',colors(afi,1:3))
    plot(time_indx_all(bin_ct*2+1:bin_ct*3),meta_corr_norm(bin_ct*2+1:bin_ct*3,afi),'Color',colors(afi,1:3))


    plot(time_indx_all(1:bin_ct),meta_corr_ave_norm(1:bin_ct,1),'k','LineWidth',3)
    plot(time_indx_all(bin_ct+1:bin_ct*2),meta_corr_ave_norm(bin_ct+1:bin_ct*2,1),'k','LineWidth',3)
    plot(time_indx_all(bin_ct*2+1:bin_ct*3),meta_corr_ave_norm(bin_ct*2+1:bin_ct*3,1),'k','LineWidth',3)
% 
%     scatter(time_indx_all(1:bin_ct),meta_corr_norm(1:bin_ct,afi),'Color',colors(afi,1:3))
%     scatter(time_indx_all(bin_ct+1:bin_ct*2),meta_corr_norm(bin_ct+1:bin_ct*2,afi),'Color',colors(afi,1:3))
%     scatter(time_indx_all(bin_ct*2+1:bin_ct*3),meta_corr_norm(bin_ct*2+1:bin_ct*3,afi),'Color',colors(afi,1:3))


%     scatter(time_indx_all(1:bin_ct),meta_corr_ave_norm(1:bin_ct,1),'k','LineWidth',3)
%     scatter(time_indx_all(bin_ct+1:bin_ct*2),meta_corr_ave_norm(bin_ct+1:bin_ct*2,1),'k','LineWidth',3)
%     scatter(time_indx_all(bin_ct*2+1:bin_ct*3),meta_corr_ave_norm(bin_ct*2+1:bin_ct*3,1),'k','LineWidth',3)

end
 
ax = gca;
ax.Layer = 'top';
ax.XLim = [0 2160];
ax.YLim = y_lim;
ax.XTick = 0:360:2160;

ax.XLabel.String = 'Experiment Time in ZT (minute)';
ax.YLabel.String = 'Normalized Correlation';
ax.FontSize = 14;
ax.LineWidth = 2;
ax.YLabel.FontSize = 16;
ax.XLabel.FontSize = 16;
title('Light/Dark/Light+CPP')
text(1700,4.5, 'CPP' ,'FontSize',16,'Color','k')

hold off

%% all three session together, raw
time_indx_all = cat(2,5:5:139*5,145*5:5:283*5,289*5:5:427*5); %in minutes
colors = parula(5);

%normalized
f = figure();
hold on
y_lim = [0 0.3];
%set(gca,'yscale','log');


rectangle('Position',[0 y_lim(1) 720 y_lim(2)-y_lim(1)],...
    'FaceColor','w','EdgeColor','none');

rectangle('Position',[720 y_lim(1) 720 y_lim(2)-y_lim(1)],...
    'FaceColor',[.85 .85 .85] ,'EdgeColor','none')

rectangle('Position',[1440 y_lim(1) 720 y_lim(2)-y_lim(1)],...
    'FaceColor','w' ,'EdgeColor','none')


for afi = 1:size(meta_corr,2)

    plot(time_indx_all(1:139),meta_corr(1:139,afi),'Color',colors(afi,1:3))
    plot(time_indx_all(140:278),meta_corr(140:278,afi),'Color',colors(afi,1:3))
    plot(time_indx_all(279:417),meta_corr(279:417,afi),'Color',colors(afi,1:3))


    plot(time_indx_all(1:139),meta_corr_ave(1:139,1),'k','LineWidth',3)
    plot(time_indx_all(140:278),meta_corr_ave(140:278,1),'k','LineWidth',3)
    plot(time_indx_all(279:417),meta_corr_ave(279:417,1),'k','LineWidth',3)

end
 
ax = gca;
ax.Layer = 'top';
ax.XLim = [0 2160];
ax.YLim = y_lim;
ax.XTick = 0:360:2160;

ax.XLabel.String = 'Experiment Time in ZT (minute)';
ax.YLabel.String = 'Correlation';
ax.FontSize = 14;
ax.LineWidth = 2;
ax.YLabel.FontSize = 16;
ax.XLabel.FontSize = 16;
title('Light/Dark/Light+CPP- Raw Correlation')
text(1700,0.25, 'CPP' ,'FontSize',16,'Color','k')

hold off