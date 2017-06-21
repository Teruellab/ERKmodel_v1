if ~exist('erk_paired','var')
    load('/Volumes/labdata/brooks/Data/erk_sim.mat')
end

%%
setcolors;
loadcolormaps;
%%
for idx = 1:size(erk_paired,3)
    erk1 = squeeze(erk_paired(:,:,idx));
    erk2 = squeeze(erk_unpaired(:,:,idx));

    % - - - HEATMAPS - - - - - - - - - - - - 
    graph_lim = prctile([erk1(:);erk2(:)],[2 99]);
    [~,ord1] = sort(sum(erk1(:,1:400),2),'descend');
    [~,ord2] = sort(sum(erk2(:,1:400),2),'descend');
    figure('Position',positionfig(1000,400),'Name',['RasGTP: ', num2str(ras_doses(idx))]);
    ha = tight_subplot(1,2,0.05 ,[0.15 .1]);
    imagesc(erk1(ord1,:),'Parent',ha(1)), set(ha(1),'CLim',graph_lim)
    imagesc(erk2(ord2,:),'Parent',ha(2)), set(ha(2),'CLim',graph_lim)
    colormap(colormaps.viridis)
    title(ha(1),'Paired ERK and MEK')
    title(ha(2),'Unpaired')
    for i = 1:length(ha)
        set(ha(i),'YTick',[],'XTick',1:3600:7201, 'XTickLabel',{'0','1','2'})
        xlabel(ha(i),'Time (hrs)')
    end
    
% - - POPULATION AVERAGES - - - - 
    figure('Position',positionfig(1000,400),'Name',['RasGTP: ', num2str(ras_doses(idx))]);
    ha = tight_subplot(1,2,0.1 ,0.15);
    plot((1:size(erk1,2))/3600, [mean(erk1)',mean(erk2)'],'Parent',ha(1),'LineWidth',2)
    plot((1:size(erk1,2))/3600, [median(erk1)',median(erk2)'],'Parent',ha(2),'LineWidth',2)
    for j = 1:2
        xlabel(ha(j),'Time (hrs)')
        ylabel(ha(j),'Phosphorylated ERK')
        set(ha(j),'ColorOrder',cell2mat(colors.peacock'))
        legend(ha(j),{'paired ERK/MEK', 'unpaired'},'Location','southeast')
    end

% - - - - - (MAX) HISTOGRAMS - - - - -
    
    ovr_erk1 = sum(erk1(:,1:end),2); 
    ovr_erk2 = sum(erk2(:,1:end),2); 
    
%     ovr_erk1 = mean(erk1(:,end-50:end),2); 
%     ovr_erk2 = mean(erk2(:,end-50:end),2); 
    graph_lim = (prctile([ovr_erk1;ovr_erk2],[1 99]));
    
    ovr_erk1(ovr_erk1>graph_lim(2)) = graph_lim(2);
    ovr_erk2(ovr_erk2>graph_lim(2)) = graph_lim(2);
    ovr_erk1(ovr_erk1<graph_lim(1)) = graph_lim(1);
    ovr_erk2(ovr_erk2<graph_lim(1)) = graph_lim(1);
    %erk1(erk1>graph_lim(2)) = graph_lim(2);
    %erk2(erk2>graph_lim(2)) = graph_lim(2);
    figure('Name',['RasGTP: ', num2str(ras_doses(idx))]),
    hold on;
    histogram(ovr_erk1,linspace(graph_lim(1),graph_lim(2),32))
    histogram(ovr_erk2,linspace(graph_lim(1),graph_lim(2),32))

    hold off
end

%%
idx = [3 9];

for k = 1:length(idx)
t = (0:7200)/3600;
traj{1} = squeeze(erk_paired(:,:,idx(k)));
traj{2} = squeeze(erk_unpaired(:,:,idx(k)));
    
for j = 1:length(traj)
% Define discrete data subsets (sorted by early activity)
[~,order] = sort(prctile(traj{j}(:,1:2000),97,2),'descend');
rez = 0.10;

figs.splitavg(j+(k-1)*2) = figure('Position',positionfig(320,200),'PaperPositionMode','auto');
ha = tight_subplot(1,1);
set(ha(1),'ColorOrder',colormaps.plasma([220:-round(rez*219):1],:))
hold on
for i = 0:rez:(1-rez)
    lo = round(i*length(order))+1;
    hi = round((i+rez)*length(order));
    plot(t, nanmedian(traj{j}(order(lo:hi),:)),'LineWidth',1.5)
end
%plot(t, nanmedian(all_norm{1}),'LineWidth',1.5,'Color','k')
hold off
grid on
set(ha(1),'XLim',[-0.05 2],'Box','on','XTick',[0 1 2],'LineWidth',1.2,'YLim',[0 prctile(traj{2}(:),99.9)*1.1])
disp(['Max x = ',num2str(prctile(traj{2}(:),99.9)*1.1)])
end
end

%%

idx = [3 9];
peak_time = cell(2,2);
peak_height = cell(2,2);
ss_val = cell(2,2);

for k = 1:length(idx)
t = (0:7200)/3600;
traj{1} = squeeze(erk_paired(:,:,idx(k)));
traj{2} = squeeze(erk_unpaired(:,:,idx(k)));
    
for j = 1:length(traj)
% Define discrete data subsets (sorted by early activity)
    peak_time{k,j} = zeros(size(traj{j},1),1);
    peak_height{k,j} = max(traj{j},[],2);
    ss_val{k,j} = mean(traj{j}(:,end-200:end),2);
    for i = 1:size(traj{j},1)
        curr = traj{j}(i,:);
        peak_time{k,j}(i) = find(curr>prctile(curr,99.9),1,'first'); 
    end
    
    
end
end


%%

idx = 16;
erk1 = squeeze(erk_paired(:,:,idx));
erk2 = squeeze(erk_unpaired(:,:,idx));

ovr_erk1 = sum(erk1(:,1:end),2); 
ovr_erk2 = sum(erk2(:,1:end),2); 

%     ovr_erk1 = mean(erk1(:,end-50:end),2); 
%     ovr_erk2 = mean(erk2(:,end-50:end),2); 
graph_lim = (prctile([ovr_erk1;ovr_erk2],[2 98]));

ovr_erk1(ovr_erk1>graph_lim(2)) = graph_lim(2);
ovr_erk2(ovr_erk2>graph_lim(2)) = graph_lim(2);
ovr_erk1(ovr_erk1<graph_lim(1)) = graph_lim(1);
ovr_erk2(ovr_erk2<graph_lim(1)) = graph_lim(1);
%erk1(erk1>graph_lim(2)) = graph_lim(2);
%erk2(erk2>graph_lim(2)) = graph_lim(2);

bins = linspace(graph_lim(1),graph_lim(2),30);

n1 = histcounts(ovr_erk2,bins); n1 = [n1 0]/sum(n1);
n2 = histcounts(ovr_erk1,bins); n2 = [n2 0]/sum(n2);

[x1, y1] = stairs(bins,n1); x1 = [x1(1); x1]; y1 = [0 ; y1];
[x2, y2] = stairs(bins,n2); x2 = [x2(1); x2]; y2 = [0 ; y2];



figs.AUChist = figure('Name',['RasGTP: ', num2str(ras_doses(idx))],'Position',positionfig(255,160));
ha = tight_subplot(1,1);
hold on;
area(x1,y1,'FaceColor',colors.lavender,'FaceAlpha',0.9,'EdgeColor',colors.lavender/2)
area(x2,y2,'FaceColor',colors.orange,'FaceAlpha',0.4,'EdgeColor',colors.lavender/2)
hold off
set(gca,'XLim',[min(bins)-range(bins)*0.03 max(bins)+range(bins)*0.03],'XTickLabel',{},'YTickLabel',{},...
'YLim',[0 0.15],'YTick',0:.05:.25)
legend({'Correlated MEK/ERK', 'Uncorrelated'})
    
    

figure, 
ha = tight_subplot(1,1);
plot((1:size(erk1,2))/3600, [mean(erk1)',mean(erk2)'],'Parent',ha(1),'LineWidth',2)
