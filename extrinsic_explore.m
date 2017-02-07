if ~exist('erk_paired','var')
    load('/Volumes/labdata/brooks/Data/erk_sim.mat')
end
setcolors;
loadcolormaps;
%%
for idx = 1:13
    erk1 = squeeze(erk_paired(:,:,idx));
    erk2 = squeeze(erk_unpaired(:,:,idx));




    % - - - HEATMAPS - - - - - - - - - - - - 
%     graph_lim = prctile([erk1(:);erk2(:)],[2 99]);
%     [~,ord1] = sort(sum(erk1(:,1:1000),2),'descend');
%     [~,ord2] = sort(sum(erk2(:,1:1000),2),'descend');
%     figure('Position',positionfig(1000,400));
%     ha = tight_subplot(1,2,0.05 ,[0.15 .1]);
%     imagesc(erk1(ord1,:),'Parent',ha(1)), set(ha(1),'CLim',graph_lim)
%     imagesc(erk2(ord2,:),'Parent',ha(2)), set(ha(2),'CLim',graph_lim)
%     colormap(colormaps.viridis)
%     title(ha(1),'Paired ERK and MEK')
%     title(ha(2),'Unpaired')
%     for i = 1:length(ha)
%         set(ha(i),'YTick',[],'XTick',1:3600:7201, 'XTickLabel',{'0','1','2'})
%         xlabel(ha(i),'Time (hrs)')
%     end
    
% - - POPULATION AVERAGES - - - - 
    figure('Position',positionfig(1000,400));
    ha = tight_subplot(1,2,0.1 ,0.15);
    plot((0:7200)/3600, [mean(erk1)',mean(erk2)'],'Parent',ha(1),'LineWidth',2)
    plot((0:7200)/3600, [median(erk1)',median(erk2)'],'Parent',ha(2),'LineWidth',2)
    for j = 1:2
        xlabel(ha(j),'Time (hrs)')
        ylabel(ha(j),'Phosphorylated ERK')
        set(ha(j),'ColorOrder',cell2mat(colors.peacock'))
        legend(ha(j),{'paired ERK/MEK', 'unpaired'},'Location','southeast')
    end

% - - - - - (MAX) HISTOGRAMS - - - - -
    graph_lim = prctile([erk1(:);erk2(:)],[2 99.9]);
    erk1(erk1>graph_lim(2)) = graph_lim(2);
    erk2(erk2>graph_lim(2)) = graph_lim(2);
    figure,histogram(max(erk1,[],2),linspace(0,graph_lim(2),32))
    hold on; histogram(max(erk2,[],2),linspace(0,graph_lim(2),32))
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

