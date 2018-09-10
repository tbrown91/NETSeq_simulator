load ../Scer/NETseq_geneNames.mat;
load ../compare_paramsAllGenes/NETseqBestFitBinned.mat;
load ../compare_paramsAllGenes/NETseqGeneParamsBinned.mat;

model_idx = 7;

%%
double_initDists = importdata('WT_NETseqDoubleInitParams_MovingStalledBacktracked.out');
double_initDists(:,1:265)=0;
half_initDists = importdata('WT_NETseqHalfInitParams_MovingStalledBacktracked.out');
half_initDists(:,1:265)=0;
double_elongDists = importdata('WT_NETseqDoubleElongParams_MovingStalledBacktracked.out');
double_elongDists(:,1:265)=0;
half_elongDists = importdata('WT_NETseqHalfElongParams_MovingStalledBacktracked.out');
half_elongDists(:,1:265)=0;

triple_initDists = importdata('WT_NETseqTripleInitParams_MovingStalledBacktracked.out');
triple_initDists(:,1:265)=0;
third_initDists = importdata('WT_NETseqThirdInitParams_MovingStalledBacktracked.out');
third_initDists(:,1:265)=0;
triple_elongDists = importdata('WT_NETseqTripleElongParams_MovingStalledBacktracked.out');
triple_elongDists(:,1:265)=0;
third_elongDists = importdata('WT_NETseqThirdElongParams_MovingStalledBacktracked.out');
third_elongDists(:,1:265)=0;

quadruple_initDists = importdata('WT_NETseqQuadrupleInitParams_MovingStalledBacktracked.out');
quadruple_initDists(:,1:265)=0;
quarter_initDists = importdata('WT_NETseqQuarterInitParams_MovingStalledBacktracked.out');
quarter_initDists(:,1:265)=0;
quadruple_elongDists = importdata('WT_NETseqQuadrupleElongParams_MovingStalledBacktracked.out');
quadruple_elongDists(:,1:265)=0;
quarter_elongDists = importdata('WT_NETseqQuarterElongParams_MovingStalledBacktracked.out');
quarter_elongDists(:,1:265)=0;

no_backtrackingDists = importdata('WT_NETseqNoBacktrackingParams_MovingStalledBacktracked.out');
no_backtrackingDists(:,1:265)=0;
no_backtrackReleaseDists = importdata('WT_NETseqNoBacktrackReleaseParams_MovingStalledBacktracked.out');
no_backtrackReleaseDists(:,1:265)=0;

%%

double_stallDists = importdata('WT_NETseqDoubleStallParams_MovingStalledBacktracked.out');
double_stallDists(:,1:265)=0;
half_stallDists = importdata('WT_NETseqHalfStallParams_MovingStalledBacktracked.out');
half_stallDists(:,1:265)=0;
double_backtrackDists = importdata('WT_NETseqDoubleElongParams_MovingStalledBacktracked.out');
double_backtrackDists(:,1:265)=0;
half_backtrackDists = importdata('WT_NETseqHalfElongParams_MovingStalledBacktracked.out');
half_backtrackDists(:,1:265)=0;

triple_stallDists = importdata('WT_NETseqTripleStallParams_MovingStalledBacktracked.out');
triple_stallDists(:,1:265)=0;
third_stallDists = importdata('WT_NETseqThirdStallParams_MovingStalledBacktracked.out');
third_stallDists(:,1:265)=0;
triple_backtrackDists = importdata('WT_NETseqTripleElongParams_MovingStalledBacktracked.out');
triple_backtrackDists(:,1:265)=0;
third_backtrackDists = importdata('WT_NETseqThirdElongParams_MovingStalledBacktracked.out');
third_backtrackDists(:,1:265)=0;

quadruple_stallDists = importdata('WT_NETseqQuadrupleStallParams_MovingStalledBacktracked.out');
quadruple_stallDists(:,1:265)=0;
quarter_stallDists = importdata('WT_NETseqQuarterStallParams_MovingStalledBacktracked.out');
quarter_stallDists(:,1:265)=0;
quadruple_backtrackDists = importdata('WT_NETseqQuadrupleElongParams_MovingStalledBacktracked.out');
quadruple_backtrackDists(:,1:265)=0;
quarter_backtrackDists = importdata('WT_NETseqQuarterElongParams_MovingStalledBacktracked.out');
quarter_backtrackDists(:,1:265)=0;

%% Normalise distributions
double_initBinned = zeros(size(double_initDists,1),100);
double_initDists=double_initDists./sum(double_initDists,2);
for g=1:size(double_initDists,1)
    for i=1:100
        double_initBinned(g,i)=sum(double_initDists(g,(i-1)*10+251:i*10+250));
    end
end
triple_initBinned = zeros(size(triple_initDists,1),100);
triple_initDists=triple_initDists./sum(triple_initDists,2);
for g=1:size(triple_initDists,1)
    for i=1:100
        triple_initBinned(g,i)=sum(triple_initDists(g,(i-1)*10+251:i*10+250));
    end
end
quadruple_initBinned = zeros(size(quadruple_initDists,1),100);
quadruple_initDists=quadruple_initDists./sum(quadruple_initDists,2);
for g=1:size(quadruple_initDists,1)
    for i=1:100
        quadruple_initBinned(g,i)=sum(quadruple_initDists(g,(i-1)*10+251:i*10+250));
    end
end

half_initBinned = zeros(size(half_initDists,1),100);
half_initDists=half_initDists./sum(half_initDists,2);
for g=1:size(half_initDists,1)
    for i=1:100
        half_initBinned(g,i)=sum(half_initDists(g,(i-1)*10+251:i*10+250));
    end
end
third_initBinned = zeros(size(third_initDists,1),100);
third_initDists=third_initDists./sum(third_initDists,2);
for g=1:size(third_initDists,1)
    for i=1:100
        third_initBinned(g,i)=sum(third_initDists(g,(i-1)*10+251:i*10+250));
    end
end
quarter_initBinned = zeros(size(quarter_initDists,1),100);
quarter_initDists=quarter_initDists./sum(quarter_initDists,2);
for g=1:size(quarter_initDists,1)
    for i=1:100
        quarter_initBinned(g,i)=sum(quarter_initDists(g,(i-1)*10+251:i*10+250));
    end
end

%%
double_elongBinned = zeros(size(double_elongDists,1),100);
double_elongDists=double_elongDists./sum(double_elongDists,2);
for g=1:size(double_elongDists,1)
    for i=1:100
        double_elongBinned(g,i)=sum(double_elongDists(g,(i-1)*10+251:i*10+250));
    end
end
triple_elongBinned = zeros(size(triple_elongDists,1),100);
triple_elongDists=triple_elongDists./sum(triple_elongDists,2);
for g=1:size(triple_elongDists,1)
    for i=1:100
        triple_elongBinned(g,i)=sum(triple_elongDists(g,(i-1)*10+251:i*10+250));
    end
end
quadruple_elongBinned = zeros(size(quadruple_elongDists,1),100);
quadruple_elongDists=quadruple_elongDists./sum(quadruple_elongDists,2);
for g=1:size(quadruple_elongDists,1)
    for i=1:100
        quadruple_elongBinned(g,i)=sum(quadruple_elongDists(g,(i-1)*10+251:i*10+250));
    end
end

half_elongBinned = zeros(size(half_elongDists,1),100);
half_elongDists=half_elongDists./sum(half_elongDists,2);
for g=1:size(half_elongDists,1)
    for i=1:100
        half_elongBinned(g,i)=sum(half_elongDists(g,(i-1)*10+251:i*10+250));
    end
end
third_elongBinned = zeros(size(third_elongDists,1),100);
third_elongDists=third_elongDists./sum(third_elongDists,2);
for g=1:size(third_elongDists,1)
    for i=1:100
        third_elongBinned(g,i)=sum(third_elongDists(g,(i-1)*10+251:i*10+250));
    end
end
quarter_elongBinned = zeros(size(quarter_elongDists,1),100);
quarter_elongDists=quarter_elongDists./sum(quarter_elongDists,2);
for g=1:size(quarter_elongDists,1)
    for i=1:100
        quarter_elongBinned(g,i)=sum(quarter_elongDists(g,(i-1)*10+251:i*10+250));
    end
end

%%
double_stallBinned = zeros(size(double_stallDists,1),100);
double_stallDists=double_stallDists./sum(double_stallDists,2);
for g=1:size(double_stallDists,1)
    for i=1:100
        double_stallBinned(g,i)=sum(double_stallDists(g,(i-1)*10+251:i*10+250));
    end
end
triple_stallBinned = zeros(size(triple_stallDists,1),100);
triple_stallDists=triple_stallDists./sum(triple_stallDists,2);
for g=1:size(triple_stallDists,1)
    for i=1:100
        triple_stallBinned(g,i)=sum(triple_stallDists(g,(i-1)*10+251:i*10+250));
    end
end
quadruple_stallBinned = zeros(size(quadruple_stallDists,1),100);
quadruple_stallDists=quadruple_stallDists./sum(quadruple_stallDists,2);
for g=1:size(quadruple_stallDists,1)
    for i=1:100
        quadruple_stallBinned(g,i)=sum(quadruple_stallDists(g,(i-1)*10+251:i*10+250));
    end
end

half_stallBinned = zeros(size(half_stallDists,1),100);
half_stallDists=half_stallDists./sum(half_stallDists,2);
for g=1:size(half_stallDists,1)
    for i=1:100
        half_stallBinned(g,i)=sum(half_stallDists(g,(i-1)*10+251:i*10+250));
    end
end
third_stallBinned = zeros(size(third_stallDists,1),100);
third_stallDists=third_stallDists./sum(third_stallDists,2);
for g=1:size(third_stallDists,1)
    for i=1:100
        third_stallBinned(g,i)=sum(third_stallDists(g,(i-1)*10+251:i*10+250));
    end
end
quarter_stallBinned = zeros(size(quarter_stallDists,1),100);
quarter_stallDists=quarter_stallDists./sum(quarter_stallDists,2);
for g=1:size(quarter_stallDists,1)
    for i=1:100
        quarter_stallBinned(g,i)=sum(quarter_stallDists(g,(i-1)*10+251:i*10+250));
    end
end

%%
double_backtrackBinned = zeros(size(double_backtrackDists,1),100);
double_backtrackDists=double_backtrackDists./sum(double_backtrackDists,2);
for g=1:size(double_backtrackDists,1)
    for i=1:100
        double_backtrackBinned(g,i)=sum(double_backtrackDists(g,(i-1)*10+251:i*10+250));
    end
end
triple_backtrackBinned = zeros(size(triple_backtrackDists,1),100);
triple_backtrackDists=triple_backtrackDists./sum(triple_backtrackDists,2);
for g=1:size(triple_backtrackDists,1)
    for i=1:100
        triple_backtrackBinned(g,i)=sum(triple_backtrackDists(g,(i-1)*10+251:i*10+250));
    end
end
quadruple_backtrackBinned = zeros(size(quadruple_backtrackDists,1),100);
quadruple_backtrackDists=quadruple_backtrackDists./sum(quadruple_backtrackDists,2);
for g=1:size(quadruple_backtrackDists,1)
    for i=1:100
        quadruple_backtrackBinned(g,i)=sum(quadruple_backtrackDists(g,(i-1)*10+251:i*10+250));
    end
end

half_backtrackBinned = zeros(size(half_backtrackDists,1),100);
half_backtrackDists=half_backtrackDists./sum(half_backtrackDists,2);
for g=1:size(half_backtrackDists,1)
    for i=1:100
        half_backtrackBinned(g,i)=sum(half_backtrackDists(g,(i-1)*10+251:i*10+250));
    end
end
third_backtrackBinned = zeros(size(third_backtrackDists,1),100);
third_backtrackDists=third_backtrackDists./sum(third_backtrackDists,2);
for g=1:size(third_backtrackDists,1)
    for i=1:100
        third_backtrackBinned(g,i)=sum(third_backtrackDists(g,(i-1)*10+251:i*10+250));
    end
end
quarter_backtrackBinned = zeros(size(quarter_backtrackDists,1),100);
quarter_backtrackDists=quarter_backtrackDists./sum(quarter_backtrackDists,2);
for g=1:size(quarter_backtrackDists,1)
    for i=1:100
        quarter_backtrackBinned(g,i)=sum(quarter_backtrackDists(g,(i-1)*10+251:i*10+250));
    end
end

%%
no_backtrackingBinned = zeros(size(no_backtrackingDists,1),100);
no_backtrackingDists=no_backtrackingDists./sum(no_backtrackingDists,2);
for g=1:size(no_backtrackingDists,1)
    for i=1:100
        no_backtrackingBinned(g,i)=sum(no_backtrackingDists(g,(i-1)*10+251:i*10+250));
    end
end

no_backtrackReleaseBinned = zeros(size(no_backtrackReleaseDists,1),100);
no_backtrackReleaseDists=no_backtrackReleaseDists./sum(no_backtrackReleaseDists,2);
for g=1:size(no_backtrackReleaseDists,1)
    for i=1:100
        no_backtrackReleaseBinned(g,i)=sum(no_backtrackReleaseDists(g,(i-1)*10+251:i*10+250));
    end
end

%% Make figures of metagenes
% 
% fig=figure('Position',[0 100 1720 1720]);
% subplot(1,2,1);
% plot(mean(NETseqBestFitBinned{model_idx}(:,26:end),1),'k','linewidth',2);
% hold on;
% plot(mean(half_initBinned,1),'color',[0 0.33 0.33],'linewidth',2);
% plot(mean(third_initBinned,1),'color',[0 0.67 0.67],'linewidth',2);
% plot(mean(quarter_initBinned,1),'color',[0 1 1],'linewidth',2);
% legend('Wild-type','\times1/2 initiation rate','\times1/3 initiation rate','\times1/4 initiation rate');
% ylim([0 0.03]);
% ylabel('Mean reads per 10bp per gene');
% ax=gca;
% ax.FontSize=16;
% ax.XTick = [1,50 100];
% ax.XTickLabel={'TSS','+500','+1000'};
% 
% subplot(1,2,2);
% plot(mean(NETseqBestFitBinned{model_idx}(:,26:end),1),'k','linewidth',2);
% hold on;
% plot(mean(double_initBinned,1),'color',[0 0.33 0.33],'linewidth',2);
% plot(mean(triple_initBinned,1),'color',[0 0.67 0.67],'linewidth',2);
% plot(mean(quadruple_initBinned,1),'color',[0 1 1],'linewidth',2);
% legend('Wild-type','\times2 initiation rate','\times3 initiation rate','\times4 initiation rate');
% ylim([0 0.03]);
% ylabel('Mean reads per 10bp per gene');
% ax=gca;
% ax.FontSize=16;
% ax.XTick = [1,50 100];
% ax.XTickLabel={'TSS','+500','+1000'};
% saveas(fig,'Vary_InitiationRate.svg');
% saveas(fig,'Vary_InitiationRate.eps');
% 
% fig=figure('Position',[0 100 1720 600]);
% subplot(1,2,1);
% plot(mean(NETseqBestFitBinned{model_idx}(:,26:end),1),'k','linewidth',2);
% hold on;
% plot(mean(half_elongBinned,1),'color',[0 0.33 0.33],'linewidth',2);
% plot(mean(third_elongBinned,1),'color',[0 0.67 0.67],'linewidth',2);
% plot(mean(quarter_elongBinned,1),'color',[0 1 1],'linewidth',2);
% legend('Wild-type','\times1/2 elongation rate','\times1/3 elongation rate','\times1/4 elongation rate');
% ylim([0 0.03]);
% ylabel('Mean reads per 10bp per gene');
% ax=gca;
% ax.FontSize=16;
% ax.XTick = [1,50 100];
% ax.XTickLabel={'TSS','+500','+1000'};
% 
% subplot(1,2,2);
% plot(mean(NETseqBestFitBinned{model_idx}(:,26:end),1),'k','linewidth',2);
% hold on;
% plot(mean(double_elongBinned,1),'color',[0 0.33 0.33],'linewidth',2);
% plot(mean(triple_elongBinned,1),'color',[0 0.67 0.67],'linewidth',2);
% plot(mean(quadruple_elongBinned,1),'color',[0 1 1],'linewidth',2);
% legend('Wild-type','\times2 elongation rate','\times3 elongation rate','\times4 elongation rate');
% ylim([0 0.03]);
% ylabel('Mean reads per 10bp per gene');
% ax=gca;
% ax.FontSize=16;
% ax.XTick = [1,50 100];
% ax.XTickLabel={'TSS','+500','+1000'};
% saveas(fig,'Vary_ElongationRate.svg');
% saveas(fig,'Vary_ElongationRate.eps');

fig=figure('Position',[0 100 1720 600]);
subplot(1,2,1);
plot(mean(NETseqBestFitBinned{model_idx}(:,26:end),1),'k','linewidth',2);
hold on;
plot(mean(half_stallBinned,1),'color',[0 0.33 0.33],'linewidth',2);
plot(mean(third_stallBinned,1),'color',[0 0.67 0.67],'linewidth',2);
plot(mean(quarter_stallBinned,1),'color',[0 1 1],'linewidth',2);
legend('Wild-type','\times1/2 stall rate','\times1/3 stall rate','\times1/4 stall rate');
ylim([0 0.03]);
ylabel('Mean reads per 10bp per gene');
ax=gca;
ax.FontSize=16;
ax.XTick = [1,50 100];
ax.XTickLabel={'TSS','+500','+1000'};

subplot(1,2,2);
plot(mean(NETseqBestFitBinned{model_idx}(:,26:end),1),'k','linewidth',2);
hold on;
plot(mean(double_stallBinned,1),'color',[0 0.33 0.33],'linewidth',2);
plot(mean(triple_stallBinned,1),'color',[0 0.67 0.67],'linewidth',2);
plot(mean(quadruple_stallBinned,1),'color',[0 1 1],'linewidth',2);
legend('Wild-type','\times2 stall rate','\times3 stall rate','\times4 stall rate');
ylim([0 0.03]);
ylabel('Mean reads per 10bp per gene');
ax=gca;
ax.FontSize=16;
ax.XTick = [1,50 100];
ax.XTickLabel={'TSS','+500','+1000'};
saveas(fig,'Vary_StallRate.svg');
saveas(fig,'Vary_StallRate.eps');

fig=figure('Position',[0 100 1720 600]);
subplot(1,2,1);
plot(mean(NETseqBestFitBinned{model_idx}(:,26:end),1),'k','linewidth',2);
hold on;
plot(mean(half_backtrackBinned,1),'color',[0 0.33 0.33],'linewidth',2);
plot(mean(third_backtrackBinned,1),'color',[0 0.67 0.67],'linewidth',2);
plot(mean(quarter_backtrackBinned,1),'color',[0 1 1],'linewidth',2);
legend('Wild-type','\times1/2 backtrack rate','\times1/3 backtrack rate','\times1/4 backtrack rate');
ylim([0 0.03]);
ylabel('Mean reads per 10bp per gene');
ax=gca;
ax.FontSize=16;
ax.XTick = [1,50 100];
ax.XTickLabel={'TSS','+500','+1000'};

subplot(1,2,2);
plot(mean(NETseqBestFitBinned{model_idx}(:,26:end),1),'k','linewidth',2);
hold on;
plot(mean(double_backtrackBinned,1),'color',[0 0.33 0.33],'linewidth',2);
plot(mean(triple_backtrackBinned,1),'color',[0 0.67 0.67],'linewidth',2);
plot(mean(quadruple_backtrackBinned,1),'color',[0 1 1],'linewidth',2);
legend('Wild-type','\times2 backtrack rate','\times3 backtrack rate','\times4 backtrack rate');
ylim([0 0.03]);
ylabel('Mean reads per 10bp per gene');
ax=gca;
ax.FontSize=16;
ax.XTick = [1,50 100];
ax.XTickLabel={'TSS','+500','+1000'};
saveas(fig,'Vary_BacktrackRate.svg');
saveas(fig,'Vary_BacktrackRate.eps');

%%
% fig=figure('Position',[200 200 640 600]);
% plot(mean(NETseqBestFitBinned{model_idx}(:,26:end),1),'k','linewidth',2);
% hold on;
% plot(mean(no_backtrackingBinned,1),'c','linewidth',2);
% plot(mean(no_backtrackReleaseBinned,1),'b','linewidth',2);
% legend('Wild-type','No backtracking','No backtrack release');
% ylim([0 0.12]);
% ylabel('Mean reads per 10bp per gene');
% title('Backtracking extremes','fontsize',24);
% ax=gca;
% ax.FontSize=16;
% ax.XTick = [1,50 100];
% ax.XTickLabel={'TSS','+500','+1000'};
% 
% saveas(fig,'backtracking_extremes.svg');
% saveas(fig,'backtracking_extremes.eps');
% close(fig);

%% Display processivity of the parameter sets

double_initCounts = importdata('WT_NETseqDoubleInitParams_EventCounts.out');
double_initCounts = double_initCounts.data;
half_initCounts = importdata('WT_NETseqHalfInitParams_EventCounts.out');
half_initCounts = half_initCounts.data;
double_elongCounts = importdata('WT_NETseqDoubleElongParams_EventCounts.out');
double_elongCounts = double_elongCounts.data;
half_elongCounts = importdata('WT_NETseqHalfElongParams_EventCounts.out');
half_elongCounts = half_elongCounts.data;

triple_initCounts = importdata('WT_NETseqTripleInitParams_EventCounts.out');
triple_initCounts = triple_initCounts.data;
third_initCounts = importdata('WT_NETseqThirdInitParams_EventCounts.out');
third_initCounts = third_initCounts.data;
triple_elongCounts = importdata('WT_NETseqTripleElongParams_EventCounts.out');
triple_elongCounts = triple_elongCounts.data;
third_elongCounts = importdata('WT_NETseqThirdElongParams_EventCounts.out');
third_elongCounts = third_elongCounts.data;

quadruple_initCounts = importdata('WT_NETseqQuadrupleInitParams_EventCounts.out');
quadruple_initCounts = quadruple_initCounts.data;
quarter_initCounts = importdata('WT_NETseqQuarterInitParams_EventCounts.out');
quarter_initCounts = quarter_initCounts.data;
quadruple_elongCounts = importdata('WT_NETseqQuadrupleElongParams_EventCounts.out');
quadruple_elongCounts = quadruple_elongCounts.data;
quarter_elongCounts = importdata('WT_NETseqQuarterElongParams_EventCounts.out');
quarter_elongCounts = quarter_elongCounts.data;
%%
double_stallCounts = importdata('WT_NETseqDoubleStallParams_EventCounts.out');
double_stallCounts = double_stallCounts.data;
half_stallCounts = importdata('WT_NETseqHalfStallParams_EventCounts.out');
half_stallCounts = half_stallCounts.data;
double_backtrackCounts = importdata('WT_NETseqDoubleBacktrackParams_EventCounts.out');
double_backtrackCounts = double_backtrackCounts.data;
half_backtrackCounts = importdata('WT_NETseqHalfBacktrackParams_EventCounts.out');
half_backtrackCounts = half_backtrackCounts.data;

triple_stallCounts = importdata('WT_NETseqTripleStallParams_EventCounts.out');
triple_stallCounts = triple_stallCounts.data;
third_stallCounts = importdata('WT_NETseqThirdStallParams_EventCounts.out');
third_stallCounts = third_stallCounts.data;
triple_backtrackCounts = importdata('WT_NETseqTripleBacktrackParams_EventCounts.out');
triple_backtrackCounts = triple_backtrackCounts.data;
third_backtrackCounts = importdata('WT_NETseqThirdBacktrackParams_EventCounts.out');
third_backtrackCounts = third_backtrackCounts.data;

quadruple_stallCounts = importdata('WT_NETseqQuadrupleStallParams_EventCounts.out');
quadruple_stallCounts = quadruple_stallCounts.data;
quarter_stallCounts = importdata('WT_NETseqQuarterStallParams_EventCounts.out');
quarter_stallCounts = quarter_stallCounts.data;
quadruple_backtrackCounts = importdata('WT_NETseqQuadrupleBacktrackParams_EventCounts.out');
quadruple_backtrackCounts = quadruple_backtrackCounts.data;
quarter_backtrackCounts = importdata('WT_NETseqQuarterBacktrackParams_EventCounts.out');
quarter_backtrackCounts = quarter_backtrackCounts.data;
%%

no_BacktrackCounts = importdata('WT_NETseqNoBacktrackingParams_EventCounts.out');
no_BacktrackCounts=no_BacktrackCounts.data;
no_BacktrackReleaseCounts = importdata('WT_NETseqNoBacktrackReleaseParams_EventCounts.out');
no_BacktrackReleaseCounts=no_BacktrackReleaseCounts.data;
%%
proc_vec=nan(size(double_initCounts,1),14);
WT_procVec = 100*NETseqGeneParamsBinned{model_idx}(:,end)./(NETseqGeneParamsBinned{model_idx}(:,end)+NETseqGeneParamsBinned{model_idx}(:,end-1));

proc_vec(:,1)=100*quarter_initCounts(:,3)./(quarter_initCounts(:,2)+quarter_initCounts(:,3));
proc_vec(:,2)=100*third_initCounts(:,3)./(third_initCounts(:,2)+third_initCounts(:,3));
proc_vec(:,3)=100*half_initCounts(:,3)./(half_initCounts(:,2)+half_initCounts(:,3));
proc_vec(:,4)=100*double_initCounts(:,3)./(double_initCounts(:,2)+double_initCounts(:,3));
proc_vec(:,5)=100*triple_initCounts(:,3)./(triple_initCounts(:,2)+triple_initCounts(:,3));
proc_vec(:,6)=100*quadruple_initCounts(:,3)./(quadruple_initCounts(:,2)+quadruple_initCounts(:,3));

proc_vec(:,7)=100*quarter_elongCounts(:,3)./(quarter_elongCounts(:,2)+quarter_elongCounts(:,3));
proc_vec(:,8)=100*third_elongCounts(:,3)./(third_elongCounts(:,2)+third_elongCounts(:,3));
proc_vec(:,9)=100*half_elongCounts(:,3)./(half_elongCounts(:,2)+half_elongCounts(:,3));
proc_vec(:,10)=100*double_elongCounts(:,3)./(double_elongCounts(:,2)+double_elongCounts(:,3));
proc_vec(:,11)=100*triple_elongCounts(:,3)./(triple_elongCounts(:,2)+triple_elongCounts(:,3));
proc_vec(:,12)=100*quadruple_elongCounts(:,3)./(quadruple_elongCounts(:,2)+quadruple_elongCounts(:,3));

proc_vec(:,13)=100*no_BacktrackReleaseCounts(:,3)./(no_BacktrackReleaseCounts(:,2)+no_BacktrackReleaseCounts(:,3));
proc_vec(:,14)=100*no_BacktrackCounts(:,3)./(no_BacktrackCounts(:,2)+no_BacktrackCounts(:,3));
%%
proc_vecStall=nan(size(double_stallCounts,1),6);
proc_vecStall(:,1)=100*quarter_stallCounts(:,3)./(quarter_stallCounts(:,2)+quarter_stallCounts(:,3));
proc_vecStall(:,2)=100*third_stallCounts(:,3)./(third_stallCounts(:,2)+third_stallCounts(:,3));
proc_vecStall(:,3)=100*half_stallCounts(:,3)./(half_stallCounts(:,2)+half_stallCounts(:,3));
proc_vecStall(:,4)=100*double_stallCounts(:,3)./(double_stallCounts(:,2)+double_stallCounts(:,3));
proc_vecStall(:,5)=100*triple_stallCounts(:,3)./(triple_stallCounts(:,2)+triple_stallCounts(:,3));
proc_vecStall(:,6)=100*quadruple_stallCounts(:,3)./(quadruple_stallCounts(:,2)+quadruple_stallCounts(:,3));
proc_vecBacktrack=nan(size(double_backtrackCounts,1),6);
proc_vecBacktrack(:,1)=100*quarter_backtrackCounts(:,3)./(quarter_backtrackCounts(:,2)+quarter_backtrackCounts(:,3));
proc_vecBacktrack(:,2)=100*third_backtrackCounts(:,3)./(third_backtrackCounts(:,2)+third_backtrackCounts(:,3));
proc_vecBacktrack(:,3)=100*half_backtrackCounts(:,3)./(half_backtrackCounts(:,2)+half_backtrackCounts(:,3));
proc_vecBacktrack(:,4)=100*double_backtrackCounts(:,3)./(double_backtrackCounts(:,2)+double_backtrackCounts(:,3));
proc_vecBacktrack(:,5)=100*triple_backtrackCounts(:,3)./(triple_backtrackCounts(:,2)+triple_backtrackCounts(:,3));
proc_vecBacktrack(:,6)=100*quadruple_backtrackCounts(:,3)./(quadruple_backtrackCounts(:,2)+quadruple_backtrackCounts(:,3));

%% Differences in processivity

proc_vecDiffs = proc_vec-WT_procVec;
proc_vecStall = proc_vecStall-WT_procVec;
proc_vecBacktrack = proc_vecBacktrack-WT_procVec;

fig=figure('Position',[1400 100 1200 1800]);
subplot(5,1,1);
violin(proc_vecDiffs(:,1:6),'facecolor',[0 0 1]);
% boxplot(proc_vecDiffs(:,1:6));
ylim([-50 50]);
grid on;
ylabel('Processivity difference');
xlabel('Initiation rate factor');
ax=gca;
ax.XTick=[1:6];
ax.XTickLabel={'\times^1/_4','\times^1/_3','\times^1/_2','\times2','\times3','\times4'};
ax.FontSize=16;
subplot(5,1,2);
violin(proc_vecDiffs(:,7:12),'facecolor',[0 0 1]);
% boxplot(proc_vecDiffs(:,7:12));
ylim([-50 50]);
grid on;
ylabel('Processivity difference');
xlabel('Elongation rate factor');
ax=gca;
ax.XTick=[1:6];
ax.XTickLabel={'\times^1/_4','\times^1/_3','\times^1/_2','\times2','\times3','\times4'};
ax.FontSize=16;

subplot(5,1,3);
violin(proc_vecStall,'facecolor',[0 0 1]);
ylim([-50 50]);
grid on;
ylabel('Processivity difference');
xlabel('Stall rate factor');
ax=gca;
ax.XTick=[1:6];
ax.XTickLabel={'\times^1/_4','\times^1/_3','\times^1/_2','\times2','\times3','\times4'};
ax.FontSize=16;
subplot(5,1,4);
violin(proc_vecBacktrack,'facecolor',[0 0 1]);
% boxplot(proc_vecDiffs(:,7:12));
ylim([-50 50]);
grid on;
ylabel('Processivity difference');
xlabel('Backtrack rate factor');
ax=gca;
ax.XTick=[1:6];
ax.XTickLabel={'\times^1/_4','\times^1/_3','\times^1/_2','\times2','\times3','\times4'};
ax.FontSize=16;

subplot(5,1,5);
violin(proc_vecDiffs(:,13:14),'facecolor',[0 0 1]);
% boxplot(proc_vecDiffs(:,13:14));
ylim([-100 100]);
grid on;
ylabel('Processivity difference');
ax=gca;
ax.XTick=[1,2];
ax.XTickLabel={'No backtrack release','No backtracking'};
ax.FontSize=16;

% saveas(fig,'pairwise_processivity.svg');
% saveas(fig,'pairwise_processivity.eps');
%% Proximal-distal Ratios

init_proxDistRatios = zeros(size(NETseqBestFitBinned{model_idx},1),7);
init_proxDistRatios(:,1) = sum(quarter_initBinned(:,1:25),2);
init_proxDistRatios(:,2) = sum(third_initBinned(:,1:25),2);
init_proxDistRatios(:,3) = sum(half_initBinned(:,1:25),2);
init_proxDistRatios(:,4) = sum(NETseqBestFitBinned{model_idx}(:,26:50),2);
init_proxDistRatios(:,5) = sum(double_initBinned(:,1:25),2);
init_proxDistRatios(:,6) = sum(triple_initBinned(:,1:25),2);
init_proxDistRatios(:,7) = sum(quadruple_initBinned(:,1:25),2);

init_proxDistRatioDiffs = zeros(size(NETseqBestFitBinned{model_idx},1),6);
for i=1:3
    init_proxDistRatioDiffs(:,i)=init_proxDistRatios(:,i)-init_proxDistRatios(:,4);
end
for i=1:3
    init_proxDistRatioDiffs(:,i+3)=init_proxDistRatios(:,i+4)-init_proxDistRatios(:,4);
end

fig = figure('position',[0 100 1720 600]);
boxplot(init_proxDistRatioDiffs)
ylabel('Difference in proximal-distal ratio');
xlabel('Initiation rate factor');
ax=gca;
ax.TickLabelInterpreter='tex';
ax.XTickLabel={'\times^1/_4','\times^1/_3','\times^1/_2','\times2','\times3','\times4'};
ax.FontSize=16;
ylim([-.3 .3])
saveas(fig,'prox_distRatioInit.svg');
saveas(fig,'prox_distRatioInit.eps');

%%
elong_proxDistRatios = zeros(size(NETseqBestFitBinned{model_idx},1),7);
elong_proxDistRatios(:,1) = sum(quarter_elongBinned(:,1:25),2);
elong_proxDistRatios(:,2) = sum(third_elongBinned(:,1:25),2);
elong_proxDistRatios(:,3) = sum(half_elongBinned(:,1:25),2);
elong_proxDistRatios(:,4) = sum(NETseqBestFitBinned{model_idx}(:,26:50),2);
elong_proxDistRatios(:,5) = sum(double_elongBinned(:,1:25),2);
elong_proxDistRatios(:,6) = sum(triple_elongBinned(:,1:25),2);
elong_proxDistRatios(:,7) = sum(quadruple_elongBinned(:,1:25),2);

elong_proxDistRatioDiffs = zeros(size(NETseqBestFitBinned{model_idx},1),6);
for i=1:3
    elong_proxDistRatioDiffs(:,i)=elong_proxDistRatios(:,i)-elong_proxDistRatios(:,4);
end
for i=1:3
    elong_proxDistRatioDiffs(:,i+3)=elong_proxDistRatios(:,i+4)-elong_proxDistRatios(:,4);
end

fig = figure('position',[0 100 1720 600]);
boxplot(elong_proxDistRatioDiffs)
ylabel('Difference in proximal-distal ratio');
xlabel('Elongation rate factor');
ax=gca;
ax.TickLabelInterpreter='tex';
ax.XTickLabel={'\times^1/_4','\times^1/_3','\times^1/_2','\times2','\times3','\times4'};
ax.FontSize=16;
ylim([-.3 .3])
saveas(fig,'prox_distRatioElong.svg');
saveas(fig,'prox_distRatioElong.eps');

%%
stall_proxDistRatios = zeros(size(NETseqBestFitBinned{model_idx},1),7);
stall_proxDistRatios(:,1) = sum(quarter_stallBinned(:,1:25),2);
stall_proxDistRatios(:,2) = sum(third_stallBinned(:,1:25),2);
stall_proxDistRatios(:,3) = sum(half_stallBinned(:,1:25),2);
stall_proxDistRatios(:,4) = sum(NETseqBestFitBinned{model_idx}(:,26:50),2);
stall_proxDistRatios(:,5) = sum(double_stallBinned(:,1:25),2);
stall_proxDistRatios(:,6) = sum(triple_stallBinned(:,1:25),2);
stall_proxDistRatios(:,7) = sum(quadruple_stallBinned(:,1:25),2);

stall_proxDistRatioDiffs = zeros(size(NETseqBestFitBinned{model_idx},1),6);
for i=1:3
    stall_proxDistRatioDiffs(:,i)=stall_proxDistRatios(:,i)-stall_proxDistRatios(:,4);
end
for i=1:3
    stall_proxDistRatioDiffs(:,i+3)=stall_proxDistRatios(:,i+4)-stall_proxDistRatios(:,4);
end

fig = figure('position',[0 100 1720 600]);
boxplot(stall_proxDistRatioDiffs)
ylabel('Difference in proximal-distal ratio');
xlabel('Stalling rate factor');
ax=gca;
ax.TickLabelInterpreter='tex';
ax.XTickLabel={'\times^1/_4','\times^1/_3','\times^1/_2','\times2','\times3','\times4'};
ax.FontSize=16;
ylim([-.3 .3])
saveas(fig,'prox_distRatioStall.svg');
saveas(fig,'prox_distRatioStall.eps');

%%
backtrack_proxDistRatios = zeros(size(NETseqBestFitBinned{model_idx},1),7);
backtrack_proxDistRatios(:,1) = sum(quarter_backtrackBinned(:,1:25),2);
backtrack_proxDistRatios(:,2) = sum(third_backtrackBinned(:,1:25),2);
backtrack_proxDistRatios(:,3) = sum(half_backtrackBinned(:,1:25),2);
backtrack_proxDistRatios(:,4) = sum(NETseqBestFitBinned{model_idx}(:,26:50),2);
backtrack_proxDistRatios(:,5) = sum(double_backtrackBinned(:,1:25),2);
backtrack_proxDistRatios(:,6) = sum(triple_backtrackBinned(:,1:25),2);
backtrack_proxDistRatios(:,7) = sum(quadruple_backtrackBinned(:,1:25),2);

backtrack_proxDistRatioDiffs = zeros(size(NETseqBestFitBinned{model_idx},1),6);
for i=1:3
    backtrack_proxDistRatioDiffs(:,i)=backtrack_proxDistRatios(:,i)-backtrack_proxDistRatios(:,4);
end
for i=1:3
    backtrack_proxDistRatioDiffs(:,i+3)=backtrack_proxDistRatios(:,i+4)-backtrack_proxDistRatios(:,4);
end

fig = figure('position',[0 100 1720 600]);
boxplot(backtrack_proxDistRatioDiffs)
ylabel('Difference in proximal-distal ratio');
xlabel('Backtracking rate factor');
ax=gca;
ax.TickLabelInterpreter='tex';
ax.XTickLabel={'\times^1/_4','\times^1/_3','\times^1/_2','\times2','\times3','\times4'};
ax.FontSize=16;
ylim([-.3 .3])
saveas(fig,'prox_distRatioBacktrack.svg');
saveas(fig,'prox_distRatioBacktrack.eps');