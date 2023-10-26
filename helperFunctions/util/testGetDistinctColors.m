% just runs through the variables for getDistinctColors
%%
colors = getDistinctColors();

figure
hold on
for i = 1:length(colors)
     plot(1:10,0.5*i*(1:10),'-','MarkerEdgeColor',colors{i},'MarkerFaceColor',colors{i},'color',colors{i})
end

%%
colors = getDistinctColors('RANDOM_ORDER',1);

figure
hold on
for i = 1:length(colors)
     plot(1:10,0.5*i*(1:10),'-','MarkerEdgeColor',colors{i},'MarkerFaceColor',colors{i},'color',colors{i})
end

%%
colors = getDistinctColors('SELECT_ORDER',2);

figure
hold on
for i = 1:length(colors)
     plot(1:10,0.5*i*(1:10),'-','MarkerEdgeColor',colors{i},'MarkerFaceColor',colors{i},'color',colors{i})
end