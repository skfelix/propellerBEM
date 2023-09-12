function plotPerf(perf,title,lgd,group,plotVelocity)



% fontBold = 'bold';
fontSize = 12;
lgdSize = 12;
solidLine = '-';
lineWidth = 1;
lineVec = {'-','--',':','-.'};
colorVec = {'k','r','b','g','m','c','#D95319','#7E2F8E'};
markerColorVec = {'k','r','b','g','m','c','#D95319','#7E2F8E'};
markerVec = {'o','s','^','d','v','+','x','*'};
markerQuantity = 10;

maxCpVec = zeros(1,length(perf));
maxCtVec = zeros(1,length(perf));
maxTVec = zeros(1,length(perf));
xMaxVec = zeros(1,length(perf));

ls = lineVec{group};

if plotVelocity
    xVec = 'V';
else
    xVec = 'J';
end

for i = 1:length(perf)
    xMaxVec(i) = max(perf{i}.(xVec));
end
xMax = max(xMaxVec);


sgtitle(title)

%% Power coefficient - Cp
subplot(2,1,2);
hold on
for i = 1:length(perf)
    plot(perf{i}.(xVec),perf{i}.Cp,...
        'LineStyle',ls, ...
        'Color',colorVec{i}, ...
        'Marker',markerVec{i}, ...
        'MarkerEdgeColor', markerColorVec{i}, ...
        'LineWidth',lineWidth, ...
        'MarkerIndices',1:round((length(perf{i}.J)/markerQuantity)):length(perf{i}.J))
    maxCpVec(i) = max(perf{i}.Cp);
end
maxCp = max(maxCpVec)*1.1;

% maxCp = 1.1*Cp1;
if plotVelocity
    xlabel('$V$ [m/s]','Interpreter','Latex','FontSize',fontSize);
else
    xlabel('$J = \frac{V}{nD}$','Interpreter','Latex','FontSize',fontSize);
end
ylabel('$C_P$','Interpreter','Latex','FontSize',14);
grid on; grid minor;
legend(lgd,'Interpreter','Latex','Location','best','FontSize',lgdSize,'AutoUpdate','off')
axis([0 xMax 0 maxCp])

%% Thrust Coefficient - Ct
subplot(2,1,1);
hold on
j = 1;
for i = 1:length(perf)
    plot(perf{i}.(xVec),perf{i}.Ct,...
        'LineStyle',ls, ...
        'Color',colorVec{i}, ...
        'Marker',markerVec{i}, ...
        'MarkerEdgeColor', markerColorVec{i}, ...
        'LineWidth',lineWidth, ...
        'MarkerIndices',1:round((length(perf{i}.J)/markerQuantity)):length(perf{i}.J))
    maxCtVec(i) = max(perf{i}.Ct);
end
if j > 1; plot(V_Cp_max,Ct_Cp_max,'--k','lineWidth',2); end
maxCt = max(maxCtVec)+0.02;
if plotVelocity
    xlabel('$V$ [m/s]','Interpreter','Latex','FontSize',fontSize);
else
    xlabel('$J = \frac{V}{nD}$','Interpreter','Latex','FontSize',fontSize);
end
ylabel('$C_T$','Interpreter','Latex','FontSize',14);
grid on; grid minor;
% legend(lgd,'Interpreter','Latex','Location','best','FontSize',lgdSize)
axis([0 xMax 0 maxCt]);
