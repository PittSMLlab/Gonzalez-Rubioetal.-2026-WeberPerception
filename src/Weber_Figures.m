% Weber_Figures.m
% Reproduces all figures from: "Weber's Law in walking: sensory scaling
% is observed in multi-sensory, dynamic tasks"
% Gonzalez-Rubio, Iturralde, & Torres-Oviedo
%
% Figure 1: Experimental protocol schematic (not generated here)
% Figure 2: Psychometric functions + JND box plots (Weber fraction - Relative & absolute)
% Figure 3: Sensitivity comparison (beta2 interaction + paired JND bars + difference JND bars)
% Figure 4: DDM reaction time fits + DDM-derived JND box plots
 
clear all;
close all;
clc;
 
% Make sure to have the code folder in your path, and that matlab current
% folder is where the data is located.
 
%% Data loading
 
% Load all participant data from datalogs
% Data=loadAllDataIntoTable;
Data=readtable('Data\Data.csv');
 
% Compute Weber Fraction and speed-related columns
% Speed mapping: 0=reference(1050), 1=slow(700), 2=comfortable(1400), 3=fast(1750)
speed_values=containers.Map({0,1,2,3},[1050,700,1400,1750]);
 
Data.meanSpeed=nan(size(Data,1),1);
for s=keys(speed_values)
    this_speed=s{1};
    idx=Data.speed==this_speed;
    Data.WF(idx)=Data.pertSize(idx)./speed_values(this_speed);
    Data.meanSpeed(idx)=speed_values(this_speed);
    Data.invSpeed(idx)=1./speed_values(this_speed);
end
Data.pertSize_normalized=Data.pertSize./1050;
 
% Remove no-response trials
Data=Data(~Data.noResponse,:);

%% Define Figure colors
 
colorRef=[0 0.4470 0.7410];        % Blue - reference speed (1.05 m/s)
colorSlow=[0.9290 0.6940 0.1250];  % Yellow - slow speed (0.7 m/s)
colorComf=[0.85 0.325 0.098];      % Orange - comfortable speed (1.4 m/s)
colorFast=[0.494 0.184 0.556];     % Purple - fast speed (1.75 m/s)
colorList=[colorSlow;colorComf;colorFast];
 
%% GLME fits
% Model 1: Absolute perturbation size (pertSize) as predictor
% Model 2: Weber Fraction (WF = pertSize/meanSpeed) as predictor
% Both models include speed as an interaction term
 
groups={'G1','G2','G3'}; % G1=Slow, G2=Comfortable, G3=Fast
 
for model=1:2
 
    if model==1
        frml='leftResponse~1+pertSize+pertSize:speed+(pertSize*speed-speed|subID)';
    elseif model==2
        frml='leftResponse~1+WF+WF:speed+(WF*speed-speed|subID)';
    end
 
    X=Data;
    X.speed=categorical(X.speed);
    X.subID=categorical(X.subID);
 
    mm{model}=fitglme(X,frml,'Distribution','binomial','Link','logit','FitMethod','Laplace');
 
    % Extract fixed and random effects
    [~,~,fe]=fixedEffects(mm{model});
    [~,~,re]=randomEffects(mm{model});
 
    % Organize random effects by coefficient name
    reEstimate=[];
    se_re=[];
    coef=mm{model}.CoefficientNames;
    names={};
    for i=1:length(coef)
        if strcmp(coef{i},'(Intercept)')
            names{i}='b0';
        elseif strcmp(coef{i},'pertSize') || strcmp(coef{i},'WF')
            names{i}='b1';
        elseif strcmp(coef{i},'pertSize:speed_1') || strcmp(coef{i},'speed_1:pertSize') || strcmp(coef{i},'speed_1:WF') || strcmp(coef{i},'WF:speed_1')
            names{i}='b2';
        elseif strcmp(coef{i},'pertSize:speed_2') || strcmp(coef{i},'speed_2:pertSize') || strcmp(coef{i},'speed_2:WF') || strcmp(coef{i},'WF:speed_2')
            names{i}='b3';
        elseif strcmp(coef{i},'pertSize:speed_3') || strcmp(coef{i},'speed_3:pertSize') || strcmp(coef{i},'speed_3:WF') || strcmp(coef{i},'WF:speed_3')
            names{i}='b4';
        end
        idx=strcmp(re.Name,coef{i});
        reEstimate=[reEstimate,re(idx,'Estimate').Estimate];
        se_re=[se_re,re(idx,'SEPred').SEPred];
    end
 
    % Individual predicted coefficients: fixed + random effects
    allBetas=fe.Estimate'+reEstimate;
    betas{model}=table(allBetas(:,1),allBetas(:,2),allBetas(:,3),allBetas(:,4),allBetas(:,5),'VariableNames',names);
 
    % Compute individual JNDs from model parameters: JND = ln(3) / slope
    % Reference uses b1; each testing condition uses b1+b2/b3/b4
    JND_All{model}=table(...
        log(3)*(betas{model}.b1).^(-1),...
        log(3)*(betas{model}.b1+betas{model}.b2).^(-1),...
        log(3)*(betas{model}.b1+betas{model}.b3).^(-1),...
        log(3)*(betas{model}.b1+betas{model}.b4).^(-1),...
        'VariableNames',{'Reference','Slow','Comfortable','Fast'});
end
 
%% JNDs by groups
% Participant indices: G2(Comfortable)=1:13, G1(Slow)=14:26, G3(Fast)=27:39
 
for model=1:2
    for g=1:length(groups)
        if strcmp(groups{g},'G1')
            JND_ref=JND_All{model}.Reference(14:26,:);
            JND_test=JND_All{model}.Slow(14:26,:);
        elseif strcmp(groups{g},'G2')
            JND_ref=JND_All{model}.Reference(1:13,:);
            JND_test=JND_All{model}.Comfortable(1:13,:);
        elseif strcmp(groups{g},'G3')
            JND_ref=JND_All{model}.Reference(27:end,:);
            JND_test=JND_All{model}.Fast(27:end,:);
        end
        JND_Choices{model}{g}=table(JND_ref,JND_test,'VariableNames',{'Reference','Testing'});
    end
end
 
%% LOAD DDM-DERIVED JNDs
% JNDs estimated from reaction times using drift-diffusion model fits
% (generated by Python code in DDM/DDM RT fit)
 
csvFiles=dir(fullfile('Data\JNDs_RT','*.csv'));
fileNames={csvFiles.name};
countPS=1;
countWF=1;
for i=1:length(fileNames)
    if contains(fileNames{i},'pertSize')
        JND_RT{1}{countPS}=readtable(['Data\JNDs_RT\' fileNames{i}]);
        countPS=countPS+1;
    else
        JND_RT{2}{countWF}=readtable(['Data\JNDs_RT\' fileNames{i}]);
        countWF=countWF+1;
    end
end
 
clearvars -except Data colorRef colorSlow colorComf colorFast colorList JND_Choices JND_All mm JND_RT betas groups
 
%% FIGURE 2
% A) Psychometric functions: relative (top) and absolute (bottom) for each group
% B) JND box plots: Weber fraction (top) and absolute JND (bottom)
 
fh=figure('Units','pixels','InnerPosition',[100 100 4*300 2*300],'Name','Figure 2');
[cmap,unsignedMap]=probeColorMap(23);
 
% Panel A: Psychometric curves
for model=1:2
    for g=1:length(groups)
 
        if strcmp(groups{g},'G1')
            b0=betas{model}.b0(14:26);
            b_ref=betas{model}.b1(14:26);
            b_test=betas{model}.b1(14:26)+betas{model}.b2(14:26);
        elseif strcmp(groups{g},'G2')
            b0=betas{model}.b0(1:13);
            b_ref=betas{model}.b1(1:13);
            b_test=betas{model}.b1(1:13)+betas{model}.b3(1:13);
        elseif strcmp(groups{g},'G3')
            b0=betas{model}.b0(27:end);
            b_ref=betas{model}.b1(27:end);
            b_test=betas{model}.b1(27:end)+betas{model}.b4(27:end);
        end
 
        sp=subplot(2,4,(model-1)*4+g);
        hold on;
 
        xx=[-350:1:350];
        if model==2
            xx=[-350:1:350]./700;
        end
 
        for s=1:length(b0)
            yRef=1./(1+exp(-(b0(s)+b_ref(s).*xx)));
            plot(xx-(-(b0(s)./b_ref(s))),yRef,'Color',colorRef,'LineWidth',1);
 
            yTest=1./(1+exp(-(b0(s)+b_test(s).*xx)));
            plot(xx-(-(b0(s)./b_test(s))),yTest,'Color',colorList(g,:),'LineWidth',1);
        end
 
        if model==1
            xlim([-400 400]);
        else
            xlim([-400 400]./700);
        end
 
        axes(sp);
        grid on;
        ylabel('proportion of left choices');
        xlabel('$\Delta V$','Interpreter','latex');
        if model==2
            xlabel('$\frac{\Delta V}{\bar{V}}$','Interpreter','latex');
        end
        uistack(sp,'bottom');
    end
end
 
% Panel B: JND box plots
% Top: Weber Fraction (%)
subplot(2,4,4);
hold on;
plotIndiv=1;
plotMean=1;
 
Boxplot(JND_Choices{2}{1}.Testing.*100,colorList(1,:),plotIndiv,plotMean,700);
Boxplot([JND_Choices{2}{1}.Reference;JND_Choices{2}{2}.Reference;JND_Choices{2}{3}.Reference].*100,colorRef,plotIndiv,plotMean,1050);
Boxplot(JND_Choices{2}{2}.Testing.*100,colorList(2,:),plotIndiv,plotMean,1400);
Boxplot(JND_Choices{2}{3}.Testing.*100,colorList(3,:),plotIndiv,plotMean,1750);
 
grid on;
xlim([400 2000]);
ylim([0 11]);
ylabel('JND (% of mean walking speed)');
xlabel('$\bar{V}$ (mm/s)','Interpreter','latex');
hold off;
 
% Bottom: Absolute JND (mm/s)
subplot(2,4,8);
hold on;
plotIndiv=1;
plotMean=1;
 
Boxplot(JND_Choices{1}{1}.Testing,colorList(1,:),plotIndiv,plotMean,700);
Boxplot([JND_Choices{1}{1}.Reference;JND_Choices{1}{2}.Reference;JND_Choices{1}{3}.Reference],colorRef,plotIndiv,plotMean,1050);
Boxplot(JND_Choices{1}{2}.Testing,colorList(2,:),plotIndiv,plotMean,1400);
Boxplot(JND_Choices{1}{3}.Testing,colorList(3,:),plotIndiv,plotMean,1750);
 
grid on;
xlim([400 2000]);
ylim([0 140]);
ylabel('JND (mm/s)');
xlabel('$\bar{V}$ (mm/s)','Interpreter','latex');
hold off;
 
set(gcf,'renderer','painters');
 
%% FIGURE 3
% A) Interaction parameter beta2,c with CIs (departure from Weber's Law)
% B) Paired JND bar plots: reference vs testing for each group
% C) Paired difference (Testing - Reference) with CI for each group
 
figure('Units','pixels','InnerPosition',[100 100 4*300 1.5*300],'Name','Figure 3');
 
for g=1:length(groups)
 
    if strcmp(groups{g},'G1')
        b0=betas{2}.b0(14:26);
        b_ref=betas{2}.b1(14:26);
        b_test=betas{2}.b1(14:26)+betas{2}.b2(14:26);
    elseif strcmp(groups{g},'G2')
        b0=betas{2}.b0(1:13);
        b_ref=betas{2}.b1(1:13);
        b_test=betas{2}.b1(1:13)+betas{2}.b3(1:13);
    elseif strcmp(groups{g},'G3')
        b0=betas{2}.b0(27:end);
        b_ref=betas{2}.b1(27:end);
        b_test=betas{2}.b1(27:end)+betas{2}.b4(27:end);
    end
 
    % Panel A: beta2,c interaction parameter
    subplot(1,3,1);
    hold on;
 
    errorbar(1,mm{2}.Coefficients.Estimate(3),mm{2}.Coefficients.Estimate(3)-mm{2}.Coefficients.Lower(3),'k');
    plot(1,mm{2}.Coefficients.Estimate(3),'_','MarkerSize',10,'Color',colorList(1,:),'LineWidth',2);
    errorbar(2,mm{2}.Coefficients.Estimate(4),mm{2}.Coefficients.Estimate(4)-mm{2}.Coefficients.Lower(4),'k');
    plot(2,mm{2}.Coefficients.Estimate(4),'_','MarkerSize',10,'Color',colorList(2,:),'LineWidth',2);
    yline(0,'Color',colorRef,'LineWidth',1);
    errorbar(3,mm{2}.Coefficients.Estimate(5),mm{2}.Coefficients.Estimate(5)-mm{2}.Coefficients.Lower(5),'k');
    plot(3,mm{2}.Coefficients.Estimate(5),'_','MarkerSize',10,'Color',colorList(3,:),'LineWidth',2);
    xlim([0 4]);
    hold off;
 
    xticks([1:3]);
    xticklabels({'Slower','Comfortable','Faster'});
    ylabel('$\bar{\beta}_{2}$','Interpreter','latex');
 
    subplot(1,3,2);
    hold on;
 
    xPos=(g-1)*3+[1 2];
    b=bar(xPos,[mean(JND_Choices{1}{g}.Reference) mean(JND_Choices{1}{g}.Testing)],'FaceColor','flat','EdgeColor','none');
    b.CData=[colorRef;colorList(g,:)];

    % Individual participants with connecting lines
    plot((xPos+[0.25 -0.25])',[JND_Choices{1}{g}.Reference JND_Choices{1}{g}.Testing]','-','Color',0.7*ones(1,3),'LineWidth',0.5);
    plot(xPos(1)+0.25,JND_Choices{1}{g}.Reference,'o','MarkerSize',6,'MarkerFaceColor',0.7*ones(1,3),'MarkerEdgeColor','none');
    plot(xPos(2)-0.25,JND_Choices{1}{g}.Testing,'o','MarkerSize',6,'MarkerFaceColor',0.7*ones(1,3),'MarkerEdgeColor','none');
        
    % % Within-subject CIs (Loftus & Masson, 1994)
    % subjMean=(JND_Choices{1}{g}.Reference+JND_Choices{1}{g}.Testing)./2;
    % grandMean=mean([JND_Choices{1}{g}.Reference;JND_Choices{1}{g}.Testing]);
    % refAdj=JND_Choices{1}{g}.Reference-subjMean+grandMean;
    % testAdj=JND_Choices{1}{g}.Testing-subjMean+grandMean;
    % n=length(JND_Choices{1}{g}.Reference);
    % ciRef=1.96*(std(refAdj)./sqrt(n));
    % ciTest=1.96*(std(testAdj)./sqrt(n));
    % ee=errorbar(xPos,[mean(JND_Choices{1}{g}.Reference),mean(JND_Choices{1}{g}.Testing)],[ciRef ciTest],'k','LineStyle','none','DisplayName','95% CI');
    
    % Standard errors
    n=length(JND_Choices{1}{g}.Reference);
    seRef=std(JND_Choices{1}{g}.Reference)./sqrt(n);
    seTest=std(JND_Choices{1}{g}.Testing)./sqrt(n);
    ee=errorbar(xPos,[mean(JND_Choices{1}{g}.Reference),mean(JND_Choices{1}{g}.Testing)],[seRef seTest],'k','LineStyle','none','DisplayName','SE');

    % Paired t-test with Bonferroni correction
    [~,pJND]=ttest(JND_Choices{1}{g}.Reference,JND_Choices{1}{g}.Testing);
    if pJND<(0.05/3)
        text(xPos(1),115,['p=' num2str(pJND)]);
        plot(xPos,[110 110],'-k','LineWidth',2);
    end
    ylabel('JND (mm/s)');
    hold off;
 
    % Panel C: Paired absolute JND difference (Testing - Reference) with CI
    subplot(1,3,3);
    hold on;
 
    % Compute within-subject difference and its CI
    diff=JND_Choices{1}{g}.Testing-JND_Choices{1}{g}.Reference;
    meanDiff=mean(diff);
    seDiff=std(diff)./sqrt(length(diff));
    % ciDiff=1.96*seDiff;

    % Individual participant differences
    plot(g*ones(length(diff),1)+0.25,diff,'o','MarkerSize',6,'MarkerFaceColor',0.7*ones(1,3),'MarkerEdgeColor','none');
  
    % Plot the difference as a bar with CI
    bDiff=bar(g,meanDiff,'FaceColor',colorList(g,:),'EdgeColor','none');
    % errorbar(g,meanDiff,ciDiff,'k','LineStyle','none');
    errorbar(g,meanDiff,seDiff,'k','LineStyle','none');
 
end
 
% Panel C: add zero line and labels
subplot(1,3,3);
yline(0,'Color','k','LineWidth',1);
xlim([0 4]);
xticks([1:3]);
xticklabels({'Slow','Comfortable','Fast'});
ylabel('$\Delta$JND (mm/s)','Interpreter','latex');
hold off;
 
set(gcf,'renderer','painters');
 
%% FIGURE 4
% A) DDM reaction time fits (generated in Python) + DDM-derived JND bar plots
% B) DDM-derived JND box plots (absolute, mm/s)
 
figure('Units','pixels','InnerPosition',[100 100 2*400 400],'Name','Figure 4');
 
% Panel A: DDM-derived JND bar plots
subplot(1,2,1);
for g=1:length(groups)
    hold on;
 
    xPos=(g-1)*3+[1 2];
    b=bar(xPos,[mean(JND_RT{1}{g}.Reference) mean(JND_RT{1}{g}.Testing)],'FaceColor','flat','EdgeColor','none');
    b.CData=[colorRef;colorList(g,:)];
 
    % Individual participants with connecting lines
    plot((xPos+[0.05 -0.05])',[JND_RT{1}{g}.Reference JND_RT{1}{g}.Testing]','-','Color',0.7*ones(1,3),'LineWidth',0.5);
    plot(xPos(1)+0.05,JND_RT{1}{g}.Reference,'o','MarkerSize',4,'MarkerFaceColor',0.7*ones(1,3),'MarkerEdgeColor','none');
    plot(xPos(2)-0.05,JND_RT{1}{g}.Testing,'o','MarkerSize',4,'MarkerFaceColor',0.7*ones(1,3),'MarkerEdgeColor','none');
 
    % Bring bars to front
    b=bar(xPos,[mean(JND_RT{1}{g}.Reference) mean(JND_RT{1}{g}.Testing)],'FaceColor','flat','EdgeColor','none','FaceAlpha',0.7);
    b.CData=[colorRef;colorList(g,:)];
 
    % Within-subject CIs (Loftus & Masson, 1994)
    subjMean=(JND_RT{1}{g}.Reference+JND_RT{1}{g}.Testing)./2;
    grandMean=mean([JND_RT{1}{g}.Reference;JND_RT{1}{g}.Testing]);
    refAdj=JND_RT{1}{g}.Reference-subjMean+grandMean;
    testAdj=JND_RT{1}{g}.Testing-subjMean+grandMean;
    n=length(JND_RT{1}{g}.Reference);
    ciRef=1.96*(std(refAdj)./sqrt(n));
    ciTest=1.96*(std(testAdj)./sqrt(n));
    ee=errorbar(xPos,[mean(JND_RT{1}{g}.Reference),mean(JND_RT{1}{g}.Testing)],[ciRef ciTest],'k','LineStyle','none','DisplayName','95% CI');
 
    % Paired t-test with Bonferroni correction
    [~,pJND]=ttest(JND_RT{1}{g}.Reference,JND_RT{1}{g}.Testing);
    if pJND<(0.05/3)
        text(xPos(1),115,['p=' num2str(pJND)]);
        plot(xPos,[110 110],'-k','LineWidth',2);
    end
end
hold off;
grid on;
ylabel('JND (mm/s)');
ylim([0 300]);
 
% Panel B: DDM-derived JND box plots
subplot(1,2,2);
hold on;
plotIndiv=1;
plotMean=1;
 
Boxplot(JND_RT{1}{1}.Testing,colorList(1,:),plotIndiv,plotMean,700);
Boxplot([JND_RT{1}{1}.Reference;JND_RT{1}{2}.Reference;JND_RT{1}{3}.Reference],colorRef,plotIndiv,plotMean,1050);
Boxplot(JND_RT{1}{2}.Testing,colorList(2,:),plotIndiv,plotMean,1400);
Boxplot(JND_RT{1}{3}.Testing,colorList(3,:),plotIndiv,plotMean,1750);
 
grid on;
hold off;
xlim([400 2000]);
ylabel('JND from RT (mm/s)');
xlabel('$\bar{V}$ (mm/s)','Interpreter','latex');
 
set(gcf,'renderer','painters');
 
