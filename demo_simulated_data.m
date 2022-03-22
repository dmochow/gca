% 03.22.22
% using the VAR system in Example 1 from:
%
% Stokes, P. A., & Purdon, P. L. (2017). 
% A study of problems encountered in Granger causality analysis from a neuroscience perspective. 
% Proceedings of the national academy of sciences, 114(34), E7063-E7072.
% 
% run GCA and determine whether the underlying structure is recovered 

clear all; close all; clc
addpath(genpath('.'));
rng(1);

%% figure options
LOAD_DATA=1;
nRuns=100;
fs=14; % font size in plots
colors=cbrewer('qual','Set2',8);
gray=[0.7 0.7 0.7];

%% generate the VAR process
opts = setDefaultOptions;

if LOAD_DATA
    % load precomputed data here
    loadMatFilename='rez_simulated_var_2pair_15-Mar-2022.mat';
    load(loadMatFilename);
else
    % data will be generated and then plotted
    for n=1:nRuns
        tic
        fprintf('Currently simulating run %d of %d ... \n',n,nRuns);
        rng('shuffle'); % randomize across runs -- mostly this changes the mixing matrix
        stats(n,1) = runSimulatedVar(opts);
        toc
    end
end

%% assign variables to make figure;
figStats=stats(nRuns); % pick which one you want to use for figure
S=figStats.S;
X=figStats.X;
yhat=figStats.yhat;
zhat=figStats.zhat;
G_src=figStats.G_src;
G_gca=figStats.G_gca;
G_sensor=figStats.G_sensor;
N=size(figStats.X,1);
nSensorPairs=size(figStats.G_sensor,1);
G_src_2D=figStats.G_src_2D;
G_sensor_2D=figStats.G_sensor_2D;
G_gca_2D=figStats.G_gca_2D;

%% for each run, figure out the performance, keeping in mind the ambiguity in the order of the recovered pairs
isSwapped=nan(nRuns,1);
for n=1:nRuns
   
    % 1-->2, 2-->3
    % scenario 1: y_1~s_1, z_1~s_2, y_2~s_2, z_2~s_3
    % scenario 2: y_2~s_1, z_2~s_2, y_1~s_2, z_1~s_3
    
    testRho1=corrcoef12(stats(n).yhat(:,1),stats(n).S(1,:)) ; % ordered
    testRho2=corrcoef12(stats(n).yhat(:,2),stats(n).S(1,:)); % swapped
    
    if abs(testRho1) > abs(testRho2) % scenario 1: ordered
    
        isSwapped(n)=0; 
        rho_y1s1(n,1) = corrcoef12(stats(n).yhat(:,1),stats(n).S(1,:));
        rho_z1s2(n,1) = corrcoef12(stats(n).zhat(:,1),stats(n).S(2,:));
        rho_y2s2(n,1) = corrcoef12(stats(n).yhat(:,2),stats(n).S(2,:));
        rho_z2s3(n,1) = corrcoef12(stats(n).zhat(:,2),stats(n).S(3,:));
        
        % estimated
        R=cov(stats(n).X);
        A1=R*stats(n).What(:,1) / (stats(n).What(:,1).'*R*stats(n).What(:,1));
        A2=R*stats(n).What(:,2) / (stats(n).What(:,2).'*R*stats(n).What(:,2));
        A3=R*stats(n).Vhat(:,2) / (stats(n).Vhat(:,2).'*R*stats(n).Vhat(:,2));
        Aest = [A1 A2 A3];
        Aest = Aest./ repmat(sqrt(sum(Aest.^2)) , size(Aest,1) , 1 ); % normalize
        Aest= abs (Aest); % rectify
                
    else % scenario 2: swapped
        
        isSwapped(n)=1; 
        rho_y1s1(n,1) = corrcoef12(stats(n).yhat(:,2),stats(n).S(1,:));
        rho_z1s2(n,1) = corrcoef12(stats(n).zhat(:,2),stats(n).S(2,:));
        rho_y2s2(n,1) = corrcoef12(stats(n).yhat(:,1),stats(n).S(2,:));
        rho_z2s3(n,1) = corrcoef12(stats(n).zhat(:,1),stats(n).S(3,:));
        
        % estimated
        R=cov(stats(n).X);
        A1=R*stats(n).What(:,2) / (stats(n).What(:,2).'*R*stats(n).What(:,2));
        A2=R*stats(n).What(:,1) / (stats(n).What(:,1).'*R*stats(n).What(:,1));
        A3=R*stats(n).Vhat(:,1) / (stats(n).Vhat(:,1).'*R*stats(n).Vhat(:,1));
        Aest = [A1 A2 A3];
        Aest = Aest./ repmat(sqrt(sum(Aest.^2)) , size(Aest,1) , 1 ); % normalize
        Aest= abs (Aest); % rectify
        
    end
    
    Atrue=stats(n).A;
    Atrue = Atrue ./ repmat(sqrt(sum(Atrue.^2)) , size(Atrue,1) , 1 ); % normalize
    Atrue = abs (Atrue ); % rectify
    
    rhoA(n,1)=corrcoef12(Atrue(:),Aest(:));
    
end

% measure r^2 between source and recoverd signals
r2_y1s1=rho_y1s1.^2;
r2_z1s2=rho_z1s2.^2;
r2_y2s2=rho_y2s2.^2;
r2_z2s3=rho_z2s3.^2; 

% average r^2 values across runs
[sems_r2_11,mus_r2_11]=nansem(r2_y1s1,1);
[sems_r2_12,mus_r2_12]=nansem(r2_z1s2,1);
[sems_r2_22,mus_r2_22]=nansem(r2_y2s2,1);
[sems_r2_23,mus_r2_23]=nansem(r2_z2s3,1);

[sems_Gsensor,mus_Gsensor]=nansem(cat(2,stats.G_sensor),2);
[sems_Gsrc,mus_Gsrc]=nansem(cat(2,stats.G_src),2);
[sems_Ggca,mus_Ggca]=nansem(cat(2,stats.G_gca),2);

[sems_rhoA,mus_rhoA]=nansem(rhoA,1);
[sems_rhoA2,mus_rhoA2]=nansem(rhoA.^2,1);

[maxval,maxind]=max(mus_Gsensor);
gsrc=cat(2,stats.G_src); gsrc=gsrc([1 3],:);
ggca=cat(2,stats.G_gca); 
gobs=cat(2,stats.G_sensor);

% measure difference between gca G and observed G
pval1=signrank(ggca(1,:),gobs(maxind,:))
pval2=signrank(ggca(2,:),gobs(maxind,:))

% measure difference between src G and observed G
pval3=signrank(gsrc(1,:),gobs(maxind,:))
pval4=signrank(gsrc(2,:),gobs(maxind,:))

% 12.23.21
% mixing matrices
% ground-truth
Atrue=figStats.A;
Atrue = Atrue ./ repmat(sqrt(sum(Atrue.^2)) , size(Atrue,1) , 1 ); % normalize
Atrue = abs (Atrue ); % rectify

% estimated
R=cov(figStats.X);
A1=R*figStats.What(:,1) / (figStats.What(:,1).'*R*figStats.What(:,1));
A2=R*figStats.What(:,2) / (figStats.What(:,2).'*R*figStats.What(:,2));
A3=R*figStats.Vhat(:,2) / (figStats.Vhat(:,2).'*R*figStats.Vhat(:,2));
Aest = [A1 A2 A3];
Aest = Aest./ repmat(sqrt(sum(Aest.^2)) , size(Aest,1) , 1 ); % normalize
Aest= abs (Aest); % rectify

%% latent signals
lw=1.5; 
hf=figure;
hs=subplot(221); hold on
St=S.';
xLimShow=[1900 2100];
plot(St(:,3),'Color',colors(1,:),'LineWidth',lw);
plot(St(:,2)+ 5*ones(N,1) ,'Color',colors(2,:) ,'LineWidth',lw);
plot(St(:,1)+ 10*ones(N,1) ,'Color',colors(3,:),'LineWidth',lw);
xlim(xLimShow);
axis off

%% observed signals
hf=figure;
hs=subplot(221); hold on
plot(X(:,4) ,  'Color',gray,'LineWidth',lw);
plot(X(:,3)+ 5*ones(N,1), 'Color',gray,'LineWidth',lw);
plot(X(:,2)+ 10*ones(N,1), 'Color',gray,'LineWidth',lw);
plot(X(:,1)+ 15*ones(N,1), 'Color',gray,'LineWidth',lw);
xlim(xLimShow);
axis off

%% recovered signals
hf=figure;
hs=subplot(221); hold on
plot(-zhat(:,2)  ,'Color',colors(1,:) ,'LineWidth',lw);
plot(yhat(:,2) + 2.5*ones(N,1) ,'Color',colors(2,:) ,'LineWidth',lw);
plot(-zhat(:,1) + 5*ones(N,1) ,'Color',colors(2,:),'LineWidth',lw);
plot(-yhat(:,1) + 7.5*ones(N,1) ,'Color',colors(3,:),'LineWidth',lw);
xlim(xLimShow);
axis off

%% G_src
hf=figure;
hs=subplot(221); %hold on
imagesc(G_src_2D); 
axis square; axis off; 
caxis([0 0.15])
colormap hot

%% G_sensor
hf=figure;
hs=subplot(221); %hold on
imagesc(G_sensor_2D); 
axis square; axis off; 
caxis([0 0.15])
colormap hot

%% G_gca
hf=figure;
hs=subplot(221); %hold on
imagesc(G_gca_2D); 
axis square; axis off; 
caxis([0 0.15])
colormap hot

%% G_gca
hf=figure;
hs=subplot(221); %hold on
himg=imagesc(G_gca_2D); 
axis square; axis off; 
hcb=colorbar('eastoutside');
title(hcb,'$\mathcal{G}$','Interpreter','Latex');
caxis([0 0.15])
delete(himg)
colormap hot

%% true forward model
hf=figure;
hs=subplot(221); %hold on
imagesc(Atrue); axis off; colormap cool; axis square
caxis([0 1])

%% estimated forward model
hf=figure;
hs=subplot(221); %hold on
imagesc(Aest); axis off; colormap cool; axis square
caxis([0 1])

%% forward model colorbar
hf=figure;
hs=subplot(221); %hold on
himg=imagesc(Aest); axis off; colormap cool; axis square
hcb=colorbar('eastoutside');
delete(himg)
caxis([0 1])

%% scatter s1 y1
ms=8;
hf=figure;
hs=subplot(221); %hold on
scatter(S(1,:),-yhat(:,1),'filled','MarkerFaceColor',[0.7 0.7 0.7],'MarkerFaceAlpha',0.5,'SizeData',ms);
xlim([-4 4]);
ylim([-4 4]);
set(gca,'Xtick',[-4 0 4]);
set(gca,'ytick',[-4 0 4]);

%% scatter s2 y2
ms=8;
hf=figure;
hs=subplot(221); %hold on
scatter(S(2,:),yhat(:,2),'filled','MarkerFaceColor',[0.7 0.7 0.7],'MarkerFaceAlpha',0.5,'SizeData',ms);
xlim([-4 4]);
ylim([-4 4]);
set(gca,'Xtick',[-4 0 4]);
set(gca,'ytick',[-4 0 4]);

%% scatter s2 z1
ms=8;
hf=figure;
hs=subplot(221); %hold on
scatter(S(2,:),-zhat(:,1),'filled','MarkerFaceColor',[0.7 0.7 0.7],'MarkerFaceAlpha',0.5,'SizeData',ms);
xlim([-4 4]);
ylim([-4 4]);
set(gca,'Xtick',[-4 0 4]);
set(gca,'ytick',[-4 0 4]);

%% scatter s3 z2
ms=8;
hf=figure;
hs=subplot(221); %hold on
scatter(S(3,:),zhat(:,2),'filled','MarkerFaceColor',[0.7 0.7 0.7],'MarkerFaceAlpha',0.5,'SizeData',ms);
xlim([-4 4]);
ylim([-4 4]);
set(gca,'Xtick',[-4 0 4]);
set(gca,'ytick',[-4 0 4]);

