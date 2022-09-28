%compute performance

clc
clear all
close all


%% ALL
load('quality_LPCSI.mat')
load('MOS.mat')
load('weights2.mat');

avg_MOS=ones(312,1);
avg_wl2=ones(312,1);
avg_wh2=ones(312,1);
avg_l2=ones(312,1);
avg_h2=ones(312,1);
avg_l2_ref=ones(312,1);
avg_h2_ref=ones(312,1);
avg_wl1=ones(312,1);
avg_wh1=ones(312,1);
avg_l1=ones(312,1);
avg_h1=ones(312,1);
avg_l1_ref=ones(312,1);
avg_h1_ref=ones(312,1);
avg_msssim1=ones(312,1);
avg_kld1=ones(312,1);
avg_msssim2=ones(312,1);
avg_kld2=ones(312,1);
avg_msssim3=ones(312,1);
avg_kld3=ones(312,1);
cnt=1;

for id1=1:length(MOS)
    for id2=1:length(wl2)    
       a = img_name{id2,1};
       if strfind(char(name{id1,1}), a)==1
           avg_MOS(cnt)=MOS(id1);

           avg_wl2(cnt)=wl2(id2);
           avg_wh2(cnt)=wh2(id2);
           avg_wl1(cnt)=wl1(id2);
           avg_wh1(cnt)=wh1(id2);

           avg_l2(cnt)=l2(id2);
           avg_h2(cnt)=h2(id2);
           avg_l1(cnt)=l1(id2);
           avg_h1(cnt)=h1(id2);
           
           avg_l2_ref(cnt)=l2_ref(id2);
           avg_h2_ref(cnt)=h2_ref(id2);
           avg_l1_ref(cnt)=l1_ref(id2);
           avg_h1_ref(cnt)=h1_ref(id2);

           avg_msssim1(cnt)=real(msssim_value1(id2));
           avg_kld1(cnt)=kld_value1(id2);           
           avg_msssim2(cnt)=msssim_value2(id2);
           avg_kld2(cnt)=kld_value2(id2);
           avg_msssim3(cnt)=msssim_value3(id2);
           avg_kld3(cnt)=kld_value3(id2);
           
           id1=id1+1;
           cnt=cnt+1;
           break;
       else
           id2=id2+1;
       end
    end
end

%initialize the parameters used by the nonlinear fitting function
beta(1) = 10;
beta(2) = 0;
beta(3) = mean(avg_msssim1);
beta(4) = 0.1;
beta(5) = 0.1;

%fitting a curve using the data
[bayta ehat,J] = nlinfit(avg_msssim1,avg_MOS,@logistic,beta);
%given a ssim value, predict the correspoing mos (ypre) using the fitted curve
[avg_msssim1_map junk] = nlpredci(@logistic,avg_msssim1,bayta,ehat,J);

%initialize the parameters used by the nonlinear fitting function
beta(1) = 10;
beta(2) = 0;
beta(3) = mean(avg_kld1);
beta(4) = 0.1;
beta(5) = 0.1;

%fitting a curve using the data
[bayta ehat,J] = nlinfit(avg_kld1,avg_MOS,@logistic,beta);
%given a ssim value, predict the correspoing mos (ypre) using the fitted curve
[avg_kld1_map junk] = nlpredci(@logistic,avg_kld1,bayta,ehat,J);

%initialize the parameters used by the nonlinear fitting function
beta(1) = 10;
beta(2) = 0;
beta(3) = mean(avg_msssim2);
beta(4) = 0.1;
beta(5) = 0.1;

%fitting a curve using the data
[bayta ehat,J] = nlinfit(avg_msssim2,avg_MOS,@logistic,beta);
%given a ssim value, predict the correspoing mos (ypre) using the fitted curve
[avg_msssim2_map junk] = nlpredci(@logistic,avg_msssim2,bayta,ehat,J);

%initialize the parameters used by the nonlinear fitting function
beta(1) = 10;
beta(2) = 0;
beta(3) = mean(avg_kld2);
beta(4) = 0.1;
beta(5) = 0.1;

%fitting a curve using the data
[bayta ehat,J] = nlinfit(avg_kld2,avg_MOS,@logistic,beta);
%given a ssim value, predict the correspoing mos (ypre) using the fitted curve
[avg_kld2_map junk] = nlpredci(@logistic,avg_kld2,bayta,ehat,J);

%initialize the parameters used by the nonlinear fitting function
beta(1) = 10;
beta(2) = 0;
beta(3) = mean(avg_msssim3);
beta(4) = 0.1;
beta(5) = 0.1;

%fitting a curve using the data
[bayta ehat,J] = nlinfit(avg_msssim3,avg_MOS,@logistic,beta);
%given a ssim value, predict the correspoing mos (ypre) using the fitted curve
[avg_msssim3_map junk] = nlpredci(@logistic,avg_msssim3,bayta,ehat,J);

%initialize the parameters used by the nonlinear fitting function
beta(1) = 10;
beta(2) = 0;
beta(3) = mean(avg_kld3);
beta(4) = 0.1;
beta(5) = 0.1;

%fitting a curve using the data
[bayta ehat,J] = nlinfit(avg_kld3,avg_MOS,@logistic,beta);
%given a ssim value, predict the correspoing mos (ypre) using the fitted curve
[avg_kld3_map junk] = nlpredci(@logistic,avg_kld3,bayta,ehat,J);

% avg_l2_ref = normalization(avg_l2_ref, min(avg_wh2), max(avg_wh2), max(avg_l2_ref), min(avg_l2_ref));
% avg_l2 = normalization(avg_l2, min(avg_wh2), max(avg_wh2), max(avg_l2), min(avg_l2));
% avg_wl2 = avg_l2./avg_l2_ref;

% avg_wl2 = normalization(avg_wl2, min(avg_wh2), max(avg_wh2), max(avg_wl2), min(avg_wl2));
% avg_wh2 = normalization(avg_wh2, min(avg_wl2), max(avg_wl2), max(avg_wh2), min(avg_wh2));

% avg_msssim2 = normalization(avg_msssim2, 0, 1, max(avg_msssim2), min(avg_msssim2));
% avg_kld2 = normalization(avg_kld2, 0, 1, max(avg_kld2), min(avg_kld2));
% 
% avg_wl2 = normalization(avg_wl2, 0, 1, max(avg_wl2), min(avg_wl2));
% avg_wh2 = normalization(avg_wh2, 0, 1, max(avg_wh2), min(avg_wh2));

% avg_score2 = avg_wl2.^(1).*avg_msssim2;
% avg_score2 = avg_wh2.^(1).*avg_kld2;

% avg_wl2 = avg_wl2./(avg_wl2+avg_wh2);
% avg_wh2 = avg_wh2./(avg_wl2+avg_wh2);

avg_score2 = (avg_wl2.^(2).*avg_msssim2_map)+(avg_wh2.^(2).*avg_kld2_map);

% weight1=weight2./(weight1+weight2);
% weight2=weight1./(weight1+weight2);

avg_score2_mean = 0.5.*avg_msssim2_map+0.5.*avg_kld2_map;

x = [avg_msssim1_map avg_msssim2_map avg_msssim3_map];
y = avg_MOS;
coef = x\y
avg_msssim_mean = coef(1)*avg_msssim1_map+coef(2)*avg_msssim2_map+coef(3)*avg_msssim3_map;
x = [avg_kld1_map avg_kld2_map avg_kld3_map];
y = avg_MOS;
coef = x\y
avg_kld_mean = coef(1)*avg_kld1_map+coef(2)*avg_kld2_map+coef(3)*avg_kld3_map;

%initialize the parameters used by the nonlinear fitting function
beta(1) = 10;
beta(2) = 0;
beta(3) = mean(avg_msssim_mean);
beta(4) = 0.1;
beta(5) = 0.1;

%fitting a curve using the data
[bayta ehat,J] = nlinfit(avg_msssim_mean,avg_MOS,@logistic,beta);
%given a ssim value, predict the correspoing mos (ypre) using the fitted curve
[avg_msssim_mean_map junk] = nlpredci(@logistic,avg_msssim_mean,bayta,ehat,J);

%initialize the parameters used by the nonlinear fitting function
beta(1) = 10;
beta(2) = 0;
beta(3) = mean(avg_kld_mean);
beta(4) = 0.1;
beta(5) = 0.1;

%fitting a curve using the data
[bayta ehat,J] = nlinfit(avg_kld_mean,avg_MOS,@logistic,beta);
%given a ssim value, predict the correspoing mos (ypre) using the fitted curve
[avg_kld_mean_map junk] = nlpredci(@logistic,avg_kld_mean,bayta,ehat,J);

avg_score_mean = 0.5.*avg_msssim_mean_map+0.5.*avg_kld_mean_map;

avg_score_weight = (weight1.*(avg_msssim_mean_map+bias1))+(weight2.*(avg_kld_mean_map+bias2));

%initialize the parameters used by the nonlinear fitting function
beta(1) = 10;
beta(2) = 0;
beta(3) = mean(avg_score_mean);
beta(4) = 0.1;
beta(5) = 0.1;

%fitting a curve using the data
[bayta ehat,J] = nlinfit(avg_score_mean,avg_MOS,@logistic,beta);
%given a ssim value, predict the correspoing mos (ypre) using the fitted curve
[avg_score_mean_map junk] = nlpredci(@logistic,avg_score_mean,bayta,ehat,J);

%initialize the parameters used by the nonlinear fitting function
beta(1) = 10;
beta(2) = 0;
beta(3) = mean(avg_score_weight);
beta(4) = 0.1;
beta(5) = 0.1;

%fitting a curve using the data
[bayta ehat,J] = nlinfit(avg_score_weight,avg_MOS,@logistic,beta);
%given a ssim value, predict the correspoing mos (ypre) using the fitted curve
[avg_score_weight_map junk] = nlpredci(@logistic,avg_score_weight,bayta,ehat,J);


[srocc,krocc,plcc,rmse]=cal_corr(avg_MOS, avg_score_weight_map)

[SROCC, KROCC, CC, RMSE] = resultevaluation(avg_score_weight_map, avg_MOS)

% [CC, RMSE, SROCC, KROCC] = resultevaluation(avg_score,avg_MOS)
% [CC, RMSE, SROCC, KROCC] = resultevaluation(avg_score2,avg_MOS)

% [SROCC,KROCC,PLCC,RMSE]=nonlinearfitting(avg_score, avg_MOS)

save('draw_scatter.mat','avg_msssim_mean_map','avg_kld_mean_map','avg_msssim1_map','avg_kld1_map','avg_msssim2_map','avg_kld2_map','avg_msssim3_map','avg_kld3_map','avg_MOS','avg_wl2','avg_wh2','avg_l2','avg_h2','avg_l2_ref','avg_h2_ref',...
'avg_wl1','avg_wh1','avg_l1','avg_h1','avg_l1_ref','avg_h1_ref','name');


