clear all
close all


load('draw_scatter.mat')
texture_ratio=avg_h1./avg_h2_ref;
sharpness_ratio=avg_l1./avg_l2_ref;
avg_msssim2_map=avg_msssim_mean_map;
avg_kld2_map=avg_kld_mean_map;

figure(1)
h1 = scatterhist(texture_ratio.^5+sharpness_ratio.^5,avg_MOS-avg_msssim2_map);
figure(2)
h2 = scatterhist(texture_ratio.^5+sharpness_ratio.^5,avg_MOS-avg_kld2_map);


data1 = h1(2).Children(1);
data2 = h2(2).Children(1);
dat1 = data1.Data;
dat2 = data2.Data;

d1 = avg_MOS-avg_msssim2_map;
d2 = avg_MOS-avg_kld2_map;

M1 = containers.Map(dat1,d1);
M2 = containers.Map(dat2,d2);

keySet1 = keys(M1);
keySet2 = keys(M2);
valueSet1 = values(M1);
valueSet2 = values(M2);

cnt=1;
for i = 1:length(keySet1)-38
    if mod(i,39)==1
        A1=[];
        A2=[];
        for j=i:i+38
            A1=[A1,valueSet1(j)];
            A2=[A2,valueSet2(j)];
        end
        sum1=0;
        sum11=0;
        for i=1:length(A1)
          sum1=sum1+cell2mat(A1(i));
          sum11=sum11+cell2mat(A2(i));
        end
        avg1(cnt)=sum1/length(A1); %the mean
        avg2(cnt)=sum11/length(A2); %the mean
        sum2=0;
        sum22=0;
        for i=1:length(A1)
            sum2=sum2+(cell2mat(A1(i))-avg1(cnt))^5;
            sum22=sum22+(cell2mat(A2(i))-avg2(cnt))^5;
        end
        V1(cnt)=sum2/length(A1); %Varaince
        V2(cnt)=sum22/length(A2); %Varaince 
        cnt=cnt+1;
    end
end

for i=1:length(V1)
    w1(i)=V2(i)/(V1(i)+V2(i));
    w2(i)=V1(i)/(V1(i)+V2(i));
end

cnt=1;
for i=1:length(keySet1)-38
    if mod(i,39)==1
        in(cnt)=(cell2mat(keySet1(i))+cell2mat(keySet1(i+1)))/2;
        cnt=cnt+1;
    end
end

in=in';
w1=w1';
w2=w2';
avg1=avg1';
avg2=avg2';
% in=in(1:8);
% w1=w1(1:8);
% w2=w2(1:8);
% avg1=avg1(1:8);
% avg2=avg2(1:8);

[curve1, goodness1, output1] = fit(in,w1,'poly3');
figure(3)
plot(curve1,in,w1)

[curve2, goodness2, output2] = fit(in,w2,'poly3');
hold on;
plot(curve2,'g',in,w2,'+')
% xlim([1.1 1.3])
xlabel('Sharpness and Texture Ratios','fontsize',12);
ylabel('Relative Weights','fontsize',12);

[curve3, goodness3, output3] = fit(in,avg1,'poly3');
figure(4)
plot(curve3,in,avg1)
[curve4, goodness4, output4] = fit(in,avg2,'poly3');
hold on;
plot(curve4,'g',in,avg2,'+')
xlabel('Sharpness and Texture Ratios','fontsize',12);
ylabel('Relative Biases','fontsize',12);

x=texture_ratio.^5+sharpness_ratio.^5;

weight1=curve1(x);
weight2=curve2(x);
bias1=curve3(x);
bias2=curve4(x);

save('weights2.mat','weight1','weight2','bias1','bias2','name');


