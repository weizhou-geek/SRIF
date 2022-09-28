clear all
close all
clc

addpath('./iwssim_iwpsnr/');
addpath('./LabForComputationalVision-matlabPyrTools-42e4602/');

window = fspecial('gaussian',7,7/6);
window = window/sum(sum(window));

%% 2
imgdis_dir=dir('../sr_image/');

cnt=1;
quality=ones(312,1);
kld_value5=ones(312,1);
kld_value4=ones(312,1);
msssim_value3=ones(312,1);
kld_value3=ones(312,1);
msssim_value2=ones(312,1);
kld_value2=ones(312,1);
msssim_value1=ones(312,1);
kld_value1=ones(312,1);
wl2=ones(312,1);
wh2=ones(312,1);
l2_ref=ones(312,1);
l2=ones(312,1);
h2_ref=ones(312,1);
h2=ones(312,1);
wl1=ones(312,1);
wh1=ones(312,1);
l1_ref=ones(312,1);
l1=ones(312,1);
h1_ref=ones(312,1);
h1=ones(312,1);
img_name=cell(312,1);

cnt=1;
for i=3:length(imgdis_dir)
    i
    imgdis_name1 = imgdis_dir(i).name;
    imgdis_dir1 = dir(['../sr_image/', imgdis_name1]);
    for j=3:length(imgdis_dir1)
        imgdis_name2 = imgdis_dir1(j).name;
        imgdis_dir2 = dir(['../sr_image/', imgdis_name1, '/', imgdis_name2]);
        for m=3:10
            imgdis_name3 = imgdis_dir2(m).name;
            imgdis_dir3 = dir(['../sr_image/', imgdis_name1, '/', imgdis_name2, '/', imgdis_name3]);
                    
            img_name(cnt) = cellstr([imgdis_name1, '/', imgdis_name2, '/', imgdis_name3]);
            ref_index = strfind(char(img_name(cnt)),'.bmp');
            a = char(img_name(cnt));
            ref_name = [a(1:ref_index-2), 'O.bmp'];
            
            dim = imread(['../sr_image/', char(img_name(cnt))]);
            imgref = imread(['../sr_image/', ref_name]);
            
            if(size(dim,3)==3)
                dim = rgb2gray(dim);
            end
            if(size(imgref,3)==3)
                imgref = rgb2gray(imgref);
            end  
            
            targetSize=size(dim);
            if targetSize==511
                oim = imcrop(imgref, [1.5 1.5 510 510]);       
            end
            
            if targetSize==509
                oim = imcrop(imgref, [2.5 2.5 508 508]);       
            end
            
            if targetSize==505
                oim = imcrop(imgref, [4.5 4.5 504 504]);     
            end
            
            oimg = double(oim);
             
            gauss_pyr_img_ori = gauss_pyramid(oimg, 6);
            lapl_pyr_img_ori = lapl_pyramid(gauss_pyr_img_ori);

            %2
            laplacianL5_ori = lapl_pyr_img_ori{5};
            laplacianL4_ori = lapl_pyr_img_ori{4};
            gaussianL3_ori = gauss_pyr_img_ori{3};
            laplacianL3_ori = lapl_pyr_img_ori{3};
            gaussianL2_ori = gauss_pyr_img_ori{2};
            laplacianL2_ori = lapl_pyr_img_ori{2};
            gaussianL1_ori = gauss_pyr_img_ori{1};
            laplacianL1_ori = lapl_pyr_img_ori{1};

            [LPC_SI_ref lpc_map] = lpc_si(gaussianL2_ori);
            J1 = rangefilt(uint8(gaussianL2_ori));
            m1 = entropy(J1);

            f2 = double(uint8(laplacianL2_ori)+128);
            mu            = filter2(window, f2, 'same');
            mu_sq         = mu.*mu;
            sigma         = sqrt(abs(filter2(window, f2.*f2, 'same') - mu_sq));
            structref     = (f2-mu)./(sigma+1);
            [LPC_SI_ref1 lpc_map1] = lpc_si(gaussianL1_ori);
            J11 = rangefilt(uint8(gaussianL1_ori));
            m11 = entropy(J11);

            f2 = double(uint8(laplacianL3_ori)+128);
            mu            = filter2(window, f2, 'same');
            mu_sq         = mu.*mu;
            sigma         = sqrt(abs(filter2(window, f2.*f2, 'same') - mu_sq));
            structref3     = (f2-mu)./(sigma+1);

            f2 = double(uint8(laplacianL1_ori)+128);
            mu            = filter2(window, f2, 'same');
            mu_sq         = mu.*mu;
            sigma         = sqrt(abs(filter2(window, f2.*f2, 'same') - mu_sq));
            structref1     = (f2-mu)./(sigma+1);

            f2 = double(uint8(laplacianL4_ori)+128);
            mu            = filter2(window, f2, 'same');
            mu_sq         = mu.*mu;
            sigma         = sqrt(abs(filter2(window, f2.*f2, 'same') - mu_sq));
            structref4     = (f2-mu)./(sigma+1);

            f2 = double(uint8(laplacianL5_ori)+128);
            mu            = filter2(window, f2, 'same');
            mu_sq         = mu.*mu;
            sigma         = sqrt(abs(filter2(window, f2.*f2, 'same') - mu_sq));
            structref5     = (f2-mu)./(sigma+1);

            
            dimg = double(dim);
            gauss_pyr_img_dis = gauss_pyramid(dimg, 6);
            lapl_pyr_img_dis = lapl_pyramid(gauss_pyr_img_dis); 

            %2
            laplacianL5_dis = gauss_pyr_img_dis{5};
            laplacianL4_dis = lapl_pyr_img_dis{4};
            gaussianL3_dis = gauss_pyr_img_dis{3};
            laplacianL3_dis = lapl_pyr_img_dis{3};
            gaussianL2_dis = gauss_pyr_img_dis{2};
            laplacianL2_dis = lapl_pyr_img_dis{2};
            gaussianL1_dis = gauss_pyr_img_dis{1};
            laplacianL1_dis = lapl_pyr_img_dis{1};

            [LPC_SI_dis lpc_map] = lpc_si(gaussianL2_dis);
            J2 = rangefilt(uint8(laplacianL2_dis)+128);
            m2 = entropy(J2);

            f2 = double(uint8(laplacianL2_dis)+128);
    %         [~,~,f2]=ZCA(f2);
            mu            = filter2(window, f2, 'same');
            mu_sq         = mu.*mu;
            sigma         = sqrt(abs(filter2(window, f2.*f2, 'same') - mu_sq));
            structdis2     = (f2-mu)./(sigma+1);

            [LPC_SI_dis1 lpc_map1] = lpc_si(gaussianL1_dis);
            J21 = rangefilt(uint8(laplacianL1_dis)+128);
            m21 = entropy(J21);

            f2 = double(uint8(laplacianL3_dis)+128);
    %         [~,~,f2]=ZCA(f2);
            mu            = filter2(window, f2, 'same');
            mu_sq         = mu.*mu;
            sigma         = sqrt(abs(filter2(window, f2.*f2, 'same') - mu_sq));
            structdis3     = (f2-mu)./(sigma+1);

            f2 = double(uint8(laplacianL1_dis)+128);
    %         [~,~,f2]=ZCA(f2);
            mu            = filter2(window, f2, 'same');
            mu_sq         = mu.*mu;
            sigma         = sqrt(abs(filter2(window, f2.*f2, 'same') - mu_sq));
            structdis1     = (f2-mu)./(sigma+1);

            f2 = double(uint8(laplacianL4_dis)+128);
    %         [~,~,f2]=ZCA(f2);
            mu            = filter2(window, f2, 'same');
            mu_sq         = mu.*mu;
            sigma         = sqrt(abs(filter2(window, f2.*f2, 'same') - mu_sq));
            structdis4     = (f2-mu)./(sigma+1);

            f2 = double(uint8(laplacianL5_dis)+128);
    %         [~,~,f2]=ZCA(f2);
            mu            = filter2(window, f2, 'same');
            mu_sq         = mu.*mu;
            sigma         = sqrt(abs(filter2(window, f2.*f2, 'same') - mu_sq));
            structdis5     = (f2-mu)./(sigma+1);

            [a1,b1]=imhist(structdis2(:),20);
            a1=a1/length(structdis2(:));
            [a2,b2]=imhist(structref(:),20);
            a2=a2/length(structref(:));
            kld_value2(cnt) = kld(a1, a2);

                a1=imhist(structdis3(:),20);
        %         a1=a1/length(structdis2(:));
                a2=imhist(structref3(:),20);
        %         a2=a2/length(structref(:));
                kld_value3(cnt) = kld(a1, a2);  

                a1=imhist(structdis1(:),20);
        %         a1=a1/length(structdis2(:));
                a2=imhist(structref1(:),20);
        %         a2=a2/length(structref(:));
                kld_value1(cnt) = kld(a1, a2); 

                a1=imhist(structdis4(:),20);
        %         a1=a1/length(structdis2(:));
                a2=imhist(structref4(:),20);
        %         a2=a2/length(structref(:));
                kld_value4(cnt) = kld(a1, a2); 

                a1=imhist(structdis5(:),20);
        %         a1=a1/length(structdis2(:));
                a2=imhist(structref5(:),20);
        %         a2=a2/length(structref(:));
                kld_value5(cnt) = kld(a1, a2); 

                [M N] = size(gaussianL2_ori);
                level = 5;
                min_img_width = floor(min(M, N)/(2^(level-1)));
                if mod(min_img_width, 2) == 0
                    min_img_width = min_img_width-1;
                end
                msssim_value2(cnt) = iwssim(uint8(gaussianL2_ori), uint8(gaussianL2_dis), 1, 5, [0.01 0.03], 255, [0.0448 0.2856 0.3001 0.2363 0.1333], fspecial('gaussian', min_img_width, 1.5));
%                     msssim_value2(cnt) = msssim(uint8(gaussianL2_ori), uint8(gaussianL2_dis), [0.01 0.03], fspecial('gaussian', min_img_width, 1.5));
%                     msssim_value2(cnt) = ssim_index(uint8(gaussianL2_ori), uint8(gaussianL2_dis));

                [M N] = size(gaussianL3_ori);
                level = 5;
                min_img_width = floor(min(M, N)/(2^(level-1)));
                if mod(min_img_width, 2) == 0
                    min_img_width = min_img_width-1;
                end
%                     msssim_value3(cnt) = msssim(uint8(gaussianL3_ori), uint8(gaussianL3_dis), [0.01 0.03], fspecial('gaussian', min_img_width, 1.5));
                msssim_value3(cnt) = iwssim(uint8(gaussianL3_ori), uint8(gaussianL3_dis), 1, 5, [0.01 0.03], 255, [0.0448 0.2856 0.3001 0.2363 0.1333], fspecial('gaussian', min_img_width, 1.5));
%                     msssim_value3(cnt) = ssim_index(uint8(gaussianL3_ori), uint8(gaussianL3_dis));

                [M N] = size(gaussianL1_ori);
                level = 5;
                min_img_width = floor(min(M, N)/(2^(level-1)));
                if mod(min_img_width, 2) == 0
                    min_img_width = min_img_width-1;
                end
%                     msssim_value1(cnt) = msssim(uint8(gaussianL1_ori), uint8(gaussianL1_dis), [0.01 0.03], fspecial('gaussian', min_img_width, 1.5));
                msssim_value1(cnt) = iwssim(uint8(gaussianL1_ori), uint8(gaussianL1_dis), 1, 5, [0.01 0.03], 255, [0.0448 0.2856 0.3001 0.2363 0.1333], fspecial('gaussian', min_img_width, 1.5));
%                     msssim_value1(cnt) = ssim_index(uint8(gaussianL1_ori), uint8(gaussianL1_dis));

%                     [M N] = size(gaussianL4_ori);
%                     level = 5;
%                     min_img_width = floor(min(M, N)/(2^(level-1)));
%                     msssim_value4(cnt) = msssim(uint8(gaussianL4_ori), uint8(gaussianL4_dis), [0.01 0.03], fspecial('gaussian', min_img_width, 1.5));
%                     msssim_value4(cnt) = iwssim(uint8(gaussianL4_ori), uint8(gaussianL4_dis), 1, 5, [0.01 0.03], 255, [0.0448 0.2856 0.3001 0.2363 0.1333], fspecial('gaussian', min_img_width, 1.5));

                wl2(cnt) = LPC_SI_dis/LPC_SI_ref;
                wh2(cnt) = m2/m1;

                l2_ref(cnt) = LPC_SI_ref;
                l2(cnt) = LPC_SI_dis;
                h2_ref(cnt) = m1;
                h2(cnt) = m2;

                wl1(cnt) = LPC_SI_dis1/LPC_SI_ref1;
                wh1(cnt) = m21/m11;

                l1_ref(cnt) = LPC_SI_ref1;
                l1(cnt) = LPC_SI_dis1;
                h1_ref(cnt) = m11;
                h1(cnt) = m21;

        %         wl1(cnt) = LPC_SI_dis1/LPC_SI_ref1;
        %         wh1(cnt) = m21/m11;

        %         avg_msssim2 = msssim_value2(cnt);
        %         avg_kld2 = -kld_value2(cnt);         
        %         avg_wl2 = wl2(cnt)/(wl2(cnt)+wh2(cnt));
        %         avg_wh2 = wh2(cnt)/(wl2(cnt)+wh2(cnt));
        %         quality(cnt) = (avg_wl2.*avg_msssim2)+(avg_wh2.*avg_kld2);

                cnt = cnt+1;
        end
    end
end

% avg_msssim1 = normalization(msssim_value1, 0, 1, max(msssim_value1), min(msssim_value1));
% avg_kld1 = 1-normalization(kld_value1, 0, 1, max(kld_value1), min(kld_value1));
% avg_msssim2 = normalization(msssim_value2, 0, 1, max(msssim_value2), min(msssim_value2));
% avg_kld2 = 1-normalization(kld_value2, 0, 1, max(kld_value2), min(kld_value2));
% quality = wl1.*avg_msssim1+wh1.*avg_kld1+wl2.*avg_msssim2+wh2.*avg_kld2;

save('quality_LPCSI.mat','kld_value5','kld_value4','msssim_value1','kld_value1','msssim_value2','kld_value2','msssim_value3','kld_value3','wl2','wh2','l2_ref','l2','h2_ref','h2','wl1','wh1','l1_ref','l1','h1_ref','h1','img_name','-v6');


