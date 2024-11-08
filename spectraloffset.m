%模拟不同中心波长，光谱分辨率，FWHM导致的光谱漂移
clc
clear
close all;
[num,txt,raw]=xlsread('tra1000-2600.xlsx');
[num1,txt1,raw1]=xlsread('MMS噪声（F积分时间）.xlsx');
load supersp.mat;
load supersun.mat;
x=num(:,1);
%%%不同天顶角下双层透过率
tra4=num.*(num.^(1/cos(30*pi/180)));%30°天顶角
L0(:,1)=x;
%%%计算理想辐亮度
for i=2:30
 L30(:,i)=supersun(:,2).*cos(30*pi/180).*tra4(:,2).*supersp(:,i).*1000./2.25./pi;%30°天顶角
end
%%%%%%将矿物辐亮度平均：
AL30=mean(L30(:,2:30),2);
% % %%%%%模拟不同中心波长，光谱分辨率，FWHM导致的光谱漂移
b1=0;
n=8; %模拟常见仪器光谱分辨率,即FWHM：2-20nm，间隔2nm
    centerb=0;
    b1=b1+1;
    bands=fix(1600/n);
    for i=2:bands
%         centerb(1,1,b1)=1000;
        centerb(1,1)=1000+n/2;
        centerb(i,1)=centerb(1,1)+n*(i-1);
    end    
    for i=1:bands
    centerb(i,2)=pchip(x(:,1),AL30(:,1),centerb(i,1)/1000);
    end
    for i=1:bands
    noise(i,2)=pchip(num1(:,1),num1(:,4),centerb(i,1)/1000);
    end
    cno1=0;
    for cno=-2.0*n:0.1*n:2.0*n%模拟中心波长漂移量:-2×采样间隔～2×采样间隔，间隔0.1×采样间隔
        cno1=cno1+1;
        f1=0;
        fwhm=n;
        for f=-1.0*fwhm:0.1*fwhm:1.0*fwhm%FWHM变化量:-1×FWHM～1×FWHM，间隔0.1×FWHM 
            f1=f1+1;
            for no=0.0:0.2:2.0%噪声等级信噪比，10个等级，最大为MMS噪声的2倍，最小为0。
            for j=1:bands
               rsp=exp(-(x*1000-centerb(j,1)-cno).^2/(2*((fwhm+f)/2.355)^2));
                mAL0(j,2:101)=repmat(sum(rsp.*AL30(:,1))/sum(rsp),1,100)+normrnd(0,noise(j,2).*no,1,100);
            end
                if find(isnan(mAL0))
                    continue
                 else
                    mAL0(1:bands,1)=centerb(:,1);
                    mAL01=mAL0(1:bands,:);
                    fname = [ 'band','_',num2str(n),'_','shift','_',num2str(cno,'%4.1f'),'_','fwhm','_',num2str(f,'%4.1f'),'_','noise','_',num2str(no,'%4.1f'),'.csv'];
                   csvwrite(['F:\spectraldata1\bands_202\',fname],mAL01); 
                end
            end
        end
    end

