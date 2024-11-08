%ģ�ⲻͬ���Ĳ��������׷ֱ��ʣ�FWHM���µĹ���Ư��
clc
clear
close all;
[num,txt,raw]=xlsread('tra1000-2600.xlsx');
[num1,txt1,raw1]=xlsread('MMS������F����ʱ�䣩.xlsx');
load supersp.mat;
load supersun.mat;
x=num(:,1);
%%%��ͬ�춥����˫��͸����
tra4=num.*(num.^(1/cos(30*pi/180)));%30���춥��
L0(:,1)=x;
%%%�������������
for i=2:30
 L30(:,i)=supersun(:,2).*cos(30*pi/180).*tra4(:,2).*supersp(:,i).*1000./2.25./pi;%30���춥��
end
%%%%%%�����������ƽ����
AL30=mean(L30(:,2:30),2);
% % %%%%%ģ�ⲻͬ���Ĳ��������׷ֱ��ʣ�FWHM���µĹ���Ư��
b1=0;
n=8; %ģ�ⳣ���������׷ֱ���,��FWHM��2-20nm�����2nm
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
    for cno=-2.0*n:0.1*n:2.0*n%ģ�����Ĳ���Ư����:-2�����������2��������������0.1���������
        cno1=cno1+1;
        f1=0;
        fwhm=n;
        for f=-1.0*fwhm:0.1*fwhm:1.0*fwhm%FWHM�仯��:-1��FWHM��1��FWHM�����0.1��FWHM 
            f1=f1+1;
            for no=0.0:0.2:2.0%�����ȼ�����ȣ�10���ȼ������ΪMMS������2������СΪ0��
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

