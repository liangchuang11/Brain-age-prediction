%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define the PLS calculate stats function with the following arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Revise from Petra Vertes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
warning off           
close all              
clear                   
clc                    
disp('Re-run PLS to get explained variance and associated stats')

%import response variables
MRIdata= xlsread('\file_path\tvalue.xlsx');

%import predictor variables
genedata=xlsread('\file_path\genedata.xlsx');

X=zscore(genedata);
Y=zscore(MRIdata);

dim=15;
[XL,YL,XS,YS,BETA,PCTVAR,MSE,stats]=plsregress(X,Y,dim);
%temp=cumsum(100*PCTVAR(2,1:dim));
for i=1:dim
    Rsquared(i,1) = 100*PCTVAR(2,i);
end
%align PLS components with desired direction%
for i=1:dim
    [R1(1,i),P1(1,i)]=corr(XS(:,i),MRIdata);
end
for i=1:dim
    if R1(1,i)<0
        XS(:,i)=-1*XS(:,i);
    end
end
% if R1(3,1)<0
%     XS(:,3)=-1*XS(:,3);
% end
% if R1(4,1)<0
%     XS(:,4)=-1*XS(:,4);
% end
% if R1(5,1)<0
%     XS(:,5)=-1*XS(:,5);
% end

%calculate correlations of PLS components with MRI variables

for i=1:dim
    [R2(1,i),P2(1,i)]=corr(XS(:,i),MRIdata);
end

% [R3,p3]=corr(XS(:,3),MRIdata);
% [R4,p4]=corr(XS(:,4),MRIdata);
% [R5,p5]=corr(XS(:,5),MRIdata);
%a=[R1',p1',R2',p2',R3',p3',R4',p4',R5',p5'];


%assess significance of PLS result

for j=1:1000
    
    order=randperm(size(Y,1));
    Yp=Y(order,:);
    [XLr,YLr,XSr,YSr,BETAr,PCTVARr,MSEr,statsr]=plsregress(X,Yp,dim);
    %        temp=cumsum(100*PCTVARr(2,1:dim));
    for i=1:dim
        Rsq(i,j) = 100*PCTVARr(2,i);
    end
    j
end

for i=1:dim
    p(1,i)=length(find(Rsq(i,:)>=Rsquared(i,1)))/j;
end
% plot histogram
% hist(Rsq,30)
% hold on
% plot(Rsquared,20,'.r','MarkerSize',15)
% set(gca,'Fontsize',14)
% xlabel('R squared','FontSize',14);
% ylabel('Permuted runs','FontSize',14);
% title('p<0.0001')

%save stats
%myStats=[PCTVAR; p, j];
%csvwrite(fullfile(output_dir,'PLS_stats.csv'),myStats);
