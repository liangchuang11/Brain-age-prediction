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

for i=1:dim
    [R2(1,i),P2(1,i)]=corr(XS(:,i),MRIdata);
end


for j=1:1000
    
    order=randperm(size(Y,1));
    Yp=Y(order,:);
    [XLr,YLr,XSr,YSr,BETAr,PCTVARr,MSEr,statsr]=plsregress(X,Yp,dim);
    for i=1:dim
        Rsq(i,j) = 100*PCTVARr(2,i);
    end
    j
end

for i=1:dim
    p(1,i)=length(find(Rsq(i,:)>=Rsquared(i,1)))/j;
end

