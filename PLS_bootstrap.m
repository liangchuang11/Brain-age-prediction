warning off             % 关闭报警信息
close all               % 关闭开启的图窗
clear                   % 清空变量
clc                     % 清空命令行
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define the PLS bootstrap function with the following arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Revise from Petra Vertes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('Running PLS')

%import response variables
MRIdata= xlsread('\file_path\tvalue.xlsx');

%import predictor variables
genedata=xlsread('\file_path\genedata.xlsx');

X=zscore(genedata);
Y=zscore(MRIdata);

%number of bootstrap iterations
bootnum=1000;

dim=15;
[XL,YL,XS,YS,BETA,PCTVAR,MSE,stats]=plsregress(X,Y,dim);

%store regions IDs and weights in descending order of weight for both
%components
%[R1,p1]=corr([XS(:,1),XS(:,2)],MRIdata);

%align PLS components with desired direction%
for i=1:dim
    [R1(1,i),P1(1,i)]=corr(XS(:,i),MRIdata);
end
for i=1:dim
    if R1(1,i)<0
        XS(:,i)=-1*XS(:,i);
    end
end
[PLS1w,x1] = sort(stats.W(:,1),'descend');%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PLS1ids=genes(x1);
geneindex1=geneindex(x1);

%define variables for storing the (ordered) weights from all bootstrap runs
PLS1weights=[];
% PLS2weights=[];

%start bootstrap
disp('  Bootstrapping - could take a while')
for i=1:bootnum
    myresample = randsample(size(X,1),size(X,1),1);
    res(i,:)=myresample; %store resampling out of interest
    Xr=X(myresample,:); % define X for resampled regions
    Yr=Y(myresample,:); % define Y for resampled regions
    [XL,YL,XS,YS,BETA,PCTVAR,MSE,stats]=plsregress(Xr,Yr,dim); %perform PLS for resampled data
    
    temp=stats.W(:,1);%extract PLS1 weights%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    newW=temp(x1); %order the newly obtained weights the same way as initial PLS 
    if corr(PLS1w,newW)<0 % the sign of PLS components is arbitrary - make sure this aligns between runs
        newW=-1*newW;
    end
    PLS1weights=[PLS1weights,newW];%store (ordered) weights from this bootstrap run
    i
%     temp=stats.W(:,2);%extract PLS2 weights
%     newW=temp(x2); %order the newly obtained weights the same way as initial PLS 
%     if corr(PLS2w,newW)<0 % the sign of PLS components is arbitrary - make sure this aligns between runs
%         newW=-1*newW;
%     end
%     PLS2weights=[PLS2weights,newW]; %store (ordered) weights from this bootstrap run    
end

%get standard deviation of weights from bootstrap runs
PLS1sw=std(PLS1weights');
% PLS2sw=std(PLS2weights');

%get bootstrap weights
temp1=PLS1w./PLS1sw';
% temp2=PLS2w./PLS2sw';

%order bootstrap weights (Z) and names of regions
[Z1 ind1]=sort(temp1,'descend');
PLS1=PLS1ids(ind1);
geneindex1=geneindex1(ind1);
% [Z2 ind2]=sort(temp2,'descend');
% PLS2=PLS2ids(ind2);
% geneindex2=geneindex2(ind2);


%print out results
% fid1 = fopen(fullfile(output_dir,'PLS1_geneWeights.csv'),'w');
% for i=1:length(genes)
%   fprintf(fid1,'%s, %d, %f\n', PLS1{i}, geneindex1(i), Z1(i));
% end
% fclose(fid1);
% 
% fid2 = fopen(fullfile(output_dir,'PLS2_geneWeights.csv'),'w');
% for i=1:length(genes)
%   fprintf(fid2,'%s, %d, %f\n', PLS2{i},geneindex2(i), Z2(i));
% end
% fclose(fid2);
