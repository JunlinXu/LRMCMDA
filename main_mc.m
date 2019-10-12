
% rank:Rank of low rank matrix recovery
%interaction:Association
% sm,sd: Representing miRNA and disease similarity, respectively
% iter Number of iterations
%lam1,lam2,lam3: parameters
%tol:Threshold of convergence
%recMatrix:obtain the prediction score matrix
load interaction;
load sd;
load sm;
norm=1;
rank=3;
iter=60;
tol = 1e-8;
lam1=0.1;
lam2=1;
lam3=1;
[m,n]=size(interaction);
Map=MAPS(interaction);
I11=getAdjMA(interaction,Map);
[X S Y dist]=OptSpaceII(I11,sd,sm,rank,iter,tol,lam1,lam2,lam3);
recMatrix = X*S*Y';
    


