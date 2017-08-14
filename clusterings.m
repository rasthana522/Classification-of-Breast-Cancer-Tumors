fig_num = 0;
filename = 'BreastData.xlsx';
xlrange = 'B2:AF19';
X = xlsread(filename, xlrange);
X1 = [X(:,1) X(:,17:30)];
rng(1);

%%% Creating Color Matrix for Plotting
green = [0.3 0.8 0.3];
red = [1 0 0];
C = cell(18, 3);
[C(1:9,1),C(1:9,2),C(1:9,3)] = deal({green(1)},{green(2)},{green(3)});
[C(10:18,1),C(10:18,2),C(10:18,3)] = deal({red(1)},{red(2)},{red(3)});
color_mat = cell2mat(C);
clear C;
clust_cols = [[0 0 0]; [0.15 0 1]; [0.45 0 0.7]; [0.6 0.3 0.3]; [0.3 0.6 1]; [1 0.4 0]; [0.65 0.5 0]];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%  K-means Clustering  %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% k = 4;
% gene_labels = string({'COX2','DBC2','MYC','CCND1','CHD1','P53','HER2'});
% length = size(gene_labels);
% for index = 1:length(2)
%     gene = gene_labels(index);
%     i = 2*index;
%     X_gain_loss = X1(:,i:i+1);
%     [idx,centroids] = kmeans(X_gain_loss,k);
%     group = string(idx);
%     
%     fig_num = fig_num + 1;
%     figure(fig_num)
%     v = 1:18;
%     sz = 30;
%     scatter(X_gain_loss(intersect(v(group == '1'),1:9),1),X_gain_loss(intersect(v(group == '1'),1:9),2),sz,green,'o');
%     hold on
%     scatter(X_gain_loss(intersect(v(group == '1'),10:18),1),X_gain_loss(intersect(v(group == '1'),10:18),2),sz,'r','o');
%     scatter(X_gain_loss(intersect(v(group == '2'),1:9),1),X_gain_loss(intersect(v(group == '2'),1:9),2),sz,green,'x');
%     scatter(X_gain_loss(intersect(v(group == '2'),10:18),1),X_gain_loss(intersect(v(group == '2'),10:18),2),sz,'r','x');
%     scatter(X_gain_loss(intersect(v(group == '3'),1:9),1),X_gain_loss(intersect(v(group == '3'),1:9),2),sz,green,'+');
%     scatter(X_gain_loss(intersect(v(group == '3'),10:18),1),X_gain_loss(intersect(v(group == '3'),10:18),2),sz,'r','+');
%     scatter(X_gain_loss(intersect(v(group == '4'),1:9),1),X_gain_loss(intersect(v(group == '4'),1:9),2),sz,green,'s');
%     scatter(X_gain_loss(intersect(v(group == '4'),10:18),1),X_gain_loss(intersect(v(group == '4'),10:18),2),sz,'r','s');
%     scatter(X_gain_loss(intersect(v(group == '5'),1:9),1),X_gain_loss(intersect(v(group == '5'),1:9),2),sz,green,'p');
%     scatter(X_gain_loss(intersect(v(group == '5'),10:18),1),X_gain_loss(intersect(v(group == '5'),10:18),2),sz,'r','p');
%     scatter(X_gain_loss(intersect(v(group == '6'),1:9),1),X_gain_loss(intersect(v(group == '6'),1:9),2),sz,green,'*');
%     scatter(X_gain_loss(intersect(v(group == '6'),10:18),1),X_gain_loss(intersect(v(group == '6'),10:18),2),sz,'r','*');
%     title(strcat(sprintf('K=%d-means Clustering on Gain/Loss ',k),{' '},gene));
%     xlabel(strcat('Gain ',{' '},gene));
%     ylabel(strcat('Loss ',{' '},gene));
%     hold off
% %     str = strcat(sprintf('k%d_means',k),{'_'},gene);
% %     str = str{1};
% %     saveas(gcf,str,'jpeg');
%     clear str i idx centroids gene X_gain_loss group v sz
% end
% clear gene_labels index length k

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Hierarchical Clustering %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

k0 = 5;
TreeWeightedCorr = linkage(X1(:,2:15), 'weighted', 'correlation');
cutoff_wei_corr = get_maxcluster_cutoff(TreeWeightedCorr, k0);
fig_num = fig_num + 1;
figure(fig_num)
[H_wei_corr,cluster_wei_corr,perm_wei_corr] = dendrogram(TreeWeightedCorr, 'ColorThreshold', cutoff_wei_corr);
H_wei_corr = recolor_dendrogram(H_wei_corr, clust_cols);
set(H_wei_corr,'LineWidth',2)
title('Hierarchical Clustering with Correlation Distance');
xlabel('Sample Number');
ylabel('Height');
% str = sprintf('hier_weighted_corr_k%d',k0);
% saveas(gcf,str,'jpeg');
clear str k0 cutoff_wei_corr cluster_wei_corr perm_wei_corr H_wei_corr

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%   Laplace Embedding (KNN)   %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

k = 5;
m = 15;
[LapKX1,LambdaLapK] = lap_embed_knn(X1(:,2:15), k, m);
clear m

%%% 3-D Plots of Laplace Embedding (KNN)
% fig_num = fig_num + 1;
% figure(fig_num)
% [i1,i2,i3] = deal(3,4,5);
% scatter3(LapKX1(:,i1),LapKX1(:,i2),LapKX1(:,i3),15,color_mat,'fill');
% %axis([-0.1 0.1 -0.1 0.1 -0.1 0.1]);
% title(sprintf('Laplace Embedding (K=%d-NN) with %d,%d,%d',k,i1,i2,i3));
% xlabel(sprintf('Eig %d',i1));
% ylabel(sprintf('Eig %d',i2));
% zlabel(sprintf('Eig %d',i3));
% clear i1 i2 i3

%%% 2-D Plots of Laplace Embedding (KNN)
% fig_num = fig_num + 1;
% figure(fig_num)
% [j1,j2] = deal(2,4);
% scatter(LapKX1(:,j1),LapKX1(:,j2),30,color_mat,'fill');
% %axis([-0.1 0.1 -0.1 0.1 -0.1 0.1]);
% title(sprintf('Laplace Embedding (K=%d-NN) with %d,%d',k,j1,j2));
% xlabel(sprintf('Eig %d',j1));
% ylabel(sprintf('Eig %d',j2));
% % str = sprintf('lap_knn_k%d',k);
% % saveas(gcf,str,'jpeg');
% clear str j1 j2

%%% Hierarchical Clustering of Embedding (KNN)
k0 = 5;
TreeLapKWard = linkage(LapKX1(:,1:k0), 'ward', 'euclidean');
cutoff_lapK_ward = get_maxcluster_cutoff(TreeLapKWard, k0);
fig_num = fig_num + 1;
figure(fig_num)
[H_lapK_ward,cluster_lapK_ward,perm_lapK_ward] = dendrogram(TreeLapKWard, 'ColorThreshold', cutoff_lapK_ward);
H_lapK_ward = recolor_dendrogram(H_lapK_ward, clust_cols);
set(H_lapK_ward,'LineWidth',2)
title(sprintf('Hierarchical Clustering (after Laplace Embedding from K=%d-NN)',k));
xlabel('Sample Number');
ylabel('Height');
% str = sprintf('hier_lapK%d_ward_k%d',k,k0);
% saveas(gcf,str,'jpeg');
clear str k0 cutoff_lapK_ward cluster_lapK_ward perm_lapK_ward H_lapK_ward

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Laplace Embedding (Correlation) %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

m = 15;
[LapCX1,LambdaLapC] = lap_embed_corr(X1(:,2:15), m);
clear m

%%% 3-D Plots of Laplace Embedding (Correlation)
% fig_num = fig_num + 1;
% figure(fig_num)
% [i1,i2,i3] = deal(2,3,4);
% scatter3(LapCX1(:,i1),LapCX1(:,i2),LapCX1(:,i3),15,color_mat,'fill');
% %axis([-0.1 0.1 -0.1 0.1 -0.1 0.1]);
% title(sprintf('Laplace Embedding (Correlation) with %d,%d,%d',i1,i2,i3));
% xlabel(sprintf('Eig %d',i1));
% ylabel(sprintf('Eig %d',i2));
% zlabel(sprintf('Eig %d',i3));
% clear i1 i2 i3

%%% 2-D Plots of Laplace Embedding (Correlation)
% fig_num = fig_num + 1;
% figure(fig_num)
% [j1,j2] = deal(2,4);
% scatter(LapCX1(:,j1),LapCX1(:,j2),30,color_mat,'fill');
% %axis([-0.1 0.1 -0.1 0.1 -0.1 0.1]);
% title(sprintf('Laplace Embedding (Correlation) with %d,%d',j1,j2));
% xlabel(sprintf('Eig %d',j1));
% ylabel(sprintf('Eig %d',j2));
% % str = 'lap_corr';
% % saveas(gcf,str,'jpeg');
% clear str j1 j2

%%% Hierarchical Clustering of Embedding (Correlation)
k0 = 5;
TreeLapCWard = linkage(LapCX1(:,1:k0), 'ward', 'euclidean');
cutoff_lapC_ward = get_maxcluster_cutoff(TreeLapCWard, k0);
fig_num = fig_num + 1;
figure(fig_num)
[H_lapC_ward,cluster_lapC_ward,perm_lapC_ward] = dendrogram(TreeLapCWard, 'ColorThreshold', cutoff_lapC_ward);
H_lapC_ward = recolor_dendrogram(H_lapC_ward, clust_cols);
set(H_lapC_ward,'LineWidth',2)
title('Hierarchical Clustering (after Laplace Embedding from Correlation)');
xlabel('Sample Number');
ylabel('Height');
% str = sprintf('hier_lapC_ward_k%d',k0);
% saveas(gcf,str,'jpeg');
clear str k0 cutoff_lapC_ward cluster_lapC_ward perm_lapC_ward H_lapC_ward

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%      Functions      %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function cutoff = get_maxcluster_cutoff(Z, k)
    cutoff = Z(end-k+2,3)-eps;
end

function S = similarity_matrix_knn(X, k)
    sz = size(X);
    n1 = sz(1);
    S = double(zeros(n1));
    IDX = knnsearch(X, X, 'K', k);
    for i=1:n1
        for j=1:n1
            sij = (sum(IDX(i,:) == j) > 0);
            sji = (sum(IDX(j,:) == i) > 0);
            if (sij || sji) S(i,j) = 1; else S(i,j) = 0; end;
        end
    end
end

function [LapX,Lambda] = lap_embed_knn(X, k, m)
    S = similarity_matrix_knn(X, k);
    n = length(S);
    d = double(zeros(1,n));
    for i=1:n
        d(i) = sum(S(i,:));
    end
    P = double(zeros(n));
    for i=1:n
        for j=1:n
            P(i,j) = S(i,j)/d(i);
        end
    end
    [V,Lam] = eig(P);
    Lambda = diag(Lam);
    [list,I] = sort(abs(Lambda),'descend');
    Lambda = Lambda(I);
    Lambda = Lambda.';
    V = V(:, I);

    LapX = V(:,1:min(m,n));
    Lambda = Lambda(1:min(m,n));
end

function S = similarity_matrix_corr(X)
    S = squareform(pdist(X,'correlation'));
end

function [LapX,Lambda] = lap_embed_corr(X, m)
    S = similarity_matrix_corr(X);
    n = length(S);
    d = double(zeros(1,n));
    for i=1:n
        d(i) = sum(S(i,:));
    end
    P = double(zeros(n));
    for i=1:n
        for j=1:n
            P(i,j) = S(i,j)/d(i);
        end
    end
    [V,Lam] = eig(P);
    Lambda = diag(Lam);
    [list,I] = sort(abs(Lambda),'descend');
    Lambda = Lambda(I);
    Lambda = Lambda.';
    V = V(:, I);

    LapX = V(:,1:min(m,n));
    Lambda = Lambda(1:min(m,n));
end

function H = recolor_dendrogram(H, clust_cols)
    n = length(H);
    H_color = [1:n; 1:n; 1:n; 1:n; 1:n];
    for i = 1:n
        H_color(2:4,i) = H(i).Color;
        H_color(5,i) = 100*H_color(2,i)+10*H_color(3,i)+H_color(4,i);
    end
    H_color_sorted = transpose(sortrows(transpose(H_color), [2 3 4 1]));
    %disp(H_color_sorted);
    index_list = unique(H_color_sorted(5,:));
    num_colors = length(index_list);
    for i = 1:num_colors
        for j = 1:n
            if H_color_sorted(5,j) == index_list(i)
                H_color_sorted(2:4,j) = deal(clust_cols(i,:));
            end
        end
    end
    H_color = transpose(sortrows(transpose(H_color_sorted), 1));
    for i = 1:n
        H(i).Color = H_color(2:4,i);
    end
end



