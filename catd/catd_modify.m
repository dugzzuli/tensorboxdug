% ===============================================================
% This code presents context-aware tensor decomposition (CATD) with element-wise solver.
% 
% Input: Tensor A (N * M * 2L): road segments-dirvers-time slots tensor
%        Matrix X (2L * P): temporal feature matrix (correlation between coarse-grained traffic conditions in a city)
%        Matrix Y (N * Q): geographic feature matrix (including road segments feature and POIs)
% 
% Output: core tensor S
%         latent factor matrix R, U, T
%         A_rec = S \times_{R} R \times_{U} U \times_{T} T
%
% Note: This code uses one auxilliary toolbox named as Tensor toolbox: version 2.5.
%           See http://csmr.ca.sandia.gov/~tgkolda/TensorToolbox/.
%
% For more details about this code, please refer to our paper:
%       Travel Time Estimation of a Path using Sparse Trajectories.
%       Yilun Wang, Yu Zheng, Yexiang Xue, Eric Chang. 
%       In KDD 2014.
%
% Copyright by the paper authors.
% If you have any question, please send email to yuzheng@microsoft.com.
% July 15th, 2014.
% ===============================================================

function [S,R,U,T] = catd(A, epsilon) 

% size of core Tensor
dimR = 5;
dimU = 5;
dimT = 5;

% step size
t0 = 10000000000000; 
t = t0;

lambdaR = 0.01;
lambdaU = 0.01;
lambdaT = 0.01;
lambdaS = 0.01;

dim1 = size(A, 1);
dim2 = size(A, 2);
dim3 = size(A, 3);

%initialize S R U T F G with small random values
R = rand(dim1, dimR);
U = rand(dim2, dimU);
T = rand(dim3, dimT);
S = tenrand(dimR, dimU, dimT);

%找到非零元素的个数和�?
[indexs, values] = find(A);
%非零元素个数的索�?
turn = 1 : length(values);

% initialize function loss
loss_t = epsilon + 1;
loss_t1 = 0;

while abs( loss_t - loss_t1 )> epsilon   
    % optimize each element in randomized sequence   
    for num = 1 : length(values) - 1     
        change = randi([num + 1, length(values)]);
        temp = turn(num);
        turn(num) = turn(change);
        turn(change) = temp;
    end
   
    for num = 1 : length(values) % for every nonzero entries in A
        if (isnan(S(1, 1, 1)))  % check for NAN
            disp nanerror;
            return;
        end
        
        tnum = turn(num);
        nita = 1 / sqrt(t);  % step size
        t = t + 1;
        i = indexs(tnum, 1);
        j = indexs(tnum, 2);
        k = indexs(tnum, 3);
        
        Ri = R(i, :)';
        Uj = U(j, :)';
        Tk = T(k, :)';
               
        SRi = double(ttv(S, {Ri}, 1));
        Fijk = Uj' * SRi * Tk;  % 计算Yijk
        
        Yijk = values(tnum); % 数组中的Aijk
        Lfy = Fijk - Yijk; % 计算 
        nitaLfy = nita * Lfy; 
        
        SUj = double(ttv(S, {Uj}, 2));
        RLfy = nitaLfy * SUj * Tk;
        ULfy = nitaLfy * SRi * Tk;
        TLfy = nitaLfy * (Uj' * SRi)';
        
        SLfy = tensor(ktensor({nitaLfy * Ri, Uj, Tk}));
		
        R(i,:) = ((1 - nita * lambdaR) * Ri - RLfy )';
        U(j,:) = ((1 - nita * lambdaU) * Uj - ULfy)';       
        T(k,:) = ((1 - nita * lambdaT) * Tk - TLfy)' ;
        S = (1 - nita * lambdaS) * S - SLfy;
        
    end
    
    % compute function loss 
    c = size(values);
    for j = 1 : length(values)
        ijk = ttv(S, {R(indexs(j, 1), :)', U(indexs(j, 2), :)', T(indexs(j, 3), :)'});
        c(j) = ijk;
    end
    loss_t = loss_t1;
    loss_t1 = norm(c-values');
    ccc=abs(loss_t - loss_t1)
end   
end
