function [w_approximate] = active_ranking(X,w_star,R)
% X : (n-by-d) matrix of n objects in d dimensions
% R : 0,1,2,... dictating how many votes are used for robust version. For
%     noiseless algorithm set R=0 (default). For guidance see paper below.
%     
% author: Kevin Jamieson
% date: 9/26/2011
% website: www.cae.wisc.edu/~jamieson/
% relavent paper: Jamieson, K., Nowak, R., "Active Ranking using Pairwise
% Comparisons." NIPS 2011.

% add by gylai
% record the number of ambiguous queries

cnt_total_pairwise_comparison = 0;
cnt_unambiguous_query = 0;
cnt_ambiguous_query = 0;

% add by gylai
% nargin: number of parameters
if nargin < 3
    R = 0; % number of objects that will be used for voting in the robust case
end

Xf = X;

n = size(X,1); % number of objects
d = size(X,2); % dimension

% build hyperplanes
% add by gylai
% Geometry: temp is the perpendicular bisector of object{i} & object{j}
% embedded into a d-dimensional Euclidean space
H = zeros(d+1,n,n); % n  (d+1, n) matrices
for j = 2:n
    for i = 1:j-1
        w = (X(j,:)-X(i,:));
        % add by gylai
        % w' -> conjugate transposition of w
        w0 = -.5*(X(j,:)+X(i,:))*w';
        const = norm(abs([w,w0]))';
        w = w/const;
        w0 = w0/const;
        temp = [ w'; w0 ];
        H(:,i,j) = temp;
        % add by gylai
        % FIXME
        H(:,j,i) = -H(:,i,j); % equal to following ...
%         w = (X(i,:)-X(j,:));
%         w0 = -.5*(X(i,:)+X(j,:))*w';
%         const = norm(abs([w,w0]))';
%         w = w/const;
%         w0 = w0/const;
%         temp = [ w'; w0 ];
%         H(:,j,i) = temp;

%         a = input('pause');
    end
end

% fprintf('<info> matrix of hyperplanes H:\n');
% for j = 2:n
%     for i = 1:j-1
%         fprintf('H(%d, %d):\n', i, j);
%         disp(H(:, i, j));
%     end
% end

% add by gylai
% nchoosek(n,2) -> C^{2}_{n}
% all is initialized to 0
error = zeros( nchoosek(n,2) , 1 );

% randomly enumerate the objects
% add by gylai
% [n*1] vector of random numbers in [0, 1]
% v -> vector of ordered values
% hi -> [n*1] vector of random sequence of 1,2,...,n
% add ends
[v hi] = sort(rand(n,1));
hi =[5;
     2;
     9;
     8;
     6;
     3;
     4;
     1;
    10;
     7]
% disp(v);
% disp(hi);
%%%%%%%%%%%

% start algorithm
Qh = zeros(n,n); % for noise-free case -> preference label of pairwise comparisons
Qk = zeros(n,n); % for noise case
known = []; % records for object index of pairwise comparisons

placeHolder = -1;
secondPhase = 0;
jj = 1;
while jj < n
    jj = jj + 1;
    % add by gylai
    % compare object {hi(1), hi(2), ... , hi(jj)}
    fprintf('jj: %d\n', jj);
    
    % add by gylai
    % number of undefined preference relations in Qh
    goodInds = ( Qh( hi(jj), hi(1:jj-1) ) )==0;
    goodInds
    
    while sum(goodInds)>0
        % add by gylai
        % {1,2,..., jj-1}
        toSort = 1:(jj-1);
        goodInds
        toSort = toSort(goodInds)
        Qhyp = Qh(hi(toSort),hi(toSort));
        Qhyp
        
        index = compare_sort(1:sum(goodInds),Qhyp);
        index
        % add by gylai
        % 
        list = toSort(index)
        
        bis = floor((length(list)+1)/2)
        ii = list(bis);
        above = list(bis+1:end)
        below = list(1:bis-1)
        a = input('in:');
        % add by gylai
        % judge the preference relation between object {i} & object {j}
        i = hi(ii);
        j = hi(jj);
        cnt_total_pairwise_comparison = cnt_total_pairwise_comparison + 1;
        X = zeros(size(known,1),d+1);
        Y = zeros(size(known,1),1);
        for kkk = 1:size(known,1)
            X(kkk,:) = H(:, known(kkk,1), known(kkk,2))';
            Y(kkk) = Qh(known(kkk,1), known(kkk,2));
        end
        
        % add by gylai
        % assume that the label for Qh(i,j) is 1
        [w,fval,flag,accuracy_p]=linear_program_solver([ X; H(:,i,j)' ],[Y; 1],0);
        % TODO: this w may be useful
        % add by gylai
        % assume that the label for Qh(i,j) is -1
        [w,fval,flag,accuracy_n]=linear_program_solver([ X; H(:,i,j)' ],[Y; -1],0);
        
        if j == placeHolder
            secondPhase = 1;
        end
        
        if (accuracy_p < 100) && (accuracy_n < 100)
            fprintf('broke - maxiter is too low for convergence. restarting.\n');
            % add by gylai
            % 
            known = known(1:end-1,:);
        % not ambiguous: 
        elseif accuracy_p(1) < 100
            cnt_unambiguous_query = cnt_unambiguous_query + 1;
            Qh(i,j) = -1;
            Qh(j,i) = -Qh(i,j);
        % not ambiguous: 
        elseif accuracy_n(1) < 100
            cnt_unambiguous_query = cnt_unambiguous_query + 1;
            Qh(i,j) = 1;
            Qh(j,i) = -Qh(i,j);
        else
            % disp([ 'accuracy_p: ' num2str(accuracy_p, '%1.3g') '    accuracy_n: ' num2str(accuracy_n, '%1.3g') ]);
            % accuracy_p == 100 && accuracy_n == 100
            % If you are here, the query is Ambiguous.
            % add by gylai
            % ask the DM to provide her or his preference
            cnt_ambiguous_query = cnt_ambiguous_query + 1;
%             fprintf('<info> the query is Ambiguous, query No.%d\n', cnt_ambiguous_query);
            if R==0 % noiseless case
                known = [ known; [i j] ];
%                 fprintf('<info> number of query requests: %d\n', length(known));
                % artificial DM
                if norm(w_star-Xf(i,:)) < norm(w_star-Xf(j,:))
                    % object i ranks higher than object j
                    Qh(i, j) = 1;
                    % FIXME
                    % add by gylai
                else
                    % object j ranks higher than object i
                    Qh(i, j) = -1;
                    % add ends
                end
%                 reply = input(['Which comes first? [1 or 2]\n 1. ',objects{i},'\n 2. ',objects{j},'\n'], 's');
%                 if strcmp(reply,'1')
%                     Qh(i,j) = 1;
%                 else
%                     Qh(i,j) = -1;
%                 end
                Qh(j,i) = -Qh(i,j);
                Qk(i,j) = 1;
                Qk(j,i) = 1;
                
                % add by gylai
                % note that the value of accuracy must be 100 with label Qh(i, j)
                [w,fval,flag,accuracy]=linear_program_solver([ X; H(:, i, j)' ],[Y; Qh(i, j)],0);
                w_approximate = w(1 : end - 1) / w(end);
                w_approximate = w_approximate'
                % add by gylai
                % disp('============== w ==============');
                % FIXME
                % disp(size(w));
                % disp(w);
                % add ends
                
                err = 0;
                for jjj = 2:n
                    for iii = 1:jjj-1
                        % add by gylai
                        % error by counting differences of the binary relationships
                        % induced by w_star and w respectively
                        % FIXME: w-Xf(iii,:)
                        
                        if ((norm(w_star-Xf(iii,:))<norm(w_star-Xf(jjj,:))) && (norm(w_approximate-Xf(iii,:))>norm(w_approximate-Xf(jjj,:))) ...
                            || (norm(w_star-Xf(iii,:))>norm(w_star-Xf(jjj,:))) && (norm(w_approximate-Xf(iii,:))<norm(w_approximate-Xf(jjj,:))) ...
                            )
                            err = err + 1;
                        end
                    end
                end
                
                fprintf('the number of pairwise errors when length(known) is %d: %d\n', ...
                            length(known), err);
                error( length(known) ) = err;
                
%                 for j=2:n
            else % add by gylai: noise case
                [ y_a,Qk,Qh,flag ] = predictLabelNewest(i,j,Qh,Qk,H,known,secondPhase,objects,R);
                if ~flag
                    Qh(i,j) = y_a;
                    Qh(j,i) = -Qh(i,j);
                    known = [ known; [i j] ];
                    Qk(i,j) = 1;
                    Qk(j,i) = 1;
                else
                    if placeHolder == -1
                        placeHolder = j;
                    end
                    
                    Qh(i,j) = 0;
                    
                    if ~isempty(known)
                        known( (known(:,1)==i)&(known(:,2)==j),: ) = [];
                    end
                    
                    hi(jj:end) = [ hi(jj+1:end); hi(jj) ];
                end
            end
        end
        
        % add by gylai
        if Qh(i,j) == 1
            % hi(jj) == j
            cnt_total_pairwise_comparison = cnt_total_pairwise_comparison + length(below);
            cnt_unambiguous_query = cnt_unambiguous_query + length(below);
            Qh( hi( below ), hi(jj) ) = 1;
            Qh( hi( jj ), hi(below) ) = -1;
        % add by gylai
        elseif Qh(i,j) == -1
            cnt_total_pairwise_comparison = cnt_total_pairwise_comparison + length(above);
            cnt_unambiguous_query = cnt_unambiguous_query + length(above);
            Qh( hi(jj), hi( above ) ) = 1;
            Qh( hi( above ), hi( jj ) ) = -1;
        end
        goodInds = (Qh( hi(jj), hi(1:jj-1) )==0);
        
    end % [ while sum(goodInds)>0 ] ends
    Qh;
    % add by gylai
    % DEBUG
    Q_gold = zeros(n,n);
    for object_x = 2 : n
        for object_y = 1 : object_x - 1
            if norm(w_star-Xf(object_x, :)) < norm(w_star-Xf(object_y, :))
                Q_gold(object_x, object_y) = 1;
            else
                Q_gold(object_x, object_y) = -1;
            end
            Q_gold(object_y, object_x)=-Q_gold(object_x, object_y); 
        end
    end
    Q_gold;
    if isequal(Qh, Q_gold)
        fprintf('<info> exact identification of the ranking is guaranteed.\n');
        fprintf('<info> details: total number of pairwise comparisons: %d\n', cnt_total_pairwise_comparison);
        fprintf('<info> details: number of unambiguous pairwise comparisons: %d\n', cnt_unambiguous_query);
        fprintf('<info> details: number of ambiguous pairwise comparisons: %d\n', cnt_ambiguous_query);
        fprintf('<info> details: %% of queries requested from the DM: %.2f\n', 100.0 * cnt_ambiguous_query / cnt_total_pairwise_comparison);
    end
    
    w_star
    w_approximate

    % DEBUG ends
%     ranking = compare_sort(1:n,Qh);
%     fprintf('\nYour ranking from most to least preferred: \n')
%     for i = 1:n
%         disp(ranking(i));
%     end
end

% uncomment by gylai
% ranking = compare_sort(1:n,Qh);
% fprintf('\nYour ranking from most to least preferred: \n')
% for i = 1:n
%     fprintf([objects{ranking(i)},'\n'])
% end


