%
% created by gylai on Mar 27, 2022
%

function best_solution=parallel_run()
    %% Problem Definition
    nObj=3;             % Number of Objective Functions
    nVar=7;             % Number of Decision Variables

    CostFunction=@(x) DTLZ1(x,nObj);      % Cost Function

    VarSize=[1 nVar];   % Size of Decision Variables Matrix

    VarMin= 0;          % Lower Bound of Variables
    VarMax= 1;          % Upper Bound of Variables

    % Number of Objective Functions
    nObj=numel(CostFunction(unifrnd(VarMin,VarMax,VarSize)));

    %% NSGA-II Parameters
    MaxIt=250;      % Maximum Number of Iterations
    nPop=100;        % Population Size

    %% Initialization
    % assumed reference point
    w_star=[ 0.5/3 0.5/3 0.5/3 ];
    interact_frequency=50;

    % approximated reference point via active ranking
    w_approximate=[];

    % individual
    empty_individual.Position=[];
    empty_individual.Cost=[];
    empty_individual.Rank=[];
    empty_individual.DominationSet=[];
    empty_individual.DominatedCount=[];
    empty_individual.CrowdingDistance=[];
    empty_individual.Preference=[];

    % linear program related
    A=[];
    b=[];
    d=nObj+1;
    f = [zeros(2*d,1);1];
    LB = [ zeros(2*d,1); -1];
    UB = inf(2*d+1,1);

    % active pairwise comparisons -> point p_active and its label y_active
    p_active=[];
    y_active=[];

    pop=repmat(empty_individual,nPop,1);
    for i=1:nPop

        pop(i).Position=unifrnd(VarMin,VarMax,VarSize);

        pop(i).Cost=CostFunction(pop(i).Position);

    end

    % Non-Dominated Sorting
    [pop, F]=NonDominatedSorting(pop);
    % 
    % % Calculate Crowding Distance
    % pop=CalcCrowdingDistance(pop,F);
    % 
    % % Sort Population
    % [pop, F]=SortPopulation(pop);

    %% PLAR Main Loop
    for it=1:MaxIt

        % add by gylai
        % active ranking
        if mod(it,interact_frequency)==0
            X=[];
            for i=1:nPop
                X=[X;pop(i).Cost'];
            end

            % IMPORTANT!
            % make sure there are no identical items
            X=unique(X,'rows');
            [x,y]=active_pairwise_comparisons(p_active,y_active,X,w_star);

            p_active=[p_active;x];
            y_active=[y_active;y];

            [A,b]=linear_constraints_accumulater(x,y,A,b);
            % disp('w_approximate by accumulater:');
            w_approximate=approximate_reference_point(A,b,f,LB,UB,d);
            % disp(w_approximate);

        end

        % Crossover
        pop_offspring=repmat(empty_individual,nPop/2,2);
        [v1 idx1] = sort(rand(nPop,1));
        [v2 idx2] = sort(rand(nPop,1));
        for k=1:nPop/2

            % FIXME
            p1=TournamentByRank(pop(idx1(2*k-1)),pop(idx1(2*k)));

            p2=TournamentByRank(pop(idx2(2*k-1)),pop(idx2(2*k)));

            [pop_offspring(k,1).Position, pop_offspring(k,2).Position]=SBX(p1.Position,p2.Position);

        end
        pop_offspring=pop_offspring(:);
        % Mutation
        for k=1:nPop

            pop_offspring(k).Position=PolynomialMutation(pop_offspring(k).Position);

            pop_offspring(k).Cost=CostFunction(pop_offspring(k).Position);

        end

        % Merge
        pop=[pop
             pop_offspring]; 

        % Non-Dominated Sorting
        [pop, F]=NonDominatedSorting(pop);

        % Calculate Crowding Distance
        % pop=CalcCrowdingDistance(pop,F);
        if it>=interact_frequency
            % Calculate Preference
            pop=CalcPreference(pop,w_approximate);
            % Sort Population via preference
            pop=SortPopulationViaPreference(pop);
        else
            % Calculate Crowding Distance
            pop=CalcCrowdingDistance(pop,F);
            % Sort Population
            [pop, F]=SortPopulation(pop);
        end


        % Truncate
        % mixed pop (parent & offspring) -> new pop
        pop=pop(1:nPop);

        % add by gylai
        % in the new population

        % Non-Dominated Sorting
        [pop, F]=NonDominatedSorting(pop);

        % Calculate Crowding Distance
        % pop=CalcCrowdingDistance(pop,F);

    %     % Calculate Preference
    %     pop=CalcPreference(pop,w_approximate);
    % 
    %     % Sort Population via preference
         % pop=SortPopulationViaPreference(pop);

        % Store F1
        % F1 -> nondominated solutions in the new population
        F1=pop(F{1});

        % Show Iteration Information
        % disp(['Iteration ' num2str(it) ': Number of F1 Members = ' num2str(numel(F1))]);

        % Plot F1 Costs
    %     figure(1);
    %     PlotCosts(F1);
    %     pause(0.01);

    end

    %% Results
    % Save F1 Costs
    % SaveCosts(F1);
    % Save the best one
    % SaveIndividual(F1(1));  
    best_solution=F1(1).Cost;
end
