function winner=TournamentByRank(individual_1,individual_2)
    if individual_1.Rank<individual_2.Rank
        winner=individual_1;
    elseif individual_1.Rank>individual_2.Rank
        winner=individual_2;
    else
        if rand(1)<0.5
            winner=individual_1;
        else
            winner=individual_2;
        end
    end
end