function pop=SortPopulationViaPreference(pop)
    % Sort Based on Preference
    [~, CDSO]=sort([pop.Preference],'ascend');
    pop=pop(CDSO);
   
    % Sort Based on Rank
    [~, RSO]=sort([pop.Rank]);
    pop=pop(RSO);
   
end