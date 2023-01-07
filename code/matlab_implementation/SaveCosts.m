function SaveCosts(pop)
    Costs=[pop.Cost]';
    [nInd, nObj]=size(Costs);
    fid=fopen('final_population.txt', 'w');
    for i=1:nInd
        for j=1:nObj
            fprintf(fid, '%f\t', Costs(i, j));
        end
        fprintf(fid, '\n');
    end
    fclose(fid);
    

end