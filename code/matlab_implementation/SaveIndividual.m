function SaveIndividual(ind)
    Cost=[ind.Cost]';
    [nInd, nObj]=size(Cost);
    fid=fopen('final_individual.txt', 'w');
    for i=1:nInd
        for j=1:nObj
            fprintf(fid, '%f\n', Cost(i, j));
        end
        fprintf(fid, '\n');
    end
    fclose(fid);
end