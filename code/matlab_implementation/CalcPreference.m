%
% created by gylai on Jan 5, 2022
%

function pop=CalcPreference(pop,reference_point)
        num_pop=numel(pop);
        dimensionality=numel(reference_point);
        for i=1:num_pop
            dis=0;
            for j=1:dimensionality
                dis=dis+power(pop(i).Cost(j)-reference_point(j),2);
            end
            dis=sqrt(dis);
            pop(i).Preference=dis;
        end
end