% do polynomial mutation on decision variables x
function y=PolynomialMutation(x)
    nVar=numel(x);
    muatation_pro=1.0/nVar;
    eta_m=20.0;
    x_lowerbound=0;
    x_upperbound=1;
    for i=1:nVar
        x_i=x(i);
        if rand(1)<muatation_pro
            delta1 = (x_i - x_lowerbound) / (x_upperbound - x_lowerbound);
            delta2 = (x_upperbound - x_i) / (x_upperbound - x_lowerbound);
            rnd = rand(1);
            mut_pow = 1.0 / (eta_m + 1.0);
            if rnd<0.5
                xy = 1.0 - delta1;
                val = 2.0 * rnd + (1.0 - 2.0 * rnd) * (power(xy, (eta_m + 1.0)));
                deltaq = power(val, mut_pow) - 1.0;
            else
                xy = 1.0 - delta2;
                val = 2.0 * (1.0 - rnd) + 2.0 * (rnd - 0.5) * (power(xy, (eta_m + 1.0)));
                deltaq = 1.0 - (power(val, mut_pow));
            end
            x_i = x_i + deltaq * (x_upperbound - x_lowerbound);
        end
        if x_i<x_lowerbound || x_i>x_upperbound
            x_i=rand(1);
        end
        y(i)=x_i;
    end
end