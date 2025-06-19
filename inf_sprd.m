function chance = inf_sprd(d,is_infected)
% Chance of infection decreases exponentially with distance (d)
% Only infected trees can spread infection
    if is_infected ==1
        chance = exp(-(d+6)/2);
    else
        chance = 0;
    end
end
