function dielectric = leaf_vwc(t_inf,dielectric_healthy,dielectric_infected, time_to_die)
% Approximating that the decrease in vwc and dielectric is linear as tree
% becomes more and more sick.

    for i = 1:length(t_inf)
        if t_inf(i) <= time_to_die
            dielectric(i) = (dielectric_infected - dielectric_healthy) / time_to_die * t_inf(i) + dielectric_healthy;
        else
            dielectric(i) = dielectric_infected;
        end
    end
end
