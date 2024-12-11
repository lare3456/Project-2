% changing initial pressure
changing_P_gauge = linspace(300000,10^6,173);
empty_maxDIS1 = zeros(size(changing_P_gauge));
empty_maxHEIGHT1 = zeros(size(changing_P_gauge));

for i = 1:length(changing_P_gauge)

    const.P_gauge = changing_P_gauge(i);
    const.pressure_i = const.P_atm + const.P_gauge;
    [t,statevector] = ode45(@(t,statevector) phase_funct(t,statevector,const),tspan,statevector_i);
    empty_maxDIS1(i) = max(statevector(:,1));
    empty_maxHEIGHT1(i) = max(statevector(:,3));

end

const = getConst();

% for marker on plot
dis_from_target1 = abs(92-empty_maxDIS1);
[min_dis1,mindex1] = min(dis_from_target1);




% changing initial volume of water
changing_volume_water = linspace(0,0.002,173);
empty_maxDIS2 = zeros(size(changing_volume_water));
empty_maxHEIGHT2 = zeros(size(changing_volume_water));

for i = 1:length(changing_volume_water)

    const.ViH2O = changing_volume_water(i);
    [t,statevector] = ode45(@(t,statevector) phase_funct(t,statevector,const),tspan,statevector_i);
    empty_maxDIS2(i) = max(statevector(:,1));
    empty_maxHEIGHT2(i) = max(statevector(:,3));

end

% for marker
dis_from_target2 = abs(92-empty_maxDIS2);
[min_dis2,mindex2] = min(dis_from_target2);

% volume of water for minimum distance from target
valVol_min_dis = changing_volume_water(1);
disp('Volume of water to hit target: ')
disp(valVol_min_dis)

valPress_min_dis = changing_P_gauge(83);
disp('Pressure to hit target: ')
disp(valPress_min_dis)
