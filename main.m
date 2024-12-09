% for marker on plot
dis_from_target1 = abs(92-empty_maxDIS1);
[min_dis1,mindex1] = min(dis_from_target1);


figure(6)
plot(changing_P_gauge, empty_maxDIS1)
title("Initial Pressure vs. Horizontal Distance")
xlabel("Initial Pressure (Pa)")
ylabel("Horizontal Distance (m)")
target = 92;
yline(target)
hold on;
scatter(changing_P_gauge(mindex1), target, "*", LineWidth=2)
