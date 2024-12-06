clear;
clc;
close all;

%load const struct
const = getConst();

%calculate initial velocity
v0x = const.v0 * cos(const.theta0);
v0z = const.v0 * sin(const.theta0);

tspan = [0, 5];
%initial conditions
%x pos, x vel, z pos, z vel, rocket mass, air volume, air mass
statevector0 = [const.x0, v0x, const.z0, v0z, const.mr0, const.vair0, const.mair0];

%variable parameters
vars = getVars();

xdistance = zeros(1, length(vars.theta));
theta0 = const.theta0;
for i = 1:length(vars.theta)
    const.theta0 = vars.theta(i);
    [t, statevector] = ode45(@(t,statevector) newState(t,statevector,const),tspan,statevector0);
    xdistance(i) = statevector(end,1);
end
const.theta0 = theta0;

figure;
plot(linspace(0, 90, 1000), xdistance);
xlabel('Launch Angle');
ylabel('X Distance');

xdistance = zeros(1, length(vars.cd));
cd = const.cd;
for i = 1:length(vars.cd)
    const.cd = vars.cd(i);
    [t, statevector] = ode45(@(t,statevector) newState(t,statevector,const),tspan,statevector0);
    xdistance(i) = statevector(end,1);
end
const.cd = cd;

figure;
plot(vars.cd, xdistance);
xlabel('Cd');
ylabel('X Distance');

function dstatevector = newState(~, statevector, const)
%get state values
x = statevector(1);
vx = statevector(2);
z = statevector(3);
vz = statevector(4);
mr = statevector(5);
vair = statevector(6);
mair = statevector(7);

%calculate air pressure using pend
pend = const.p0 * (const.v0air / vair) ^ const.gamma;
pair = pend * (mair / const.mair0) ^ const.gamma;

%heading vector calculation
vvec = [vx, vz];
vnorm = sqrt(vx^2 + vz^2);
distance = sqrt((x-const.x0)^2 + (z-const.z0)^2);
%if the rocket is still on the stand, define heading manually using launch
%angle. Otherwise define dynamically using current velocity
if distance < const.lstand
    heading = [cos(const.theta0), sin(const.theta0)];
else
    heading = vvec ./ vnorm;
end

%drag is the same in all phases
drag = 0.5 * const.cd * const.areab * const.rhoair * vnorm^2 .*heading;

if vair < const.vb %phase 1: water propulsion
    pair = const.p0 * (const.v0air / vair) ^ const.gamma;
    vexit = sqrt(2 * (pair - const.pa) / const.rhowater);
    dmwater = const.cdis * const.rhowater * const.areat * vexit;
    dmtotal = -dmwater;
    thrust = dmwater * vexit;
    dvair = const.cdis * const.areat * vexit;
    dmair = 0;
elseif pair > const.pa %phase 2: air propulsion
    pcrit = pair * (2 / (const.gamma + 1)) ^ (const.gamma / (const.gamma - 1));
    currdensity = mair / const.vb;
    tair = pair / (currdensity * const.rair);
    %define exit conditions based on choked vs non-choked
    if pcrit > const.pa %choked flow
        pexit = pcrit;
        texit = tair * (2 / (const.gamma + 1));
        vexit = sqrt(const.gamma * const.rair * texit);
        rhoexit = pexit / (const.rair * texit);
    else %non choked flow
        pexit = const.pa;
        machexit = sqrt(((pair / const.pa) ^ ((const.gamma - 1) / const.gamma) - 1) / ((const.gamma - 1) / 2));
        texit = tair / (1 + ((const.gamma - 1) / 2) * machexit^2);
        rhoexit = pexit / (const.rair * texit);
        vexit = machexit * sqrt(const.gamma * const.rair * texit);
    end
    %use exit conditions to define statevector result
    dmair = const.cdis * rhoexit * const.areat * vexit;
    dmtotal = -dmair;
    thrust = dmair * vexit + (pexit - const.pa) * const.areat;
    dvair = 0;
else %phase 3: no propulsion/ballistic
    %if the rocket hit the ground, no change in anything
    if z <= 0
        dstatevector = [0; 0; 0; 0; 0; 0; 0];
        return;
    end
    thrust = 0;
    dmtotal = 0;
    dvair = 0;
    dmair = 0;
end

%convert thrust to a vector
thrust = thrust .* heading;

%air is flowing out, not in
dmair = -dmair;

Fnet = thrust - drag - mr.*[0, const.g];
anet = Fnet ./ mr;
dstatevector = [vx; anet(1); vz; anet(2); dmtotal; dvair; dmair];
end

function const = getConst()
    const.g = 9.81; %m/s2
    const.cdis = 0.78;
    const.rhoair = 0.961; %kg/m3
    const.vb = 0.002; %m3
    const.pa = 12.1 * 6894.76; %psia to Pa
    const.gamma = 1.4;
    const.rhowater = 1000; %kg/m3
    const.de = 2.1 / 100; %cm to m
    const.db = 10.5 / 100; %cm to m
    const.rair = 287; %J/kgK
    const.mb = 0.15; %kg
    const.cd = 0.425;
    const.p0 = 48 * 6894.76 + const.pa; %psig to Pa
    const.v0water = 0.0005; %m3
    const.v0air = const.vb - const.v0water; %m3
    const.t0 = 310; %K
    const.v0 = 0; %m/s
    const.theta0 = deg2rad(40); %deg to radians
    const.x0 = 0; %m
    const.z0 = 0.25; %m
    const.lstand = 0.5; %m
    const.mairideal = const.vb * const.rhoair; %kg
    const.areab = pi * (const.db / 2)^2; %m2
    const.areat = pi * (const.de / 2)^2; %m2
    const.vair0 = const.vb - const.v0water; %m3
    const.mair0 = (const.p0 * const.vair0) / (const.rair * const.t0); %kg
    const.mwater0 = const.v0water * const.rhowater; %kg
    const.mr0 = const.mb + const.mair0 + const.mwater0; %kg
end

function vars = getVars()
    vars.theta = deg2rad(linspace(0, 90, 1000));
    vars.cd = linspace(0.1, 1.5, 100);
end