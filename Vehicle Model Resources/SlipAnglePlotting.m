clear; close all; clc;

Fz_range = linspace(800, 2000, 50);
SA_range = linspace(-15, 15, 50);
SR = 0;
Pressure    = 70;
inc = linspace(0, 5, 6);
idx = 1;
V = 30;
load("Hoosier_R20_16(18)x75(60)-10x8(7).mat")
Model = struct( 'Pure', 'Pacejka', 'Combined', 'MNC' );

[Fz, SA] = meshgrid(Fz_range, SA_range);

tiledlayout(2, 3)

for i=1:6
nexttile

[Fx, Fy, Mz, Mx, My] = ContactPatchLoads(Tire, SA, SR, Fz, ...
    Pressure, inc(i), V, idx, Model);

surf(Fz, SA, Fy)
xlabel("Fz (N)")
ylabel("Slip angle (deg)")
zlabel("Fy (N)")
title("inclination: " + inc(i) + "deg")

end