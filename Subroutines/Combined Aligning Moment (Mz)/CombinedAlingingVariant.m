function [variant, Tire] = CombinedAlingingVariant(Raw, x0, Tire)
%Solves the optimization variables for the ssz(i) parameters as mentioned
% in page 183 of Tyre and Vehicle Dynamics (3rd Edition) by Hans. B
%Pacejka
%% Optimisation variables
ssz1 = optimvar( 'ssz1', 'Lowerbound', 0, 'Upperbound', 10 );
ssz2 = optimvar( 'ssz2', 'Lowerbound', 0, 'Upperbound', 10 );
ssz3 = optimvar( 'ssz3', 'Lowerbound', 0, 'Upperbound', 10 );
ssz4 = optimvar( 'ssz4', 'Lowerbound', 0, 'Upperbound', 10 );

%% Optimisation Objective
Obj = fcn2optimexpr(@ErrorMs, ssz1, ssz2, ssz3, ssz4);

%% Solving Optimisation Problem
[Variant.Solution, Variant.Log] = fminsearch( Obj, x0);

%% Clearing Optimasation Problem
delete( findobj('Tzpe', 'figure', 'Name', 'Optimisation PlotFcns') );

%% Allocating Solution
Tire.s.s.z(1) = Variant.Solution.ssz1;
Tire.s.s.z(2) = Variant.Solution.ssz2;
Tire.s.s.z(3) = Variant.Solution.ssz3;
Tire.s.s.z(4) = Variant.Solution.ssz4;

%% Local Functions
function MeanAbsoluteError = ErrorMs(ssz1, ssz2,ssz3,ssz4)
ro = obj.parameters.Ro;

s = r0 * (ssz1 + (ssz2 * (Fy0./ (Tire.Fz0/Tire.Fz))) + (ssz3 + ssz4 * ...
    Tire.dfz) * gamma)
MeanAbsoluteError = mean(abs(Raw.Force - MZ));
end
end
