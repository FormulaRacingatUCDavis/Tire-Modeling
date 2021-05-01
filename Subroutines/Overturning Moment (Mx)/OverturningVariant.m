function Tire = OverturningVariant(Raw, Tire)
%% Overturning Variant - Fits p,q Parameters for Overturning Equations
% Performs Fitting Process on q and p parameters for Overturning Moment
% referenced in Pacejka's "Tire and Vehicle Dynamics" 3rd Edition in
% Section 4.3.2. This will use a contrained fitting process in finding a
% local minima using the fmincon function (although no contraints are used)
%
% Inputs:
%   Raw     - Allocated Data
%   Tire    - Tire Model
%
% Outputs:
%   Variant - Parameters for Overturning Equation
%   Tire    - Updated Tire Model
%
% Authors:
% Carlos Lopez       (calopez@ucdavis.edu)        [Dec 2020 - June 2022]
% Blake Christierson (bechristierson@ucdavis.edu) [May 2021 - June 2021]
%
% Last Updated: 01-May-2021

%% Optimization Variables
qsx1  = optimvar( 'qsx1' , 'Lowerbound',  0   , 'Upperbound',  0    );
qsx2  = optimvar( 'qsx2' , 'Lowerbound',- 5   , 'Upperbound',  5    );
qsx3  = optimvar( 'qsx3' , 'Lowerbound',  0   , 'Upperbound',  1    );
qsx4  = optimvar( 'qsx4' , 'Lowerbound',- 5   , 'Upperbound',  5    );
qsx5  = optimvar( 'qsx5' , 'Lowerbound',  0   , 'Upperbound',  5    );
qsx6  = optimvar( 'qsx6' , 'Lowerbound',  0   , 'Upperbound',  5    );
qsx7  = optimvar( 'qsx7' , 'Lowerbound',- 5   , 'Upperbound',  5    );
qsx8  = optimvar( 'qsx8' , 'Lowerbound',- 5   , 'Upperbound',  5    );
qsx9  = optimvar( 'qsx9' , 'Lowerbound',  0   , 'Upperbound',  5    );
qsx10 = optimvar( 'qsx10', 'Lowerbound',- 5   , 'Upperbound',  5    );
qsx11 = optimvar( 'qsx11', 'Lowerbound',  0   , 'Upperbound',  5    );

ppmx1 = optimvar( 'ppmx1', 'Lowerbound',- 5   , 'Upperbound',  5    );

%% Optimization Initialization
x0.qsx1  = 0;
x0.qsx2  = 0;
x0.qsx3  = 0.2;
x0.qsx4  = 0;
x0.qsx5  = 0.01;
x0.qsx6  = 0.01;
x0.qsx7  = 0;
x0.qsx8  = 0.01;
x0.qsx9  = 0.01;
x0.qsx10 = 0;
x0.qsx11 = 0.01;

x0.ppmx1 = 0;

%% Optimization Objective 
Obj = fcn2optimexpr(@ErrorMx, qsx1,qsx2,qsx3,qsx4,qsx5,qsx6,qsx7, ...
    qsx8, qsx9, qsx10, qsx11, ppmx1);

%% Optimization Constraint
Constr = [];

%% Solving Optimization Problem
[Variant.Solution, Variant.Log] = Runfmincon( Obj, x0, Constr, 1);

%% Clearing Optimization Figure
delete( findobj('Type', 'figure', 'Name', 'optimization PlotFcns') );

%% Allocating Solution
Tire.Pacejka.q.s.x(1)  = Variant.Solution.qsx1  ;
Tire.Pacejka.q.s.x(2)  = Variant.Solution.qsx2  ;
Tire.Pacejka.q.s.x(3)  = Variant.Solution.qsx3  ;
Tire.Pacejka.q.s.x(4)  = Variant.Solution.qsx4  ;
Tire.Pacejka.q.s.x(5)  = Variant.Solution.qsx5  ;
Tire.Pacejka.q.s.x(6)  = Variant.Solution.qsx6  ;
Tire.Pacejka.q.s.x(7)  = Variant.Solution.qsx7  ;
Tire.Pacejka.q.s.x(8)  = Variant.Solution.qsx8  ;
Tire.Pacejka.q.s.x(9)  = Variant.Solution.qsx9  ;
Tire.Pacejka.q.s.x(10) = Variant.Solution.qsx10 ;
Tire.Pacejka.q.s.x(11) = Variant.Solution.qsx11 ;

Tire.Pacejka.p.P.Mx(1) = Variant.Solution.ppmx1 ;

%% Local Functions 
function RMSE = ErrorMx(qsx1, qsx2, qsx3, qsx4, qsx5, qsx6, qsx7, ...
        qsx8, qsx9, qsx10, qsx11, ppmx1)

    [~, Fy, ~, ~, ~] = ContactPatchLoads(Tire, rad2deg([Raw.Alpha]), [Raw.Kappa], ...
        [Raw.Load], [Raw.Pressure], [Raw.Inclination], 10, 1, ...
        struct('Pure', 'Pacejka', 'Combined', 'MNC'));

    Mx = Tire.Pacejka.Ro .* [Raw.Load] .* (qsx1 - qsx2 .* [Raw.Inclination] .* ...
        (1 + ppmx1 .* [Raw.dPi]) + qsx3 .* (Fy ./ Tire.Pacejka.Fzo) + ...
        qsx4 .* cos( qsx5 .* (atan( qsx6 .* ([Raw.Load] ./ Tire.Pacejka.Fzo) )).^2 ) .* ...
        sin( qsx7 .* [Raw.Inclination] + qsx8 .* atan( qsx9 .* (Fy ./ Tire.Pacejka.Fzo) ) ) + ...
        qsx10 .* atan( qsx11 .* ([Raw.Load] ./ Tire.Pacejka.Fzo) ) .* [Raw.Inclination] );

    RMSE = sqrt( mean( ([Raw.Moment] - Mx).^2) );
end

end