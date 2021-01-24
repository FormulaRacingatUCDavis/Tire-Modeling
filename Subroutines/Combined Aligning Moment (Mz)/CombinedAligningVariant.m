function [variant, Tire] = CombinedAligningVariant(Raw, x0, Tire)
%Solves the optimization variables for the ssz(i) parameters as mentioned
% in page 183 of Tyre and Vehicle Dynamics (3rd Edition) by Hans. B
%Pacejka

%% Solving Lateral Optimization Variables
rvy6 = optimvar( 'rvy6', 'Lowerbound', 0, 'Upperbound', 5 );
rvy5 = optimvar( 'rvy5', 'Lowerbound', 0, 'Upperbound', 5 );
rvy4 = optimvar( 'rvy4', 'Lowerbound', 0, 'Upperbound', 5 );
rvy3 = optimvar( 'rvy3', 'Lowerbound', 0, 'Upperbound', 5 );
rvy2 = optimvar( 'rvy2', 'Lowerbound', 0, 'Upperbound', 5 );
rvy1 = optimvar( 'rvy1', 'Lowerbound', 0, 'Upperbound', 5 );
rhy1 = optimvar( 'rhy1', 'Lowerbound', 0, 'Upperbound', 5 );
rhy2 = optimvar( 'rhy2', 'Lowerbound', 0, 'Upperbound', 5 );
rey1 = optimvar( 'rey1', 'Lowerbound', 0, 'Upperbound', 5 );
rey2 = optimvar( 'rey2', 'Lowerbound', 0, 'Upperbound', 5 );
rcy1 = optimvar( 'rcy1', 'Lowerbound', 0, 'Upperbound', 5 );
rby1= optimvar( 'rby1', 'Lowerbound', 0, 'Upperbound', 5 );
rby2= optimvar( 'rby2', 'Lowerbound', 0, 'Upperbound', 5 );
rby3 = optimvar( 'rby3', 'Lowerbound', 0, 'Upperbound', 5 );
rby4 = optimvar( 'rby4', 'Lowerbound', 0, 'Upperbound', 5 );

Obj = fcn2optimexpr(@ErrorMZ, rvy6, rvy5, rvy4, rvy3, rvy2, rvy1, rhy1, rhy2, ...
    rey1, rey2, rcy1, rby1, rby2, rby3, rby4);

Constr = optimineq(1);
Constr(1) = (rby1 + (rby4 * gamma_star^2)) * cos(atan(rby2(alpha_star - rby3 ...
    ))) > 0;
[Variant.Solution, Variant.Log] = Runfmincon(Obj, x0, Constr, 3);
Tire.r.v.y(6) = Variant.Solution.rvy6;
Tire.r.v.y(5) = Variant.Solution.rvy5;
Tire.r.v.y(4) = Variant.Solution.rvy4;
Tire.r.v.y(3) = Variant.Solution.rvy3;
Tire.r.v.y(2) = Variant.Solution.rvy2;
Tire.r.v.y(1) = Variant.Solution.rvy1;
Tire.r.h.y(1) = Variant.Solution.rhy1;
Tire.r.h.y(2) = Variant.Solution.rhy2;
Tire.r.e.y(1) = Variant.Solution.rey1;

%% Optimisation variables
ssz1 = optimvar( 'ssz1', 'Lowerbound', 0, 'Upperbound', 10 );
ssz2 = optimvar( 'ssz2', 'Lowerbound', 0, 'Upperbound', 10 );
ssz3 = optimvar( 'ssz3', 'Lowerbound', 0, 'Upperbound', 10 );
ssz4 = optimvar( 'ssz4', 'Lowerbound', 0, 'Upperbound', 10 );

%% Optimization Objective
Obj = fcn2optimexpr(@ErrorMs, ssz1, ssz2, ssz3, ssz4);

%% Solving Optimisation Problem
[Variant.Solution, Variant.Log] = fminsearch( Obj, x0);

%% Clearing Optimazation Problem
delete( findobj('Type', 'figure', 'Name', 'Optimization PlotFcns') );

%% Allocating Solution
Tire.s.s.z(1) = Variant.Solution.ssz1;
Tire.s.s.z(2) = Variant.Solution.ssz2;
Tire.s.s.z(3) = Variant.Solution.ssz3;
Tire.s.s.z(4) = Variant.Solution.ssz4;

%% Local Functions
function MeanAbsoluteError = ErrorMs(ssz1, ssz2,ssz3,ssz4)
ro = obj.parameters.Ro;

s = r0 * (ssz1 + (ssz2 * (Fy0./ (Tire.Fz0/Tire.Fz))) + (ssz3 + ssz4 * ...
    Tire.dfz) * gamma);
MeanAbsoluteError = mean(abs(Raw.Force - MZ));
end
    function [Solution, Log] = Runfmincon( Obj, x0, Constr, n )
        % Initializing Logs
        Log(n).x = [];
        Log(n).fval = [];
        
        % Creating Array of Initial Vectors (Latin Hypercube Sampling)
        if nargin > 3 && n > 1
            Jitter = lhsdesign( n-1, numel(fieldnames(x0)) ) - 0.5;
            
            Initial = (struct2array(x0)/2) .* Jitter + struct2array(x0);
            Initial(:, find( ~struct2array(x0) ) ) = ...
                Jitter( :, find( ~struct2array(x0) ) ) .* 0.2; %#ok<FNDSB>
            
            Variables = fieldnames( x0 );
            for ii = 1 : n-1
                for jj = 1 : numel( Variables )
                    x0(ii+1).(Variables{jj}) = Initial(ii, jj);
                end
            end
        end
        
        % Optimization Problem
        if isempty(Constr)
            Prob = optimproblem( 'Objective', Obj );
        else
            Prob = optimproblem( 'Objective', Obj, 'Constraints', Constr );
        end
        
        % Optimization Options
        Opts = optimoptions( 'fmincon', ...
            'Algorithm', 'sqp', ...
            'MaxFunctionEvaluations', 10000, ...
            'MaxIterations', 10000, ...
            'Display', 'off', ...
            'PlotFcn', {@optimplotx, @optimplotfval}, ...
            'OutputFcn', @OptimLogging );
        
        % Solving Optimization Problem(s)
        for ii = 1 : n
            try
                [Solution(ii), Feval(ii)] = solve( Prob, x0(ii), ...
                    'solver', 'fmincon', 'options', Opts ); %#ok<AGROW>
            catch
                continue
            end
        end

        % Selecting Optimal Solution
        [~, MinIdx] = min( Feval );
        Solution = Solution(MinIdx);
        
        function stop = OptimLogging( x, optimValues, state )
            stop = false;
            
            if strcmp(state, 'iter')
                Log(ii).x = [Log(ii).x x];
                Log(ii).fval = [Log(ii).fval optimValues.fval];
            end
        end
    end
    function MeanAbsoluteError = ErrorMZC(ssz1, ssz2,ssz3, ssz4, Dt, Ct, ...
            Bt, Dr, Cr, Br, R0, G_yk, Fy0, alpha_ r, alpha_t, Kxk, Kyalphaprime, kappa)
        alphateq = sqrt(alpha_t.^2 + ((Kxk./Kyalphaprime)^2 * kappa))*sgn(alpha_t);
        alphareq = sqrt(alpha_r.^2 + ((Kxk./Kyalphaprime)^2 * kappa))*sgn(alpha_r);
        
        s = R0 * (ssz1 + (ssz2*(Tire.Fy0./Tire.Fz0)) + (ssz3 + (ssz4 * dfz) * sin(gamma));
        b_r = Tire.Br * alphareq;
        b_t = Tire.Br * alphateq;
        
        Mzr = Dr * cos(Cr * atan(b_r));
        Fyp = Tire.G_yk * Fy0;
        t = Dt * cos(Ct * atan(b_t - Et *(b_t - atan(b_t)))) * (vcx./(Vc + 0.1));
        Mzp = -t. * Fyp;
        Mz = Mzp + Mzr + (s * Fx);
            
end
