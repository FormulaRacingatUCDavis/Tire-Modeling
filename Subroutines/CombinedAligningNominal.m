function [Nominal] = CombinedAligningNominal(Raw)

%% Defining Initial Point
[x0, lb, ub] = MZInitialBounds(Raw.Slip, Raw.Force, Raw.Kappa);
s = optimvar("s", "LowerBound", lb.C, "UpperBound", ub.s);
%% Defining Optimization Objective
Obj = fcn2optimexpr( @ErrorMZNom, s, [Raw.Slip; Raw.Force; Raw.Kappa]);

%% Defining Optimization Problem
Prob = optimproblem("Objective", Obj);

%% Defining Optimization Options
Opts = optimoptions("fminsearch", "MaxFunctionEvaluations", 2000, ...
    "MaxIterations", 2000, "Display", "off");

%% Solving Optimization Problem
Solution = solve(Prob, x0, "solver", "fminsearch", "options", Opts);

%% Setting up for Lateral Coefficient Solving
alpha_star = tan(gamma) * sgn(Vcx);
dfz = (Fz - Tire.Fz0)./Tire.Fz0;
D_vyk = mu_y * Fz * (rvy1 + (rvy2 * dfz)) + (rvy3 * sin(gamma)) * ...
    cos(atan(rvy4 * alpha_star));
S_vyk = D_vyk.* sin(rvy5 * atan(rvy6 * kappa));
S_hyk = rhy1 + (rhy2 * dfz);
E_yk = rey1 + (rey2 * dfz);
C_yk = rcy1;
B_yk = (rby1 + (rby4 * sin(gamma).^2)) * cos(atan(rby2 * (alpha_star - rby3)));

K_s = kappa + S_hyk;
B_k = B_yk * K_s;
B_s = B_yk * S_hyk;

G_yk0 = cos(C_yk * atan(B_s - E_yk(B_s - atan(B_s))));
G_yk = cos(C_yk * atan(B_k - E_yk * (B_k - atan(B_k))))./G_yk0;
Fy = (G_yk.Fy0) + S_vyk;

%% Setting up for Longitudinal Coefficient Solving
alpha_s = alpha_star + S_hxalpha;
b_k = B_xalpha * alpha_s;
b_s = B_xalpha * S_hxalpha;
G_xalpha0 = cos(C_xalpha * atan(b_s - E_xalpha*(b_s - atan(b_s))));
G_xalpha = cos(c_xalpha * atan(b_k - E_xalpha*(b_k - atan(b_k))))./G_xalpha0;
Fx = G_xalpha * Fx0;



%% Getting Coefficients for Combined Lateral
    function [Fy_mnc, Fy_mnc_eq] = Comstock_Lateral(Kappa, Slip)
        Fy_mnc = abs( Fx0 .* Fy0 ./ sqrt( (Kappa - Kappa0).^2 .* Fy0.^2 + Fx0.^2 .*...
            tan(Alpha - Alpha0).^2 ) .* sqrt( (1 - abs(Kappa - Kappa0)).^2 ...
            .* cos(Alpha - Alpha0).^2 .* Fy0.^2 +...
            sin(Alpha - Alpha0).^2 .* Kxk.^2 ) ./ ...
            ( Kxk .* cos(Alpha - Alpha0)) ) .* sign(Fy0);
        Fy_mnc_eq = [];
    end

Lat_Constr = @Comstock_lateral;
Constr = optimineq(3);
Constr(1) = G_yk > 0;
Constr(2) = B_yk > 0;
Constr(3) = Fy == Fy_mnc;
Obj = fcn2optimexpr(D_yvk, S_vyk, S_hyk, E_yk, C_yk, B_yk, G_yk0, G_yk, Fy);
[Solution] = Runfmincon(Obj, Constr, 3, Lat_Constr);

%% Getting Coefficients for Combined Longitudinal
    function [Fx_mnc, Fx_mnc_eq] = Comstock_Longitudinal(Kappa, Slip)
        Fx_mnc =  abs(Fx0 .* Fy0 ./ sqrt( (Kappa - Kappa0).^2 .* Fy0.^2 + ...
            Fx0.^2 .* tan(Alpha - Alpha0).^2 ) .*sqrt( (Kappa - Kappa0).^2 ...
            .* Kya.^2 + (1 - abs(Kappa - Kappa0)).^2 .* cos(Alpha - Alpha0).^2 ...
            .* Fx0.^2 ) ./Kya) .* sign(Fx0);
        Fx_mnc_eq = [];
    end
Long_Constr = @Comstock_Longitudinal;
Constr_x = optimineq(4);
Constr_x(1) = G_xalpha > 0;
Constr_x(2) = B_xalpha > 0;
Constr_x(3) = E_xalpha <= 0.99;
Constr_x(4) = Fx_mnc == Fx;
Obj_x = fcn2optimexpr(G_xalpha, G_xalpha0, C_xalpha, E_xalpha, B_xalpha, ...
    S_hxalpha, Fx);
[Solution] = Runfmincon(Obj_x, Constr_x, 3, Long_Constr);
%% Allocating Solution
Nominal.D_yvk = Solution.D_vyk;
Nominal.S_vyk = Solution.S_vyk;
Nominal.S_hyk = Solution.S_hyk;
Nominal.E_yk = Solution.E_yk;
Nominal.C_yk = Solution.C_yk;
Nominal.B_yk = Solution.B_yk;
Nominal.G_yk0 = Solution.G_yk0;
Nominal.G_yk = Solution.G_yk;

Nominal.G_xalpha = Solution.G_xalpha;
Nominal.G_xalpha0 = Solution.G_xalpha0;
Nominal.C_xalpha = Solution.C_xalpha;
Nominal.E_xalpha = Solution.E_xalpha;
Nominal.B_xalpha = Solution.B_xalpha;
Nominal.S_hxalpha = Solution.S_hxalpha;


Nominal.s0 = x0.s;

Nominal.s = Solution.s;

Nominal.Residual = MZNomResidual(Nominal.s, [Raw.Slip; Raw.Force] );

%% Local Functions
    %function [x0, lb, ub] = MZInitialBounds(SlipAngle, SlipRatio, Aligning Moment)
        %[~, MaxIdx] = maxk(AligningMoment, 5);
        %[~, MinIdx] = mink(AligningMoment, 5);
        
        %MaxIdx = round(mean(MaxIdx));
       % MinIdx = round(mean(MinIdx));
        
       % LinearFit = fitlm(SlipAngle( (SlipAngle > -1) & (SlipAngle<1)), ...
            %AligningMoment( (SlipAngle >-1) & (SlipAngle<1)));
        %% Bounding
       % lb.s = 0.1;
        
       % ub.s = inf;
    end
    function MeanSquareError = ErrorMZNom(s, Data)
        %Magic Formula Evaluation
    end
end
