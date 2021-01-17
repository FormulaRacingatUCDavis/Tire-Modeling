function [] = Debugging(Tire)
%% Parameters 
Pi = 11;
Slip = 0;
gamma = 2;
%% Dt and Et Graphing
% Code does not affect fit of tire. Only Estimating Dt/Et parameters
dFz = @(Fz)(Fz - Tire.Fzo) ./ Tire.Fzo;
dPi = @(Pi)(Pi - Tire.Pio) ./ Tire.Pio;
figure
subplot(2,1,1)
Dt = @(Fz)(Tire.Ro .*(Fz./Tire.Fzo)) .* (Tire.q.D.z(1) + Tire.q.D.z(2).*dFz(Fz) .*...
    (1 - Tire.p.P.z(1) .* dPi(Pi))) .* (1 + Tire.q.D.z(3) .* abs(gamma) + ...
    Tire.q.D.z(4) .* gamma.^2);
 fplot(Dt,[0,2500])
 title('Dt vs Load')
 Ct = 0.8;
 Bt = 2.5;
 subplot(2,1,2)
 Et = @(Fz)(Tire.q.E.z(1) + Tire.q.E.z(2).*dFz(Fz) + Tire.q.E.z(3) .* dFz(Fz).^2)...
     .*(1 + (Tire.q.E.z(4) + Tire.q.E.z(5).* gamma).*(2/pi).*atan(Bt .* Ct .* Slip));
 fplot(Et, [0,2500])
 title('Et vs Load')
 %}
 %% Plotting Bt, Br, Dr
 Bt = @(Fz) (Tire.q.B.z(1) + Tire.q.B.z(2) + Tire.q.B.z(3).*dFz(Fz).^2).*(1 + ...
     Tire.q.B.z(5).*abs(gamma) + Tire.q.B.z(6).*(gamma).^2); 
 
 Cy = Tire.p.C.y(1);
 
 Dy = @(Fz)(Tire.p.D.y(1) + Tire.p.D.y(2).*dFz(Fz)).*(1 + Tire.p.P.y(3).*dFz(Fz) + ...
     Tire.p.P.y(4).*dPi(Pi).^2).*(1 - Tire.p.D.y(3).* gamma.^2).* Fz;
 
 Kya = @(Fz) (Tire.p.K.y(1) .* Tire.Fzo) .* (1 + Tire.p.P.y(1).*dPi(Pi)) .* (1 - ...
     Tire.p.K.y(3) .* abs(gamma)) .* sin(Tire.p.K.y(4) .* atan((Fz./Tire.Fzo)./...
     ((Tire.p.K.y(2) + Tire.p.K.y(5).*gamma.^2).*(1 + Tire.p.P.y(2).*dPi(Pi)))));
 
 
By = @(Fz) Kya(Fz) ./ ( Cy.*Dy(Fz) );
 
Br = @(Fz) Tire.q.B.z(10) .* By(Fz) .* Cy;

Dr = @(Fz) Tire.Ro .* Fz .* ( ( Tire.q.D.z(6) + Tire.q.D.z(7).*dFz(Fz) ) + ...
     ( ( Tire.q.D.z(8) + Tire.q.D.z(9).*dFz(Fz) ) .* ( 1 + Tire.p.P.z(2).*dPi(Pi) ) + ...
     ( Tire.q.D.z(10) + Tire.q.D.z(11).*dFz(Fz) ) .* abs(gamma) ).*gamma );
 
figure
subplot(3,1,1)
fplot(Br, [0,2500])
title('Br vs Load')
subplot(3,1,2)
fplot(Dr, [0,2500])
title('Dr vs Load')
subplot(3,1,3)
fplot(Bt, [0,2500])
title('Bt vs Load')
 
 
 %% Plotting Trail function
 
Cy = Tire.p.C.y(1);
        
Dy = @(Fz) (Tire.p.D.y(1) + Tire.p.D.y(2).*dFz(Fz)) .* ...
            (1 + Tire.p.P.y(3).*dPi(Pi) + Tire.p.P.y(4).*dPi(Pi).^2) .* ...
            (1 - Tire.p.D.y(3).*gamma.^2).*Fz;
        
Kya = @(Fz) Tire.p.K.y(1) .* Tire.Fzo .* ( 1 + Tire.p.P.y(1).*dPi(Pi) ) .* ...
            ( 1 - Tire.p.K.y(3).*abs(gamma) ) .* sin( Tire.p.K.y(4) .* ...
            atan( (Fz./Tire.Fzo) ./ ...
            ( ( Tire.p.K.y(2) + Tire.p.K.y(5).*gamma.^2 ) .* ( 1 + Tire.p.P.y(2).*dPi(Pi) ) ) ) );
        
Kyg0 = @(Fz) Fz.*(Tire.p.K.y(6) + Tire.p.K.y(7).*dFz(Fz) ) .* (1 + Tire.p.P.y(5).*dPi(Pi) );
        
By = @(Fz) Kya(Fz) ./ ( Cy.*Dy(Fz) );
        
Vyg = @(Fz) Fz.*(Tire.p.V.y(3) + Tire.p.V.y(4).*dFz(Fz) ).*gamma;
        
Vy = @(Fz) Fz.*(Tire.p.V.y(1) + Tire.p.V.y(2).*dFz(Fz) ) + Vyg(Fz);
        
Hy = @(Fz) (Tire.p.H.y(1) + Tire.p.H.y(2).*dFz(Fz) ) .* ...
            (Kyg0(Fz).*gamma - Vyg(Fz) ) ./ Kya(Fz);
        
Ey = @(Fz, Hy) ( Tire.p.E.y(1) + Tire.p.E.y(2).*dFz(Fz) ) .* ...
            ( 1 + Tire.p.E.y(5).*gamma.^2 - ...
            ( Tire.p.E.y(3) + Tire.p.E.y(4).*gamma ).*sign(Slip + Hy) );
        
Ht = @(Fz) Tire.q.H.z(1) + Tire.q.H.z(2).*dFz(Fz) + ...
    (Tire.q.H.z(3) + Tire.q.H.z(4).*dFz(Fz)).*gamma;

Fyo = @(Fz) Dy(Fz) .* ...
    sin( Cy .* atan( (1-Ey(Fz, Hy(Fz) )) .* ...
    By(Fz).*(Slip + Hy(Fz) ) + ...
    Ey(Fz, Hy(Fz) ).*atan( ...
    By(Fz).*(Slip + Hy(Fz) ) ) ) ) + Vy(Fz);

t0 = @(Fz) Dt(Fz) .* cos( Ct .* ...
    atan( (1-Et(Fz)).*(Bt(Fz).*(Slip+Ht(Fz)) + ...
    Et(Fz).*atan( Bt(Fz).*(Slip+Ht(Fz)) ) ) ) ) .* cosd( Slip );



Mzop = @(Fz) -t0(Fz) .* Fyo(Fz);

figure
subplot(2,1,1)
fplot(Mzop, [0,2500])
title("Trail vs Load")

%%Plotting residual function
%Br and Dr defined above
        
Cr = 1;
               
Hf = @(Fz) Hy(Fz) + Vy(Fz) ./ Kya(Fz);



Mzro = @(Fz) Dr(Fz) .* cos( Cr .* ...
    atan( Br(Fz).*(Slip+Hf(Fz)) ) ) .* cosd( Slip );

subplot(2,1,2)
fplot(Mzro, [0,2500])
title("Residual vs Load")
end