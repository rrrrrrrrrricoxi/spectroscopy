clear all
clc
close all

omega_IR = 2700:3900; % cm^-1
omega_VIS = 12500;

[n_saph_fcn, k_saph_fcn, n_gold_fcn, k_gold_fcn, n_water_fcn, k_water_fcn] = get_refractive_indices();

G = [22, 43, 60, 75];

for jj = 1:length(G)

    % incidence angles onto substrate (travelling in air)
    gamma = deg2rad(G(jj));
    % theta0_IR  = 45;
    % theta0_VIS = 45;
    % theta0_SFG = 45;

    d = 5; % film thickness, nm

    Lxx = zeros(1,length(omega_IR));
    for j = 1:length(omega_IR)

        w_IR = omega_IR(j);   lambdaIR_um  = lambd(w_IR) / 1000;
        w_VIS = 12500;        lambdaVIS_um = lambd(w_VIS) / 1000;
        w_SFG = w_IR + w_VIS; lambdaSFG_um = lambd(w_SFG) / 1000;

        n_saph_IR  = n_saph_fcn(lambdaIR_um);
        n_saph_VIS = n_saph_fcn(lambdaVIS_um);
        n_saph_SFG = n_saph_fcn(lambdaSFG_um);
        k_saph_IR  = 0.0;
        k_saph_VIS = 0.0;
        k_saph_SFG = 0.0;

        n_gold_IR  = n_gold_fcn(lambdaIR_um);
        k_gold_IR  = k_gold_fcn(lambdaIR_um);
        n_gold_VIS = n_gold_fcn(lambdaVIS_um);
        k_gold_VIS = k_gold_fcn(lambdaVIS_um);
        n_gold_SFG = n_gold_fcn(lambdaSFG_um);
        k_gold_SFG = k_gold_fcn(lambdaSFG_um);

        n_water_IR  = n_water_fcn(w_IR);
        k_water_IR  = k_water_fcn(w_IR);
        n_water_VIS = n_water_fcn(w_VIS);
        k_water_VIS = k_water_fcn(w_VIS);
        n_water_SFG = n_water_fcn(w_SFG);
        k_water_SFG = k_water_fcn(w_SFG);

        n_IR  = [1.0, n_saph_IR  + 1i*k_saph_IR,   n_gold_IR  + 1i*k_gold_IR,  n_water_IR  + 1i*k_water_IR,  1.0];
        n_VIS = [1.0, n_saph_VIS + 1i*k_gold_VIS,  n_gold_VIS + 1i*k_gold_VIS, n_water_VIS + 1i*k_water_VIS, 1.0];
        n_SFG = [1.0, n_saph_SFG + 1i*k_water_SFG, n_gold_SFG + 1i*k_gold_SFG, n_water_SFG + 1i*k_water_SFG, 1.0];

        n_IR = n_IR(2:4);
        n_VIS = n_VIS(2:4);
        n_SFG = n_SFG(2:4);

        theta_IR = zeros(1,3); theta_VIS = zeros(1,3); theta_SFG = zeros(1,3);
        sin_IR = zeros(1,3);   sin_VIS = zeros(1,3);   sin_SFG = zeros(1,3);
        cos_IR = zeros(1,3);   cos_VIS = zeros(1,3);   cos_SFG = zeros(1,3);
        tan_IR = zeros(1,3);   tan_VIS = zeros(1,3);   tan_SFG = zeros(1,3);

        [theta_IR(1),sin_IR(1),cos_IR(1),tan_IR(1)] = get_theta_1(n_IR(1),gamma);
        [theta_IR(2),sin_IR(2),cos_IR(2),tan_IR(2)] = snell(n_IR(1),n_IR(2),theta_IR(1));
        [theta_IR(3),sin_IR(3),cos_IR(3),tan_IR(3)] = snell(n_IR(2),n_IR(3),theta_IR(2));

        [theta_VIS(1),sin_VIS(1),cos_VIS(1),tan_VIS(1)] = get_theta_1(n_VIS(1),gamma);
        [theta_VIS(2),sin_VIS(2),cos_VIS(2),tan_VIS(2)] = snell(n_VIS(1),n_VIS(2),theta_VIS(1));
        [theta_VIS(3),sin_VIS(3),cos_VIS(3),tan_VIS(3)] = snell(n_VIS(2),n_VIS(3),theta_VIS(2));

        [theta_SFG(2),sin_SFG(2),cos_SFG(2),tan_SFG(2)] = get_sfg_angle(w_IR,w_VIS,theta_IR(1),theta_VIS(1));
        [theta_SFG(3),sin_SFG(3),cos_SFG(3),tan_SFG(3)] = snell(n_SFG(2),n_SFG(3),theta_SFG(2));

        %[theta_IR, sin_IR, cos_IR, tan_IR, n_IR]      = get_snell(n_IR,  theta0_IR);
        %[theta_VIS, sin_VIS, cos_VIS, tan_VIS, n_VIS] = get_snell(n_VIS, theta0_VIS);
        %[theta_SFG, sin_SFG, cos_SFG, tan_SFG, n_SFG] = get_snell(n_SFG, theta0_SFG);

        % Fresnel reflection/transmission coefficients (p polarization)
        % for p polarization, n2/n1*t = r + 1
        r12p = (n_IR(2)*cos_IR(1) - n_IR(1)*cos_IR(2))/(n_IR(1)*cos_IR(2) + n_IR(2)*cos_IR(1));
        t12p = (2*n_IR(1)*cos_IR(1))/(n_IR(1)*cos_IR(2) + n_IR(2)*cos_IR(1));
        r23p = (n_IR(3)*cos_IR(2) - n_IR(2)*cos_IR(3))/(n_IR(2)*cos_IR(3) + n_IR(3)*cos_IR(2));
        t23p = (2*n_IR(2)*cos_IR(2))/(n_IR(2)*cos_IR(3) + n_IR(3)*cos_IR(2));

        delta_IR = @(x) 2*pi*n_IR(2)*d/(lambd(x)*cos_IR(2)) - ...
                        2*pi*n_IR(1)*d/lambd(omega_VIS) * (tan_IR(2) + tan_SFG(2)) * sin_IR(1);

        beta = @(x) 2*pi/lambd(x) * d * n_IR(2) * cos_IR(2);

        delta = delta_IR(w_IR);
        b = beta(w_IR);
        Lxx(j) = exp(1i*delta) * (t12p/(1+r12p*r23p*exp(2i*b))) * (1-r23p) * cos_IR(2) / cos_IR(1);

    end

    subplot(2,2,jj)
    plot(omega_IR,abs(Lxx).^2,'Linewidth',2)
    xlabel('\omega_{IR} / cm^{-1}')
    ylabel('|L_{xx}^{II}|^2')
    axis([-inf,inf,0,4])
    title(sprintf('gamma = %4.1f',G(jj)))
    set(gca,'Fontsize',18,'Linewidth',2,'Box','off')
    grid on

end


function [theta,s,c,t,n] = get_snell(n,theta0)

    theta = zeros(1,length(n));
    s = zeros(1,length(n));
    c = zeros(1,length(n));
    t = zeros(1,length(n));

    theta(1) = deg2rad(theta0);
    c(1) = cos(theta(1));
    s(1) = sin(theta(1));
    t(1) = s(1)/c(1);

    for i = 1:length(n)-1
        [theta_refr, s2, c2] = snell(n(i),n(i+1),theta(i));
        theta(i+1) = theta_refr;
        s(i+1) = s2;
        c(i+1) = c2;
        t(i+1) = s2/c2;
    end
    
    theta = theta(2:4);
    c = c(2:4);
    s = s(2:4);
    t = t(2:4);
    n = n(2:4);

end

function [lambda] = lambd(omega)
    lambda = 1.0./omega * 1e7;
end

function [omega] = towvnm(lambda)
    omega = 1./(lambda * 1e-7);
end

function [theta2,s2,c2,t2] = snell(n1,n2,theta1)
    s1 = 1/(2*1i)*(exp(1i*theta1) - exp(-1i*theta1));
    s2 = n1/n2 * s1;
    %s2 = n2/n1 * s1;
    %c2 = sqrt(1-s2*conj(s2));
    c2 = sqrt(1-s2^2);
    theta2 = acos(c2);
    t2 = s2/c2;
end

function [theta_SFG,s,c,t] = get_sfg_angle(omega_IR, omega_VIS, theta_IR, theta_VIS)
    theta_SFG = atan( (omega_IR * sin(theta_IR) + omega_VIS * sin(theta_VIS))/...
                      (omega_IR * cos(theta_IR) + omega_VIS * sin(theta_VIS)) );

    s = sin(theta_SFG);
    c = cos(theta_SFG);
    t = s/c;
end

function [th1,s1,c1,t1] = get_theta_1(n,gamma)
    th1 = pi/3 - asin(1/n * sin(pi/3 - gamma));
    s1 = sin(th1);
    c1 = cos(th1);
    t1 = s1/c1;
end