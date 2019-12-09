function calc_AB(alpha::Number, beta::Number)
    #=
    Calculate A & B from alpha & beta (Dorren 1987)
    =#
    if (beta - alpha) > (pi / 2.0)# spot out of view
        return 0.0, 0.0
    end
    cosalpha = cos(alpha)
    sinalpha = sin(alpha)
    cosbeta = cos(beta)
    sinbeta = sin(beta)
    tanbeta = sinbeta / cosbeta
    if (beta + alpha) <= (pi / 2.0)#spot fully visble
        delta = 0.0
        sindelta = 0.0
        cosdelta = 1.0
        zeta = 0.0
        sinzeta = 0.0
        coszeta = 1.0
    else
        cosdelta = 1.0 / tan(alpha) / tan(beta)
        delta = acos(cosdelta)
        sindelta = sin(delta)
        sinzeta = sindelta * sinalpha
        zeta = asin(sinzeta)
    end
    if beta <= (pi / 2.0)
        T = atan(sinzeta * tanbeta)
    else
        T = pi - atan(-sinzeta * tanbeta)
    end
    biga = zeta + (pi - delta) * cosbeta * sinalpha^2 -
    sinzeta * sinbeta * cosalpha
    bigb = (1.0/3.0) * (pi - delta) *
        (-2 * cosalpha^3 - 3 * sinbeta^2 * cosalpha * sinalpha^2) +
        (2.0/3.0) * (pi - T) + (1.0/6.0) * sinzeta * sin(2 * beta) *
        (2 - 3 * cosalpha^2)
    return biga, bigb
end

function calc_ab(ustar::Number, uspot::Number, fratio::Number)
    #=
    Calculate a & b from u_star & u_spot & F_spot/F_star (Dorren
    1987)
    =#
    littlea = (1 - ustar) - (1 - uspot) * fratio
    littleb = ustar - uspot * fratio
    return littlea, littleb
end

function dorren_F(ustar::Number, uspot::Number, fratio::Number, alpha::Number, beta::Number)
    #=
    Calculate F (fraction of stellar disk hidden by spot) from
    u_star, u_spot, F_spot/F_star, alpha & beta, following Dorren (1987)
    =#
    biga, bigb = calc_AB(alpha, beta)
    littlea, littleb = calc_ab(ustar, uspot, fratio)
    F = (littlea * biga + littleb * bigb) / pi / (1 - ustar / 3.0)
    if F <0
        F = 0
    end
    return -F
end
