module StommelModel


    export Model, F, dFdx, dFdp

    const stiffness = 1.0

    mutable struct Model

        c :: Float64
        μ :: Float64
        ν :: Float64
        p :: Float64
        ξ :: Float64
        x :: AbstractArray{Float64, 1}

        function Model(;
                c,
                μ,
                ν,
                p,
                ξ,
        )

            return new(
                c,
                μ,
                ν,
                p,
                ξ,
                zeros(Float64, 2),
            )
        end
    end

    function cal_ψ(
        m :: Model;
        x :: AbstractArray = m.x,
        p :: Float64 = m.p,

    )

        return m.μ * (x[1] - x[2]) - m.ν * p * m.ξ 
        #return m.μ^0.5 * (x[1] - x[2]) - m.ν * p * m.ξ 
    
    end


    function cal_M(
        m :: Model;
        x :: AbstractArray = m.x,
        p :: Float64 = m.p,
    )
        
        return 1.0 + _abs( cal_ψ(m; x=x, p=p),  stiffness)
        #return 1.0 + cal_ψ(m; x=x, p=p)^2
    
    end

    function cal_dMdx(
        m :: Model;
        x :: AbstractArray = m.x,
        p :: Float64 = m.p,
    )

        ψ = cal_ψ(m; x=x, p=p)
        dψdx1 =   m.μ
        dψdx2 = - m.μ

        dabsψ = _dabs(ψ, stiffness)

        return [ dabsψ * dψdx1 ; dabsψ * dψdx2 ]
 

        #dψdx1 =   m.μ^0.5
        #dψdx2 = - m.μ^0.5
        #return [ 2 * ψ * dψdx1 ; 2 * ψ * dψdx2 ]
        
    end

    function cal_dMdp(
        m :: Model;
        x :: AbstractArray = m.x,
        p :: Float64 = m.p,
    )

        ψ = cal_ψ(m; x=x, p=p)
        dabsψ = _dabs(ψ, stiffness)
        
        return dabsψ * ( - m.ν * m.ξ )
        #return 2 * ψ * ( - m.ν * m.ξ )

       
    end

    function F(
        m :: Model;
        x :: AbstractArray = m.x,
        p :: Float64 = m.p,
    )

        M = cal_M(m; x=x, p=p)

        return [
            (- m.c * (x[1] - 1.0) - x[1] * M) ;
            (p - x[2] * M) ; 
        ]

    end

    function dFdx(
        m :: Model;
        x :: AbstractArray = m.x,
        p :: Float64 = m.p,
    )

        M    = cal_M(m; x=x, p=p)
        dMdx = cal_dMdx(m; x=x, p=p)

        return [
            (- m.c - M - x[1] * dMdx[1])     (- x[1] * dMdx[2])     ;
            (- x[2] * dMdx[1])               (- M - x[2] * dMdx[2]) ;
        ]

    end

    function dFdp(
        m :: Model;
        x :: AbstractArray = m.x,
        p :: Float64 = m.p,
    )

        dMdp = cal_dMdp(m; x=x, p=p)

        return [
                - x[1] * dMdp ; 
            1.0 - x[2] * dMdp ;
        ]

    end



    function _abs(
        x :: Float64,
        k :: Float64,
    )

        if isinf(k) && k > 0
            return abs(x)
        elseif k*x > 50
            return abs(x)
        else
            return 2 / k * ( log( 1 + exp(k * x) ) - log(2) ) - x
        end
    end

    function _dabs(
        x :: Float64,
        k :: Float64,
    )

        if isinf(k) && k > 0
            return sign(x)
        elseif k*x > 50 
            return sign(x)
        else
            return 2 * ( exp(k*x) / (1 + exp(k*x)) ) - 1.0
        end
    end


end
