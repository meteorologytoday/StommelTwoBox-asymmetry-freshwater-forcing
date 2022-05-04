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

        use_Θ :: Bool

        function Model(;
                c,
                μ,
                ν,
                p,
                ξ,
                use_Θ = false,
        )

            return new(
                c,
                μ,
                ν,
                p,
                ξ,
                zeros(Float64, 2),
                use_Θ,
            )
        end
    end


    function pos(x)
        return (x + _abs(x, stiffness)) / 2
    end

    function neg(x)
        return (x - _abs(x, stiffness)) / 2
    end

    w = 0.1
    δ = 0.5
    c = 0.9

    function cal_Θ(m::Model; p:: Float64=m.p)
        if m.use_Θ
            return ( 1 - δ * neg(m.ξ) * 0.5 * (1 + tanh( (p - c) / w )) )
        else
            return 1.0
        end 
    end

    function cal_dΘ(m::Model; p:: Float64=m.p)
        d = 1e-7 
        dΘ =  (cal_Θ(m;p=p+d) - cal_Θ(m;p=p-d)) / (2*d)
        if ! isfinite(dΘ)
            throw(ErrorException())
        end
        return dΘ
        #return δ * m.ξ * 0.5 * (1 - ( tanh( (p - c) / w ) )^2 ) * 1/w
    end


    function cal_ψ(
        m :: Model;
        x :: AbstractArray = m.x,
        p :: Float64 = m.p,

    )
        return m.μ * (x[1] - x[2]) * cal_Θ(m; p=p) - m.ν * p * m.ξ
        #return m.μ * (x[1] - x[2]) - m.ν * p * m.ξ  
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
        
        return dabsψ * ( m.μ * (x[1] - x[2]) * cal_dΘ(m; p=p) -  m.ν * m.ξ )
        #return - dabsψ *  m.ν * m.ξ * ( step(m;p=p) + p * dstep(m; p=p) )

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
