module StommelModel


    export Model, F, dFdx


    mutable struct Θ_Func
        Θ    :: Function
        dΘdp :: Function
 
        function Θ_Func(
                Θ    :: Union{Function, Integer} = 0,
                dΘdp :: Union{Function, Integer} = 0,
        )

            if typeof(Θ) <: Integer

                if Θ == 0
                    Θ    = (p, ) -> 0.0
                    dΘdp = (p, ) -> 0.0
                elseif Θ == 1
                    fac = 0.3
                    Θ    = (p, ) -> fac * p
                    dΘdp = (p, ) -> fac
                elseif Θ == 2

                    fac = 1.0
                    cen = 1.0
                    wid = 0.3
                    
                    Θ    = (p, ) -> ((1.0 + tanh( (p - cen)/wid ))/2) * fac - ((1.0 + tanh( (-cen)/wid ))/2) * fac
                    dΘdp = (p, ) -> (((sech( (p - cen)/wid ))^2)/2 / wid) * fac

                else

                    throw(ErrorException("P can only be 0, 1, or 2. Now it is $P")) 

                end

            end

            return new(Θ, dΘdp)
        end  
    end

    mutable struct Model
        α :: Float64
        μ :: Float64
        p :: Float64
        ξ :: Float64
        x :: AbstractArray{Float64, 1}


        Θ_func :: Θ_Func

        function Model(;
                α,
                μ,
                p,
                ξ,
                Θ_func,
        )

            return new(
                α,
                μ,
                p,
                ξ,
                zeros(Float64, 2),
                Θ_func,
            )
        end
    end

    function MOC(
        x :: AbstractArray{Float64, 1},
        α :: Float64,
        μ :: Float64,
        p :: Float64,
        ξ :: Float64,
        Θ_func :: Θ_Func,
    )
        Δ = x[1] - x[2]
        
        return μ * Δ * sqrt(1 - ξ * Θ_func.Θ(p))
    
    end


    function MIXING(
        x :: AbstractArray{Float64, 1},
        α :: Float64,
        μ :: Float64,
        p :: Float64,
        ξ :: Float64,
        Θ_func :: Θ_Func,
    )
        Δ = x[1] - x[2]
        
        return 1 + μ^2 * Δ^2 * (1 - ξ * Θ_func.Θ(p))
    
    end

    function dMIXINGdx(
        x :: AbstractArray{Float64, 1},
        α :: Float64,
        μ :: Float64,
        p :: Float64,
        ξ :: Float64, 
        Θ_func :: Θ_Func,
    )
        Δ = x[1] - x[2]
        return 2 * μ^2 * Δ * (1 - ξ * Θ_func.Θ(p))
        
    end

    function dMIXINGdy(
        x :: AbstractArray{Float64, 1},
        α :: Float64,
        μ :: Float64,
        p :: Float64,
        ξ :: Float64, 
        Θ_func :: Θ_Func,
    )
        Δ = x[1] - x[2]
        return - 2 * μ^2 * Δ * (1 - ξ * Θ_func.Θ(p))
        
    end

    function dMIXINGdp(
        x :: AbstractArray{Float64, 1},
        α :: Float64,
        μ :: Float64,
        p :: Float64,
        ξ :: Float64, 
        Θ_func :: Θ_Func,
    )
        Δ = x[1] - x[2]
        
        return - μ^2 * Δ^2 * ξ * Θ_func.dΘdp(p)
    
    end



    function F(
        x :: AbstractArray{Float64, 1},
        α :: Float64,
        μ :: Float64,
        p :: Float64,
        ξ :: Float64, 
        Θ_func :: Θ_Func,
    )

        C = MIXING(x, α, μ, p, ξ, Θ_func)

        return [
            (- α * (x[1] - 1.0) - x[1] * C) ;
            (p - x[2] * C) ; 
        ]
    end

    function dFdx(
        x :: AbstractArray{Float64, 1},
        α :: Float64,
        μ :: Float64,
        p :: Float64,
        ξ :: Float64,
        Θ_func :: Θ_Func,
    )

        C    =    MIXING(x, α, μ, p, ξ, Θ_func)
        dCdx = dMIXINGdx(x, α, μ, p, ξ, Θ_func)
        dCdy = dMIXINGdy(x, α, μ, p, ξ, Θ_func)

        return [
            (- α - C - x[1] * dCdx)       (- x[1] * dCdy)     ;
            (- x[2] * dCdx)               (- C - x[2] * dCdy) ;
        ]

    end

    function dFdp(
        x :: AbstractArray{Float64, 1},
        α :: Float64,
        μ :: Float64,
        p :: Float64,
        ξ :: Float64,
        Θ_func :: Θ_Func,
    )
        
        dCdp = dMIXINGdp(x, α, μ, p, ξ, Θ_func)

        return [
                - x[1] * dCdp ; 
            1.0 - x[2] * dCdp ;
        ]
    end

    function dFdξ(
        x :: AbstractArray{Float64, 1},
        α :: Float64,
        μ :: Float64,
        p :: Float64,
        ξ :: Float64,
    )
        return [
            x[1] ;
            x[2] ;
        ]
    end

end
