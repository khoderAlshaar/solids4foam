#ifndef VanAlbadaLimiter_H
#define VanAlbadaLimiter_H

namespace Foam
{

class VanAlbadaLimiter
{
public:
    VanAlbadaLimiter() {}

    // For scalar fields
    inline void limiter
    (
        scalar& lim,
        const scalar& cellVolume,
        const scalar& deltaOneMax,
        const scalar& deltaOneMin,
        const scalar& extrapolate
    )
    {
        if (mag(extrapolate) < SMALL)
        {
            return;
        }
        // Numerator: 2 * delta
        // Denominator: (delta^2 + SMALL to prevent divide by zero)
        scalar r = (extrapolate) / (mag(extrapolate) + SMALL); // Gradient ratio
        lim = min(lim, (2.0*r)/(r*r+1.0+SMALL));
    }

    // For vector fields
    template<class Type>
    inline void limiter
    (
        Type& lim,
        const scalar& cellVolume,
        const Type& deltaOneMax,
        const Type& deltaOneMin,
        const Type& extrapolate
    )
    {
        for (direction cmpt = 0; cmpt < Type::nComponents; cmpt++)
        {
            limiter
            (
                lim.component(cmpt),
                cellVolume,
                deltaOneMax.component(cmpt),
                deltaOneMin.component(cmpt),
                extrapolate.component(cmpt)
            );
        }
    }
};

} // End namespace Foam

#endif
