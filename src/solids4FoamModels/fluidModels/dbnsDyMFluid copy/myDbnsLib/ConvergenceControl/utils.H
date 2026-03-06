// utils.H - small helper macros and safe lookup helpers

#ifndef UTILS_H
#define UTILS_H

// Safe dictionary lookup helper template
template<typename T>
T lookupOrDefault(const dictionary& dict, const word& name, const T& def)
{
    if (dict.found(name))
    {
        T tmp;
        dict.lookup(name) >> tmp;
        return tmp;
    }
    return def;
}

#endif // UTILS_H
