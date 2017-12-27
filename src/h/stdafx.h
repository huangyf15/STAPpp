#pragma once

//	Clear an array
template <class type> inline void clear(type* a, unsigned int N)
{
    for (unsigned int i = 0; i < N; i++)
        a[i] = 0;
}
