/*
Misc. util functions.
*/

#ifndef __is_util_h__
#define __is_util_h__

inline int
index_wrap(const int i, const int n) {
    int r = i - n * (int)( (double)i / (double)n );
    if (r < 0) {
        r += n;
    }
    return r;
}

#endif
