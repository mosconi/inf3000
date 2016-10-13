#ifndef __BALANCE_H__
#define __BALANCE_H__

#ifdef __cplusplus
extern "C" {
#endif

    typedef struct balance_t balance_t;

    balance_t *
    balance_new();

    void
    balance_destroy(balance_t **);

    int64_t
    balance_lookup(balance_t *, int64_t, int64_t);

#ifdef _cplusplus
}
#endif
    
#endif
