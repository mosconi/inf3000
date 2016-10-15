#ifndef __BALANCE_H__
#define __BALANCE_H__

#ifdef __cplusplus
extern "C" {
#endif

    typedef struct balance_t balance_t;

    balance_t *
    balance_new(char *);

    void
    balance_destroy(balance_t **);

    int
    balance_cmp(balance_t *, uint64_t, uint64_t);

    uint64_t
    balance_resource1(balance_t *);
    
    uint64_t
    balance_resource2(balance_t *);
    
    uint64_t
    balance_target(balance_t *);

    uint64_t
    balance_weightcost(balance_t *);

    void
    balance_set_cost(balance_t *, char *);

    
    void
    balance_test(bool);

#ifdef _cplusplus
}
#endif
    
#endif
