#ifndef __INSTANCE_H__
#define __INSTANCE_H__

#ifdef __cplpusplus
extern "C" {
#endif

    typedef struct instance_t instance_t;

    instance_t *
    instance_new(char *);

    instance_t *
    instance_new_string(char *);

    void
    instance_destroy(instance_t **);

    uint64_t
    instance_nres(instance_t *);

    uint64_t
    instance_nmach(instance_t *);

    uint64_t
    instance_nserv(instance_t *);

    uint64_t
    instance_nproc(instance_t *);

    uint64_t
    instance_nbalance(instance_t *);

    uint64_t
    instance_wpmc(instance_t *);

    uint64_t
    instance_wsmc(instance_t *);

    uint64_t
    instance_wmmc(instance_t *);


    void
    instance_test(bool);

#ifdef __cplpusplus
}
#endif

#endif
