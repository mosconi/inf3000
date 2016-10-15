#ifndef __MODEL_H__
#define __MODEL_H__

#ifdef __cplpusplus
extern "C" {
#endif

    typedef struct model_t model_t;

    model_t *
    model_new(char *);

    model_t *
    model_new_string(char *);

    void
    model_destroy(model_t **);

    uint64_t
    model_nres(model_t *);

    uint64_t
    model_nmach(model_t *);

    uint64_t
    model_nserv(model_t *);

    uint64_t
    model_nproc(model_t *);

    uint64_t
    model_nbalance(model_t *);

    uint64_t
    model_wpmc(model_t *);

    uint64_t
    model_wsmc(model_t *);

    uint64_t
    model_wmmc(model_t *);


    void
    model_test(bool);

#ifdef __cplpusplus
}
#endif

#endif
