#ifndef __MODEL_H__
#define __MODEL_H__

#ifdef __cplpusplus
extern "C" {
#endif

    typedef struct model_t model_t;

    model_t *
    model_new(char *);

    model_t *
    model_load(char *);

    void
    model_destroy(model_t **);

    int64_t
    model_nres(model_t *);

    int64_t
    model_nmach(model_t *);

    int64_t
    model_nserv(model_t *);

    int64_t
    model_nproc(model_t *);

    int64_t
    model_nbalance(model_t *);

    int64_t
    model_wpmc(model_t *);

    int64_t
    model_wsmc(model_t *);

    int64_t
    model_wmmc(model_t *);

    resource_t *
    model_resource(model_t*, int64_t);

    machine_t *
    model_machine(model_t*, int64_t);

    process_t *
    model_process(model_t*, int64_t);

    service_t *
    model_service(model_t*, int64_t);

    balance_t *
    model_balance(model_t*, int64_t);

    bool
    model_validate(model_t *, int64_t , int64_t *);

    int64_t
    model_calculate(model_t *, int64_t , int64_t *, int64_t, int64_t *);

    
    void
    model_test(bool);

#ifdef __cplpusplus
}
#endif

#endif
