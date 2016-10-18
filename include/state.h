#ifndef __STATE_H__
#define __STATE_H__

#ifdef __cplusplus
extern "C" {
#endif

    typedef struct state_t state_t;

    state_t *
    state_new(model_t *);

    void
    state_destroy(state_t **);

    void
    state_load(state_t *, char *);

    void
    state_move(state_t *, int64_t, int64_t);

    void
    state_swap(state_t *, int64_t, int64_t);

    int64_t
    state_machine_curr(state_t *, int64_t);

    int64_t
    state_machine_begin(state_t *, int64_t);

    bool
    state_machine_changed(state_t *, int64_t);

    int64_t
    state_obj_resource(state_t *, int64_t);

    int64_t
    state_obj_balance(state_t *, int64_t);

    int64_t
    state_obj_pmc(state_t *);

    int64_t
    state_obj_smc(state_t *);
    
    int64_t
    state_obj_mmc(state_t *);

    bool
    state_validate(state_t *);

    void
    state_test(bool);
    
#ifdef __cplusplus
}
#endif

#endif
