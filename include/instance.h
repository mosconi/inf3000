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

#ifdef __cplpusplus
}
#endif

#endif
