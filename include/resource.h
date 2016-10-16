#ifndef __RESOURCE_H__
#define __RESOURCE_H__

#ifdef __cplusplus
extern "C" {
#endif

    typedef struct resource_t resource_t;

    resource_t *
    resource_new(char *);

    void
    resource_destroy(resource_t **);

    void
    resource_test(bool );

    bool
    resource_transient(resource_t *);

    int64_t
    resource_loadcost(resource_t *);
    
#ifdef __cplusplus
}
#endif

#endif
