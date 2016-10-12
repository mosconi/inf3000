#ifndef __RESOURCE_H__
#define __RESOURCE_H__

#ifdef __cplusplus
extern "C" {
#endif

    typedef struct resource_t resource_t;

    resource_t *
    resource_new();

    void
    resource_destroy(resouce_t **);
    
#ifdef __cplusplus
}
#endif

#endif
