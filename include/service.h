#ifndef __SERVICE_H__
#define __SERVICE_H__

#ifdef __cplusplus
extern "C" {
#endif

    typedef struct service_t service_t;

    service_t *
    service_new(char *);

    void
    service_destroy(service_t **);
    
#ifdef __cplusplus
}
#endif
#endif
