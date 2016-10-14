#ifndef __PROCESS_H__
#define __PROCESS_H__

#ifdef __cplusplus
extern "C" {
#endif

    typedef struct process_t process_t;

    process_t *
    process_new(size_t, char *);

    void
    process_destroy(process_t **);
    
#ifdef __cplusplus
}
#endif
#endif
