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

    uint64_t
    process_service(process_t *);

    uint64_t
    process_movecost(process_t *);

    uint64_t
    process_requirement(process_t *, size_t );

    void
    process_test(bool);
    

#ifdef __cplusplus
}
#endif
#endif
