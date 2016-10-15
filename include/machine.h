#ifndef __MACHINE_H__
#define __MACHINE_H__

#ifdef __cplusplus
extern "C" {
#endif
    
    typedef struct machine_t machine_t;

    machine_t *
    machine_new(size_t , size_t , char *);
    
    void
    machine_destroy(machine_t**);

    uint64_t
    machine_location(machine_t*);

    uint64_t
    machine_neigh(machine_t*);
    
    uint64_t
    machine_cap(machine_t*, uint64_t );
    
    uint64_t
    machine_safecap(machine_t*, uint64_t );
    
    uint64_t
    machine_mvcost(machine_t*, uint64_t );

    void
    machine_test(bool);
    
#ifdef __cplusplus
}
#endif

#endif
