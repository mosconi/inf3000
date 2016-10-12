#ifndef __MACHINE_H__
#define __MACHINE_H__

typedef struct machine_t machine_t;

    machine_t *
    machine_new(size_t nres, size_t nmach, char *line);
    
    void
    machine_destroy(machine_t **);

    int64_t
    machine_location(machine_t *);

    int64_t
    machine_neigh(machine_t *);
    
    int64_t
    machine_cap(machine_t *, int64_t );
    
    int64_t
    machine_safecap(machine_t *, int64_t);
    
    int64_t
    machine_mvcost(machine_t *, int64_t);

#endif
