#include "roadef.h"

#ifdef __cplusplus
extern "C" {
#endif

    typedef 
    struct balance_t {
	int64_t r1;
	int64_t r2;
	int64_t target;
	int64_t wbc;
    } balance_t;
    
    struct instance_t {
	size_t nof_resources;
	bool *resouce_transient;
	int64_t *wlc;
	size_t nof_mach;
	int64_t *mach_location;
	int64_t *mach_neighborhood;
	int64_t *mach_location;
	int64_t *mach_cap;
	
	size_t nof_proc;
	
	size_t nof_serv;
	size_t nof_locations;
	size_t nof_neighborhood;
	size_t nof_balance;
	balance_t *balance;
	int64_t wpmc;
	int64_t wsmc;
	int64_t wmmc;
    } ;

    
#ifdef __cplusplus
};
#endif	
