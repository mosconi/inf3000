#include "roadef.h"
    
    struct instance_t {
	size_t nof_resources;
	bool *resouce_transient;
	int64_t *wlc;
	size_t nof_mach;
	machine_t *machines;
	
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

    
