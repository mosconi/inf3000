#include "roadef.h"

int
main(int argc, char **argv){

    bool verbose=true;
    if (argc >1 && streq("-q",argv[1]))
	verbose=false;

    resource_test(verbose);
    machine_test(verbose);
    service_test(verbose);
    process_test(verbose);
    balance_test(verbose);
    model_test(verbose);
    state_test(verbose);
    
}
