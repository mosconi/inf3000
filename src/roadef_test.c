#include "roadef.h"


static void
roadef_test(bool verbose){
    char string = strdup("2\n1 100\n0 10\n4\n0 0 30 400 16 80 0 1 4 5\n0 0 10 240 8 160 1 0 3 4\n1 1 15 100 12 80 4 3 0 2\n1 2 10 100 8 80 5 4 2 0\n2\n2 0 \n1 1 0 \n3\n0 12 10 100\n0 10 20 100\n1 6 200 1\n1 \n0 1 20\n10\n1 10 100\n");
    instance_t *instance = instance_new_string(string);
    instance_destroy(&instance);
}

int
main(int argc, char **argv){

    bool verbose=false;
    process_test(verbose);
    roadef_test(verbose);

    
}
