#include "roadef.h"


int
main(int argc, char **argv) {

    int opt;

    char *modelfile=NULL;
    char *assigment=NULL;

    while ((opt = getopt(argc, argv, "a:m:")) != -1 ){
	switch (opt) {
	case 'a':
	    if (assigment) free(assigment);
	    assigment = strdup(optarg);
	    break;
	case 'm':
	    if (modelfile) free(modelfile);
	    modelfile = strdup(optarg);
	    break;
	default:
	    if (assigment) free(assigment);
	    if (modelfile) free(modelfile);
	    exit(-1);
	}
    }

    if (!modelfile || streq("",modelfile))
	exit(-2);

    if (!assigment || streq("",assigment))
	exit(-3);

    model_t *model = model_load(modelfile);
    state_t *state = state_new(model);
    state_load(state,assigment);
    

    printf("    memoria (current): %zu B\n", getCurrentRSS() );
    printf("    memoria (Peak):    %zu B\n", getPeakRSS() );

    state_destroy(&state);
    model_destroy(&model);
    if (assigment) free(assigment);
    if (modelfile) free(modelfile);
    
}
