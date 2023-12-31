//
// Created by edwa0468 on 31/12/2023.
//

#include <stdio.h>
#include <stdlib.h>
#include "colours.h"
#include "error.h"

int error_and_return(char * msg) {
    fprintf(stderr, "%sERROR: %s%s\n", RED, msg, ENDC);
    return 2;
}

void error_and_exit(char * msg) {
    fprintf(stderr, "%sERROR: %s%s\n", RED, msg, ENDC);
    exit(2);
}