#ifndef __ROADEF_H__
#define __ROADEF_H__

#define ROADEF_VERSION_MAJOR 1
#define ROADEF_VERSION_MINOR 0
#define ROADEF_VERSION_PATCH 0

#define ROADEF_MAKE_VERSION(major, minor, patch) \
    ((major) * 10000 + (minor) * 100 + (patch))

#define ROADEF_VERSION \
    ROADEF_MAKE_VERSION(ROADEF_VERSION_MAJOR, \
                        ROADEF_VERSION_MINOR, \
                        ROADEF_VERSION_PATCH)

#include "prelude.h"
#include "resource.h"
#include "machine.h"
#include "service.h"
#include "process.h"
#include "balance.h"
#include "model.h"
#include "state.h"

#endif
