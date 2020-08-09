#pragma once

#include <math.h>

namespace adrc {

typedef struct
{
    float v1;
    float v2;
    float r2;
    float h2;
}TD_Controller;

typedef struct
{
    float v1;
    float v2;
    float r0;
    float h0;
}TD;

typedef struct
{
    float beta1;
    float beta2;
    float alpha;
    float delta;
    float u;
    float b0;
    /* ESO */
    float z1;
    float z2;
}ESO;

typedef struct
{
    float beta1;
    float beta2;
    float u;
    float b0;
    /* LESO */
    float z1;
    float z2;
}LESO;	/* Linear ESO */

typedef struct
{
    float h1;
    float r1;
    float c;
}NLSEF;

}
