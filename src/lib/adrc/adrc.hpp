#pragma once
#include <mathlib/math/Functions.hpp>

using namespace math;

namespace adrc {

struct TD_Controller
{
    float v1;
    float v2;
    float r2;
    float h2;
};

struct TD
{
    float v1;
    float v2;
    float r0;
    float h0;
};

struct ESO
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
};

struct LESO
{
    float beta1;
    float beta2;
    float u;
    float b0;
    /* LESO */
    float z1;
    float z2;
};	/* Linear ESO */

struct NLSEF
{
    float h1;
    float r1;
    float c;
};

inline float fhan(float v1, float v2, float r0, float h0)
{
    float d = h0 * h0 * r0;
    float a0 = h0 * v2;
    float y = v1 + a0;
    float a1 = sqrt(d*(d + 8.0f*fabsf(y)));
    float a2 = a0 + signNoZero(y)*(a1-d)*0.5f;
    float sy = (signNoZero(y+d) - signNoZero(y-d))*0.5f;
    float a = (a0 + y - a2)*sy + a2;
    float sa = (signNoZero(a+d) - signNoZero(a-d))*0.5f;

    return -r0*(a/d - signNoZero(a))*sa - r0*signNoZero(a);
}

inline float fal(float e, float alpha, float delta)
{
    if(fabsf(e) <= delta){
        return e / (pow(delta, 1.0f-alpha));
    }else{
        return pow(fabsf(e), alpha) * signNoZero(e);
    }
}

inline void td_init(TD* td_t, float r0, float h0)
{
    td_t->r0 = r0;
    td_t->h0 = h0;
    td_t->v1 = td_t->v2 = 0.0f;
}
inline void td(TD* td_t, float v, float dt)
{
    float fv = fhan(td_t->v1 - v, td_t->v2, td_t->r0, td_t->h0);

    td_t->v1 += dt * td_t->v2;
    td_t->v2 += dt * fv;
}
inline void td_control_init(TD_Controller* td_controller, float r2, float h2)
{
    td_controller->v1 = 0.0f;
    td_controller->v2 = 0.0f;
    td_controller->r2 = r2;
    td_controller->h2 = h2;
}
inline float td_control(TD_Controller* td_controller, float err, float dt)
{
    float fv = fhan(-err, td_controller->v2, td_controller->r2, td_controller->h2);
    td_controller->v1 += dt * td_controller->v2;
    td_controller->v2 += dt * fv;

    return td_controller->v2;
}
inline void eso_init(ESO* eso_t, float beta1, float beta2, float alpha, float delta, float b0)
{
    eso_t->beta1 = beta1;
    eso_t->beta2 = beta2;
    eso_t->u = 0.0f;
    eso_t->alpha = alpha;
    eso_t->delta = delta;
    eso_t->b0 = b0;

    eso_t->z1 = eso_t->z2 = 0.0f;
}
inline void eso(ESO* eso_t, float y, float dt)
{
    float e = eso_t->z1 - y;
    float fe = fal(e, eso_t->alpha, eso_t->delta);

    eso_t->z1 += dt*(eso_t->z2 + eso_t->b0*eso_t->u - eso_t->beta1*e);
    eso_t->z2 -= dt*eso_t->beta2*fe;
}
inline void leso_init(LESO* leso_t, float w, float b0)
{
    // (s + w)^2 = s^2 + beta_1 * s + beta_2
    leso_t->beta1 = 2.0f*w;
    leso_t->beta2 = w*w;
    leso_t->u = 0.0f;
    leso_t->b0 = b0;

    leso_t->z1 = leso_t->z2 = 0.0f;
}
inline void leso(LESO* leso_t, float y, float dt)
{
    float e = leso_t->z1 - y;

    leso_t->z1 += dt*(leso_t->z2 + leso_t->b0*leso_t->u - leso_t->beta1*e);
    leso_t->z2 -= dt*leso_t->beta2*e;
}
inline void nlsef_init(NLSEF* nlsef_t, float r1, float h1, float c)
{
    nlsef_t->h1 = h1;
    nlsef_t->r1 = r1;
    nlsef_t->c = c;
}
inline float nlsef(NLSEF* nlsef_t, float e1, float e2)
{
    float u0 = -fhan(e1, nlsef_t->c*e2, nlsef_t->r1, nlsef_t->h1);

    return u0;
}

}
