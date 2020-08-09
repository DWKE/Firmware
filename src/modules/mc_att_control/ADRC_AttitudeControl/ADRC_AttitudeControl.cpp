/**
 * @file ADRC_AttitudeControl.cpp
 */

#include <ADRC_AttitudeControl.hpp>

#include <mathlib/math/Limits.hpp>
#include <mathlib/math/Functions.hpp>

#include <px4_log.h>
#include <px4_app.h>
#include <px4_tasks.h>

using namespace matrix;
using namespace adrc;

ADRC_AttitudeControl::ADRC_AttitudeControl()
{
    _td_controller.push_back(TD_Controller());
    _td_controller.push_back(TD_Controller());
    _td.push_back(TD());
    _td.push_back(TD());
    _nlsef.push_back(NLSEF());
    _nlsef.push_back(NLSEF());
    _leso.push_back(LESO());
    _leso.push_back(LESO());
}

void ADRC_AttitudeControl::setGains(float td_control_r2, float td_control_h2, float td_r0, float leso_w,
                                    float nlsef_r1, float nlsef_h1, float nlsef_c, float gamma, float nlsef_ki)
{
    _td_control_r2 = td_control_r2;
    _td_control_h2 = td_control_h2;
    _td_r0 = td_r0;
    _leso_w = leso_w;
    _nlsef_r1 = nlsef_r1;
    _nlsef_h1 = nlsef_h1;
    _nlsef_c = nlsef_c;
    _gamma = gamma;
    _nlsef_ki = nlsef_ki;
}

void ADRC_AttitudeControl::att_init(float h)
{
    td_control_init(&_td_controller[0], _td_control_r2, _td_control_h2*h);
    td_control_init(&_td_controller[1], _td_control_r2, _td_control_h2*h);
    td_init(&_td[0], _td_r0, h);
    td_init(&_td[1], _td_r0, h);
    leso_init(&_leso[0], _leso_w, 400);
    leso_init(&_leso[1], _leso_w, 400);
    nlsef_init(&_nlsef[0], _nlsef_r1, _nlsef_h1*h, _nlsef_c);
    nlsef_init(&_nlsef[1], _nlsef_r1, _nlsef_h1*h, _nlsef_c);

    int_i[0] = int_i[1] = 0.0f;
}
void ADRC_AttitudeControl::att_reset(float h)
{

}

float ADRC_AttitudeControl::sign(float val)
{
    if(val >= 0.0f)
        return 1.0f;
    else
        return -1.0f;
}

float ADRC_AttitudeControl::fhan(float v1, float v2, float r0, float h0)
{
    float d = h0 * h0 * r0;
    float a0 = h0 * v2;
    float y = v1 + a0;
    float a1 = sqrt(d*(d + 8.0f*fabsf(y)));
    float a2 = a0 + sign(y)*(a1-d)*0.5f;
    float sy = (sign(y+d) - sign(y-d))*0.5f;
    float a = (a0 + y - a2)*sy + a2;
    float sa = (sign(a+d) - sign(a-d))*0.5f;

    return -r0*(a/d - sign(a))*sa - r0*sign(a);
}

float ADRC_AttitudeControl::fal(float e, float alpha, float delta)
{
    if(fabsf(e) <= delta){
        return e / (pow(delta, 1.0f-alpha));
    }else{
        return pow(fabsf(e), alpha) * sign(e);
    }
}

void ADRC_AttitudeControl::td_init(TD* td_t, float r0, float h0)
{
    td_t->r0 = r0;
    td_t->h0 = h0;
    td_t->v1 = td_t->v2 = 0.0f;
}
void ADRC_AttitudeControl::td(TD* td_t, float v, float dt)
{
    float fv = fhan(td_t->v1 - v, td_t->v2, td_t->r0, td_t->h0);

    td_t->v1 += dt * td_t->v2;
    td_t->v2 += dt * fv;
}
void ADRC_AttitudeControl::td_control_init(TD_Controller* td_controller, float r2, float h2)
{
    td_controller->v1 = 0.0f;
    td_controller->v2 = 0.0f;
    td_controller->r2 = r2;
    td_controller->h2 = h2;
}
float ADRC_AttitudeControl::td_control(TD_Controller* td_controller, float err, float dt)
{
    float fv = fhan(-err, td_controller->v2, td_controller->r2, td_controller->h2);
    td_controller->v1 += dt * td_controller->v2;
    td_controller->v2 += dt * fv;

    return td_controller->v2;
}
void ADRC_AttitudeControl::eso_init(ESO* eso_t, float beta1, float beta2, float alpha, float delta, float b0)
{
    eso_t->beta1 = beta1;
    eso_t->beta2 = beta2;
    eso_t->u = 0.0f;
    eso_t->alpha = alpha;
    eso_t->delta = delta;
    eso_t->b0 = b0;

    eso_t->z1 = eso_t->z2 = 0.0f;
}
void ADRC_AttitudeControl::eso(ESO* eso_t, float y, float dt)
{
    float e = eso_t->z1 - y;
    float fe = fal(e, eso_t->alpha, eso_t->delta);

    eso_t->z1 += dt*(eso_t->z2 + eso_t->b0*eso_t->u - eso_t->beta1*e);
    eso_t->z2 -= dt*eso_t->beta2*fe;
}
void ADRC_AttitudeControl::leso_init(LESO* leso_t, float w, float b0)
{
    // (s + w)^2 = s^2 + beta_1 * s + beta_2
    leso_t->beta1 = 2.0f*w;
    leso_t->beta2 = w*w;
    leso_t->u = 0.0f;
    leso_t->b0 = b0;

    leso_t->z1 = leso_t->z2 = 0.0f;
}
void ADRC_AttitudeControl::leso(LESO* leso_t, float y, float dt)
{
    float e = leso_t->z1 - y;

    leso_t->z1 += dt*(leso_t->z2 + leso_t->b0*leso_t->u - leso_t->beta1*e);
    leso_t->z2 -= dt*leso_t->beta2*e;
}
void ADRC_AttitudeControl::nlsef_init(NLSEF* nlsef_t, float r1, float h1, float c)
{
    nlsef_t->h1 = h1;
    nlsef_t->r1 = r1;
    nlsef_t->c = c;
}
float ADRC_AttitudeControl::nlsef(NLSEF* nlsef_t, float e1, float e2)
{
    float u0 = -fhan(e1, nlsef_t->c*e2, nlsef_t->r1, nlsef_t->h1);

    return u0;
}

Vector3f ADRC_AttitudeControl::att_dis_comp(Vector2f in)
{
    Vector3f out;
    out(0) = in(0) = _gamma * _leso[0].z2 / _leso[0].b0;
    out(1) = in(1) = _gamma * _leso[1].z2 / _leso[1].b0;

    _leso[0].u = out(0);
    _leso[1].u = out(1);

    return out;
}
Vector3f ADRC_AttitudeControl::att_control(Vector3f err, const Vector3f gyr, float dt)
{
    Vector3f sp_rate;
    float rate_err[3];
    Vector2f u0;

    sp_rate(0) = td_control(&_td_controller[0], err(0), dt);
    sp_rate(1) = td_control(&_td_controller[1], err(1), dt);
    rate_err[0] = sp_rate(0) - gyr(0);
    rate_err[1] = sp_rate(1) - gyr(1);

    td(&_td[0], rate_err[0], dt);
    td(&_td[1], rate_err[1], dt);

    u0(0) = nlsef(&_nlsef[0], rate_err[0], _td[0].v2)/_leso[0].b0;
    u0(1) = nlsef(&_nlsef[1], rate_err[1], _td[1].v2)/_leso[1].b0;

    int_i[0] += rate_err[0] * _nlsef_ki * _nlsef[0].h;
    int_i[1] += rate_err[1] * _nlsef_ki * _nlsef[1].h;



    return sp_rate;
}
void ADRC_AttitudeControl::add_observer_update(const float gyr[3], float bth)
{

}

matrix::Vector3f ADRC_AttitudeControl::update(matrix::Quatf q, matrix::Quatf qd, const matrix::Vector3f rate, float dt, const float yawspeed_feedforward)
{
    // ensure input quaternions are exactly normalized because acosf(1.00001) == NaN
    q.normalize();
    qd.normalize();

    // calculate reduced desired attitude neglecting vehicle's yaw to prioritize roll and pitch
    const Vector3f e_z = q.dcm_z();
    const Vector3f e_z_d = qd.dcm_z();
    Quatf qd_red(e_z, e_z_d);

    if (fabs(qd_red(1)) > (1.f - 1e-5f) || fabs(qd_red(2)) > (1.f - 1e-5f)) {
        // In the infinitesimal corner case where the vehicle and thrust have the completely opposite direction,
        // full attitude control anyways generates no yaw input and directly takes the combination of
        // roll and pitch leading to the correct desired yaw. Ignoring this case would still be totally safe and stable.
        qd_red = qd;

    } else {
        // transform rotation from current to desired thrust vector into a world frame reduced desired attitude
        qd_red *= q;
    }

    // mix full and reduced desired attitude
    Quatf q_mix = qd_red.inversed() * qd;
    q_mix *= math::signNoZero(q_mix(0));
    // catch numerical problems with the domain of acosf and asinf
    q_mix(0) = math::constrain(q_mix(0), -1.f, 1.f);
    q_mix(3) = math::constrain(q_mix(3), -1.f, 1.f);
    qd = qd_red * Quatf(cosf(_yaw_w * acosf(q_mix(0))), 0, 0, sinf(_yaw_w * asinf(q_mix(3))));

    // quaternion attitude control law, qe is rotation from q to qd
    const Quatf qe = q.inversed() * qd;

    // using sin(alpha/2) scaled rotation axis as attitude error (see quaternion definition by axis angle)
    // also taking care of the antipodal unit quaternion ambiguity
    const Vector3f eq = 2.f * math::signNoZero(qe(0)) * qe.imag();
    //PX4_INFO("%f, %f, %f", (double)eq(0), (double)eq(1), (double)eq(2));

    // calculate angular rates setpoint
//    matrix::Vector3f rate_setpoint = eq.emult(_proportional_gain);
    matrix::Vector3f rate_setpoint = att_control(eq, rate, dt);
    PX4_INFO("%f, %f, %f", (double)rate_setpoint(0), (double)rate_setpoint(1), (double)rate_setpoint(2));


    // Feed forward the yaw setpoint rate.
    // yaw_sp_move_rate is the feed forward commanded rotation around the world z-axis,
    // but we need to apply it in the body frame (because _rates_sp is expressed in the body frame).
    // Therefore we infer the world z-axis (expressed in the body frame) by taking the last column of R.transposed (== q.inversed)
    // and multiply it by the yaw setpoint rate (yaw_sp_move_rate).
    // This yields a vector representing the commanded rotatation around the world z-axis expressed in the body frame
    // such that it can be added to the rates setpoint.
    rate_setpoint += q.inversed().dcm_z() * yawspeed_feedforward;

    // limit rates
    for (int i = 0; i < 3; i++) {
        rate_setpoint(i) = math::constrain(rate_setpoint(i), -_rate_limit(i), _rate_limit(i));
    }

    return rate_setpoint;
}


/****************************************************************************************************************/

//void ADRC_AttitudeControl::setProportionalGain(const matrix::Vector3f &proportional_gain)
//{
//	_proportional_gain = proportional_gain;

//	// prepare yaw weight from the ratio between roll/pitch and yaw gains
//	const float roll_pitch_gain = (proportional_gain(0) + proportional_gain(1)) / 2.f;
//	_yaw_w = math::constrain(proportional_gain(2) / roll_pitch_gain, 0.f, 1.f);

//	_proportional_gain(2) = roll_pitch_gain;
//}

//matrix::Vector3f ADRC_AttitudeControl::update(matrix::Quatf q, matrix::Quatf qd, const float yawspeed_feedforward)
//{
//	// ensure input quaternions are exactly normalized because acosf(1.00001) == NaN
//	q.normalize();
//	qd.normalize();

//	// calculate reduced desired attitude neglecting vehicle's yaw to prioritize roll and pitch
//	const Vector3f e_z = q.dcm_z();
//	const Vector3f e_z_d = qd.dcm_z();
//	Quatf qd_red(e_z, e_z_d);

//    if (fabs(qd_red(1)) > (1.f - 1e-5f) || fabs(qd_red(2)) > (1.f - 1e-5f)) {
//		// In the infinitesimal corner case where the vehicle and thrust have the completely opposite direction,
//		// full attitude control anyways generates no yaw input and directly takes the combination of
//		// roll and pitch leading to the correct desired yaw. Ignoring this case would still be totally safe and stable.
//		qd_red = qd;

//	} else {
//		// transform rotation from current to desired thrust vector into a world frame reduced desired attitude
//		qd_red *= q;
//	}

//	// mix full and reduced desired attitude
//	Quatf q_mix = qd_red.inversed() * qd;
//	q_mix *= math::signNoZero(q_mix(0));
//	// catch numerical problems with the domain of acosf and asinf
//	q_mix(0) = math::constrain(q_mix(0), -1.f, 1.f);
//	q_mix(3) = math::constrain(q_mix(3), -1.f, 1.f);
//	qd = qd_red * Quatf(cosf(_yaw_w * acosf(q_mix(0))), 0, 0, sinf(_yaw_w * asinf(q_mix(3))));

//	// quaternion attitude control law, qe is rotation from q to qd
//	const Quatf qe = q.inversed() * qd;

//	// using sin(alpha/2) scaled rotation axis as attitude error (see quaternion definition by axis angle)
//	// also taking care of the antipodal unit quaternion ambiguity
//	const Vector3f eq = 2.f * math::signNoZero(qe(0)) * qe.imag();

//	// calculate angular rates setpoint
//	matrix::Vector3f rate_setpoint = eq.emult(_proportional_gain);

//	// Feed forward the yaw setpoint rate.
//	// yaw_sp_move_rate is the feed forward commanded rotation around the world z-axis,
//	// but we need to apply it in the body frame (because _rates_sp is expressed in the body frame).
//	// Therefore we infer the world z-axis (expressed in the body frame) by taking the last column of R.transposed (== q.inversed)
//	// and multiply it by the yaw setpoint rate (yaw_sp_move_rate).
//	// This yields a vector representing the commanded rotatation around the world z-axis expressed in the body frame
//	// such that it can be added to the rates setpoint.
//	rate_setpoint += q.inversed().dcm_z() * yawspeed_feedforward;

//	// limit rates
//	for (int i = 0; i < 3; i++) {
//		rate_setpoint(i) = math::constrain(rate_setpoint(i), -_rate_limit(i), _rate_limit(i));
//	}

//	return rate_setpoint;
//}
