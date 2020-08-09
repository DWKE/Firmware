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
}

void ADRC_AttitudeControl::setGains(float td_control_r2, float td_control_h2)
{
    _td_control_r2 = td_control_r2;
    _td_control_h2 = td_control_h2;
}

void ADRC_AttitudeControl::att_init(float h)
{
    td_control_init(&_td_controller[0], _td_control_r2, _td_control_h2*h);
    td_control_init(&_td_controller[1], _td_control_r2, _td_control_h2*h);
}

Vector3f ADRC_AttitudeControl::att_control(Vector3f err, float dt)
{
    Vector3f sp_rate;

    sp_rate(0) = td_control(&_td_controller[0], err(0), dt);
    sp_rate(1) = td_control(&_td_controller[1], err(1), dt);

    return sp_rate;
}

matrix::Vector3f ADRC_AttitudeControl::update(matrix::Quatf q, matrix::Quatf qd, const float yawspeed_feedforward, float dt)
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

    // calculate angular rates setpoint
//    matrix::Vector3f rate_setpoint = eq.emult(_proportional_gain);
    matrix::Vector3f rate_setpoint = att_control(eq, dt);
    rate_setpoint(2) = eq(2) * _proportional_gain(2);
//    PX4_INFO("%f, %f, %f", (double)rate_setpoint(0), (double)rate_setpoint(1), (double)rate_setpoint(2));


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

void ADRC_AttitudeControl::setProportionalGain(const matrix::Vector3f &proportional_gain)
{
    _proportional_gain = proportional_gain;

    // prepare yaw weight from the ratio between roll/pitch and yaw gains
    const float roll_pitch_gain = (proportional_gain(0) + proportional_gain(1)) / 2.f;
    _yaw_w = math::constrain(proportional_gain(2) / roll_pitch_gain, 0.f, 1.f);

    _proportional_gain(2) = roll_pitch_gain;
}

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
