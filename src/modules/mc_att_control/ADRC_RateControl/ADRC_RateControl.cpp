/****************************************************************************
 *
 *   Copyright (c) 2019 PX4 Development Team. All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in
 *    the documentation and/or other materials provided with the
 *    distribution.
 * 3. Neither the name PX4 nor the names of its contributors may be
 *    used to endorse or promote products derived from this software
 *    without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
 * FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
 * COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
 * INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
 * BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS
 * OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED
 * AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 * LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
 * ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 *
 ****************************************************************************/

/**
 * @file ADRC_RateControl.cpp
 */

#include <ADRC_RateControl.hpp>
#include <px4_defines.h>

#include <mathlib/math/Limits.hpp>

#include <px4_log.h>
#include <px4_app.h>
#include <px4_tasks.h>

using namespace matrix;
using namespace adrc;

ADRC_RateControl::ADRC_RateControl()
{
    _td.push_back(TD());
    _td.push_back(TD());
    _leso.push_back(LESO());
    _leso.push_back(LESO());
    _nlsef.push_back(NLSEF());
    _nlsef.push_back(NLSEF());
}

void ADRC_RateControl::setGains(const Vector3f &P, const Vector3f &I, const Vector3f &D,
                                float td_r0, float leso_w, float nlsef_r1, float nlsef_h1, float nlsef_c, float gamma, float nlsef_ki)
{
	_gain_p = P;
	_gain_i = I;
	_gain_d = D;
    _td_r0 = td_r0;
    _leso_w = leso_w;
    _nlsef_r1 = nlsef_r1;
    _nlsef_h1 = nlsef_h1;
    _nlsef_c = nlsef_c;
    _gamma = gamma;
    _nlsef_ki = nlsef_ki;
}

void ADRC_RateControl::rate_init(float h)
{
    td_init(&_td[0], _td_r0, h);
    td_init(&_td[1], _td_r0, h);
    leso_init(&_leso[0], _leso_w, 400);
    leso_init(&_leso[1], _leso_w, 400);
    nlsef_init(&_nlsef[0], _nlsef_r1, _nlsef_h1*h, _nlsef_c);
    nlsef_init(&_nlsef[1], _nlsef_r1, _nlsef_h1*h, _nlsef_c);
}

Vector3f ADRC_RateControl::att_dis_comp(Vector2f in)
{
    Vector3f out;
    out(0) = in(0) - _gamma * _leso[0].z2 / _leso[0].b0;
    out(1) = in(1) - _gamma * _leso[1].z2 / _leso[1].b0;

    out(0) = math::constrain(out(0), -0.5f, 0.5f);
    out(1) = math::constrain(out(1), -0.5f, 0.5f);

    _leso[0].u = out(0);
    _leso[1].u = out(1);

    return out;
}

Vector3f ADRC_RateControl::rate_control(const Vector3f rate, const Vector3f rate_sp, const float dt)
{
    Vector3f rate_err;
    Vector2f u0;
    Vector3f torque;

    rate_err = rate_sp - rate;

    td(&_td[0], rate_err(0), dt);
    td(&_td[1], rate_err(1), dt);

    u0(0) = nlsef(&_nlsef[0], rate_err(0), _td[0].v2) / _leso[0].b0;
    u0(1) = nlsef(&_nlsef[1], rate_err(1), _td[1].v2) / _leso[1].b0;

    u0(0) = math::constrain(u0(0), -0.5f, 0.5f);
    u0(1) = math::constrain(u0(1), -0.5f, 0.5f);

    torque = att_dis_comp(u0);

    return torque;
}

void ADRC_RateControl::setDTermCutoff(const float loop_rate, const float cutoff, const bool force)
{
	// only do expensive filter update if the cutoff changed
	if (force || fabsf(_lp_filters_d.get_cutoff_freq() - cutoff) > 0.01f) {
		_lp_filters_d.set_cutoff_frequency(loop_rate, cutoff);
		_lp_filters_d.reset(_rate_prev);
	}
}

void ADRC_RateControl::setSaturationStatus(const MultirotorMixer::saturation_status &status)
{
	_mixer_saturation_positive[0] = status.flags.roll_pos;
	_mixer_saturation_positive[1] = status.flags.pitch_pos;
	_mixer_saturation_positive[2] = status.flags.yaw_pos;
	_mixer_saturation_negative[0] = status.flags.roll_neg;
	_mixer_saturation_negative[1] = status.flags.pitch_neg;
	_mixer_saturation_negative[2] = status.flags.yaw_neg;
}

Vector3f ADRC_RateControl::update(const Vector3f rate, const Vector3f rate_sp, const float dt, const bool landed)
{
	// angular rates error
	Vector3f rate_error = rate_sp - rate;

	// prepare D-term based on low-pass filtered rates
	Vector3f rate_filtered(_lp_filters_d.apply(rate));
	Vector3f rate_d;

	if (dt > FLT_EPSILON) {
		rate_d = (rate_filtered - _rate_prev_filtered) / dt;
	}

	// PID control with feed forward
    Vector3f torque = _gain_p.emult(rate_error) + _rate_int - _gain_d.emult(rate_d) + _gain_ff.emult(rate_sp);
    torque += rate_control(rate, rate_sp, dt);
//    torque(2) = _gain_p(2)*rate_error(2) + _rate_int(2) - _gain_d(2)*rate_d(2) + _gain_ff(2)*rate_sp(2);
    PX4_INFO("%f, %f, %f", (double)torque(0), (double)torque(1), (double)torque(2));

	_rate_prev = rate;
	_rate_prev_filtered = rate_filtered;

	// update integral only if we are not landed
	if (!landed) {
		updateIntegral(rate_error, dt);
	}

	return torque;
}

void ADRC_RateControl::updateIntegral(Vector3f &rate_error, const float dt)
{
	for (int i = 0; i < 3; i++) {
		// prevent further positive control saturation
		if (_mixer_saturation_positive[i]) {
			rate_error(i) = math::min(rate_error(i), 0.f);
		}

		// prevent further negative control saturation
		if (_mixer_saturation_negative[i]) {
			rate_error(i) = math::max(rate_error(i), 0.f);
		}

		// I term factor: reduce the I gain with increasing rate error.
		// This counteracts a non-linear effect where the integral builds up quickly upon a large setpoint
		// change (noticeable in a bounce-back effect after a flip).
		// The formula leads to a gradual decrease w/o steps, while only affecting the cases where it should:
		// with the parameter set to 400 degrees, up to 100 deg rate error, i_factor is almost 1 (having no effect),
		// and up to 200 deg error leads to <25% reduction of I.
		float i_factor = rate_error(i) / math::radians(400.f);
		i_factor = math::max(0.0f, 1.f - i_factor * i_factor);

		// Perform the integration using a first order method
		float rate_i = _rate_int(i) + i_factor * _gain_i(i) * rate_error(i) * dt;

		// do not propagate the result if out of range or invalid
		if (PX4_ISFINITE(rate_i)) {
			_rate_int(i) = math::constrain(rate_i, -_lim_int(i), _lim_int(i));
		}
	}
}

void ADRC_RateControl::getRateControlStatus(rate_ctrl_status_s &rate_ctrl_status)
{
	rate_ctrl_status.rollspeed_integ = _rate_int(0);
	rate_ctrl_status.pitchspeed_integ = _rate_int(1);
	rate_ctrl_status.yawspeed_integ = _rate_int(2);
}
