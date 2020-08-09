/**
 * @file ADRC_AttitudeControl.hpp
 * @author Matthias Grob	<maetugr@gmail.com>
 */

#pragma once

#include <matrix/matrix/math.hpp>
#include <adrc/adrc.hpp>
#include <iostream>
#include <vector>

using namespace adrc;

class ADRC_AttitudeControl
{
public:
    ADRC_AttitudeControl();
    ~ADRC_AttitudeControl() = default;

	/**
	 * Set proportional attitude control gain
	 * @param proportional_gain 3D vector containing gains for roll, pitch, yaw
	 */
	void setProportionalGain(const matrix::Vector3f &proportional_gain);

	/**
	 * Set hard limit for output rate setpoints
	 * @param rate_limit [rad/s] 3D vector containing limits for roll, pitch, yaw
	 */
	void setRateLimit(const matrix::Vector3f &rate_limit) { _rate_limit = rate_limit; }

	/**
	 * Run one control loop cycle calculation
	 * @param q estimation of the current vehicle attitude unit quaternion
	 * @param qd desired vehicle attitude setpoint
	 * @param yawspeed_feedforward [rad/s] yaw feed forward angular rate in world frame
	 * @return [rad/s] body frame 3D angular rate setpoint vector to be executed by the rate controller
	 */
    matrix::Vector3f update(matrix::Quatf q, matrix::Quatf qd, float yawspeed_feedforward, float dt);

    /*****************************************************************************************/

    void setGains(float td_control_r2, float td_control_h2, float td_r0, float leso_w, float nlsef_r1, float nlsef_h1, float nlsef_c, float gamma, float nlsef_ki);
    void att_init(float dt);
    void att_reset(float h);
    matrix::Vector3f att_dis_comp(matrix::Vector2f in);
    matrix::Vector3f att_control(matrix::Vector3f err, float dt);
    void add_observer_update(const float gyr[3], float bth);

    /***************************************************************************************/

private:
	matrix::Vector3f _proportional_gain;
	matrix::Vector3f _rate_limit;
    float _yaw_w = 0.0f; /**< yaw weight [0,1] to prioritize roll and pitch */

    /*****************************************************************************************/

    float _td_control_r2;
    float _td_control_h2;
    float _td_r0;
    float _leso_w;
    float _nlsef_r1;
    float _nlsef_h1;
    float _nlsef_c;
    float _gamma;
    float _nlsef_ki;

    std::vector<TD_Controller> _td_controller;
    std::vector<TD> _td;
    std::vector<NLSEF> _nlsef;
    std::vector<LESO> _leso;

    float int_i[2] = {0.0f, 0.0f};

    float sign(float val);
    float fhan(float v1, float v2, float r0, float h0);
    float fal(float e, float alpha, float delta);
    void td_init(TD* td_t, float r0, float h0);
    void td(TD* td, float v, float dt);
    void td_control_init(TD_Controller* td_controller, float r2, float h2);
    float td_control(TD_Controller* td_controller, float err, float dt);
    void eso_init(ESO* eso_t, float beta1, float beta2, float alpha, float delta, float b0);
    void eso(ESO* eso_t, float y, float dt);
    void leso_init(LESO* leso_t, float w, float b0);
    void leso(LESO* leso_t, float y, float dt);
    void nlsef_init(NLSEF* nlsef_t, float r1, float h1, float c);
    float nlsef(NLSEF* nlsef_t, float e1, float e2);

    /*****************************************************************************************/
};
