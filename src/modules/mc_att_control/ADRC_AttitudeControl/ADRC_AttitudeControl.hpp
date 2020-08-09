/**
 * @file ADRC_AttitudeControl.hpp
 * @author Matthias Grob	<maetugr@gmail.com>
 */

#pragma once

#include <matrix/matrix/math.hpp>
#include <adrc/adrc.hpp>
#include <iostream>
#include <vector>

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

    void setGains(float td_control_r2, float td_control_h2);
    void att_init(float dt);
    matrix::Vector3f att_control(matrix::Vector3f err, float dt);

    /***************************************************************************************/

private:
	matrix::Vector3f _proportional_gain;
	matrix::Vector3f _rate_limit;
    float _yaw_w = 0.0f; /**< yaw weight [0,1] to prioritize roll and pitch */

    /*****************************************************************************************/

    float _td_control_r2;
    float _td_control_h2;

    std::vector<adrc::TD_Controller> _td_controller;

    float int_i[2] = {0.0f, 0.0f};

    /*****************************************************************************************/
};
