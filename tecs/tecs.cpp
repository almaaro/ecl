/****************************************************************************
 *
 *   Copyright (c) 2017 Estimation and Control Library (ECL). All rights reserved.
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
 * 3. Neither the name ECL nor the names of its contributors may be
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

#include "tecs.h"

#include <ecl.h>
#include <geo/geo.h>

using math::constrain;
using math::max;
using math::min;

static constexpr float DT_MIN = 0.001f;	///< minimum allowed value of _dt (sec)
static constexpr float DT_MAX = 1.0f;	///< max value of _dt allowed before a filter state reset is performed (sec)

/**
 * @file tecs.cpp
 *
 * @author Paul Riseborough
 */

/*
 * This function implements a complementary filter to estimate the climb rate when
 * inertial nav data is not available. It also calculates a true airspeed derivative
 * which is used by the airspeed complimentary filter.
 */
void TECS::update_vehicle_state_estimates(float airspeed, const matrix::Dcmf &rotMat,
		const matrix::Vector3f &accel_body, bool altitude_lock, bool in_air,
		float altitude, float vz)
{
	// calculate the time lapsed since the last update
	uint64_t now = ecl_absolute_time();
	float dt = constrain((now - _state_update_timestamp) * 1.0e-6f, DT_MIN, DT_MAX);

	bool reset_altitude = false;

	if (_state_update_timestamp == 0 || dt > DT_MAX) {
		dt = DT_DEFAULT;
		reset_altitude = true;
	}

	if (!altitude_lock || !in_air) {
		reset_altitude = true;
	}

	if (reset_altitude) {
		_states_initialized = false;
	}

	_state_update_timestamp = now;
	_EAS = airspeed;

	_in_air = in_air;

	// Set the velocity and position state to the the INS data
	_vert_vel_state = -vz;
	_vert_pos_state = altitude;

	// Update and average speed rate of change if airspeed is being measured
	if (ISFINITE(airspeed) && airspeed_sensor_enabled()) {
		// Assuming the vehicle is flying X axis forward, use the X axis measured acceleration
		// compensated for gravity to estimate the rate of change of speed
		float speed_deriv_raw = rotMat(2, 0) * CONSTANTS_ONE_G + accel_body(0);

		// Apply some noise filtering
		_speed_derivative = 0.95f * _speed_derivative + 0.05f * speed_deriv_raw;

	} else {
		_speed_derivative = 0.0f;
	}

	if (!_in_air) {
		_states_initialized = false;
	}

}

void TECS::_update_speed_states(float airspeed_setpoint, float indicated_airspeed, float EAS2TAS)
{
	// Calculate the time in seconds since the last update and use the default time step value if out of bounds
	uint64_t now = ecl_absolute_time();
	const float dt = constrain((now - _speed_update_timestamp) * 1.0e-6f, DT_MIN, DT_MAX);

	// Convert equivalent airspeed quantities to true airspeed
	_EAS_setpoint = airspeed_setpoint;
	_TAS_setpoint  = _EAS_setpoint * EAS2TAS;
	_TAS_max   = _indicated_airspeed_max * EAS2TAS;
	_TAS_trim  = _indicated_airspeed_trim * EAS2TAS;
	_TAS_min   = _indicated_airspeed_min * EAS2TAS;

	// If airspeed measurements are not being used, fix the airspeed estimate to halfway between
	// min and max limits
	if (!ISFINITE(indicated_airspeed) || !airspeed_sensor_enabled()) {
		_EAS = 0.5f * (_indicated_airspeed_min + _indicated_airspeed_max);

	} else {
		_EAS = indicated_airspeed;
	}

	// If first time through or not flying, reset airspeed states
	if (_speed_update_timestamp == 0 || !_in_air) {
		_tas_rate_state = 0.0f;
		_tas_state = (_EAS * EAS2TAS);
	}

	// Obtain a smoothed airspeed estimate using a second order complementary filter

	// Update TAS rate state
	float tas_error = (_EAS * EAS2TAS) - _tas_state;
	float tas_rate_state_input = tas_error * _tas_estimate_freq * _tas_estimate_freq;

	// limit integrator input to prevent windup
	if (_tas_state < 3.1f) {
		tas_rate_state_input = max(tas_rate_state_input, 0.0f);

	}

	// Get smoothed _EAS from smoothed _tas_state
	_EAS = _tas_state / EAS2TAS;

	// Update TAS state
	_tas_rate_state = _tas_rate_state + tas_rate_state_input * dt;
	float tas_state_input = _tas_rate_state + _speed_derivative + tas_error * _tas_estimate_freq * 1.4142f;
	_tas_state = _tas_state + tas_state_input * dt;

	// Limit the airspeed state to a minimum of 3 m/s
	_tas_state = max(_tas_state, 3.0f);
	_speed_update_timestamp = now;

}

void TECS::_update_speed_setpoint()
{
	// Set the airspeed demand to the minimum value if an
	// uncontrolled descent condition exists to maximise climb rate
	if (_uncommanded_descent_recovery){
		_TAS_setpoint = 1.1f * _TAS_min;
	}

	// Calculate limits for the demanded rate of change of speed based on physical performance limits
	// with a 50% margin to allow the total energy controller to correct for errors.
	float velRateMax = _STE_rate_max / _tas_state;
	float velRateMin = _STE_rate_min / _tas_state;

	_TAS_setpoint_adj = constrain(_TAS_setpoint, _TAS_min, _TAS_max);

	// calculate the demanded rate of change of speed proportional to speed error
	// and apply performance limits
	_TAS_rate_setpoint = constrain((_TAS_setpoint_adj - _tas_state) * _speed_error_gain, velRateMin, velRateMax);

}

void TECS::_update_height_setpoint(float desired, float state)
{
	// Detect first time through and initialize previous value to demand
	if (ISFINITE(desired) && fabsf(_hgt_setpoint_in_prev) < 0.1f) {
		_hgt_setpoint_in_prev = desired;
	}

	// Apply a 2 point moving average to demanded height to reduce
	// intersampling noise effects.
	if (ISFINITE(desired)) {
		_hgt_setpoint = 0.5f * (desired + _hgt_setpoint_in_prev);

	} else {
		_hgt_setpoint = _hgt_setpoint_in_prev;
	}

	_hgt_setpoint_in_prev = _hgt_setpoint;

	// Apply a rate limit to respect vehicle performance limitations
	if ((_hgt_setpoint - _hgt_setpoint_prev) > (_max_climb_rate * _dt)) {
		_hgt_setpoint = _hgt_setpoint_prev + _max_climb_rate * _dt;

	} else if ((_hgt_setpoint - _hgt_setpoint_prev) < (-_max_sink_rate * _dt)) {
		_hgt_setpoint = _hgt_setpoint_prev - _max_sink_rate * _dt;
	}

	_hgt_setpoint_prev = _hgt_setpoint;

	// Apply a first order noise filter
	_hgt_setpoint_adj = 0.1f * _hgt_setpoint + 0.9f * _hgt_setpoint_adj_prev;

	// Calculate the demanded climb rate proportional to height error plus a feedforward term to provide
	// tight tracking during steady climb and descent manoeuvres.
	_hgt_rate_setpoint = (_hgt_setpoint_adj - state) * _height_error_gain + _height_setpoint_gain_ff *
				 (_hgt_setpoint_adj - _hgt_setpoint_adj_prev) / _dt;
	_hgt_setpoint_adj_prev = _hgt_setpoint_adj;

	// if externally demanded height rate (=flare)
	if(_use_position_control_hgt_rate) {
		_hgt_rate_setpoint = _position_control_hgt_rate;
	}

	// Limit the rate of change of height demand to respect vehicle performance limits
	if (_hgt_rate_setpoint > _max_climb_rate) {
		_hgt_rate_setpoint = _max_climb_rate;

	} else if (_hgt_rate_setpoint < -_max_sink_rate) {
		_hgt_rate_setpoint = -_max_sink_rate;
	}

	// Limit the rate of change of height demand to respect vehicle performance limits
	_hgt_rate_setpoint = constrain(_hgt_rate_setpoint, -_max_sink_rate, _max_climb_rate);
}

void TECS::_detect_underspeed()
{
	if (!_detect_underspeed_enabled) {
		_underspeed_detected = false;
		return;
	}

	if (_tas_state < 0.9f * _TAS_min || ((_vert_pos_state < _hgt_setpoint_adj || _tas_state < (0.5f * (_TAS_setpoint + _TAS_min))) && _underspeed_detected)) {

		_underspeed_detected = true;

	} else {
		_underspeed_detected = false;
	}
}

void TECS::_update_energy_estimates()
{
	// Calculate specific energy demands in units of (m**2/sec**2)
	_SPE_setpoint = _hgt_setpoint_adj * CONSTANTS_ONE_G; // potential energy
	_SKE_setpoint = 0.5f * _TAS_setpoint_adj * _TAS_setpoint_adj; // kinetic energy

	// Calculate specific energy rate demands in units of (m**2/sec**3)
	_SPE_rate_setpoint = _hgt_rate_setpoint * CONSTANTS_ONE_G; // potential energy rate of change
	_SKE_rate_setpoint = _tas_state * _TAS_rate_setpoint; // kinetic energy rate of change

	// Calculate specific energies in units of (m**2/sec**2)
	_SPE_estimate = _vert_pos_state * CONSTANTS_ONE_G; // potential energy
	_SKE_estimate = 0.5f * _tas_state * _tas_state; // kinetic energy

	// Calculate specific energy rates in units of (m**2/sec**3)
	_SPE_rate = _vert_vel_state * CONSTANTS_ONE_G; // potential energy rate of change
	_SKE_rate = _tas_state * _speed_derivative;// kinetic energy rate of change
}

void TECS::_update_throttle_setpoint(const float throttle_cruise, const matrix::Dcmf &rotMat)
{
	// Calculate total energy error
	_STE_error = _SPE_setpoint - _SPE_estimate + _SKE_setpoint - _SKE_estimate;

	// The STE rate setpoint must now be constrained so that we will never end up in a situation where
	// the total energy is OK, but the airspeed is dangerously low because we have too much potential energy.
	// This is prevented by first bleeding off any excess potential energy into kinetic energy and then reducing the excess
	// airspeed. However, if the pitch is controlling also the airspeed, the risk of this underspeeding is reduced.
	float STE_rate_min_adj = (0.5f * _pitch_speed_weight) * _STE_rate_min + (1.0f - 0.5f * _pitch_speed_weight) * _SKE_rate_setpoint;

	//Also adjust the safety margin to the current airspeed. Full effect at min airspeed, no effect at trim airspeed.
	float scaler = constrain((_EAS - _indicated_airspeed_min) / (max(0.1f, _indicated_airspeed_trim - _indicated_airspeed_min)), 0.0f, 1.0f);
	STE_rate_min_adj = (1.0f - scaler) * STE_rate_min_adj + scaler * _STE_rate_min;

	// Calculate demanded rate of change of total energy, respecting vehicle limits
	float STE_rate_setpoint = constrain((_SPE_rate_setpoint + _SKE_rate_setpoint), STE_rate_min_adj, _STE_rate_max);

	// Calculate the total energy rate error, applying a first order IIR filter
	// to reduce the effect of accelerometer noise
	_STE_rate_error = 0.2f * (STE_rate_setpoint - _SPE_rate - _SKE_rate) + 0.8f * _STE_rate_error;

	// Change to feed forward and proportional control
	STE_rate_setpoint = STE_rate_setpoint + _STE_rate_error * _throttle_damping_gain;

	// Calculate the throttle demand

	if (_throttle_integrator_gain > 0.0f) {
			// Calculate throttle integrator state upper and lower limits with allowance for
			// 10% throttle saturation to accommodate noise on the demand.
			float integrator_margin = 0.1f * (_STE_rate_max - STE_rate_min_adj);
			float integ_state_max = _STE_rate_max - STE_rate_setpoint + integrator_margin;
			float integ_state_min = STE_rate_min_adj - STE_rate_setpoint  - integrator_margin;

			if (_climbout_mode_active) {
			// During climbout, set the integrator to maximum throttle to prevent transient throttle drop
			// at end of climbout when we transition to closed loop throttle control
			_STE_integ_state = integ_state_max;

			} else {
					// If the throttle would be over max/min then decay the integrator to settle the compensated value at max/min throttle
					// (maximum effect to _STE_integ_state is +-0.1). Same method as with the pitch integrator. The reason is to
					// prevent _STE_integ_state from remaining at a very low value for example during descends that require zero throttle
					// and drained _STE_integ_state at the start of the descend. This will reduce airspeed undershoot at the end of the
					// descend.

					// Else add to the integrator value normally.

					float ste_integ = STE_rate_setpoint + _STE_integ_state;
					if (ste_integ > _STE_rate_max) {
							_STE_integ_state = _STE_integ_state - (ste_integ - _STE_rate_max) * _dt;

					} else if (ste_integ < STE_rate_min_adj) {
							_STE_integ_state = _STE_integ_state - (ste_integ - STE_rate_min_adj) * _dt;

					} else {
							_STE_integ_state = _STE_integ_state + (_STE_rate_error * _throttle_integrator_gain) * _dt;
					}

					// Respect integrator limits during closed loop operation.
					_STE_integ_state = constrain(_STE_integ_state, integ_state_min, integ_state_max);
			}

	} else {
			_STE_integ_state = 0.0f;
	}

	STE_rate_setpoint = constrain(STE_rate_setpoint + _STE_integ_state, STE_rate_min_adj, _STE_rate_max);

	// Calculate a predicted throttle from the demanded rate of change of energy, using the cruise throttle
	// as the starting point. Assume:
	// Specific total energy rate = _STE_rate_max is achieved when throttle is set to _throttle_setpoint_max
	// Specific total energy rate = 0 at cruise throttle
	// Specific total energy rate = _STE_rate_min is achieved when throttle is set to _throttle_setpoint_min
	float throttle_predicted = 0.0f;

	if (_advanced_thr_calc_initialized) {
			const float as_squared = _EAS * _EAS;
			/* Throttle calculations */

			// This is (delta v at max throttle at  _adv_thr_calc_airspeed) / (v2 at max throttle at _indicated_airspeed_trim). Used to scale the required delta v for throttle.
			// The adjusted delt v is calculated from thrust, which is assumed to have a linear relationship to airspeed.
			const float max_delta_v_airspeed_coefficient = (sqrtf(as_squared + (_max_thrust_as_coefficient *
																			   (_EAS - _indicated_airspeed_trim) + _thrust_trim_as_max_climb)
															/ _thrust_coefficient) - _EAS) / _delta_v_trim_as_max_climb;

			// required thrust to overcome air drag, climb rate and acceleration. Also de-normalizimg this from _auw.
			const float required_thrust = (_cd_i_specific / as_squared + _cd_o_specific * as_squared + STE_rate_setpoint / _tas_state) * _auw;

			// The calculated delta v to produce the required thrust at the current airspeed
			_required_delta_v = sqrtf(max(0.001f, required_thrust / _thrust_coefficient + as_squared)) - _EAS;
			_required_as_elev = _required_delta_v * _motor_airstream_at_elevator_scaler + _EAS;

			// Adjusting the delta v to match the new maximum delta v at the current airspeed
			const float delta_v_trim_as_level_adj = _delta_v_trim_as_level * max_delta_v_airspeed_coefficient;
			const float delta_v_trim_as_max_climb_adj = _delta_v_trim_as_max_climb * max_delta_v_airspeed_coefficient;

			// Getting the required throttle by interpolating and finding the required delta v.
			if (_required_delta_v > delta_v_trim_as_level_adj) {
					throttle_predicted = (_required_delta_v - delta_v_trim_as_level_adj) / (delta_v_trim_as_max_climb_adj - delta_v_trim_as_level_adj) *
											 (_throttle_setpoint_max - throttle_cruise) + throttle_cruise;
			} else {
					throttle_predicted = (_required_delta_v / delta_v_trim_as_level_adj) * (throttle_cruise - _throttle_setpoint_min) + _throttle_setpoint_min;
			}
	}
	else{
			// Adjust the demanded total energy rate to compensate for induced drag rise in turns.
			// Assume induced drag scales linearly with normal load factor.
			// The additional normal load factor is given by (1/cos(bank angle) - 1)
			float cosPhi = sqrtf((rotMat(0, 1) * rotMat(0, 1)) + (rotMat(1, 1) * rotMat(1, 1)));
			STE_rate_setpoint = STE_rate_setpoint + _load_factor_correction * (1.0f / constrain(cosPhi, 0.1f, 1.0f) - 1.0f);

			if (STE_rate_setpoint >= 0) {
					// throttle is between cruise and maximum
					throttle_predicted = throttle_cruise + STE_rate_setpoint / _STE_rate_max * (_throttle_setpoint_max - throttle_cruise);

			} else {
					// throttle is between cruise and minimum
					throttle_predicted = throttle_cruise + STE_rate_setpoint / _STE_rate_min * (_throttle_setpoint_min - throttle_cruise);

			}
	}

	// Constrain to throttle limits
	_throttle_setpoint = throttle_predicted;
	_throttle_setpoint = constrain(_throttle_setpoint, _throttle_setpoint_min, _throttle_setpoint_max);

	// Rate limit the throttle demand
	if (fabsf(_throttle_slewrate) > 0.01f) {
		float throttle_increment_limit = _dt * (_throttle_setpoint_max - _throttle_setpoint_min) * _throttle_slewrate;
		_throttle_setpoint = constrain(_throttle_setpoint, _last_throttle_setpoint - throttle_increment_limit,
						   _last_throttle_setpoint + throttle_increment_limit);
	}

	_last_throttle_setpoint = _throttle_setpoint;

}

void TECS::_detect_uncommanded_descent()
{
	/*
	 * This function detects a condition that can occur when the demanded airspeed is greater than the
	 * aircraft can achieve in level flight. When this occurs, the vehicle will continue to reduce height
	 * while attempting to maintain speed.
	*/

	// Calculate rate of change of total specific energy
	float STE_rate = _SPE_rate + _SKE_rate;

	// If total energy is very low and reducing, throttle is high, and we are not in an underspeed condition, then enter uncommanded descent recovery mode
	bool enter_mode = !_uncommanded_descent_recovery && !_underspeed_detected && (_STE_error > 200.0f) && (STE_rate < 0.0f)
			  && (_throttle_setpoint >= _throttle_setpoint_max * 0.9f);

	// If we enter an underspeed condition or recover the required total energy, then exit uncommanded descent recovery mode
	bool exit_mode = _uncommanded_descent_recovery && (_underspeed_detected || (_STE_error < 0.0f));

	if (enter_mode) {
		_uncommanded_descent_recovery = true;

	} else if (exit_mode) {
		_uncommanded_descent_recovery = false;

	}
}

void TECS::_update_pitch_setpoint()
{
	/*
	 * The SKE_weighting variable controls how speed and height control are prioritised by the pitch demand calculation.
	 * A weighting of 1 givea equal speed and height priority
	 * A weighting of 0 gives 100% priority to height control and must be used when no airspeed measurement is available.
	 * A weighting of 2 provides 100% priority to speed control and is used when:
	 * a) an underspeed condition is detected.
	 * b) during climbout where a minimum pitch angle has been set to ensure height is gained. If the airspeed
	 * rises above the demanded value, the pitch angle demand is increased by the TECS controller to prevent the vehicle overspeeding.
	 * The weighting can be adjusted between 0 and 2 depending on speed and height accuracy requirements.
	*/

	// Calculate the weighting applied to control of specific kinetic energy error
	float SKE_weighting = constrain(_pitch_speed_weight, 0.0f, 2.0f);

	if ((_underspeed_detected || _climbout_mode_active) && airspeed_sensor_enabled()) {
		SKE_weighting = 2.0f;

	} else if (!airspeed_sensor_enabled()) {
		SKE_weighting = 0.0f;
	}

	// Calculate the weighting applied to control of specific potential energy error
	float SPE_weighting = 2.0f - SKE_weighting;

	SPE_weighting = min(SPE_weighting, 1.0f);
	SKE_weighting = min(SKE_weighting, 1.0f);

	// Calculate the specific energy balance demand which specifies how the available total
	// energy should be allocated to speed (kinetic energy) and height (potential energy)
	float SEB_setpoint = _SPE_setpoint * SPE_weighting - _SKE_setpoint * SKE_weighting;

	// Because the airspeed is more critical for an airplane than the altitude, limit the _SPE_rate_setpoint
	// so that the demanded rate in airspeed can be achieved.
	float SPE_rate_setpoint_adj = min(_SPE_rate_setpoint, _STE_rate_max - _SKE_rate_setpoint);

	// Calculate the specific energy balance rate demand
	float SEB_rate_setpoint = SPE_rate_setpoint_adj * SPE_weighting - _SKE_rate_setpoint * SKE_weighting;

	// Calculate the specific energy balance and balance rate error
	_SEB_error = SEB_setpoint - (_SPE_estimate * SPE_weighting - _SKE_estimate * SKE_weighting);
	_SEB_rate_error = SEB_rate_setpoint - (_SPE_rate * SPE_weighting - _SKE_rate * SKE_weighting);

	// Calculate derivative from change in climb angle to rate of change of specific energy balance
	float climb_angle_to_SEB_rate = _tas_state * CONSTANTS_ONE_G;

		if (_pitch_integrator_gain > 0.0f) {
				// Prevent the integrator changing in a direction that will increase pitch demand saturation
				// Decay the integrator if the pitch demand from the previous time step is saturated
				if (_pitch_setpoint_unc > _pitch_setpoint_max) {
						_pitch_integ_state = _pitch_integ_state - (_pitch_setpoint_unc - _pitch_setpoint_max) * climb_angle_to_SEB_rate * _dt;

				} else if (_pitch_setpoint_unc < _pitch_setpoint_min) {
						_pitch_integ_state = _pitch_integ_state - (_pitch_setpoint_unc - _pitch_setpoint_min) * climb_angle_to_SEB_rate * _dt;

				} else {
						_pitch_integ_state = _pitch_integ_state + _SEB_rate_error * _pitch_integrator_gain * _dt;
				}

		} else {
				_pitch_integ_state = 0.0f;
		}

	// Calculate a specific energy correction that doesn't include the integrator contribution
	float SEB_correction = _SEB_rate_error * _pitch_damping_gain + SEB_rate_setpoint;

	// During climbout, bias the demanded pitch angle so that a zero speed error produces a pitch angle
	// demand equal to the minimum pitch angle set by the mission plan. This prevents the integrator
	// having to catch up before the nose can be raised to reduce excess speed during climbout.
	if (_climbout_mode_active) {
		SEB_correction += _pitch_setpoint_min * climb_angle_to_SEB_rate;
	}

	// Sum the correction terms and convert to a pitch angle demand. This calculation assumes:
	// a) The climb angle follows pitch angle with a lag that is small enough not to destabilise the control loop.
	// b) The offset between climb angle and pitch angle (angle of attack) is constant, excluding the effect of
	// pitch transients due to control action or turbulence.
	_pitch_setpoint_unc = (SEB_correction + _pitch_integ_state) / climb_angle_to_SEB_rate;

	/* Calculate and add the pitch offsets */

	// never calculate the offset for airspeeds lower than the minimum. This will avoid increasing the AOA in stalls.
	float EAS_adj = max(_indicated_airspeed_min, _EAS);
	float cl = _cl_coefficient / (EAS_adj * EAS_adj);

	//Then calculate the needed pitch. Take the flap setting into account.
	float offset = (1.0f - _flaps_applied) * _pitchsp_offset_rad + _flaps_applied * _pitchsp_offset_flaps_rad;
	float psp_offset_adj = offset + _cl_to_alpha_rad_slope * (cl - _cl_cruise_trim_as);

	_pitch_setpoint_unc += psp_offset_adj;

	_pitch_setpoint = constrain(_pitch_setpoint_unc, _pitch_setpoint_min, _pitch_setpoint_max);

	// Comply with the specified vertical acceleration limit by applying a pitch rate limit
	float ptchRateIncr = _dt * _vert_accel_limit / _tas_state;

	if ((_pitch_setpoint - _last_pitch_setpoint) > ptchRateIncr) {
		_pitch_setpoint = _last_pitch_setpoint + ptchRateIncr;

	} else if ((_pitch_setpoint - _last_pitch_setpoint) < -ptchRateIncr) {
		_pitch_setpoint = _last_pitch_setpoint - ptchRateIncr;
	}

	_last_pitch_setpoint = _pitch_setpoint;

}

void TECS::_initialize_states(float pitch, float throttle_cruise, float baro_altitude, float pitch_min_climbout,
				  float EAS2TAS)
{
	if (_pitch_update_timestamp == 0 || _dt > DT_MAX || !_in_air || !_states_initialized) {
		// On first time through or when not using TECS of if there has been a large time slip,
		// states must be reset to allow filters to a clean start
		_vert_vel_state = 0.0f;
		_vert_pos_state = baro_altitude;
		_tas_rate_state = 0.0f;
		_tas_state = _EAS * EAS2TAS;
		_STE_integ_state =  0.0f;
		_pitch_integ_state = 0.0f;
		_last_throttle_setpoint = (_in_air ? throttle_cruise : 0.0f);;
		_last_pitch_setpoint = constrain(pitch, _pitch_setpoint_min, _pitch_setpoint_max);
		_pitch_setpoint_unc = _last_pitch_setpoint;
		_hgt_setpoint_adj_prev = baro_altitude;
		_hgt_setpoint_adj = _hgt_setpoint_adj_prev;
		_hgt_setpoint_prev = _hgt_setpoint_adj_prev;
		_hgt_setpoint_in_prev = _hgt_setpoint_adj_prev;
		_TAS_setpoint_last = _EAS * EAS2TAS;
		_TAS_setpoint_adj = _TAS_setpoint_last;
		_underspeed_detected = false;
		_uncommanded_descent_recovery = false;
		_STE_rate_error = 0.0f;

		if (_dt > DT_MAX || _dt < DT_MIN) {
			_dt = DT_DEFAULT;
		}

	} else if (_climbout_mode_active) {
		// During climbout use the lower pitch angle limit specified by the
		// calling controller
		_pitch_setpoint_min	   = pitch_min_climbout;

		// throttle lower limit is set to a value that prevents throttle reduction
		_throttle_setpoint_min  = _throttle_setpoint_max - 0.01f;

		// height demand and associated states are set to track the measured height
		_hgt_setpoint_adj_prev  = baro_altitude;
		_hgt_setpoint_adj       = _hgt_setpoint_adj_prev;
		_hgt_setpoint_prev      = _hgt_setpoint_adj_prev;

		// airspeed demand states are set to track the measured airspeed
		_TAS_setpoint_last      = _EAS * EAS2TAS;
		_TAS_setpoint_adj       = _EAS * EAS2TAS;

		// disable speed and decent error condition checks
		_underspeed_detected = false;
		_uncommanded_descent_recovery = false;
	}

	_states_initialized = true;
}

void TECS::_update_STE_rate_lim(float throttle_cruise, const matrix::Dcmf &rotMat)
{
	// Calculate the specific total energy upper rate limits from the max throttle climb rate
	const float rate_max = _max_climb_rate * CONSTANTS_ONE_G;

	// Calculate the specific total energy lower rate limits from the min throttle sink rate
	const float rate_min = - _min_sink_rate * CONSTANTS_ONE_G;

	const float rate_min_flaps = - (_flaps_applied * _min_sink_rate_flaps + (1.0f - _flaps_applied) * _min_sink_rate) * CONSTANTS_ONE_G;

	if (_use_advanced_thr_calculation) {

		// Some error checks
		if (_auw < 0.01f || _wingspan < 0.01f || _indicated_airspeed_trim > 0.95f * _indicated_airspeed_max ||
				throttle_cruise > 0.95f * _throttle_setpoint_max){
			goto throttle_calculation_default;
		}

		/* Airspeed-dependent drag coefficients */
		// All of these calculations assume that the air density is constant (sea level 15C).
		// This also assumes that in the conditions where the parameters such as CLIMB_MAX, SINK_MIN and ASPD_MAX were measured,
		// EAS2TAS was 1 meaning that all _indicated_airspeed_[min|trim|max] were equal to their TAS values.

		// The additional normal load factor is given by 1/cos(bank angle)
		float cosPhi = sqrtf((rotMat(0, 1) * rotMat(0, 1)) + (rotMat(1, 1) * rotMat(1, 1)));
		const float lift = _auw * CONSTANTS_ONE_G / constrain(cosPhi, 0.1f, 1.0f);

		// _Cd_i_specific: Vehicle specific induced drag coefficient, which equals to 1/2*S*rho*Cd_i
		// Cd_i_specific = ... assuming planar wing. Efficiency factor of 0.85 for a starting point (should possibly be a parameter).
		_cd_i_specific = lift * lift / (0.5f * M_PI_F * CONSTANTS_AIR_DENSITY_SEA_LEVEL_15C * _wingspan * _wingspan * _auw * 0.85f);

		// _Cd_o_specific: Vehicle specific parasitic drag coefficient, which equals to 1/2*A*rho*Cd_o
		// Cd_o_specific: subtracting induced drag from total drag at a known airspeed to calculate parasitic drag
		_cd_o_specific = (-rate_min_flaps - _cd_i_specific / _indicated_airspeed_trim) /
				(_indicated_airspeed_trim * _indicated_airspeed_trim * _indicated_airspeed_trim);

		/* Throttle calculation */
		if (!_advanced_thr_calc_initialized){
			// Calculated thrust at different throttle settings and airspeeds.
			_thrust_trim_as_max_climb = (rate_max - rate_min) / _indicated_airspeed_trim * _auw;
			const float thrust_trim_as_level = -rate_min / _indicated_airspeed_trim * _auw;
			const float thrust_max_as_level = (_cd_i_specific / _indicated_airspeed_max / _indicated_airspeed_max +
							   _cd_o_specific * _indicated_airspeed_max * _indicated_airspeed_max) * _auw;

			// This tells how much the maximum thrust changes as airspeed changes. Assuming a linear relationship.
			// Sources: https://www.electricrcaircraftguy.com/2014/04/propeller-static-dynamic-thrust-equation-background.html (part "Application & Conjecturing")
			// But! Not that simple. The relationship is not linear if the airspeed decreases too much, ie. not for static thrust. Probably because the
			// propeller stalling. No source for this stall hypothesis. Picture: https://www.researchgate.net/figure/Typical-propeller-thrust-curves-as-a-function-of-advance-ratio-J_fig4_281946347
			_max_thrust_as_coefficient = (thrust_max_as_level - _thrust_trim_as_max_climb) / (_indicated_airspeed_max - _indicated_airspeed_trim);

			// v2 is the speed of airflow after the propeller. Source: https://www.grc.nasa.gov/www/k-12/airplane/propth.html
			// Source states that F = 0.5*PI*(d/2)^2*rho*(v2^2 - airspeed^2) where d is the propeller diameter and rho is the air density.
			_thrust_coefficient = 0.125f * M_PI_F * (_propeller_diameter) * (_propeller_diameter) * CONSTANTS_AIR_DENSITY_SEA_LEVEL_15C;
			_delta_v_trim_as_level = sqrtf(_indicated_airspeed_trim * _indicated_airspeed_trim + thrust_trim_as_level
								 / _thrust_coefficient) - _indicated_airspeed_trim;
			_delta_v_trim_as_max_climb = sqrtf(_indicated_airspeed_trim * _indicated_airspeed_trim + _thrust_trim_as_max_climb
								 / _thrust_coefficient) - _indicated_airspeed_trim;

			// Some more error checks
			if (_thrust_coefficient > 0.001f && _delta_v_trim_as_max_climb > 0.1f){

				_as_elev_trim_as_level_sq = _delta_v_trim_as_level * _motor_airstream_at_elevator_scaler + _indicated_airspeed_trim;
				_as_elev_trim_as_level_sq = _as_elev_trim_as_level_sq * _as_elev_trim_as_level_sq;

				_as_elev_max_as_level_sq = (sqrtf(thrust_max_as_level / _thrust_coefficient + _indicated_airspeed_max * _indicated_airspeed_max)
											- _indicated_airspeed_max) * _motor_airstream_at_elevator_scaler + _indicated_airspeed_max;
				_as_elev_max_as_level_sq = _as_elev_max_as_level_sq * _as_elev_max_as_level_sq;

				_advanced_thr_calc_initialized = true;
			} else {
				goto throttle_calculation_default;
			}
		}

		if (_EAS > 1.0f) {
			// _STE_rate_min equals to the sum of parasitic and induced drag power.
			// Drag force = _Cd_i / _EAS /_EAS + _Cd_o_specific * _EAS *_EAS;
			// Drag power = Drag force * _tas_state

			float EAS_sq = _EAS * _EAS;
			_STE_rate_min = - (_cd_i_specific / EAS_sq + _cd_o_specific * EAS_sq) * _tas_state;

			_STE_rate_max = ((_EAS - _indicated_airspeed_trim) * _max_thrust_as_coefficient + _thrust_trim_as_max_climb) * _tas_state / _auw + _STE_rate_min;

		} else {
			goto throttle_calculation_default;
		}
	}
	else {
throttle_calculation_default:

		_advanced_thr_calc_initialized = false;

		_as_elev_trim_as_level_sq = -1.0f;

		_STE_rate_max = rate_max;
		_STE_rate_min = rate_min;
	}
}

void TECS::update_pitch_throttle(const matrix::Dcmf &rotMat, float pitch, float baro_altitude, float hgt_setpoint,
				 float EAS_setpoint, float indicated_airspeed, float eas_to_tas, bool climb_out_setpoint, float pitch_min_climbout,
				 float throttle_min, float throttle_max, float throttle_cruise, float pitch_limit_min, float pitch_limit_max)
{
	// Calculate the time since last update (seconds)
	uint64_t now = ecl_absolute_time();
	_dt = constrain((now - _pitch_update_timestamp) * 1e-6f, DT_MIN, DT_MAX);

	// Set class variables from inputs
	_throttle_setpoint_max = throttle_max;
	_throttle_setpoint_min = throttle_min;
	_pitch_setpoint_max = pitch_limit_max;
	_pitch_setpoint_min = pitch_limit_min;
	_climbout_mode_active = climb_out_setpoint;

	// Initialize selected states and variables as required
	_initialize_states(pitch, throttle_cruise, baro_altitude, pitch_min_climbout, eas_to_tas);

	// Don't run TECS control algorithms when not in flight
	if (!_in_air) {
		return;
	}

	// Update the true airspeed state estimate
	_update_speed_states(EAS_setpoint, indicated_airspeed, eas_to_tas);

	// Calculate rate limits for specific total energy
	_update_STE_rate_lim(throttle_cruise, rotMat);

	// Detect an underspeed condition
	_detect_underspeed();

	// Detect an uncommanded descent caused by an unachievable airspeed demand
	_detect_uncommanded_descent();

	// Calculate the demanded true airspeed
	_update_speed_setpoint();

	// Calculate the demanded height
	_update_height_setpoint(hgt_setpoint, baro_altitude);

	// Calculate the specific energy values required by the control loop
	_update_energy_estimates();

	//initialize pitch setpoint offsets if needed
	_initialize_pitchsp_offset();

	// Calculate the pitch demand
	_update_pitch_setpoint();

	// Calculate the throttle demand
	_update_throttle_setpoint(throttle_cruise, rotMat);

	// Update time stamps
	_pitch_update_timestamp = now;

	// Set TECS mode for next frame
	if (_underspeed_detected) {
		_tecs_mode = ECL_TECS_MODE_UNDERSPEED;

	} else if (_uncommanded_descent_recovery) {
		_tecs_mode = ECL_TECS_MODE_BAD_DESCENT;

	} else if (_climbout_mode_active) {
		_tecs_mode = ECL_TECS_MODE_CLIMBOUT;

	} else {
		// This is the default operation mode
		_tecs_mode = ECL_TECS_MODE_NORMAL;
	}

}

void TECS::_initialize_pitchsp_offset() {

	if (!_pitchsp_offset_initialized) {
		// sanity check
		if (_wing_area < 0.01f || _cl_to_alpha_rad_slope < 0.01f) {
			_pitchsp_offset_initialized = false;
			return;
		}

		//first calculating the lift coefficient at cruise flight at trim airspeed.
		// cl = 2*L/(rho*A*V^2)
		_cl_cruise_trim_as = 2.0f * CONSTANTS_ONE_G * _auw / (CONSTANTS_AIR_DENSITY_SEA_LEVEL_15C * _wing_area * _indicated_airspeed_trim * _indicated_airspeed_trim);
		_cl_offset_clean_cruise_trim_as = _cl_cruise_trim_as - _pitchsp_offset_rad / _cl_to_alpha_rad_slope;

		//Setting the flaps should ideally only change the offset, not the slope angle.
		_cl_offset_flaps_cruise_trim_as = _cl_cruise_trim_as - _pitchsp_offset_flaps_rad / _cl_to_alpha_rad_slope;

		//the required cl for any airspeed can be calculated by dividing _cl_coefficient by airspeed squared.
		_cl_coefficient = _cl_cruise_trim_as * _indicated_airspeed_trim * _indicated_airspeed_trim;

		_pitchsp_offset_initialized = true;
	}
}
