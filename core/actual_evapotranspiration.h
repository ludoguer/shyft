///	Copyright 2012 Statkraft Energi A/S
///
///	This file is part of Shyft.
///
///	Shyft is free software: you can redistribute it and/or modify it under the terms of
/// the GNU Lesser General Public License as published by the Free Software Foundation,
/// either version 3 of the License, or (at your option) any later version.
///
///	Shyft is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
/// without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
/// PURPOSE.  See the	GNU Lesser General Public License for more details.
///
///	You should have received a copy of the GNU Lesser General Public License along with
/// Shyft, usually located under the Shyft root directory in two files named COPYING.txt
/// and COPYING_LESSER.txt.	If not, see <http://www.gnu.org/licenses/>.
///
/// Adapted from early enki method programmed by Kolbjørn Engeland and Sjur Kolberg
///
#pragma once

#include "utctime_utilities.h"
/*!\tfile
contains the actual evatransporation parameters and algorithm
*/
namespace shyft {
    namespace core {
		namespace actual_evapotranspiration {
			/**<  keeps the parameters (potential calibration/localization) for the AE */
			struct parameter {
				double ae_scale_factor = 1.5; ///< [mm] default value is 1.5
				parameter(double ae_scale_factor=1.5) : ae_scale_factor(ae_scale_factor) {}
			};

			/**<  keeps the formal response of the AE method */
			struct response {
				double ae = 0.0;
			};


			/** \brief actual_evapotranspiration calculates actual evapotranspiration
			 * based on supplied parameters.
			 *
			 * If part of the area you are calculating ae for is snow covered,
			 * make sure that you multiply p.t.evap*(1-sca).
			 *
			 *  As the water level in the 'ground-layer' approaches zero, the actual evatranspiration goes to zero (reasonable!)
			 *  If the water level in the 'ground-layer' is > ae_scale_factor the actual evatranspiration ~ potential (reasonable!)
			 *  An interesting article can be found here: http://onlinelibrary.wiley.com/doi/10.1029/2006JD008351/full
			 *
			 * \param water_level [mm]
			 * \param potential_evapotranspiration
			 * \param scale_factor typically 1.5
			 * \return calculated actual evapotranspiration
			 */
			inline double calculate_step(const double water_level,const double potential_evapotranspiration,const double scale_factor, const utctimespan dt) {
                const double dh=double(dt)/calendar::HOUR;
				const double scale = scale_factor*dh/3.0;
                return potential_evapotranspiration*(1.0-exp(-water_level/scale))/dh;
			}


		};
	};
};
