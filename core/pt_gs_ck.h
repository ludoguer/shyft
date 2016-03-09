#pragma once

#include "priestley_taylor.h"
#include "kirchner.h"
#include "gamma_snow.h"
#include "actual_evapotranspiration.h"
#include "precipitation_correction.h"
#include "pt_gs_k.h"

namespace shyft {
  namespace core {
    namespace pt_gs_ck {
        using namespace std;

        /** \brief the pt_gs_ck parameter is exactly the same as pt_gs_k parameter
        * except that the k-part is evaluated at catchment level
        */
        typedef parameter shyft::core::pt_gs_k::parameter;

        /** \brief Simple state struct for the PTGSCK method stack
         * This struct contains the states of the methods used in the PTGSK assembly.
         * \sa gamma_snow \sa kirchner \sa pt_gs_k \sa priestley_taylor
         */
        struct cell_state {
            state() {}
            state(gamma_snow::state gs) : gs(gs) {}
            gamma_snow::state gs;
        };
        struct catchment_state {
            catchment_state() {}
            catchment_state(kirchner::state kirchner):kirchner(kirchner){}
            kirchner::state kirchner;
        }


        /** \brief Simple response struct for the PTGSK method stack
         *
         * This struct contains the responses of the methods used in the PTGSK assembly.
         */
        struct cell_response {
            priestley_taylor::response pt;
            gamma_snow::response gs;
            actual_evapotranspiration::response ae;
        };
        struct catchment_response {
            kirchner::response kirchner;
            // Stack response
            double total_discharge;
        };

        /** \brief Calculation Model using assembly of PriestleyTaylor, GammaSnow and Catchment level Kirchner
         *
         * This model first uses PriestleyTaylor for calculating the potential
         * evapotranspiration based on time series data for temperature,
         * radiation and relative humidity. Then it uses the GammaSnow method
         * to calculate the snow/ice adjusted runoff using time series data for
         * precipitation and wind speed in addition to the time series used in
         * the PriestleyTaylor method. The PriestleyTaylor potential evaporation is
         * is used to calculate the actual evapotranspiration.
         * At catchment level, we calculate actual evapotranspiration taking into
         * account the water-content in kirchner state.
         *
         * Kirchner is run with the gamma snow output, as well as the calculated
         * catchment-level actual evaporation response from the two methods above to
         * calculate the discharge.
         *
         *
         * \tparam TS Time series type that implements:
         *    - TS::source_accessor_type --> Type of accessor used to retrieve time series data.
         *    - TS.accessor(const T& time_axis) const --> TS::source_accessor_type. Accessor object
         *      for this time series.
         *
         * \tparam TS::source_accessor_type Time series source accessor type that implements:
         *    - TS::source_accessor_type(const TS& time_series, const T& time_axis) --> Construct
         *      accessor for the given time series and time axis.
         *    - TS::source_accessor_type.value(size_t i) --> double, -value of the time series for
         *      period i in the time axis.
         * \tparam T Time axis type that implements:
         *    - T.size() const --> Number of periods in the time axis.
         *    - T(size_t i) const --> shyft::core::utcperiod, - Start and end as shyft::core::utctime
         *      of the i'th period.
         * \tparam S State type that implements:
         *    - S::gs_state_type --> State type for the GammaSnow method.
         *    - S::kirchner_state_type --> State type for the Kirchner method.
         *    - S.gs --> S::gs_state_type, - State variables for the GammaSnow method
         *    - S.kirchner --> S::kirchner_state_type, - State variables for the Kirchner method
         * \tparam R Response type that implements:
         *    - R::gs_response_type --> Response type for the GammaSnow routine.
         *    - R.gs --> R::gs_response_type, -Response object passed to the GammaSnow routine.
         * \tparam P Parameter type that implements:
         *    - P::pt_parameter_type --> Parameter type for the PriestleyTaylor method
         *    - P::gs_parameter_type --> Parameter type for the GammaSnow method.
         *    - P::ae_parameter_type --> Parameter type for the ActualEvapotranspiration method.
         *    - P::kirchner_parameter_type --> Parameter type for the Kirchner method.
         *    - P.pt --> P::pt_parameter_type --> Parameters for the PriestleyTaylor method.
         *    - P.gs --> P::gs_parameter_type --> Parameters for thge GammaSnow method.
         *    - P.ae --> P::ae_parameter_type --> Parameters for thge ActualEvapotranspiration method.
         *    - P.kirchner --> P::kirchner_parameter_type --> Parameters for the Kirchner method.
         * \tparam SC State collector type that implements:
         *    - SC.collect(utctime t, const S& state) --> Possibly save some states at time t.
         * \tparam RC Response collector type that implements:
         *    - RC.collect(utctime t, const R& response) --> Possibly save some responses at time t.
         */
#ifndef SWIG
        template<template <typename, typename> class A, class R, class T_TS, class P_TS, class WS_TS, class RH_TS, class RAD_TS, class T,
        class S, class GCD, class P, class SC, class RC>
        void run_pt_gs_ck(const GCD& geo_cell_data,
            const P& parameter,
            const T& time_axis,
            const T_TS& temp,
            const P_TS& prec,
            const WS_TS& wind_speed,
            const RH_TS& rel_hum,
            const RAD_TS& rad,
            S& state,
            SC& state_collector,
            RC& response_collector
            ) {
            // Access time series input data through accessors of template A (typically a direct accessor).
            using temp_accessor_t = A<T_TS, T>;
            using prec_accessor_t = A<P_TS, T>;
            using wind_speed_accessor_t = A<WS_TS, T>;
            using rel_hum_accessor_t = A<RH_TS, T>;
            using rad_accessor_t = A<RAD_TS, T>;

            auto temp_accessor = temp_accessor_t(temp, time_axis);
            auto prec_accessor = prec_accessor_t(prec, time_axis);
            auto wind_speed_accessor = wind_speed_accessor_t(wind_speed, time_axis);
            auto rel_hum_accessor = rel_hum_accessor_t(rel_hum, time_axis);
            auto rad_accessor = rad_accessor_t(rad, time_axis);

            // Initialize the method stack
            precipitation_correction::calculator p_corr(parameter.p_corr.scale_factor);
            priestley_taylor::calculator pt(parameter.pt.albedo, parameter.pt.alpha);
            gamma_snow::calculator<typename P::gs_parameter_t, typename S::gs_state_t, typename R::gs_response_t> gs;
            kirchner::calculator<kirchner::trapezoidal_average, typename P::kirchner_parameter_t> kirchner(parameter.kirchner);
            //
            gs.set_glacier_fraction(geo_cell_data.land_type_fractions_info().glacier());
            // Get the initial states
            auto &gs_state = state.gs;
            double q = state.kirchner.q;
            R response;
            const double forest_fraction=geo_cell_data.land_type_fractions_info().forest();
            const double altitude= geo_cell_data.mid_point().z;
            // Step through times in axis
            for (size_t i = 0; i < time_axis.size(); ++i) {
                utcperiod period = time_axis(i);
                double temp = temp_accessor.value(i);
                double rad = rad_accessor.value(i);
                double rel_hum = rel_hum_accessor.value(i);
                double prec = p_corr.calc(prec_accessor.value(i));
                state_collector.collect(i, state);///< \note collect the state at the beginning of each period (the end state is saved anyway)
                //
                // Land response:
                //

                // PriestleyTaylor (scale by timespan since it delivers results in mm/s)
                double pot_evap = pt.potential_evapotranspiration(temp, rad, rel_hum)*period.timespan();
                response.pt.pot_evapotranspiration=pot_evap;

                // Gamma Snow
                gs.step(gs_state, response.gs, period.start, period.timespan(), parameter.gs,
                        temp, rad, prec, wind_speed_accessor.value(i), rel_hum,forest_fraction,altitude);

                // Actual Evapotranspiration
                double act_evap = actual_evapotranspiration::calculate_step(q,
                                                    (1-response.gs.sca)*pot_evap,
                                                    parameter.ae.ae_scale_factor,
                                                    period.timespan());
                response.ae.ae = act_evap;

                response_collector.collect(i, response);///< \note collect the response valid for the i'th period (current state is now at the end of period)
            }
            response_collector.set_end_response(response);
        }

        /** \brief catchment routing provides catchment level
         * routing/handling of data
         * SiH Comments:
         * maybe easier to take it the opposite way around.
         *  look at the catchment as the 'super-cell',
         *   and delegate the details of calculations to the
         *   super-cell.
         *  why ?:
         *   - easy to do cell-level parallell first
         *     then do catchment level after that.
         *   - easy to pass data from cell to catchment
         *   - we can still have a routing model.
         *  who not ?:
         *   - if few cells, we do not gain full effect of parallel core (e.g. one catchment have 10 cells, the other 90).
         *  conclusion:
         *   seems like the best choice to have integrated pt_gs_ck, or top-layer cell, then catchment kirchner
         */
        struct catchment_routing {
            size_t cid;///< All cells that have this catchment id belongs to this catchment_routing
            catchment_routing(size_t cid=0):cid(cid) {}

            #ifndef SWIG
            /** \brief calculate avg discharge of cells, currently no routing
            */
            template <class C, class P, class S>
            void calculate( const vector<C>& cells, const P& parameter,S& state) {
                if(cells.size()==0)
                    return;
                const auto& time_axis= cells[0].rc.avg_discharge.time_axis;
                pts_t precip(time_axis,0);
                pts_t evap(time_axis,0);
                pts_t direct_response(time_axis,0);

                // (1) collect precipitation and actual evapotranspiration from the cells, also compute and land areatypes from all cells..
                double kirchner_area=0.0;
                double area=0.0;
                for(const C& c: cells) { // could be parallell, slice in time or cells
                    if(c.geo.catchment_id()==cid) {
                        precip.add(c.rc.precipitation);
                        evap.add(c.rc.actual_evapotranspiration);
                        direct_response(c.rc.direct_response);
                        area += c.geo.area(); // needed,but could be pre-computed and stashed away on catchment-structure
                        kirchner_area += c.area()*(c.geo.land_type_fractions().forest() + c.geo.land_type_fractions().unspecified());
                    }
                }

                // (2) feed total output from cells into kirchner and compute final response
                kirchner::calculator<kirchner::trapezoidal_average, typename P::kirchner_parameter_t> kirchner(parameter.kirchner);
                avg_discharge=pts_t(time_axis,0.0);
                for(size_t i=0;i<time_axis.size();++i) {
                    utcperiod period = time_axis(i);
                    double q_avg;
                    kirchner.step(period.start, period.end, state.kirchner.q, q_avg, precip.value(i), evap.value(i));
                    avg_discharge.set(i,direct_response.value(i) + q_avg);
                }
            }
            #endif
            pts_t avg_discharge;
        };

#endif
    } // pt_gs_ck
  } // core
} // shyft
