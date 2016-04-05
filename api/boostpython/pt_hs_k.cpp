#include "boostpython_pch.h"

#include "core/utctime_utilities.h"
#include "core/priestley_taylor.h"
#include "core/actual_evapotranspiration.h"
#include "core/precipitation_correction.h"
#include "core/hbv_snow.h"
#include "core/kirchner.h"
#include "core/pt_hs_k.h"
#include "api/api.h"
#include "api/pt_hs_k.h"
#include "core/pt_hs_k_cell_model.h"
#include "core/region_model.h"
#include "core/model_calibration.h"
#include "expose_statistics.h"
#include "expose.h"

static char const* version() {
   return "v1.0";
}

namespace expose {
    namespace pt_hs_k {
        using namespace boost::python;
        using namespace shyft::core;
        using namespace shyft::core::pt_hs_k;

        static void parameter_state_response() {
            class_<parameter,bases<>,std::shared_ptr<parameter>>("PTHSKParameter",
                              "Contains the parameters to the methods used in the PTHSK assembly\n"
                              "priestley_taylor,hbv_snow,actual_evapotranspiration,precipitation_correction,kirchner\n"
                )
                .def(init<const priestley_taylor::parameter&,const hbv_snow::parameter&,const actual_evapotranspiration::parameter&,const kirchner::parameter&,const precipitation_correction::parameter&>(args("pt","snow","ae","k","p_corr"),"create object with specified parameters"))
                .def(init<const parameter&>(args("p"),"clone a parameter"))
                .def_readwrite("pt",&parameter::pt,"priestley_taylor parameter")
				.def_readwrite("ae", &parameter::ae, "actual evapotranspiration parameter")
                .def_readwrite("snow",&parameter::snow,"hbv-snow parameter")
                .def_readwrite("kirchner",&parameter::kirchner,"kirchner parameter")
                .def_readwrite("p_corr",&parameter::p_corr,"precipitation correction parameter")
                .def("size",&parameter::size,"returns total number of calibration parameters")
                .def("set",&parameter::set,args("p"),"set parameters from vector/list of float, ordered as by get_name(i)")
                .def("get",&parameter::get,args("i"),"return the value of the i'th parameter, name given by .get_name(i)")
                .def("get_name",&parameter::get_name,args("i"),"returns the i'th parameter name, see also .get()/.set() and .size()")
                ;
            register_ptr_to_python<std::shared_ptr<parameter> >();

            typedef std::map<int,parameter> PTHSKParameterMap;
            class_<PTHSKParameterMap>("PTHSKParameterMap","dict (int,parameter)  where the int is 0-based catchment_id")
                .def(map_indexing_suite<PTHSKParameterMap>())
            ;

            class_<state>("PTHSKState")
                .def(init<hbv_snow::state,kirchner::state>(args("snow","k"),"initializes state with hbv-snow gs and kirchner k"))
                .def_readwrite("snow",&state::snow,"hbv-snow state")
                .def_readwrite("kirchner",&state::kirchner,"kirchner state")
                ;

            typedef std::vector<state> PTHSKStateVector;
            class_<PTHSKStateVector,bases<>,std::shared_ptr<PTHSKStateVector> >("PTHSKStateVector")
                .def(vector_indexing_suite<PTHSKStateVector>())
                ;
            class_<response>("PTHSKResponse","This struct contains the responses of the methods used in the PTHSK assembly")
                .def_readwrite("pt",&response::pt,"priestley_taylor response")
                .def_readwrite("snow",&response::snow,"hbc-snow response")
                .def_readwrite("ae",&response::ae,"actual evapotranspiration response")
                .def_readwrite("kirchner",&response::kirchner,"kirchner response")
                .def_readwrite("total_discharge",&response::total_discharge,"total stack response")
                ;
        }

        static void
        collectors() {
            typedef shyft::core::pt_hs_k::all_response_collector PTHSKAllCollector;
            class_<PTHSKAllCollector>("PTHSKAllCollector", "collect all cell response from a run")
                .def_readonly("destination_area",&PTHSKAllCollector::destination_area,"a copy of cell area [m2]")
                .def_readonly("avg_discharge",&PTHSKAllCollector::avg_discharge,"Kirchner Discharge given in [m³/s] for the timestep")
                //.def_readonly("snow_sca",&PTHSKAllCollector::snow_sca," hbv snow covered area fraction, sca.. 0..1 - at the end of timestep (state)")
                //.def_readonly("snow_swe",&PTHSKAllCollector::snow_swe,"hbv snow swe, [mm] over the cell sca.. area, - at the end of timestep")
                .def_readonly("snow_outflow",&PTHSKAllCollector::snow_outflow," hbv snow output [m³/s] for the timestep")
                .def_readonly("ae_output",&PTHSKAllCollector::ae_output,"actual evap mm/h")
                .def_readonly("pe_output",&PTHSKAllCollector::pe_output,"pot evap mm/h")
                .def_readonly("end_reponse",&PTHSKAllCollector::end_reponse,"end_response, at the end of collected")
            ;

            typedef shyft::core::pt_hs_k::discharge_collector PTHSKDischargeCollector;
            class_<PTHSKDischargeCollector>("PTHSKDischargeCollector", "collect all cell response from a run")
                .def_readonly("destination_area",&PTHSKDischargeCollector::destination_area,"a copy of cell area [m2]")
                .def_readonly("avg_discharge",&PTHSKDischargeCollector::avg_discharge,"Kirchner Discharge given in [m³/s] for the timestep")
                .def_readonly("snow_sca",&PTHSKDischargeCollector::snow_sca," hbv snow covered area fraction, sca.. 0..1 - at the end of timestep (state)")
                .def_readonly("snow_swe",&PTHSKDischargeCollector::snow_swe,"hbv snow swe, [mm] over the cell sca.. area, - at the end of timestep")
                .def_readonly("end_reponse",&PTHSKDischargeCollector::end_response,"end_response, at the end of collected")
                .def_readwrite("collect_snow",&PTHSKDischargeCollector::collect_snow,"controls collection of snow routine")
                ;
            typedef shyft::core::pt_hs_k::null_collector PTHSKNullCollector;
            class_<PTHSKNullCollector>("PTHSKNullCollector","collector that does not collect anything, useful during calibration to minimize memory&maximize speed")
                ;

            typedef shyft::core::pt_hs_k::state_collector PTHSKStateCollector;
            class_<PTHSKStateCollector>("PTHSKStateCollector","collects state, if collect_state flag is set to true")
                .def_readwrite("collect_state",&PTHSKStateCollector::collect_state,"if true, collect state, otherwise ignore (and the state of time-series are undefined/zero)")
                .def_readonly("kirchner_discharge",&PTHSKStateCollector::kirchner_discharge,"Kirchner state instant Discharge given in m^3/s")
                .def_readonly("snow_swe",&PTHSKStateCollector::snow_swe,"")
                .def_readonly("snow_sca",&PTHSKStateCollector::snow_sca,"")
            ;

        }

        static void
        cells() {
              typedef shyft::core::cell<parameter, environment_t, state, state_collector, all_response_collector> PTHSKCellAll;
              typedef shyft::core::cell<parameter, environment_t, state, null_collector, discharge_collector> PTHSKCellOpt;
              expose::cell<PTHSKCellAll>("PTHSKCellAll","tbd: PTHSKCellAll doc");
              expose::cell<PTHSKCellOpt>("PTHSKCellOpt","tbd: PTHSKCellOpt doc");
              expose::statistics::hbv_snow<PTHSKCellAll>("PTHSKCell");//it only gives meaning to expose the *All collect cell-type
              expose::statistics::actual_evapotranspiration<PTHSKCellAll>("PTHSKCell");
              expose::statistics::priestley_taylor<PTHSKCellAll>("PTHSKCell");
              expose::statistics::kirchner<PTHSKCellAll>("PTHSKCell");
        }

        static void
        models() {
            typedef shyft::core::region_model<pt_hs_k::cell_discharge_response_t> PTHSKOptModel;
            typedef shyft::core::region_model<pt_hs_k::cell_complete_response_t> PTHSKModel;
            expose::model<PTHSKModel>("PTHSKModel","PTHSK");
            expose::model<PTHSKOptModel>("PTHSKOptModel","PTHSK");
        }

        static void
        state_io() {
            expose::state_io<shyft::api::pt_hs_k_state_io,shyft::core::pt_hs_k::state>("PTHSKStateIo");
        }


        static void
        model_calibrator() {
            expose::model_calibrator<shyft::core::region_model<pt_hs_k::cell_discharge_response_t>>("PTHSKOptimizer");
        }
    }
}

BOOST_PYTHON_MODULE(_pt_hs_k)
{

    boost::python::scope().attr("__doc__")="SHyFT python api for the pt_hs_k model";
    boost::python::def("version", version);
    expose::pt_hs_k::state_io();
    expose::pt_hs_k::parameter_state_response();
    expose::pt_hs_k::cells();
    expose::pt_hs_k::models();
    expose::pt_hs_k::collectors();
    expose::pt_hs_k::model_calibrator();
}

