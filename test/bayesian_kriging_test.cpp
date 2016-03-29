#include "test_pch.h"
#include "bayesian_kriging_test.h"
#include "mocks.h"
#include "core/inverse_distance.h"
#include "core/bayesian_kriging.h"
#include "core/timeseries.h"
#include "core/geo_point.h"

#include <armadillo>
#include <ctime>
#include <random>

using namespace std;
namespace shyfttest {
	const double EPS = 1.0e-8;
	using namespace shyft::core;
	using namespace shyft::timeseries;
	using shyft::core::geo_point;

	namespace btk_structs {

		//! Simple source class for testing Bayesian Temperature Kriging.
		class Source {
		public:
			typedef xpts_t ts_source;
			typedef shyft::timeseries::average_accessor<ts_source, shyfttest::point_timeaxis> source_accessor;
		private:
			geo_point coord;
			const ts_source temperature_ts;
		public:
			Source(const geo_point& coord, const ts_source& temperatures)
				: coord(coord), temperature_ts(std::move(temperatures)) {
				// Do nothing
			}
			const geo_point& mid_point() const { return coord; }
			const ts_source& temperatures() const { return temperature_ts; }
			source_accessor temperature_accessor(const shyfttest::point_timeaxis& time_axis) const { return source_accessor(temperature_ts, time_axis); }
		};

		//! Simple destination class for use in the IWD algorithms to represent a cell that can have temperature,
		//! radiation and precipitation set.
		class MCell {
		private:
			geo_point coord;
		public:
			arma::vec temperatures;
			template<typename T> explicit MCell(const geo_point& coord, const T& time_axis)
				: coord(coord),
				temperatures((arma::uword)time_axis.size()) {
				temperatures.fill(std::nan(""));
			}
			const geo_point& mid_point() const { return coord; }
			void set_temperature(const size_t time_idx, const double temperature) {
				temperatures[(arma::uword)time_idx] = temperature;
			}
			const double temperature(const size_t time_idx) const { return temperatures[(arma::uword)time_idx]; }

		};


		//! Simple parameter class
		class Parameter {
		private:
			double gradient = -0.006; // Prior expectation of temperature gradient [C/m]
			double gradient_sd = 0.0025; // Prior standard deviation of temperature gradient in [C/m]
			double sill_value = 25.0; // Value of semivariogram at range
			double nug_value = 0.5; // Nugget magnitude
			double range_value = 200000.0; // Point where semivariogram flattens out
			double  zscale_value = 20.0; // Height scale used during distance computations.
		public:
			Parameter() {}
			Parameter(double temperature_gradient, double temperature_gradient_sd)
				: gradient(temperature_gradient / 100),
				gradient_sd(temperature_gradient_sd / 100) {}

			const double temperature_gradient(shyft::core::utcperiod period) const {
				return gradient;
			}

			const double& temperature_gradient_sd() const {
				return gradient_sd;
			}

			double sill() const { return sill_value; }
			double nug() const { return nug_value; }
			double range() const { return range_value; }
			double zscale() const { return zscale_value; }
		};


		//! TODO: Accessor class
		class AccessorMocker {};


		// Make some aliases for the tests below
		typedef std::vector<Source> SourceList;
		typedef std::vector<MCell> DestinationList;
	}; // Namespace btk_structs
}; // Namespace test

using std::begin;
using std::end;
using namespace shyfttest::btk_structs;
using shyft::timeseries::point_timeaxis;
using namespace shyft::core;

typedef shyfttest::xpts_t xpts_t;
std::vector<utctime> ctimes{ 0, 3600 };
typedef shyft::timeseries::point_timeaxis point_timeaxis;
typedef std::vector<shyft::timeseries::point> point_vector_t;

void build_sources_and_dests(const size_t num_sources_x, const size_t num_sources_y,
	const size_t num_dests_x, const size_t num_dests_y,
	const size_t ts_size, const shyft::timeseries::utctimespan dt,
	const point_timeaxis& time_axis, bool insert_nans, SourceList& sources, DestinationList& dests) {
	const double x_min = 0.0; // [m]
	const double x_max = 100000.0; // [m]
	const double y_min = 0.0; // [m]
	const double y_max = 1000000.0; // [m]

	sources.reserve(num_sources_x*num_sources_y);
	dests.reserve(num_dests_x*num_dests_y);
	geo_point pt;
	double lower_bound = 0.0;
	double upper_bound = 10.0;
	std::uniform_real_distribution<double> unif(lower_bound, upper_bound);
	std::default_random_engine re;
    vector<utctime> times;times.reserve(ts_size);
    for (size_t l = 0; l < ts_size; ++l)
        times.emplace_back(l*dt);
    times.emplace_back(shyft::core::max_utctime);
    point_timeaxis dta(times);
	for (size_t i = 0; i < num_sources_x; ++i) {
		pt.x = x_min + i*(x_max - x_min) / (num_sources_x - 1);
		for (size_t j = 0; j < num_sources_y; ++j) {
			pt.y = y_min + j*(y_max - y_min) / (num_sources_y - 1);
			pt.z = 500 * std::sin(pt.x / x_max) + std::sin(pt.y / y_max) / 2;
			vector<double> pts; pts.reserve(ts_size);
			double b_t = unif(re);
			//std::cout << "Base temp at pos (i,j) = " << i << ", " << j << ") = " << b_t << std::endl;
			for (size_t l = 0; l < ts_size; ++l)
				pts.emplace_back( b_t + pt.z*(0.6 / 100));
			sources.emplace_back(pt, xpts_t(dta,pts));
		}
	}
	for (size_t i = 0; i < num_dests_x; ++i) {
		pt.x = x_min + i*(x_max - x_min) / (num_dests_x - 1);
		for (size_t j = 0; j < num_dests_y; ++j) {
			pt.y = y_min + j*(y_max - y_min) / (num_dests_y - 1);
			pt.z = 500 * (std::sin(pt.x / x_max) + std::sin(pt.y / y_max)) / 2;
			dests.emplace_back(pt, time_axis);
		}
	}
}
using namespace shyft::core::bayesian_kriging;

void bayesian_kriging_test::test_covariance_calculation() {
	Parameter params;

	const arma::uword n = 100;
	arma::mat C(n, n, arma::fill::eye);
	arma::vec dqs(n*(n - 1) / 2); // Strict upper part of distance matrix, as vector
	arma::vec covs(n*(n - 1) / 2); // Strict upper part of cov matrix, as vector
	for (arma::uword i = 0; i < n*(n - 1) / 2; ++i) dqs.at(i) = static_cast<double>((i + 1)*(i + 1)); // NB: Not valid dists
	// Start timer
	//const std::clock_t start = std::clock();
	utils::cov(dqs, covs, params);
	arma::uword v_idx = 0;
	for (arma::uword i = 0; i < n; ++i)
		for (arma::uword j = i + 1; j < n; ++j)
			C.at(i, j) = covs.at(v_idx++);
	C.diag() *= utils::zero_dist_cov(params);
	//const std::clock_t total = std::clock() - start;
	//std::cout << "Computing upper covariance matrix with nnz " << n*(n + 1) / 2 << " took: " << 1000 * (total) / (double)(CLOCKS_PER_SEC) << " ms" << std::endl;

	// Compute source-target cov
	const arma::uword m = 200 * 200;
	arma::vec dqms(n*m), covs2(n*m);
	v_idx = 0;
	for (arma::uword i = 0; i < n; ++i)
		for (arma::uword j = 0; j < m; ++j)
			dqms.at(i*m + j) = static_cast<double>((i + 1)*(j + 1));
	const std::clock_t start2 = std::clock();
	utils::cov(dqms, covs2, params);
	arma::mat CC(covs2);
	CC.reshape(n, m);
	const std::clock_t total2 = std::clock() - start2;
	std::cout << "Computing full covariance matrix with nnz " << n*m << " took: " << 1000 * (total2) / (double)(CLOCKS_PER_SEC) << " ms" << std::endl;
}

void bayesian_kriging_test::test_build_covariance_matrices() {
	return;
	Parameter params;
	const point_timeaxis time_axis(ctimes);
	SourceList sources;
	DestinationList destinations;
	build_sources_and_dests(3, 3, 15, 15, 2, 10, time_axis, false, sources, destinations);
	arma::mat K, k;
	utils::build_covariance_matrices(begin(sources), end(sources), begin(destinations), end(destinations), params, K, k);
	TS_ASSERT_EQUALS(K.n_rows, (size_t)9);
	TS_ASSERT_EQUALS(K.n_cols, (size_t)9);
	TS_ASSERT_EQUALS(k.n_rows, (size_t)9);
	TS_ASSERT_EQUALS(k.n_cols, (size_t)15 * 15);
}

void bayesian_kriging_test::test_build_elevation_matrices() {
	return;
	Parameter params;
	const point_timeaxis time_axis(ctimes);
	SourceList sources;
	DestinationList destinations;
	build_sources_and_dests(3, 3, 15, 15, 2, 10, time_axis, false, sources, destinations);
	arma::mat S_e, D_e;
	utils::build_elevation_matrices(begin(sources), end(sources), begin(destinations), end(destinations), S_e, D_e);
	TS_ASSERT_EQUALS(S_e.n_cols, (size_t)2);
	TS_ASSERT_EQUALS(S_e.n_rows, (size_t)9);
	TS_ASSERT_EQUALS(D_e.n_cols, (size_t)15 * 15);
	TS_ASSERT_EQUALS(D_e.n_rows, (size_t)2);
}

void bayesian_kriging_test::test_interpolation() {
	Parameter params;
	SourceList sources;
	DestinationList destinations;
	using namespace shyft::timeseries;
	using namespace shyft::core;
	using namespace shyfttest;
	size_t n_s = 2;
	size_t n_d = 4;
	size_t n_times = 2; // Daily interpolations for three years
	shyft::timeseries::utctime dt = 10;
	vector<utctime> times; times.reserve(n_times);
	for (size_t i = 0; i < n_times; ++i)
		times.emplace_back(dt*i);
	const point_timeaxis time_axis(times);
	build_sources_and_dests(n_s, n_s, n_d, n_d, n_times, dt, time_axis, true, sources, destinations);
	const std::clock_t start = std::clock();
	btk_interpolation<average_accessor<shyfttest::xpts_t, point_timeaxis>>(begin(sources), end(sources), begin(destinations), end(destinations), time_axis, params);
	const std::clock_t total = std::clock() - start;
	std::cout << "Calling compute with n_sources, n_dests, and n_times = " << n_s*n_s << ", " << n_d*n_d << ", " << n_times << " took: " << 1000 * (total) / (double)(CLOCKS_PER_SEC) << " ms" << std::endl;
	MCell d = destinations[destinations.size() - 1];
	std::cout << "Temp at altitude " << d.mid_point().z << " is " << d.temperature(0) << std::endl;
	d = destinations[0];
	if (getenv("SHYFT_BTK_VERBOSE")) {
		std::cout << "Temp at altitude " << d.mid_point().z << " is " << d.temperature(0) << std::endl;
		for (auto d : destinations) {
			std::cout << d.mid_point().z << std::endl;
			std::cout << "Max/Min: " << *std::max_element(d.temperatures.begin(), d.temperatures.end()) <<
				", " << *std::min_element(d.temperatures.begin(), d.temperatures.end()) << std::endl;
		}
	}
}
