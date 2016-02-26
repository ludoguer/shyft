#include "test_pch.h"
#include "core/utctime_utilities.h"
#include "utctime_utilities_test.h"
using namespace std;
using namespace shyft;
using namespace shyft::core;

void utctime_utilities_test::test_utctime() {
    calendar c(0);
    YMDhms unixEra(1970,01,01,00,00,00);
    YMDhms y_null;

    TS_ASSERT(y_null.is_null());
    TS_ASSERT(c.time(y_null)==no_utctime);
    TS_ASSERT(c.time(YMDhms::max())==max_utctime);
    TS_ASSERT(c.time(YMDhms::min())==min_utctime);

    TS_ASSERT_EQUALS(0L,c.time(unixEra));
    YMDhms r=c.calendar_units(utctime(0L));
    TS_ASSERT_EQUALS(r,unixEra);
}
void utctime_utilities_test::test_utcperiod() {
    calendar utc;
    utctime t0=utc.time(2015,1,1);
    utctime t1=utc.time(2015,2,1);
    utctime t2=t1-deltahours(1);
    utctime t3=t1+deltahours(1);
    TS_ASSERT(utcperiod(t0,t1).valid());
    TS_ASSERT(is_valid(utcperiod())==false);
    TS_ASSERT(utcperiod(t1,t0).valid()==false);
    TS_ASSERT_EQUALS(utcperiod(t1,t2),utcperiod(t1,t2));
    TS_ASSERT(utcperiod(t1,t2)!=utcperiod(t0,t1));
    TS_ASSERT(utcperiod(t0,t1).contains(t0));
    TS_ASSERT(utcperiod(t0,t1).contains(t1)==false);
    TS_ASSERT(utcperiod(t0,t1).contains(no_utctime)==false);
    TS_ASSERT(utcperiod().contains(t0)==false);
    TS_ASSERT(utcperiod(t0,t1).overlaps(utcperiod(t1,t3))==false);
    TS_ASSERT(utcperiod(t1,t3).overlaps(utcperiod(t0,t1))==false);
    TS_ASSERT(utcperiod(t0,t1).overlaps(utcperiod(t2,t3))==true);
    TS_ASSERT(utcperiod(t2,t3).overlaps(utcperiod(t0,t1))==true);



}
void utctime_utilities_test::test_calendar_trim() {
    // simple trim test
    calendar cet(deltahours(1));
    utctime t=cet.time(YMDhms(2012,3,8,12,16,44));
    YMDhms t_y  (2012, 1, 1, 0, 0, 0);
    YMDhms t_m  (2012, 3, 1, 0, 0, 0);
    YMDhms t_w  (2012, 3, 5, 0, 0, 0);
    YMDhms t_d  (2012, 3, 8, 0, 0, 0);
    YMDhms t_h  (2012, 3, 8,12, 0, 0);
    YMDhms t_15m(2012, 3, 8,12,15, 0);
    TS_ASSERT_EQUALS(cet.calendar_units( cet.trim(t,calendar::YEAR)),t_y);
    TS_ASSERT_EQUALS(cet.calendar_units( cet.trim(t,calendar::MONTH)),t_m);
    TS_ASSERT_EQUALS(cet.calendar_units( cet.trim(t,calendar::WEEK)), t_w);
    TS_ASSERT_EQUALS(cet.calendar_units( cet.trim(t,calendar::DAY)),t_d);
    TS_ASSERT_EQUALS(cet.calendar_units( cet.trim(t,calendar::HOUR)),t_h);
    TS_ASSERT_EQUALS(cet.calendar_units( cet.trim(t,deltaminutes(15))),t_15m);
}


void utctime_utilities_test::test_calendar_timezone() {
    YMDhms unixEra(1970,01,01,00,00,00);
    calendar cet(deltahours(1));
    TS_ASSERT_EQUALS(deltahours(-1),cet.time(unixEra));
}

void utctime_utilities_test::test_calendar_to_string() {
    calendar utc;
    calendar osl(deltahours(1));
    calendar cet("Europe/Oslo");
    calendar xxx("America/St_Johns");// -03.30
    utctime t=utc.time(YMDhms(2012,5,8,12,16,44));
    TS_ASSERT_EQUALS(utc.to_string(t),string("2012-05-08T12:16:44Z"));
    TS_ASSERT_EQUALS(cet.to_string(t),string("2012-05-08T14:16:44+02"));
    TS_ASSERT_EQUALS(osl.to_string(t),string("2012-05-08T13:16:44+01"));
    TS_ASSERT_EQUALS(xxx.to_string(t),string("2012-05-08T09:46:44-02:30"));
}

void utctime_utilities_test::test_calendar_day_of_year() {
    calendar cet(deltahours(1));
    TS_ASSERT_EQUALS(1,cet.day_of_year(cet.time(YMDhms(2012,1,1,10,11,12))));
    TS_ASSERT_EQUALS(2,cet.day_of_year(cet.time(YMDhms(2012,1,2,0,0,0))));
    TS_ASSERT_EQUALS(366,cet.day_of_year(cet.time(YMDhms(2012,12,31,12,0,0))));


}
void utctime_utilities_test::test_calendar_month() {
    calendar cet(deltahours(1));
    TS_ASSERT_EQUALS( 1,cet.month(cet.time(YMDhms(2012,1,1,10,11,12))));
    TS_ASSERT_EQUALS( 2,cet.month(cet.time(YMDhms(2012,2,2,0,0,0))));
    TS_ASSERT_EQUALS(12,cet.month(cet.time(YMDhms(2012,12,31,23,59,59))));

}

void utctime_utilities_test::test_YMDhms_reasonable_calendar_coordinates() {
	TS_ASSERT_THROWS_ANYTHING(YMDhms(10000, 1, 1, 0, 0, 0));
    TS_ASSERT_THROWS_ANYTHING(YMDhms(-10000,1,1,0,0,0));
    TS_ASSERT_THROWS_ANYTHING(YMDhms(2000,0,1,0,0,0));
    TS_ASSERT_THROWS_ANYTHING(YMDhms(2000,13,1,0,0,0));
    TS_ASSERT_THROWS_ANYTHING(YMDhms(2000,1,0,0,0,0));
    TS_ASSERT_THROWS_ANYTHING(YMDhms(2000,1,32,0,0,0));
    TS_ASSERT_THROWS_ANYTHING(YMDhms(2000,1,1,-1,0,0));
    TS_ASSERT_THROWS_ANYTHING(YMDhms(2000,1,1,24,0,0));
    TS_ASSERT_THROWS_ANYTHING(YMDhms(2000,1,1,0,-1,0));
    TS_ASSERT_THROWS_ANYTHING(YMDhms(2000,1,1,0,60,0));
    TS_ASSERT_THROWS_ANYTHING(YMDhms(2000,1,1,0,0,-1));
    TS_ASSERT_THROWS_ANYTHING(YMDhms(2000,1,1,0,0,60));
}
void utctime_utilities_test::test_calendar_add_and_diff_units() {
	calendar cet(deltahours(1));
	int n_units = 3;
	utctimespan dts[] = { calendar::HOUR, calendar::DAY, calendar::WEEK, calendar::MONTH,calendar::YEAR };
	auto t1 = cet.time(YMDhms(2012, 1, 1));
	for (auto dt : dts) {
		auto t2 = cet.add(t1, dt, n_units);
		auto t3 = cet.add(t2, dt, -n_units);
		utctimespan rem;
		auto ndiff = cet.diff_units(t1, t2, dt, rem);//verify we get difference ok
		TS_ASSERT_EQUALS(ndiff, n_units);
		TS_ASSERT_EQUALS(t3, t1);// verify subtraction gives back original result
	}
}
void utctime_utilities_test::test_calendar_day_of_week() {
	calendar cet(deltahours(1));
	for (int i = 0; i < 7;++i)
		TS_ASSERT_EQUALS(i, cet.day_of_week(cet.time(YMDhms(2015, 8, 9+i, 10, 11, 12))));

}

void utctime_utilities_test::test_tz_info_db() {
    using namespace shyft::core;
    using namespace shyft::core::time_zone;
    tz_info_database tz_info_db;
    TS_ASSERT_THROWS_ANYTHING(tz_info_db.tz_info_from_region("Europe/Oslo"));
    // create from iso zone description:"Europe/Oslo","CET","CET","CEST","CEST","+01:00:00","+01:00:00","-1;0;3","+02:00:00","-1;0;10","+03:00:00"
    tz_info_db.add_tz_info("Europe/Oslo","CET+01CEST+01,M3.5.0/02:00,M10.5.0/03:00");
    auto eu_osl=tz_info_db.tz_info_from_region("Europe/Oslo");
    TS_ASSERT_EQUALS(eu_osl->name(),string("Europe/Oslo"));
    calendar utc;
    TS_ASSERT_EQUALS(eu_osl->is_dst(utc.time(YMDhms(2016, 3,27, 0,59,59))),false);// second before..
    TS_ASSERT_EQUALS(eu_osl->is_dst(utc.time(YMDhms(2016, 3,27, 1, 0, 0))),true); // exactly at shift into summer
    TS_ASSERT_EQUALS(eu_osl->is_dst(utc.time(YMDhms(2016,10,30, 0,59,59))),true); // second before..
    TS_ASSERT_EQUALS(eu_osl->is_dst(utc.time(YMDhms(2016,10,30, 1, 0, 0))),false);// exactly at shift into winter
    tz_info_database tz_iso_db;// do some minor testing of the internal iso tz db.
    tz_iso_db.load_from_iso_db();
    auto eu_osl_iso=tz_iso_db.tz_info_from_region("Europe/Oslo");
    TS_ASSERT_EQUALS(eu_osl_iso->name(),string("Europe/Oslo"));
    TS_ASSERT_EQUALS(eu_osl_iso->is_dst(utc.time(YMDhms(2016, 3,27, 0,59,59))),false);// second before..
    TS_ASSERT_EQUALS(eu_osl_iso->is_dst(utc.time(YMDhms(2016, 3,27, 1, 0, 0))),true); // exactly at shift into summer
    TS_ASSERT_EQUALS(eu_osl_iso->is_dst(utc.time(YMDhms(2016,10,30, 0,59,59))),true); // second before..
    TS_ASSERT_EQUALS(eu_osl_iso->is_dst(utc.time(YMDhms(2016,10,30, 1, 0, 0))),false);// exactly at shift into winter

}
void utctime_utilities_test::test_add_over_dst_transitions() {
    using namespace shyft::core;
    using namespace shyft::core::time_zone;
    tz_info_database tz_info_db;
    tz_info_db.load_from_iso_db();
    /// case 1: 23hour day, winter to summer shift.
    auto cet_info=tz_info_db.tz_info_from_region("Europe/Oslo");
    calendar cet(cet_info);
    calendar utc;
    utctime t0=cet.time(YMDhms(2016, 3,27,0,0,0));
    utctime t1=cet.add(t0,calendar::DAY,1);
    utctime t2=cet.add(t1,calendar::DAY,-1);
    TS_ASSERT_EQUALS(t1,cet.time(YMDhms(2016,3,28)));
    utctimespan rem(0);
    auto n_units=cet.diff_units(t0,t1,calendar::DAY,rem);
    TS_ASSERT_EQUALS(1,n_units);//its a local day diffeerence!
    TS_ASSERT_EQUALS(rem,utctimespan(0));// no remainder
    TS_ASSERT_EQUALS(deltahours(23),(t1-t0));// the miracle, it's 23 hours.
    TS_ASSERT_EQUALS(t0,t2);
    // case 1: into the details
    utctime t1_1=cet.add(t0,deltahours(1),1);
    utctime t1_2=cet.add(t0,deltahours(2),1);
    utctime t1_3=cet.add(t0,deltahours(3),1);
    TS_ASSERT_EQUALS(t1_1,t0+deltahours(1));
    TS_ASSERT_EQUALS(t1_2,t0+deltahours(1));
    TS_ASSERT_EQUALS(YMDhms(2016,3,27,1),cet.calendar_units(t1_2));
    TS_ASSERT_EQUALS(t1_3,t0+deltahours(2));

    /// case 2: 25 hour, summer->winter
    t0 = cet.time(YMDhms(2016,10,30));
    t1 = cet.add(t0,calendar::DAY,1);
    t2 = cet.add(t1,calendar::DAY,-1);
    n_units = cet.diff_units(t0,t1,calendar::DAY,rem);
    TS_ASSERT_EQUALS(1,n_units);//its a local day diffeerence!
    TS_ASSERT_EQUALS(rem,utctimespan(0));// no remainder
    TS_ASSERT_EQUALS(deltahours(25),(t1-t0));// the miracle, it's 25 hours.
    TS_ASSERT_EQUALS(t0,t2);
    /// case 2: into the details
    t1_1=cet.add(t0,deltahours(1),1);
    t1_2=cet.add(t0,deltahours(2),1);
    t1_3=cet.add(t0,deltahours(3),1);
    utctime
    t1_4=cet.add(t0,deltahours(4),1);
    TS_ASSERT_EQUALS(t1_1,t0+deltahours(1));TS_ASSERT_EQUALS(YMDhms(2016,10,30,1),cet.calendar_units(t1_1));
    TS_ASSERT_EQUALS(t1_2,t0+deltahours(2));TS_ASSERT_EQUALS(YMDhms(2016,10,30,2),cet.calendar_units(t1_2));
    TS_ASSERT_EQUALS(t1_3,t0+deltahours(4));TS_ASSERT_EQUALS(YMDhms(2016,10,30,3),cet.calendar_units(t1_3));
    TS_ASSERT_EQUALS(t1_4,t0+deltahours(5));
    auto xx=cet.calendar_units(t1_4);
    TS_ASSERT_EQUALS(YMDhms(2016,10,30,4),xx);


}

void utctime_utilities_test::test_calendar_trim_with_dst() {
    using namespace shyft::core;
    using namespace shyft::core::time_zone;
    tz_info_database tz_info_db;
    tz_info_db.load_from_iso_db();
    /// case 1: 23hour day, winter to summer shift.
    auto cet_info=tz_info_db.tz_info_from_region("Europe/Oslo");
    calendar cet(cet_info);// simple trim test
    utctime t=cet.time(YMDhms(2016,3,27, 0,0,0));// day of dst winter->summer
    for(int h=0;h<23;++h) {
        TS_ASSERT_EQUALS(cet.to_string(t),cet.to_string(cet.trim(t+deltahours(h),calendar::DAY)));// 23 hours should go down here.
        TS_ASSERT_EQUALS(cet.time(YMDhms(2016,3,21)),cet.trim(t+deltahours(h),calendar::WEEK));// should all end on monday 21st.
    }
    TS_ASSERT_EQUALS(cet.time(YMDhms(2016,3,28)),cet.trim(t+deltahours(23),calendar::DAY));
    TS_ASSERT_EQUALS(cet.time(YMDhms(2016,3,28)),cet.trim(t+deltahours(23),calendar::WEEK));
    /// case 2: 25 hour day, summer to winter
    t=cet.time(YMDhms(2016,10,30, 0,0,0));// day of dst summer->winter
    for(int h=0;h<25;++h) {
        TS_ASSERT_EQUALS(cet.to_string(t),cet.to_string(cet.trim(t+deltahours(h),calendar::DAY)));// 25 hours should go down here.
        TS_ASSERT_EQUALS(cet.time(YMDhms(2016,10,24)),cet.trim(t+deltahours(h),calendar::WEEK));// should all end on monday 21st.
    }
    TS_ASSERT_EQUALS(cet.time(YMDhms(2016,10,31)),cet.trim(t+deltahours(25),calendar::DAY));
    TS_ASSERT_EQUALS(cet.time(YMDhms(2016,10,31)),cet.trim(t+deltahours(25),calendar::WEEK));


}
void utctime_utilities_test::test_add_months() {
    using namespace shyft::core;
    using namespace shyft::core::time_zone;
    tz_info_database tz_info_db;
    tz_info_db.load_from_iso_db();
    for(auto& rtz:tz_info_db.region_tz_map) {
        auto cet_info=tz_info_db.tz_info_from_region(rtz.first);
        //cout<<"testing :"<<rtz.first<<endl;
        calendar cet(cet_info);
        calendar cal;
        /// case 1: trival just add increasing number of months, and insist on other time-parameters are constant.
        for(int month=1;month<13;++month) {
            for(int day=1;day<29;day+=3) {
                YMDhms dt0(2015,month,day,15,45,35);
                utctime t0=cet.time(dt0);
                for(int m=0;m<12*5;m+=7) {
                    auto t1=cet.add(t0,calendar::MONTH, m);
                    utctimespan remainder;
                    auto n_m=cet.diff_units(t0,t1,calendar::MONTH,remainder);
                    TS_ASSERT_EQUALS(remainder,utctimespan(0));// should match exactly
                    TS_ASSERT_EQUALS(n_m,m);
                    n_m= cet.diff_units(t1,t0,calendar::MONTH,remainder);
                    TS_ASSERT_EQUALS(remainder,utctimespan(0));// should match exactly
                    TS_ASSERT_EQUALS(n_m,-m);
                    auto t2=cet.add(t1,calendar::MONTH,-m);
                    auto dt1=cet.calendar_units(t1);
                    if( t0 != t2 ) {
                        // cout << "case "<< cet.to_string(t0)<<" vs "<<cet.to_string(t2)<<" error at m="<<m<<"d="<<day<<" diff is "<<t0-t2<<endl;
                        TS_ASSERT_EQUALS(t0,t2);
                        return;
                    }
                    TS_ASSERT_EQUALS(t0,t2);
                    dt1.month=dt0.month;//just put back this one to ease comparison, all other should be equal
                    dt1.year=dt0.year;
                    TS_ASSERT_EQUALS(cet.to_string(cet.time(dt1)),cet.to_string(cet.time(dt0)));
                }
            }
        }
    }

}
