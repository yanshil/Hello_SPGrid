//!#####################################################################
//! \file Heat_Transfer_Driver.h
//!#####################################################################
// Class Heat_Transfer_Driver
//######################################################################
#ifndef __Heat_Transfer_Driver__
#define __Heat_Transfer_Driver__

#include <nova/Tools/Utilities/Driver.h>
#include "Heat_Transfer_Example.h"

namespace Nova{
template<class T,int d>
class Heat_Transfer_Driver: public Driver<T,d>
{
    using TV                = Vector<T,d>;
    using Base              = Driver<T,d>;

    using Base::time;
    using Base::Compute_Dt;using Base::Write_Output_Files;

  public:
    // int substep_counter;
    // T density_advection_rt;
    // T velocity_advection_rt;
    // T source_modification_rf;
    // T projection_rt;
    // T total_rt;
    Heat_Transfer_Example<T,d>& example;

    Heat_Transfer_Driver(Heat_Transfer_Example<T,d>& example_input);
    ~Heat_Transfer_Driver() {}

//######################################################################
    void Initialize() override;
    void Advance_One_Time_Step_Explicitly(const T dt,const T time);
    // void Advance_One_Time_Step_Implicitly(const T dt,const T time);
    void Advance_To_Target_Time(const T target_time) override;
    void Simulate_To_Frame(const int frame) override;
//######################################################################
};
}
#endif