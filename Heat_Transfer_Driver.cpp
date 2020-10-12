//!#####################################################################
//! \file Heat_Transfer_Driver.cpp
//!#####################################################################
#include <chrono>
#include "Heat_Transfer_Driver.h"
using namespace std::chrono;
using namespace Nova;
//######################################################################
// Constructor
//######################################################################
template<class T,int d> Heat_Transfer_Driver<T,d>::
Heat_Transfer_Driver(Heat_Transfer_Example<T,d>& example_input)
    :Base(example_input),example(example_input)
{}
//######################################################################
// Initialize
//######################################################################
template<class T,int d> void Heat_Transfer_Driver<T,d>::
Initialize()
{
    time=(T)0.;
    substep_counter=0; 
    total_rt=(T)0.;
    // density_advection_rt=(T)0.;
    // velocity_advection_rt=(T)0.; source_modification_rf=(T)0.; projection_rt=(T)0.;
    Base::Initialize();

    example.Log_Parameters();
    example.Initialize_Sources(example.test_number);

    if(!example.restart) example.Initialize();
    else example.Read_Output_Files(example.restart_frame);

    // example.Initialize_Velocity_Field();
    // // divergence free
    // if(!example.uvf) example.Project();
}
//######################################################################
// Advance_One_Time_Step_Explicitly
//######################################################################
template<class T,int d> void Heat_Transfer_Driver<T,d>::
Advance_One_Time_Step_Explicitly(const T dt,const T time)
{
    example.Advect_Density(dt);
    // ...

    example.Advect_Face_Velocities(dt);
}
// //######################################################################
// // Advance_One_Time_Step_Implicitly
// //######################################################################
// template<class T,int d> void Heat_Transfer_Driver<T,d>::
// Advance_One_Time_Step_Implicitly(const T dt,const T time)
// {
//     high_resolution_clock::time_point tb=high_resolution_clock::now();
//     if(!example.uvf) example.Project();
//     high_resolution_clock::time_point te=high_resolution_clock::now();
//     projection_rt+=duration_cast<duration<T>>(te-tb).count();
// }
//######################################################################
// Advance_To_Target_Time
//######################################################################
template<class T,int d> void Heat_Transfer_Driver<T,d>::
Advance_To_Target_Time(const T target_time)
{
    bool done=false;
    for(int substep=1;!done;substep++){
        // high_resolution_clock::time_point tb=high_resolution_clock::now();
        Log::Scope scope("SUBSTEP","substep "+std::to_string(substep));
        substep_counter++;
        T dt=Compute_Dt(time,target_time);

        Example<T,d>::Clamp_Time_Step_With_Target_Time(time,target_time,dt,done);

        Advance_One_Time_Step_Explicitly(dt,time);
        // Advance_One_Time_Step_Implicitly(dt,time);

        Log::cout<<"dt: "<<dt<<std::endl;
        
        if(!done) example.Write_Substep("END Substep",substep,0);
        time+=dt;
        // high_resolution_clock::time_point te=high_resolution_clock::now();
        // total_rt+=duration_cast<duration<T>>(te-tb).count();
    }
}
//######################################################################
// Simulate_To_Frame
//######################################################################
template<class T,int d> void Heat_Transfer_Driver<T,d>::
Simulate_To_Frame(const int target_frame)
{
    example.frame_title="Frame "+std::to_string(example.current_frame);
    if(!example.restart) Write_Output_Files(example.current_frame);

    while(example.current_frame<target_frame){
        Log::Scope scope("FRAME","Frame "+std::to_string(++example.current_frame));

        Advance_To_Target_Time(example.Time_At_Frame(example.current_frame));
        
        example.frame_title="Frame "+std::to_string(example.current_frame);
        Write_Output_Files(++example.output_number);        
        *(example.output)<<"TIME = "<<time<<std::endl;
        int substeps=substep_counter;
        Log::cout<<"Average: "<<std::endl;
        Log::cout<<"Total substeps: "<<substeps<<std::endl;
        Log::cout<<"full timestep: "<<total_rt/substeps<<std::endl;
    }
}
//######################################################################
template class Nova::Heat_Transfer_Driver<float,2>;
template class Nova::Heat_Transfer_Driver<float,3>;
#ifdef COMPILE_WITH_DOUBLE_SUPPORT
template class Nova::Heat_Transfer_Driver<double,2>;
template class Nova::Heat_Transfer_Driver<double,3>;
#endif