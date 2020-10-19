//!#####################################################################
//! \file Heat_Transfer_Example.cpp
//!#####################################################################
#include <nova/Dynamics/Hierarchy/Grid_Hierarchy_Initializer.h>
#include <nova/Dynamics/Hierarchy/Grid_Hierarchy_Iterator.h>
#include <nova/Dynamics/Hierarchy/Advection/Grid_Hierarchy_Advection.h>
#include <nova/Dynamics/Hierarchy/Interpolation/Grid_Hierarchy_Averaging.h>
#include <nova/Tools/Grids/Grid_Iterator_Face.h>
#include <nova/Tools/Krylov_Solvers/Conjugate_Gradient.h>
#include <nova/Tools/Utilities/File_Utilities.h>

#include "Compute_Time_Step.h"
#include "Density_Modifier.h"
#include "Heat_Transfer_Example.h"
#include "Initialize_Dirichlet_Cells.h"
#include "Write_To_File_Helper.h"
#include <omp.h>
#include <chrono>
using namespace std::chrono;
using namespace Nova;
namespace Nova{
extern int number_of_threads;
}
//######################################################################
// Constructor
//######################################################################
template<class T,int d> Heat_Transfer_Example<T,d>::
Heat_Transfer_Example()
    :Base(),hierarchy(nullptr),rasterizer(nullptr)
{
    face_velocity_channels(0)           = &Struct_type::ch0;
    face_velocity_channels(1)           = &Struct_type::ch1;
    if(d==3) face_velocity_channels(2)  = &Struct_type::ch2;
    pressure_channel                    = &Struct_type::ch3;
    density_channel                     = &Struct_type::ch4;
}
//######################################################################
// Initialize
//######################################################################
template<class T,int d> void Heat_Transfer_Example<T,d>::
Initialize()
{
    Initialize_SPGrid();
    Initialize_Fluid_State(test_number);
}
//######################################################################
// Initialize_SPGrid
//######################################################################
template<class T,int d> void Heat_Transfer_Example<T,d>::
Initialize_SPGrid()
{
    Log::Scope scope("Initialize_SPGrid");
    Initialize_Rasterizer(test_number);
    for(Grid_Hierarchy_Iterator<d,Hierarchy_Rasterizer> iterator(hierarchy->Lattice(levels-1).Cell_Indices(),levels-1,*rasterizer);iterator.Valid();iterator.Next());
    Grid_Hierarchy_Initializer<Struct_type,T,d>::Flag_Ghost_Cells(*hierarchy);
    Grid_Hierarchy_Initializer<Struct_type,T,d>::Flag_Valid_Faces(*hierarchy);
    Grid_Hierarchy_Initializer<Struct_type,T,d>::Flag_Active_Faces(*hierarchy);
    Grid_Hierarchy_Initializer<Struct_type,T,d>::Flag_Active_Nodes(*hierarchy);
    Grid_Hierarchy_Initializer<Struct_type,T,d>::Flag_Shared_Nodes(*hierarchy);
    Grid_Hierarchy_Initializer<Struct_type,T,d>::Flag_Ghost_Nodes(*hierarchy);
    Grid_Hierarchy_Initializer<Struct_type,T,d>::Flag_T_Junction_Nodes(*hierarchy);
    Initialize_Dirichlet_Cells<Struct_type,T,d>(*hierarchy,domain_walls);
    //Set_Neumann_Faces_Inside_Sources();
    hierarchy->Update_Block_Offsets();
    hierarchy->Initialize_Red_Black_Partition(2*number_of_threads);
    
}
//######################################################################
// Limit_Dt
//######################################################################
template<class T,int d> void Heat_Transfer_Example<T,d>::
Limit_Dt(T& dt,const T time)
{
    T dt_convection=(T)0.;

    Vector<uint64_t,d> other_face_offsets;
    for(int axis=0;axis<d;++axis) other_face_offsets(axis)=Topology_Helper::Axis_Vector_Offset(axis);

    for(int level=0;level<levels;++level)
        Compute_Time_Step<Struct_type,T,d>(*hierarchy,hierarchy->Blocks(level),face_velocity_channels,
                                           other_face_offsets,level,dt_convection);

    if(dt_convection>(T)1e-5) dt=cfl/dt_convection;
    Log::cout<<"Time Step: "<<dt<<std::endl;
}
//######################################################################
// Advect_Density
//######################################################################
template    <class T,int d> void Heat_Transfer_Example<T,d>::
Advect_Density(const T dt)
{
    Channel_Vector node_velocity_channels;
    node_velocity_channels(0)               = &Struct_type::ch5;
    node_velocity_channels(1)               = &Struct_type::ch6;
    if(d==3) node_velocity_channels(2)      = &Struct_type::ch7;
    T Struct_type::* node_density_channel   = &Struct_type::ch8;
    T Struct_type::* temp_channel           = &Struct_type::ch9;
    Grid_Hierarchy_Averaging<Struct_type,T,d>::Average_Cell_Density_To_Nodes(*hierarchy,density_channel,node_density_channel,temp_channel);
    Grid_Hierarchy_Averaging<Struct_type,T,d>::Average_Face_Velocities_To_Nodes(*hierarchy,face_velocity_channels,node_velocity_channels,temp_channel);
    Grid_Hierarchy_Advection<Struct_type,T,d>::Advect_Density(*hierarchy,node_velocity_channels,density_channel,node_density_channel,temp_channel,dt);
}
//######################################################################
// Modify_Density_With_Sources
//######################################################################
template<class T,int d> void Heat_Transfer_Example<T,d>::
Modify_Density_With_Sources()
{
    for(int level=0;level<levels;++level)
        Density_Modifier<Struct_type,T,d>(*hierarchy,hierarchy->Blocks(level),density_channel,density_sources,level);
}
//######################################################################
// Advect_Face_Velocities
//######################################################################
template<class T,int d> void Heat_Transfer_Example<T,d>::
Advect_Face_Velocities(const T dt)
{
    Channel_Vector node_velocity_channels;
    node_velocity_channels(0)               = &Struct_type::ch5;
    node_velocity_channels(1)               = &Struct_type::ch6;
    if(d==3) node_velocity_channels(2)      = &Struct_type::ch7;
    T Struct_type::* temp_channel           = &Struct_type::ch9;
    Grid_Hierarchy_Advection<Struct_type,T,d>::Advect_Face_Velocities(*hierarchy,face_velocity_channels,node_velocity_channels,temp_channel,dt);
}
//######################################################################
// Register_Options
//######################################################################
template<class T,int d> void Heat_Transfer_Example<T,d>::
Register_Options()
{
    Base::Register_Options();

    // parse_args->Add_Double_Argument("-cfl",(T).5,"CFL number.");
    parse_args->Add_Integer_Argument("-levels",1,"Number of levels in the SPGrid hierarchy.");
    parse_args->Add_Integer_Argument("-test_number",1,"Test number.");
    parse_args->Add_Integer_Argument("-mg_levels",1,"Number of levels in the Multigrid hierarchy.");
    parse_args->Add_Integer_Argument("-threads",1,"Number of threads for OpenMP to use");
    if(d==2) parse_args->Add_Vector_2D_Argument("-size",Vector<double,2>(64.),"n","Grid resolution");
    else if(d==3) parse_args->Add_Vector_3D_Argument("-size",Vector<double,3>(64.),"n","Grid resolution");

    // for CG
    parse_args->Add_Integer_Argument("-cg_iterations",100,"Number of CG iterations.");
    parse_args->Add_Integer_Argument("-cg_restart_iterations",40,"Number of CG restart iterations.");
    parse_args->Add_Double_Argument("-cg_tolerance",1e-4,"CG tolerance");
}
//######################################################################
// Parse_Options
//######################################################################
template<class T,int d> void Heat_Transfer_Example<T,d>::
Parse_Options()
{
    Base::Parse_Options();

    // cfl=(T)parse_args->Get_Double_Value("-cfl");
    levels=parse_args->Get_Integer_Value("-levels");
    test_number=parse_args->Get_Integer_Value("-test_number");
    mg_levels=parse_args->Get_Integer_Value("-mg_levels");
    number_of_threads=parse_args->Get_Integer_Value("-threads");
    omp_set_num_threads(number_of_threads);
    if(d==2){auto cell_counts_2d=parse_args->Get_Vector_2D_Value("-size");for(int v=0;v<d;++v) counts(v)=cell_counts_2d(v);}
    else{auto cell_counts_3d=parse_args->Get_Vector_3D_Value("-size");for(int v=0;v<d;++v) counts(v)=cell_counts_3d(v);}

    cg_iterations=parse_args->Get_Integer_Value("-cg_iterations");
    cg_restart_iterations=parse_args->Get_Integer_Value("-cg_restart_iterations");
    cg_tolerance=(T)parse_args->Get_Double_Value("-cg_tolerance");

    // switch (test_number){
    // case 1:{const_density_source=false;const_density_value=(T)0.;uvf=false;}break;
    // case 2:{const_density_source=true;const_density_value=(T)1.;uvf=false;}break;
    // case 3:{const_density_source=false;const_density_value=(T)0.;uvf=false;}break;
    // case 4:{const_density_source=true;const_density_value=(T)1.;uvf=false;}break;
    // case 5:{const_density_source=false;const_density_value=(T)0.;uvf=true;}break;
    // case 6:{const_density_source=true;const_density_value=(T)1.;uvf=true;}break;
    // case 7:{const_density_source=true;const_density_value=(T)1.;uvf=false;}break;
    // case 8:{const_density_source=true;const_density_value=(T)1.;uvf=false;}break;}
}
//######################################################################
// Write_Output_Files
//######################################################################
template<class T,int d> void Heat_Transfer_Example<T,d>::
Write_Output_Files(const int frame) const
{
    File_Utilities::Create_Directory(output_directory+"/"+std::to_string(frame));
    File_Utilities::Write_To_Text_File(output_directory+"/"+std::to_string(frame)+"/frame_title",frame_title);

    // write hierarchy
    File_Utilities::Write_To_Text_File(output_directory+"/"+std::to_string(frame)+"/levels",levels);
    hierarchy->Write_Hierarchy(output_directory,frame);
    hierarchy->template Write_Channel<T>(output_directory+"/"+std::to_string(frame)+"/spgrid_density",density_channel);
}
//######################################################################
// Read_Output_Files
//######################################################################
template<class T,int d> void Heat_Transfer_Example<T,d>::
Read_Output_Files(const int frame)
{
}
//######################################################################
template class Nova::Heat_Transfer_Example<float,2>;
template class Nova::Heat_Transfer_Example<float,3>;
#ifdef COMPILE_WITH_DOUBLE_SUPPORT
template class Nova::Heat_Transfer_Example<double,2>;
template class Nova::Heat_Transfer_Example<double,3>;
#endif