//!#####################################################################
//! \file Heat_Diffusion_Example.cpp
//!#####################################################################
#include <nova/Dynamics/Hierarchy/Grid_Hierarchy_Initializer.h>
#include <nova/Dynamics/Hierarchy/Grid_Hierarchy_Iterator.h>
#include <nova/Dynamics/Hierarchy/Advection/Grid_Hierarchy_Advection.h>
#include <nova/Tools/Grids/Grid_Iterator_Face.h>
#include <nova/Tools/Krylov_Solvers/Conjugate_Gradient.h>
#include <nova/Tools/Utilities/File_Utilities.h>

#include "Heat_Diffusion_Example.h"
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
template<class T,int d> Heat_Diffusion_Example<T,d>::
Heat_Diffusion_Example()
    :Base(),hierarchy(nullptr),rasterizer(nullptr)
{
    // face_velocity_channels(0)           = &Struct_type::ch0;
    // face_velocity_channels(1)           = &Struct_type::ch1;
    // if(d==3) face_velocity_channels(2)  = &Struct_type::ch2;
    // face_qc_channels(0)                 = &Struct_type::ch3;
    // face_qc_channels(1)                 = &Struct_type::ch4;
    // if(d==3) face_qc_channels(2)        = &Struct_type::ch5;
    density_channel                     = &Struct_type::ch6;
    density_backup_channel              = &Struct_type::ch7;
}
//######################################################################
// Initialize
//######################################################################
template<class T,int d> void Heat_Diffusion_Example<T,d>::
Initialize()
{
    diffusion_rt=(T)0.; qc_advection_rt=(T)0.; 
    Initialize_SPGrid();
    Initialize_Fluid_State(test_number);
}
//######################################################################
// Initialize_SPGrid
//######################################################################
template<class T,int d> void Heat_Diffusion_Example<T,d>::
Initialize_SPGrid()
{
    Log::Scope scope("Initialize_SPGrid");
    // Initialize_Rasterizer(test_number);
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
template<class T,int d> void Heat_Diffusion_Example<T,d>::
Limit_Dt(T& dt,const T time)
{
}
//######################################################################
// Advect_Density
//######################################################################
template    <class T,int d> void Heat_Diffusion_Example<T,d>::
Advect_Density(const T dt)
{
    Channel_Vector cell_velocity_channels;
    cell_velocity_channels(0)               = &Struct_type::ch8;
    cell_velocity_channels(1)               = &Struct_type::ch9;
    if(d==3) cell_velocity_channels(2)      = &Struct_type::ch10;
    T Struct_type::* temp_channel           = &Struct_type::ch11;
    Vector<uint64_t,d> other_face_offsets;
    for(int axis=0;axis<d;++axis) other_face_offsets(axis)=Topology_Helper::Axis_Vector_Offset(axis);
    Uniform_Grid_Averaging_Helper<Struct_type,T,d>::Uniform_Grid_Average_Face_Velocities_To_Cells(*hierarchy,hierarchy->Allocator(0),hierarchy->Blocks(0),face_velocity_channels,cell_velocity_channels,other_face_offsets);
    Uniform_Grid_Advection_Helper<Struct_type,T,d>::Uniform_Grid_Advect_Density(*hierarchy,cell_velocity_channels,density_channel,temp_channel,dt);
}
//######################################################################
// Register_Options
//######################################################################
template<class T,int d> void Heat_Diffusion_Example<T,d>::
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
    // parse_args->Add_Double_Argument("-diff_coeff",(T)1e-3,"diffusion coefficient.");
    // parse_args->Add_Double_Argument("-fc",(T)0.,"fc.");
    // parse_args->Add_Double_Argument("-bv",(T)1.,"Background velocity(along y axis).");
    // parse_args->Add_Double_Argument("-sr",(T)1.,"Source rate");
    // parse_args->Add_Double_Argument("-tau",(T)1.,"tau.");
    // parse_args->Add_Option_Argument("-ficks","Fick's diffusion.");
    // parse_args->Add_Option_Argument("-nd","Turn off diffusion");
    // parse_args->Add_Option_Argument("-ed","Explicit diffusion");
    // parse_args->Add_Option_Argument("-uvf","Uniform velocity field");
    // for CG
    parse_args->Add_Integer_Argument("-cg_iterations",100,"Number of CG iterations.");
    parse_args->Add_Integer_Argument("-cg_restart_iterations",40,"Number of CG restart iterations.");
    parse_args->Add_Double_Argument("-cg_tolerance",1e-4,"CG tolerance");
}
//######################################################################
// Parse_Options
//######################################################################
template<class T,int d> void Heat_Diffusion_Example<T,d>::
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
    // diff_coeff=parse_args->Get_Double_Value("-diff_coeff");
    // bv=parse_args->Get_Double_Value("-bv");
    // source_rate=parse_args->Get_Double_Value("-sr");
    // Fc=parse_args->Get_Double_Value("-fc");
    // tau=parse_args->Get_Double_Value("-tau");
    // FICKS=parse_args->Get_Option_Value("-ficks");
    // uvf=parse_args->Get_Option_Value("-uvf");
    // nd=parse_args->Get_Option_Value("-nd");
    // explicit_diffusion=parse_args->Get_Option_Value("-ed");
    cg_iterations=parse_args->Get_Integer_Value("-cg_iterations");
    cg_restart_iterations=parse_args->Get_Integer_Value("-cg_restart_iterations");
    cg_tolerance=(T)parse_args->Get_Double_Value("-cg_tolerance");
    // if(nd){diff_coeff=(T)0.;tau=(T)1.;Fc=(T)0.;}
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
template<class T,int d> void Heat_Diffusion_Example<T,d>::
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
template<class T,int d> void Heat_Diffusion_Example<T,d>::
Read_Output_Files(const int frame)
{
}
//######################################################################
template class Nova::Heat_Diffusion_Example<float,2>;
template class Nova::Heat_Diffusion_Example<float,3>;
#ifdef COMPILE_WITH_DOUBLE_SUPPORT
template class Nova::Heat_Diffusion_Example<double,2>;
template class Nova::Heat_Diffusion_Example<double,3>;
#endif