//!#####################################################################
//! \file Heat_Transfer_Example.h
//!#####################################################################
// Class Heat_Transfer_Example
//######################################################################
#ifndef __Heat_Transfer_Example__
#define __Heat_Transfer_Example__

#include <nova/Dynamics/Hierarchy/Grid_Hierarchy.h>
#include <nova/Dynamics/Hierarchy/Grid_Topology_Helper.h>
#include <nova/Dynamics/Hierarchy/Rasterizers/Hierarchical_Rasterizer.h>
#include <nova/Geometry/Implicit_Objects/Implicit_Object.h>
#include <nova/Tools/Utilities/Example.h>
#include "Poisson_Data.h"

namespace Nova{
template<class T,int d>
class Heat_Transfer_Example: public Example<T,d>
{
    using TV                        = Vector<T,d>;
    using Base                      = Example<T,d>;
    using T_INDEX                   = Vector<int,d>;
    using Struct_type               = Poisson_Data<T>;
    using Channel_Vector            = Vector<T Struct_type::*,d>;
    using Hierarchy                 = Grid_Hierarchy<Struct_type,T,d>;
    using Hierarchy_Rasterizer      = Hierarchical_Rasterizer<Struct_type,T,d>;
    using Allocator_type            = SPGrid::SPGrid_Allocator<Struct_type,d>;
    using Flag_array_mask           = typename Allocator_type::template Array_mask<unsigned>;
    using Topology_Helper           = Grid_Topology_Helper<Flag_array_mask>;

  public:
    using Base::frame_title;using Base::output_directory;using Base::parse_args;using Base::first_frame;
    int test_number;
    T_INDEX counts;
    int levels,mg_levels,cg_iterations,cg_restart_iterations;
    T cfl;
    T cg_tolerance;
    
    Hierarchy *hierarchy;
    Hierarchy_Rasterizer *rasterizer;

    T Struct_type::* pressure_channel;
    T Struct_type::* density_channel;
    // T Struct_type::* density_backup_channel;
    
    Vector<T Struct_type::*,d> face_velocity_channels;
    // Vector<T Struct_type::*,d> face_qc_channels;
    Vector<Vector<bool,2>,d> domain_walls;

    Array<Implicit_Object<T,d>*> velocity_sources;
    Array<Implicit_Object<T,d>*> density_sources;

    Heat_Transfer_Example();

    ~Heat_Transfer_Example()
    {if(hierarchy!=nullptr) delete hierarchy;}

//######################################################################
// --- Will be implemented in Standard_test.h
    virtual void Initialize_Rasterizer(const int test_number)=0;
    virtual void Initialize_Fluid_State(const int test_number)=0;
    virtual void Initialize_Sources(const int test_number)=0;
    virtual void Set_Boundary(const int test_number)=0;
//######################################################################
    void Initialize();
    void Initialize_SPGrid();
    void Limit_Dt(T& dt,const T time) override;
    void Advect_Density(const T dt);
    // void Advect_Face_Vector(const T dt);
    // void Diffuse_Density(const T dt);
    // void Backup_Density();
    // void Ficks_Diffusion(const T dt);
    // void Non_Ficks_Diffusion(const T dt);
    void Modify_Density_With_Sources();
    // void Add_Source(const T dt);
    void Advect_Face_Velocities(const T dt);
    // void Set_Neumann_Faces_Inside_Sources();
    // void Initialize_Velocity_Field();
    // void Project();
    void Register_Options() override;
    void Parse_Options() override;
    void Read_Output_Files(const int frame);
    void Write_Output_Files(const int frame) const override;
    // void Wrtie_To_File(const int frame);
//######################################################################
};
}


#endif