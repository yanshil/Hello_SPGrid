//!#####################################################################
//! \file Standard_Tests.h
//!#####################################################################
// Class Standard_Tests
//######################################################################
#ifndef __Standard_Tests__
#define __Standard_Tests__

#include <nova/Geometry/Implicit_Objects/Box_Implicit_Object.h>
#include <nova/SPGrid/Tools/SPGrid_Clear.h>
#include <nova/Tools/Utilities/Range_Iterator.h>
#include "../../Poisson_Data.h"
#include "../../Heat_Transfer_Example.h"
#include "../../Rasterizers/Adaptive_Sphere_Rasterizer.h"
#include "../../Rasterizers/Randomized_Rasterizer.h"


namespace Nova{
template<class T,int d>
class Standard_Tests: public Heat_Transfer_Example<T,d>
{
    
    using TV                        = Vector<T,d>;
    using T_INDEX                   = Vector<int,d>;
    using Struct_type               = Poisson_Data<T>;
    using Hierarchy                 = Grid_Hierarchy<Struct_type,T,d>;
    using Flags_type                = typename Struct_type::Flags_type;
    using Base                      = Heat_Transfer_Example<T,d>;
    using Allocator_type            = SPGrid::SPGrid_Allocator<Struct_type,d>;
    using Flag_array_mask           = typename Allocator_type::template Array_mask<unsigned>;

  public:
    using Base::output_directory; using Base::test_number;using Base::counts;
    using Base::levels;
    using Base::domain_walls;
    using Base::hierarchy;
    using Base::rasterizer;
    using Base::cfl;
    using Base::velocity_sources;   using Base::density_sources;    
    using Base::density_channel;


    Standard_Tests()
        :Base()
    {}

//######################################################################
    void Parse_Options() override
    {
        Base::Parse_Options();
        output_directory="HeatDiffusion_"+std::to_string(d)+"d_"+"case_"+std::to_string(test_number)+"_Resolution_"+std::to_string(counts(0))+"x"+std::to_string(counts(1));

        for(int axis=0;axis<d;++axis) for(int side=0;side<2;++side) domain_walls(axis)(side)=false;
        
        // TV min_corner,max_corner=TV(4.);
        // max_corner(1)=(T)8.;
        TV min_corner,max_corner=TV(1);
        hierarchy=new Hierarchy(counts,Range<T,d>(min_corner,max_corner),levels);
    }
//######################################################################
    void Initialize_Rasterizer(const int test_number) override
    {
        // rasterizer=new Randomized_Rasterizer<Struct_type,T,d>(*hierarchy);
        rasterizer=new Adaptive_Sphere_Rasterizer<Struct_type,T,d>(*hierarchy,TV(.5),(T).1);
    }
//######################################################################
    void Initialize_Fluid_State(const int test_number) override
    {
        // clear density channel
        for(int level=0;level<levels;++level)
            SPGrid::Clear<Struct_type,T,d>(hierarchy->Allocator(level),hierarchy->Blocks(level),density_channel);

        for(int level=0;level<levels;++level){auto blocks=hierarchy->Blocks(level);
            auto block_size=hierarchy->Allocator(level).Block_Size();
            auto data=hierarchy->Allocator(level).template Get_Array<Struct_type,T>(density_channel);
            auto flags=hierarchy->Allocator(level).template Get_Const_Array<Struct_type,unsigned>(&Struct_type::flags);

            for(unsigned block=0;block<blocks.second;++block){uint64_t offset=blocks.first[block];
                Range_Iterator<d> range_iterator(T_INDEX(),*reinterpret_cast<T_INDEX*>(&block_size)-1);
                T_INDEX base_index(Flag_array_mask::LinearToCoord(offset));

                for(int e=0;e<Flag_array_mask::elements_per_block;++e,offset+=sizeof(Flags_type)){
                    const T_INDEX index=base_index+range_iterator.Index();
                    if(flags(offset)&Cell_Type_Interior && density_sources(0)->Inside(hierarchy->Lattice(level).Center(index))) data(offset)=(T)1.;
                    range_iterator.Next();}}}
    }
//######################################################################
    void Initialize_Sources(const int test_number) override
    {
        TV min_corner=TV({.45,0}),max_corner=TV({.55,.05});	
        Implicit_Object<T,d>* obj=new Box_Implicit_Object<T,d>(min_corner,max_corner);
        density_sources.Append(obj);
        // const T cell_width=(T)4./counts(0);
        // switch (test_number)
        // {
        // // test case 1: density&velocity source near the bottom 
        // case 1:{

        // }break;
        // case 2:{
        //     TV density_min_corner=TV({(T)1.8,(T)0.}),density_max_corner=TV({(T)2.2,(T)2.*cell_width});
        //     Implicit_Object<T,d>* density_obj=new Box_Implicit_Object<T,d>(density_min_corner,density_max_corner);
        //     density_sources.Append(density_obj);

        //     TV velocity_min_corner=TV({(T)1.8,(T)0.}),velocity_max_corner=TV({(T)2.2,(T)2.*cell_width});
        //     Implicit_Object<T,d>* velocity_obj=new Box_Implicit_Object<T,d>(velocity_min_corner,velocity_max_corner);
        //     velocity_sources.Append(velocity_obj);}break;
        // case 3:
        // case 4:{
        //     TV density_min_corner=TV({(T)2.-cell_width,(T)2.-cell_width}),density_max_corner=TV({(T)2.+cell_width,(T)2.+cell_width});
        //     Implicit_Object<T,d>* density_obj=new Box_Implicit_Object<T,d>(density_min_corner,density_max_corner);
        //     density_sources.Append(density_obj);
        //     TV velocity_min_corner=TV({(T)1.8,(T)0.}),velocity_max_corner=TV({(T)2.2,(T)2.*cell_width});
        //     Implicit_Object<T,d>* velocity_obj=new Box_Implicit_Object<T,d>(velocity_min_corner,velocity_max_corner);
        //     velocity_sources.Append(velocity_obj);
        // }break;
        // case 5:
        // case 6:
        // case 7:
        // case 8:{
        //     TV density_min_corner=TV({(T)2.-cell_width,(T)2.-cell_width}),density_max_corner=TV({(T)2.+cell_width,(T)2.+cell_width});
        //     Implicit_Object<T,d>* density_obj=new Box_Implicit_Object<T,d>(density_min_corner,density_max_corner);
        //     density_sources.Append(density_obj);
        //     TV velocity_min_corner=TV({(T)2.-cell_width,(T)2.-cell_width}),velocity_max_corner=TV({(T)2.+cell_width,(T)2.+cell_width});
        //     Implicit_Object<T,d>* velocity_obj=new Box_Implicit_Object<T,d>(velocity_min_corner,velocity_max_corner);
        //     velocity_sources.Append(velocity_obj);}break;}


    }
//######################################################################
    void Set_Boundary(const int test_number) override
    {

    }
//######################################################################
};
}

#endif