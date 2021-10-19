/******************************************************************************
*                 SOFA, Simulation Open-Framework Architecture                *
*                    (c) 2006 INRIA, USTL, UJF, CNRS, MGH                     *
*                                                                             *
* This program is free software; you can redistribute it and/or modify it     *
* under the terms of the GNU Lesser General Public License as published by    *
* the Free Software Foundation; either version 2.1 of the License, or (at     *
* your option) any later version.                                             *
*                                                                             *
* This program is distributed in the hope that it will be useful, but WITHOUT *
* ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or       *
* FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License *
* for more details.                                                           *
*                                                                             *
* You should have received a copy of the GNU Lesser General Public License    *
* along with this program. If not, see <http://www.gnu.org/licenses/>.        *
*******************************************************************************
* Authors: The SOFA Team and external contributors (see Authors.txt)          *
*                                                                             *
* Contact information: contact@sofa-framework.org                             *
******************************************************************************/
#pragma once
#include <SofaOrigami/config.h>

//#include <sofa/core/SofaDeformable>


//#include <SofaDeformable/config.h>

#include <SofaDeformable/SpringForceField.h>
#include <sofa/type/Mat.h>
#include <SofaBaseTopology/TopologySubsetIndices.h>

namespace sofa::component::interactionforcefield
{

/** SpringForceField able to evaluate and apply its stiffness.
This allows to perform implicit integration.
Stiffness is evaluated and stored by the addForce method.
When explicit integration is used, SpringForceField is slightly more efficient.
*/
class OrigamiEdge
{
public:
    int  m1, m2, m3, m4, foldType;            ///< the two extremities of the spring: masses m1 and m2, and the two faces opposite nodes and type of folde: montain +1, valley -1 or none 0
    double ks;                ///< spring stiffness
//    Real kcrease;            /// face angular stiffness
//    Real kface;            /// face stiffness
    double kd;                ///< damping factor
    double initpos;           ///< rest length of the spring
    bool elongationOnly;    ///< only forbid elongation, not compression
    bool enabled;           ///< false to disable this spring (i.e. broken)

    OrigamiEdge(int m1=0, int m2=0, int m3=0, int m4=0, double ks=0.0, double kd=0.0, double initpos=0.0, bool noCompression=false, bool enabled=true)
        : m1(m1), m2(m2), m3(m3), m4(m4), foldType(foldType), ks(ks), kd(kd), initpos(initpos), elongationOnly(noCompression), enabled(enabled)
    {
    }

    inline friend std::istream& operator >> ( std::istream& in, OrigamiEdge& s )
    {
        in>>s.m1>>s.m2>>s.m3>>s.m4>>s.foldType>>s.ks>>s.kd>>s.initpos;
        return in;
    }

    inline friend std::ostream& operator << ( std::ostream& out, const OrigamiEdge& s )
    {
        out<<s.m1<<" "<<s.m2<<" "<<s.m3<<" "<<s.m4<<" "<<s.foldType<<" "<<s.ks<<" "<<s.kd<<" "<<s.initpos<<"\n";
        return out;
    }

};


template<class DataTypes>
class OrigamiForceField : public sofa::component::interactionforcefield::SpringForceField<DataTypes>
{
public:
    SOFA_CLASS(SOFA_TEMPLATE(OrigamiForceField,DataTypes), SOFA_TEMPLATE(SpringForceField,DataTypes));

    typedef SpringForceField<DataTypes> Inherit;
    typedef typename DataTypes::VecCoord VecCoord;
    typedef typename DataTypes::VecDeriv VecDeriv;
    typedef typename DataTypes::Coord Coord;
    typedef typename DataTypes::Deriv Deriv;
    typedef typename Coord::value_type Real;

    typedef core::objectmodel::Data<VecDeriv>    DataVecDeriv;
    typedef core::objectmodel::Data<VecCoord>    DataVecCoord;
    typedef type::vector<sofa::Index> SetIndexArray;
    typedef sofa::component::topology::TopologySubsetIndices SetIndex;


    typedef typename Inherit::Spring Spring;

    typedef core::behavior::MechanicalState<DataTypes> MechanicalState;
    enum { N=DataTypes::spatial_dimensions };
    typedef type::Mat<N,N,Real> Mat;

    SetIndex d_indices1; ///< Indices of the source points on the first model
    SetIndex d_indices2; ///< Indices of the fixed points on the second model
    core::objectmodel::Data<double> angleTarget;
    core::objectmodel::Data<SReal> d_length;
    core::objectmodel::Data<type::vector<OrigamiEdge> > origamiEdges; ///< pairs of indices and face opposite nodes, stiffness, damping, rest length
    std::vector<std::vector<int> >  Triangles;
protected:
    sofa::type::vector<Mat>  dfdx;

    void addFaceConstraints(Real& potentialEnergy, VecDeriv& f1, const  VecCoord& p1, const VecDeriv& v1, VecDeriv& f2, const  VecCoord& p2, const  VecDeriv& v2, int i, std::vector<int> triangleInOrigami);

    /// Accumulate the spring force and compute and store its stiffness
    void addSpringForce(Real& potentialEnergy, VecDeriv& f1,const  VecCoord& p1,const VecDeriv& v1, VecDeriv& f2,const  VecCoord& p2,const  VecDeriv& v2, sofa::Index i, const OrigamiEdge& spring);

    /// Apply the stiffness, i.e. accumulate df given dx
    virtual void addSpringDForce(VecDeriv& df1,const  VecDeriv& dx1, VecDeriv& df2,const  VecDeriv& dx2, sofa::Index i, const Spring& spring, double kFactor, double bFactor);

    OrigamiForceField(double ks=100.0, double kd=5.0);
    OrigamiForceField(MechanicalState* object1, MechanicalState* object2, double ks=100.0, double kd=5.0);

    /// Will create the set of springs using \sa d_indices1 and \sa d_indices2 with \sa d_length
    void createSpringsFromInputs();

public:
    void init() override;

    /// Accumulate f corresponding to x,v
    void addForce(const sofa::core::MechanicalParams* mparams, DataVecDeriv& data_f1, DataVecDeriv& data_f2, const DataVecCoord& data_x1, const DataVecCoord& data_x2, const DataVecDeriv& data_v1, const DataVecDeriv& data_v2 ) override;
    /// Accumulate df corresponding to dx
    void addDForce(const core::MechanicalParams* mparams, DataVecDeriv& data_df1, DataVecDeriv& data_df2, const DataVecDeriv& data_dx1, const DataVecDeriv& data_dx2) override;
    using Inherit::addKToMatrix;
    void addKToMatrix(const sofa::core::MechanicalParams* mparams, const sofa::core::behavior::MultiMatrixAccessor* matrix) override;
    void draw(const core::visual::VisualParams* vparams) override;

};

#if !defined(SOFA_COMPONENT_FORCEFIELD_ORIGAMIFORCEFIELD_CPP)
extern template class SOFA_SOFAORIGAMI_API OrigamiForceField<defaulttype::Vec3Types>;
extern template class SOFA_SOFAORIGAMI_API OrigamiForceField<defaulttype::Rigid3Types>;
#endif

} // namespace sofa::component::interactionforcefield
