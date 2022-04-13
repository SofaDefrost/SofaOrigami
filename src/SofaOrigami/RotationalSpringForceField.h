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

#include <sofa/core/behavior/ForceField.h>
#include <sofa/core/topology/BaseMeshTopology.h>
#include <sofa/type/Vec.h>
#include <sofa/type/Mat.h>


// corotational triangle from
// @InProceedings{NPF05,
//   author       = "Nesme, Matthieu and Payan, Yohan and Faure, Fran\c{c}ois",
//   title        = "Efficient, Physically Plausible Finite Elements",
//   booktitle    = "Eurographics (short papers)",
//   month        = "august",
//   year         = "2005",
//   editor       = "J. Dingliana and F. Ganovelli",
//   keywords     = "animation, physical model, elasticity, finite elements",
//   url          = "http://www-evasion.imag.fr/Publications/2005/NPF05"
// }



namespace sofa::component::forcefield
{

/** Triangle FEM force field using the QR decomposition of the deformation gradient, inspired from http://www-evasion.imag.fr/Publications/2005/NPF05 , to handle large displacements.
  The material properties are uniform across the domain.
  Two methods are proposed, one for small displacements and one for large displacements.
  The method for small displacements has not been validated and we suspect that it is broke. Use it very carefully, and compare with the method for large displacements.
  */
template<class DataTypes>
class RotationalSpringForceField : public core::behavior::ForceField<DataTypes>
{
public:
    SOFA_CLASS(SOFA_TEMPLATE(RotationalSpringForceField, DataTypes), SOFA_TEMPLATE(core::behavior::ForceField, DataTypes));

    typedef core::behavior::ForceField<DataTypes> Inherited;
    typedef typename DataTypes::VecCoord VecCoord;
    typedef typename DataTypes::VecDeriv VecDeriv;
    typedef typename DataTypes::Coord    Coord   ;
    typedef typename DataTypes::Deriv    Deriv   ;
    typedef typename Coord::value_type   Real    ;

    typedef core::objectmodel::Data<VecCoord> DataVecCoord;
    typedef core::objectmodel::Data<VecDeriv> DataVecDeriv;

    typedef type::Vec<6, Real> Displacement;								///< the displacement vector
    typedef type::Mat<3, 3, Real> MaterialStiffness;						///< the matrix of material stiffness
    typedef sofa::type::vector<MaterialStiffness> VecMaterialStiffness;    ///< a vector of material stiffness matrices
    typedef type::Mat<6, 3, Real> StrainDisplacement;						///< the strain-displacement matrix (the transpose, actually)
    typedef sofa::type::vector<StrainDisplacement> VecStrainDisplacement;	///< a vector of strain-displacement matrices
    typedef type::Mat<3, 3, Real > Transformation;						///< matrix for rigid transformations like rotations
    /// Stiffness matrix ( = RJKJtRt  with K the Material stiffness matrix, J the strain-displacement matrix, and R the transformation matrix if any )
    typedef type::Mat<9, 9, Real> StiffnessMatrix;


    typedef sofa::core::topology::BaseMeshTopology::Index Index;
    typedef sofa::core::topology::BaseMeshTopology::Triangle Element;
    typedef sofa::core::topology::BaseMeshTopology::SeqTriangles VecElement;

    static const int SMALL = 1;										///< Symbol of small displacements triangle solver
    static const int LARGE = 0;										///< Symbol of large displacements triangle solver

protected:
    VecMaterialStiffness _materialsStiffnesses;						///< the material stiffness matrices vector
    VecStrainDisplacement _strainDisplacements;						///< the strain-displacement matrices vector

    const VecElement *_indexedElements;
    Data< VecCoord > _initialPoints; ///< the intial positions of the points

    RotationalSpringForceField();
    virtual ~RotationalSpringForceField();

    sofa::core::topology::BaseMeshTopology* m_topology;

public:

    void init() override;
    void addForce(const core::MechanicalParams* mparams, DataVecDeriv& f, const DataVecCoord& x, const DataVecDeriv& v) override;
    void addDForce(const core::MechanicalParams* mparams, DataVecDeriv& df, const DataVecDeriv& dx) override;
    SReal getPotentialEnergy(const core::MechanicalParams* /*mparams*/, const DataVecCoord&  /* x */) const override
    {
        msg_warning() << "Method getPotentialEnergy not implemented yet.";
        return 0.0;
    }

    void draw(const core::visual::VisualParams* vparams) override;

    int method;
    type::vector<int> roundCount;
    type::vector<double> signTemp;
    type::vector<type::vector<int>> adjacentVertices;
    Data<std::string> f_method; ///< Choice of method: 0 for small, 1 for large displacements
    Data<Real> f_poisson;       ///< Poisson ratio of the material
    Data<Real> f_kcrease;         ///< Young modulus of the material
    Data<Real> f_thickness;     ///< Thickness of the elements
//    Data<Real> f_damping;       ///< Damping coefficient of the material, currently unused
    Data<bool> f_planeStrain; ///< compute material stiffness corresponding to the plane strain assumption, or to the plane stress otherwise.
    Data<type::vector<std::string>> d_edgesAssignment;
    Data<Real> d_angleTarget;
    type::vector<Real> edgeInitialLength;
    type::vector< type::Mat<12, 12, Real> > edgeStiffness;
    type::vector< VecDeriv > edgeStiffness2;
    type::vector< type::vector<Index> > nodes;
    type::vector<Real> forceCrease;

    void getRotationalStiffness(const Real theta, const int edgeNum, Real &M, Real &k);
    Real getPoisson() { return f_poisson.getValue(); }
    void setPoisson(Real val);
    Real getYoung() { return f_kcrease.getValue(); }
//    Real getDamping() { return f_damping.getValue(); }
//    void setDamping(Real val) { f_damping.setValue(val); }
    void setYoung(Real val);
    int  getMethod() { return method; }
    void setMethod(int val);
    void setMethod(std::string val);

    /// Public methods to access FEM information per element. Those method should not be used internally as they add check on element id.
    const type::fixed_array <Coord, 3>& getRotatedInitialElement(Index elemId);
    const Transformation& getRotationMatrix(Index elemId);
    const MaterialStiffness& getMaterialStiffness(Index elemId);
    const StrainDisplacement& getStrainDisplacements(Index elemId);

    /// Link to be set to the topology container in the component graph.
    SingleLink<RotationalSpringForceField<DataTypes>, sofa::core::topology::BaseMeshTopology, BaseLink::FLAG_STOREPATH | BaseLink::FLAG_STRONGLINK> l_topology;

protected :

    ////////////// large displacements method
    sofa::type::vector< type::fixed_array <Coord, 3> > _rotatedInitialElements;   ///< The initials positions in its frame
    sofa::type::vector< Transformation > _rotations;
    void addKToMatrix(sofa::linearalgebra::BaseMatrix *mat, SReal k, unsigned int &offset) override; // compute and add all the element stiffnesses to the global stiffness matrix
    type::Mat3x3 diamond(Coord a, Coord b);
    type::Mat<3, 3, Real> InvalidTransform;
    type::fixed_array <Coord, 3> InvalidCoords;
    StrainDisplacement InvalidStrainDisplacement;
    /// Pointer to the utils class which store methods common to RotationalSpringForceField
};


#if  !defined(SOFA_COMPONENT_FORCEFIELD_RotationalSpringForceField_CPP)

extern template class SOFA_SOFAORIGAMI_API RotationalSpringForceField<sofa::defaulttype::Vec3Types>;


#endif //  !defined(SOFA_COMPONENT_FORCEFIELD_RotationalSpringForceField_CPP)

} // namespace sofa::component::forcefield
