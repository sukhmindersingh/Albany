#pragma once

#include <utility>
#include <vector>

#include "Albany_AbstractDiscretization.hpp"
#include "Albany_DataTypes.hpp"

// dealii
#include <deal.II/grid/tria.h>


namespace Albany {

class DealiiDiscretization : public AbstractDiscretization {

public:
  //! Constructor
  DealiiDiscretization(const Teuchos::RCP<Teuchos::ParameterList> &discParams,
                       const Teuchos::RCP<const Teuchos_Comm> &comm);

  //! Get node vector space (owned and overlapped)
  Teuchos::RCP<const Thyra_VectorSpace> getNodeVectorSpace() const override {
    return m_node_vs;
  }
  Teuchos::RCP<const Thyra_VectorSpace>
  getOverlapNodeVectorSpace() const override {
    return m_overlap_node_vs;
  }

  //! Get solution DOF vector space (owned and overlapped).
  Teuchos::RCP<const Thyra_VectorSpace> getVectorSpace() const override {
    return m_vs;
  }
  Teuchos::RCP<const Thyra_VectorSpace> getOverlapVectorSpace() const override {
    return m_overlap_vs;
  }

  //! Get Field node vector space (owned and overlapped)
  Teuchos::RCP<const Thyra_VectorSpace>
  getNodeVectorSpace(const std::string &/*field_name*/) const override {
    TEUCHOS_TEST_FOR_EXCEPTION(
        true, std::logic_error,
        "Albany::DealiiDiscretization: getNodeVectorSpace(field_name) not "
        "implemented yet");
    TEUCHOS_UNREACHABLE_RETURN(Teuchos::null);
  }
  Teuchos::RCP<const Thyra_VectorSpace>
  getOverlapNodeVectorSpace(const std::string &/*field_name*/) const override {

    TEUCHOS_TEST_FOR_EXCEPTION(
        true, std::logic_error,
        "Albany::DealiiDiscretization: getOverlapNodeVectorSpace(field_name) "
        "not implemented yet");
    TEUCHOS_UNREACHABLE_RETURN(Teuchos::null);
  }

  //! Get Field vector space (owned and overlapped)
  Teuchos::RCP<const Thyra_VectorSpace>
  getVectorSpace(const std::string &/*field_name*/) const override {
    TEUCHOS_TEST_FOR_EXCEPTION(
        true, std::logic_error,
        "Albany::DealiiDiscretization: getVectorSpace(field_name) not implemented "
        "yet");
    TEUCHOS_UNREACHABLE_RETURN(Teuchos::null);
  }
  Teuchos::RCP<const Thyra_VectorSpace>
  getOverlapVectorSpace(const std::string & /*field_name*/) const override {
    TEUCHOS_TEST_FOR_EXCEPTION(
        true, std::logic_error,
        "Albany::DealiiDiscretization: getOverlapVectorSpace(field_name) not "
        "implemented "
        "yet");
    TEUCHOS_UNREACHABLE_RETURN(Teuchos::null);
  }

  //! Create a Jacobian operator (owned and overlapped)
  Teuchos::RCP<Thyra_LinearOp> createJacobianOp() const override {
    return m_jac_factory->createOp();
  }
  Teuchos::RCP<Thyra_LinearOp> createOverlapJacobianOp() const override {
    return m_overlap_jac_factory->createOp();

  }
  
  //! Returns boolean telling code whether explicit scheme is used (needed for
  //! Aeras problems only)
  bool isExplicitScheme() const override {
    return false;
  }

  //! Get Node set lists
  const NodeSetList &getNodeSets() const override {
    TEUCHOS_TEST_FOR_EXCEPTION(
        true, std::logic_error,
        "Albany::DealiiDiscretization: getNodeSets() not "
        "implemented "
        "yet");
  }

  const NodeSetGIDsList &getNodeSetGIDs() const override {
    TEUCHOS_TEST_FOR_EXCEPTION(
        true, std::logic_error,
        "Albany::DealiiDiscretization: getNodeSetGIDs() not "
        "implemented "
        "yet");
  }

  const NodeSetCoordList &getNodeSetCoords() const override {
    TEUCHOS_TEST_FOR_EXCEPTION(
        true, std::logic_error,
        "Albany::DealiiDiscretization: getNodeSetCoords() not "
        "implemented "
        "yet");
  }

  //! Get Side set lists
  const SideSetList &getSideSets(const int ws) const override {
    TEUCHOS_TEST_FOR_EXCEPTION(
        true, std::logic_error,
        "Albany::DealiiDiscretization: getNodeSetCoords(ws) not "
        "implemented "
        "yet");
  }

  //! Get map from (Ws, El, Local Node, Eq) -> unkLID
  const Conn &getWsElNodeEqID() const override {
    TEUCHOS_TEST_FOR_EXCEPTION(
        true, std::logic_error,
        "Albany::DealiiDiscretization: getWsElNodeEqID() not "
        "implemented "
        "yet");
  }

  //! Get map from (Ws, El, Local Node) -> unkGID
  const WorksetArray<Teuchos::ArrayRCP<Teuchos::ArrayRCP<GO>>>::type &
    getWsElNodeID() const override {
    TEUCHOS_TEST_FOR_EXCEPTION(
        true, std::logic_error,
        "Albany::DealiiDiscretization: getWsElNodeID() not "
        "implemented "
        "yet");

    }

  //! Get IDArray for (Ws, Local Node, nComps) -> (local) NodeLID, works for
  //! both scalar and vector fields
  const std::vector<IDArray> &
    getElNodeEqID(const std::string &field_name) const override {
    TEUCHOS_TEST_FOR_EXCEPTION(
        true, std::logic_error,
        "Albany::DealiiDiscretization: getElNodeEqID(field_name) not "
        "implemented "
        "yet");
    }

  //! Get Dof Manager of field field_name
    const NodalDOFManager &
    getDOFManager(const std::string &field_name) const override {
      TEUCHOS_TEST_FOR_EXCEPTION(
          true, std::logic_error,
          "Albany::DealiiDiscretization: NodalDOFManager(field_name) not "
          "implemented "
          "yet");
    }

  //! Get Dof Manager of field field_name
  const NodalDOFManager &
    getOverlapDOFManager(const std::string &field_name) const override {
    TEUCHOS_TEST_FOR_EXCEPTION(
        true, std::logic_error,
        "Albany::DealiiDiscretization: getOverlapDOFManager(field_name) not "
        "implemented "
        "yet");
    }

  //! Retrieve coodinate ptr_field (ws, el, node)
  const WorksetArray<Teuchos::ArrayRCP<Teuchos::ArrayRCP<double *>>>::type &
    getCoords() const override {
    TEUCHOS_TEST_FOR_EXCEPTION(
        true, std::logic_error,
        "Albany::DealiiDiscretization: getCoords() not "
        "implemented "
        "yet");
    }

  //! Get coordinates (overlap map).
  const Teuchos::ArrayRCP<double> &getCoordinates() const override {
    TEUCHOS_TEST_FOR_EXCEPTION(
        true, std::logic_error,
        "Albany::DealiiDiscretization: getCoordinates() not "
        "implemented "
        "yet");
  }

  //! Set coordinates (overlap map) for mesh adaptation.
  void setCoordinates(const Teuchos::ArrayRCP<const double> &c) override {
    TEUCHOS_TEST_FOR_EXCEPTION(
        true, std::logic_error,
        "Albany::DealiiDiscretization: setCoordinates(c) not "
        "implemented "
        "yet");
  }

  //! The reference configuration manager handles updating the reference
  //! configuration. This is only relevant, and also only optional, in the
  //! case of mesh adaptation.
  void setReferenceConfigurationManager(
      const Teuchos::RCP<AAdapt::rc::Manager> &rcm) override {
    TEUCHOS_TEST_FOR_EXCEPTION(
        true, std::logic_error,
        "Albany::DealiiDiscretization: setReferenceConfigurationManager() not "
        "implemented "
        "yet");
  }
  //#ifdef ALBANY_CONTACT
  //! Get the contact manager
  // virtual Teuchos::RCP<const ContactManager>
  // getContactManager() const = 0;
  //#endif

  const WorksetArray<Teuchos::ArrayRCP<double>>::type &
  getSphereVolume() const override {
    TEUCHOS_TEST_FOR_EXCEPTION(
        true, std::logic_error,
        "Albany::DealiiDiscretization: getSphereVolume() not "
        "implemented "
        "yet");
  }

  const WorksetArray<Teuchos::ArrayRCP<double *>>::type &
  getLatticeOrientation() const override {
    TEUCHOS_TEST_FOR_EXCEPTION(
        true, std::logic_error,
        "Albany::DealiiDiscretization: getLatticeOrientation() not "
        "implemented "
        "yet");
  }

  //! Print the coords for mesh debugging
  void printCoords() const override {
    TEUCHOS_TEST_FOR_EXCEPTION(
        true, std::logic_error,
        "Albany::DealiiDiscretization: printCoords() not "
        "implemented "
        "yet");
  }

  //! Get sideSet discretizations map
  const SideSetDiscretizationsType &getSideSetDiscretizations() const override {
    TEUCHOS_TEST_FOR_EXCEPTION(
        true, std::logic_error,
        "Albany::DealiiDiscretization: getSideSetDiscretizations() not "
        "implemented "
        "yet");
  }

  //! Get the map side_id->side_set_elem_id
  const std::map<std::string, std::map<GO, GO>> &
  getSideToSideSetCellMap() const override {
    TEUCHOS_TEST_FOR_EXCEPTION(
        true, std::logic_error,
        "Albany::DealiiDiscretization: getSideToSideSetCellMap() not "
        "implemented "
        "yet");
  }

  //! Get the map side_node_id->side_set_cell_node_id
  const std::map<std::string, std::map<GO, std::vector<int>>> &
  getSideNodeNumerationMap() const override {
    TEUCHOS_TEST_FOR_EXCEPTION(
        true, std::logic_error,
        "Albany::DealiiDiscretization: getSideNodeNumerationMap() not "
        "implemented "
        "yet");
  }

  //! Get MeshStruct
  Teuchos::RCP<AbstractMeshStruct> getMeshStruct() const override {
    TEUCHOS_TEST_FOR_EXCEPTION(
        true, std::logic_error,
        "Albany::DealiiDiscretization: getMeshStruct() not "
        "implemented "
        "yet");
    TEUCHOS_UNREACHABLE_RETURN(Teuchos::null);
  }

  //! Set stateArrays
  void setStateArrays(StateArrays &sa) override {
    stateArrays = sa;
    return;
  }

  //! Get stateArrays
  StateArrays &getStateArrays() override {
    return stateArrays;
  }

  //! Get nodal parameters state info struct
  const StateInfoStruct &getNodalParameterSIS() const override {}

  //! Retrieve Vector (length num worksets) of element block names
  const WorksetArray<std::string>::type &getWsEBNames() const override {}

  //! Retrieve Vector (length num worksets) of Physics Index
  const WorksetArray<int>::type &getWsPhysIndex() const override {}

  //! Retrieve connectivity map from elementGID to workset
  WsLIDList &getElemGIDws() override {}
  const WsLIDList &getElemGIDws() const override {}

  //! Flag if solution has a restart values -- used in Init Cond
  bool hasRestartSolution() const override {}

  //! File time of restart solution
  double
  restartDataTime() const override {}

  //! Get number of spatial dimensions
  int getNumDim() const override {}

  //! Get number of total DOFs per node
  int getNumEq() const override {}
  
  //! Get Numbering for layered mesh (mesh structred in one direction)
  Teuchos::RCP<LayeredMeshNumbering<LO>>
  getLayeredMeshNumbering() const override {}

  // --- Methods to write solution in the output file --- //
  Teuchos::RCP<Thyra_Vector>
  getSolutionField(bool overlapped = false) const override {}
  Teuchos::RCP<Thyra_MultiVector>
  getSolutionMV(bool overlapped = false) const override {}
#if defined(ALBANY_LCM)
  void setResidualField(const Thyra_Vector &residual) override {}
#endif
  void getField(Thyra_Vector &field_vector,
                const std::string &field_name) const override {}
  void setField(const Thyra_Vector &field_vector, const std::string &field_name,
                bool overlapped) override {}

  //! Write the solution to the output file. Calls next two together.
  void writeSolution(const Thyra_Vector &solution, const double time,
                     const bool overlapped = false) override {}
  void writeSolution(const Thyra_Vector &solution,
                     const Thyra_Vector &solution_dot, const double time,
                     const bool overlapped = false) override {}
  void writeSolution(const Thyra_Vector &solution,
                     const Thyra_Vector &solution_dot,
                     const Thyra_Vector &solution_dotdot, const double time,
                     const bool overlapped = false) override {}
  void writeSolutionMV(const Thyra_MultiVector &solution, const double time,
                       const bool overlapped = false) override {}

  //! Write the solution to the mesh database.
  void writeSolutionToMeshDatabase(const Thyra_Vector &solution,
                                   const double time,
                                   const bool overlapped = false) override {}
  void writeSolutionToMeshDatabase(const Thyra_Vector &solution,
                                   const Thyra_Vector &solution_dot,
                                   const double time,
                                   const bool overlapped = false) override {}
  void writeSolutionToMeshDatabase(const Thyra_Vector &solution,
                                   const Thyra_Vector &solution_dot,
                                   const Thyra_Vector &solution_dotdot,
                                   const double time,
                                   const bool overlapped = false) override {}
  void writeSolutionMVToMeshDatabase(const Thyra_MultiVector &solution,
                                     const double time,
                                     const bool overlapped = false) override {}

  //! Write the solution to file. Must call writeSolution first.
  void writeSolutionToFile(const Thyra_Vector &solution, const double time,
                           const bool overlapped = false) override {}
  void writeSolutionMVToFile(const Thyra_MultiVector &solution,
                             const double time,
                             const bool overlapped = false) override {}

protected:

  //! Teuchos communicator
  Teuchos::RCP<const Teuchos_Comm>        comm;

  //! Unknown map and node map
  Teuchos::RCP<const Thyra_VectorSpace>   m_vs;
  Teuchos::RCP<const Thyra_VectorSpace>   m_node_vs;

  //! Overlapped unknown map and node map
  Teuchos::RCP<const Thyra_VectorSpace>   m_overlap_vs;
  Teuchos::RCP<const Thyra_VectorSpace>   m_overlap_node_vs;

  //! Jacobian matrix graph proxy (owned and overlap)
  Teuchos::RCP<ThyraCrsMatrixFactory> m_jac_factory;
  Teuchos::RCP<ThyraCrsMatrixFactory> m_overlap_jac_factory;

  StateArrays stateArrays;
};

} // namespace Albany
