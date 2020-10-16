#include "Albany_DealiiMeshStruct.hpp"

Albany::DealiiMeshStruct::DealiiMeshStruct(
    const Teuchos::RCP<Teuchos::ParameterList> &params,
    const Teuchos::RCP<const Teuchos_Comm> &commT) {
  meshSpecs.resize(1);
  meshSpecs[0] = Teuchos::rcp(new Albany::MeshSpecsStruct());
}

Albany::DealiiMeshStruct::~DealiiMeshStruct() {}

void Albany::DealiiMeshStruct::setFieldAndBulkData(
    const Teuchos::RCP<const Teuchos_Comm> &commT,
    const Teuchos::RCP<Teuchos::ParameterList> &params, const unsigned int neq_,
    const AbstractFieldContainer::FieldContainerRequirements &req,
    const Teuchos::RCP<Albany::StateInfoStruct> &sis,
    const unsigned int worksetSize,
    const std::map<std::string, Teuchos::RCP<Albany::StateInfoStruct>>
        &side_set_sis,
    const std::map<std::string,
                   AbstractFieldContainer::FieldContainerRequirements>
        &side_set_req) {}

Teuchos::ArrayRCP<Teuchos::RCP<Albany::MeshSpecsStruct>> &
Albany::DealiiMeshStruct::getMeshSpecs() {
  TEUCHOS_TEST_FOR_EXCEPTION(meshSpecs==Teuchos::null,
      std::logic_error,
      "meshSpecs accessed, but it has not been constructed" << std::endl);
  return meshSpecs;
}
const Teuchos::ArrayRCP<Teuchos::RCP<Albany::MeshSpecsStruct>> &
Albany::DealiiMeshStruct::getMeshSpecs() const {
  TEUCHOS_TEST_FOR_EXCEPTION(meshSpecs==Teuchos::null,
      std::logic_error,
      "meshSpecs accessed, but it has not been constructed" << std::endl);
  return meshSpecs;
}

Albany::AbstractMeshStruct::msType Albany::DealiiMeshStruct::meshSpecsType() {
  return Dealii_MS;
}
