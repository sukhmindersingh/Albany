#pragma once

#include "Albany_AbstractMeshStruct.hpp"

#include "Teuchos_ScalarTraits.hpp"

namespace Albany {

struct DealiiMeshStruct : public AbstractMeshStruct {
public:
  DealiiMeshStruct(const Teuchos::RCP<Teuchos::ParameterList> &params,
                   const Teuchos::RCP<const Teuchos_Comm> &commT);
  ~DealiiMeshStruct();

  void setFieldAndBulkData(
      const Teuchos::RCP<const Teuchos_Comm> &commT,
      const Teuchos::RCP<Teuchos::ParameterList> &params,
      const unsigned int neq_,
      const AbstractFieldContainer::FieldContainerRequirements &req,
      const Teuchos::RCP<Albany::StateInfoStruct> &sis,
      const unsigned int worksetSize,
      const std::map<std::string, Teuchos::RCP<Albany::StateInfoStruct>>
          &side_set_sis,
      const std::map<std::string,
                     AbstractFieldContainer::FieldContainerRequirements>
          &side_set_req);

  Teuchos::ArrayRCP<Teuchos::RCP<Albany::MeshSpecsStruct>> &getMeshSpecs();
  const Teuchos::ArrayRCP<Teuchos::RCP<Albany::MeshSpecsStruct>> &
  getMeshSpecs() const;

  msType meshSpecsType();

protected:
  Teuchos::ArrayRCP<Teuchos::RCP<Albany::MeshSpecsStruct>> meshSpecs;

};

} // namespace Albany
