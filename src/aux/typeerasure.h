#pragma once
#include <memory>

namespace Faunus {
/**
 * @namespace Type Erasure
 *
 * Class templates to assist implementing the type erasure pattern. The pairwise composition of concepts is also
 * supported.
 *
 * Type erasure pattern in general
 * https://www.modernescpp.com/index.php/c-core-guidelines-type-erasure-with-templates
 * This particular implementation was inspired by
 * https://aherrmann.github.io/programming/2014/10/19/type-erasure-with-merged-concepts/
 *
 * Basic agents are
 * - Concept – an interface
 * - Model – a wrapper implementing Concept
 * - Holder – a simple storage of an instance of a concrete implementation of Concept
 * - Container – a Model stored in an appropriate Holder
 * - Specification – a structure providing Concept, Model and Interface provided by the programmer
 * - TypeErasure – a template directing the type erasure composition from a specification
 * - MergeSpecifications – a template to combine two specifications into one
 *
 * @paragraph Example
 *
 * @code
 * struct PrintableSpecification {
 *    struct Concept {
 *        virtual ~Concept() = default;
 *        virtual void print() const = 0;
 *    };
 *
 *    template <class Holder>
 *    struct Model : public Holder, public virtual Concept {
 *        using Holder::Holder;
 *        void print() const override {
 *            Holder::get().print();
 *        };
 *    };
 *
 *    template <class Container>
 *    struct Interface : public Container {
 *        using Container::Container;
 *        void print() const {
 *            Container::get().print();
 *        };
 *    };
 * };
 * using PrintableType = TypeErasure::TypeErasure<PrintableSpecification>
 * using PrintableExportableType =
 *       TypeErasure::TypeErasure<TypeErasure::MergeSpecifications<PrintableSpecification, ExportableSpecification>>
 * @endcode
 *
 * @brief Class templates to assist implementing the type erasure pattern.
 */
namespace TypeErasure {
/**
 * @brief Internal templates and class templates.
 */
namespace Internal {
/**
 * @brief A universal storage of the concrete implementation with the move semantic in the constructor.
 * @tparam TImplementation  a concrete implementation of the required interface
 */
template <class TImplementation> class Holder {
  public:
    using Implementation = TImplementation;

    Holder(TImplementation implementation)
        : implementation(std::move(implementation))
    {
    }
    virtual ~Holder() = default;
    TImplementation& get() { return implementation; }
    const TImplementation& get() const { return implementation; }

  private:
    TImplementation implementation;
};

/**
 * @brief A universal storage of a reference to the concrete implementation.
 * @tparam TImplementation  a concrete implementation of the required interface
 */
template <class TImplementation> class ReferenceHolder {
  public:
    using Implementation = TImplementation;

    ReferenceHolder(TImplementation& implementation)
        : implementation(implementation)
    {
    }
    virtual ~ReferenceHolder() = default;
    TImplementation& get() { return implementation; }
    const TImplementation& get() const { return implementation; }

  private:
    TImplementation& implementation;
};

/**
 * @brief A container storing a model representing an particular concept in an appropriate holder.
 * @tparam TConcept  concept ecapsulated by the container
 * @tparam TModel  model template
 */
template <class TConcept, template <class> class TModel> // template template sytax for TModel<...>
class Container {
  public:
    using Concept = TConcept;

    /**
     * @brief Passes the implementation instance to the storage as an r-value reference using the move semantic.
     * @tparam TImplementation  a concrete implementation of the required interface
     * @param implementation  an implementation as an r-value, e.g., std::move(implementation)
     */
    template <class TImplementation>
    Container(TImplementation&& implementation)
        : self_ptr(std::make_shared<TModel<Holder<TImplementation>>>(std::move(implementation)))
    {
    }

    /**
     * @brief Passes the implementation instance to the storage as an l-value reference.
     * @tparam TImplementation  a concrete implementation of the required interface
     * @param implementation  an implementation as an l-value
     */
    template <class TImplementation>
    Container(TImplementation& implementation)
        : self_ptr(std::make_shared<TModel<ReferenceHolder<TImplementation>>>(implementation))
    {
    }

    const Concept& get() const { return *self_ptr; }
    Concept& get() { return *self_ptr; }

  private:
    std::shared_ptr<Concept> self_ptr;
    // std::shared_ptr<const Concept> self_ptr;
};

/**
 * @paragraph Helpers for extracting specification parts.
 */
//! @brief Extracts the concept type from the specification.
template <class TSpecification> using ConceptOf = typename TSpecification::Concept;
//! @brief Extracts the model type (including its holder type) from the specification.
template <class TSpecification, class THolder> using ModelOf = typename TSpecification::template Model<THolder>;
//! @brief Extracts the interface type from the specification.
template <class TSpecification, class TContainer>
using InterfaceOf = typename TSpecification::template Interface<TContainer>;
//! @brief Extracts the cotainer type from the specification.
template <class TSpecification>
using ContainerOf = Container<typename TSpecification::Concept, TSpecification::template Model>;
} // namespace Internal

/**
 * @brief Creates a type erasure type from the specification.
 *
 * @tparam TSpecification  specification of the concept, model and interface
 */
template <class TSpecification>
class TypeErasure : public Internal::InterfaceOf<TSpecification, Internal::ContainerOf<TSpecification>> {
    using Base = Internal::InterfaceOf<TSpecification, Internal::ContainerOf<TSpecification>>;

  public:
    using Base::Base; // include parent constructor
    using Specification = TSpecification;
};

/**
 * @brief Merges two specification together.
 * @tparam TSpecificationA the first specification
 * @tparam TSpecificationB the other specification
 */
template <class TSpecificationA, class TSpecificationB> struct MergeSpecifications {
    struct Concept : public virtual Internal::ConceptOf<TSpecificationA>,
                     public virtual Internal::ConceptOf<TSpecificationB> {};

    template <class THolder>
    struct Model : public Internal::ModelOf<TSpecificationA, Internal::ModelOf<TSpecificationB, THolder>>,
                   public virtual Concept {
        using Base = Internal::ModelOf<TSpecificationA, Internal::ModelOf<TSpecificationB, THolder>>;
        using Base::Base; // include parent constructor
    };

    template <class Container>
    struct Interface
        : public Internal::InterfaceOf<TSpecificationA, Internal::InterfaceOf<TSpecificationB, Container>> {

        using Base = Internal::InterfaceOf<TSpecificationA, Internal::InterfaceOf<TSpecificationB, Container>>;
        using Base::Base; // include parent constructor
    };
};
} // namespace TypeErasure
} // namespace Faunus
