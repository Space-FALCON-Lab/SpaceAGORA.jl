//////////////////////////////////////////////////////////////////////////
// The Global Reference Atmospheric Model (GRAM) Framework
//
// No recipient of this code should forward copies outside of the United 
// States without explicit approval by NASA Marshall Space Flight Center.
// 
// Module: GRAM common
//////////////////////////////////////////////////////////////////////////

#pragma once

#include "gram.h"

namespace GRAM {

//! \brief This is a utility class for storing variant data types with one class.
//!
//! This class gives the AtmosphereState the ability to store planet specific data yet still work within 
//! the common framework which is planet agnostic.  The Many class is designed store either an EarthAtmosphereState
//! or a MarsAtmosphereState or any other planet specific data.  Template functions are used to set and get
//! the appropriate data type.
//! \ingroup CommonGRAM Cpp_Venus Cpp_Earth Cpp_Mars Cpp_Jupiter Cpp_Saturn Cpp_Uranus Cpp_Neptune Cpp_Titan
class Many {
public:
  Many();
  Many(const Many& orig);
  ~Many();
  Many& operator=(const Many &value);

  template<typename T>
  void set(T &value);

  template<typename T>
  const T& get() const;

  template<typename T>
  T& get();

private:
  void* variantMemory = nullptr;  //!< A type-less pointer to the data.
  size_t size = 0;                //!< The size of the variant memory in bytes.
  bool dynamic = true;            //!< True when \p variantMemory has been dynamically allocated.  Otherwise false.
};

//! \brief Set the value of the variant memory.
//!
//! This method sets the variant memory to point at the \p value supplied as an arguement.  No memory is copied
//! by this method.  It only saves the reference to \p value.  If \p value is modified external to this object, then
//! it will also appear modified within this object.  If the \p value is deleted externally, then there could be 
//! some trouble.  So don't do that.
//! \param value The source value to set.
//! \tparam T The type of the source value.  This type does not need to be specified since it can be implied from the
//!           type of the value arguement.
template<typename T>
void Many::set(T &value)
{
  // Set the size.
  size = sizeof(value);

  // Store a reference to value as a type-less pointer.
  variantMemory = static_cast<void*>(&value);

  // Mark this memory as having not been dynamically allocated. (So that we don't accidentally delete it.)
  dynamic = false;
}

//! \brief Get a constant reference to the value of the variant memory.
//!
//! This method returns a type appropriate reference to the variant memory.  Since the Many class does not know the type
//! of the variant memory, the get template must specify this type (e.g. \p get<EarthAtmosphereState>()).
//! \code
//! Many manyObject;
//! EarthAtmosphereState eState;
//! manyObject.set(eState);
//! // Now get the value from the variant memory.
//! const EarthAtmosphereState& eAtmos = manyObject.get<EarthAtmosphereState>();
//! \endcode
//! The result of the code above is that eState and eAtmos are equivalent.
//! Note that the reference returned by this method is constant (non-modifiable).
//! \tparam T The type of the variant memory.
template<typename T>
const T& Many::get() const
{
  assert(variantMemory != nullptr);
  return *static_cast<const T*>(variantMemory);
}

//! \brief Get a reference to the value of the variant memory.
//!
//! This method returns a type appropriate reference to the variant memory.  Since the Many class does not know the type
//! of the variant memory, the get template must specify this type (e.g. \p get<EarthAtmosphereState>()).
//! \code
//! Many manyObject;
//! EarthAtmosphereState eState;
//! manyObject.set(eState);
//! // Now get the value from the variant memory.
//! EarthAtmosphereState& eAtmos = manyObject.get<EarthAtmosphereState>();
//! \endcode
//! The result of the code above is that eState and eAtmos are equivalent.
//! Note that the reference returned by this method is non-constant (modifiable).
//! \tparam T The type of the variant memory.
template<typename T>
T& Many::get()
{
  assert(variantMemory != nullptr);
  return *static_cast<T*>(variantMemory);
}


} // namespace