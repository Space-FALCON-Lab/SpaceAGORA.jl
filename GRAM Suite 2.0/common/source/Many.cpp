//////////////////////////////////////////////////////////////////////////
// The Global Reference Atmospheric Model (GRAM) Framework
//
// No recipient of this code should forward copies outside of the United 
// States without explicit approval by NASA Marshall Space Flight Center.
// 
// Module: GRAM common
//////////////////////////////////////////////////////////////////////////

#include <cstring>
#include "Many.h"

namespace GRAM {

//! \copydoc Atmosphere::Atmosphere()
Many::Many()
{
}

//! \copydoc Atmosphere::Atmosphere(const Atmosphere& orig)
Many::Many(const Many& orig) 
{
  // Make certain that the incoming value has content.
  if (orig.size > 0) {

    // Set the size.
    size = orig.size;

    // Allocate memory of that size.
    variantMemory = new char[size];

    // Mark this memory as having been dynamically allocated.
    dynamic = true;

    // Copy the contents of the source variant memory to this variant memory.
    memcpy(variantMemory, orig.variantMemory, size);
  }
}

//! \copydoc Atmosphere::~Atmosphere()
Many::~Many()
{
  // If the variant memory was dynamically allocated, then release that memory.
  if (dynamic && variantMemory != nullptr) {
    delete static_cast<char*>(variantMemory);
  }
}

//! \brief This defines the assignment operator.
//!
//! The assignment operator allows one to copy the contents of one Many object to another.  So if one codes:
//! \code
//! EarthAtmosphereState eAtmos;
//! Many A;
//! A.set(eAtmos);
//! Many B;
//! B = A;
//! \endcode
//! The result above will be that B has a copy of the eAtmos instance with the same values.
//! \param value  The source object to be copied.
//! \returns A reference to this object after the assignment.
Many& Many::operator=(const Many &value)
{
  // Make certain that the incoming value has content.
  if (value.size != 0) {

    // If this object doesn't have the memory same size as the source, then allocate the appropriate size.
    if (size != value.size) {

      // If the variant memory was dynamically allocated, then release that memory.
      if (dynamic && variantMemory != nullptr) {
        delete static_cast<char*>(variantMemory);
      }

      // Set the new size.
      size = value.size;

      // Allocate memory of that size.
      variantMemory = new char[size];

      // Mark this memory as having been dynamically allocated.
      dynamic = true;
    }

    // Copy the contents of the source variant memory to this variant memory.
    memcpy(variantMemory, value.variantMemory, size);
  }

  // Assignment operators have a return value so that we can chain them (A = B = C;).
  // So return a reference to this object.
  return *this;
}


} // namespace