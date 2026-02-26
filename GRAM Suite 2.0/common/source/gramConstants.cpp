//////////////////////////////////////////////////////////////////////////
// The Global Reference Atmospheric Model (GRAM) Framework
//
// No recipient of this code should forward copies outside of the United 
// States without explicit approval by NASA Marshall Space Flight Center.
// 
// Module: GRAM common
//////////////////////////////////////////////////////////////////////////

#include "error_strings.h"
#include "gram.h"

namespace GRAM {

bool CONSOLE_OUTPUT = true;

const std::string FILE_OPEN_ERROR_MESSAGE =
"Error: Unable to open file.\n"
"       Please check data paths.\n"
"       This is an unrecoverable error.\n"
"       Attempted to open: ";  // append path

const std::string FILE_PARSE_ERROR_MESSAGE =
"Error: Unable to parse data file."
"       This is an unrecoverable error.\n"
"       Attempting to parse: ";  // append path

const std::string FILE_READ_BINARY_ERROR_MESSAGE =
"Error: Unable to read binary data file."
"       This is an unrecoverable error.\n"
"       Attempting to read: ";  // append path

const std::string FILE_READ_TEXT_ERROR_MESSAGE =
"Error: Unable to read text data file."
"       This is an unrecoverable error.\n"
"       Attempting to read: ";  // append path

const std::string MEMORY_ALLOCATION_ERROR_MESSAGE =
"Error: An error occurred during memory allocation."
"       This is an unrecoverable error.\n       ";  // append details

const std::string SPICE_ERROR_MESSAGE =
"Error: A Spice error occurred.\n       ";  // append details

greal wrapDegrees(greal angle)
{
  if (angle < 0.0 || angle >= 360.0) {
    greal normalized = angle / 360.0;
    greal fraction = normalized - (int)(normalized)+1.0;
    return 360.0 * (fraction - (int)(fraction));
  }
  return angle;
}

greal wrapDegrees180(greal angle)
{
  greal wrapped = wrapDegrees(angle);
  if (wrapped > 180.0) {
    wrapped -= 360.0;
  }
  return wrapped;
}

greal wrapRadians(greal angle)
{
  if (angle > TWO_PI || angle < 0) {
    greal normalized = angle / TWO_PI;
    greal fraction = normalized - (int)(normalized)+1.0;
    return TWO_PI * (fraction - (int)(fraction));
  }
  else {
    return angle;
  }
}

greal wrapRadiansPi(greal angle)
{
  greal wrapped = wrapRadians(angle);
  if (wrapped > PI) {
    wrapped -= TWO_PI;
  }
  return wrapped;
}

} // namespace
