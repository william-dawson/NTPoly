#ifndef Logging_h
#define Logging_h

#include <string>

////////////////////////////////////////////////////////////////////////////////
namespace NTPoly {
////////////////////////////////////////////////////////////////////////////////
//! Activate the logger.
//! \param start_document if this is a new document we can write the start 
//!       document marker.
void ActivateLogger(bool start_document = false);

////////////////////////////////////////////////////////////////////////////////
//! Activate the logger with a specific file.
//! \param file_name file to print to.
//! \param start_document if this is a new document we can write the start 
//! document marker.
void ActivateLogger(const std::string file_name, bool start_document = false);

////////////////////////////////////////////////////////////////////////////////
//! Deactivate the logger.
void DeactivateLogger();

} // namespace NTPoly
#endif
