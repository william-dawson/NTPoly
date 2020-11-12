#ifndef Logging_h
#define Logging_h

#include <string>

////////////////////////////////////////////////////////////////////////////////
//! Activate the logger.
void ActivateLogger();

////////////////////////////////////////////////////////////////////////////////
//! Activate the logger with a specific file.
//!\param file_name file to print to.
void ActivateLoggerFile(const std::string file_name);

////////////////////////////////////////////////////////////////////////////////
//! Deactivate the logger.
void DeactivateLogger();

#endif