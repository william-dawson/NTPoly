#ifndef Logging_h
#define Logging_h

#include <string>

////////////////////////////////////////////////////////////////////////////////
//! Activate the logger.
void ActivateLogger(bool start_document = false);

////////////////////////////////////////////////////////////////////////////////
//! Activate the logger with a specific file.
//!\param file_name file to print to.
void ActivateLoggerFile(const std::string file_name,
                        bool start_document = false);

////////////////////////////////////////////////////////////////////////////////
//! Deactivate the logger.
void DeactivateLogger();

#endif