#include "Logging.h"
using std::string;

////////////////////////////////////////////////////////////////////////////////
extern "C" {
#include "Logging_c.h"
}

////////////////////////////////////////////////////////////////////////////////
void NTPoly::ActivateLogger(bool start_document) {
  ActivateLogger_wrp(&start_document);
}

////////////////////////////////////////////////////////////////////////////////
void NTPoly::ActivateLoggerFile(const string file_name, bool start_document) {
  int string_length = file_name.length();
  string temp = file_name;
  ActivateLoggerFile_wrp(&start_document, &temp.c_str()[0], &string_length);
}

////////////////////////////////////////////////////////////////////////////////
void NTPoly::DeactivateLogger() { DeactivateLogger_wrp(); }
