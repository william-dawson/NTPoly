#ifndef Logging_ch
#define Logging_ch

void ActivateLogger_wrp(const bool *start_document);
void ActivateLoggerFile_wrp(const bool *start_document, const char *file_name,
                            const int *name_size);
void DeactivateLogger_wrp();

#endif
