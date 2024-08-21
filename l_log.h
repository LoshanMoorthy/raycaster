#pragma once

#include <stdio.h>
#include <stdarg.h>

#include "m_macros.h"
#include "t_types.h"

int long_fmtap(
    char* buf,
    usize n,
    const char* file,
    int line,
    const char* function,
    const char* prefix,
    const char* fmt,
    va_list ap
);

// returns number of written characters excluding null terminator
int log_fmt(
    char* buf,
    usize n,
    const char* file,
    int line,
    const char* function,
    const char* prefix,
    const char* fmt,
    ...
);

void log_write(
    const char* file,
    int line,
    const char* function,
    FILE* f,
    const char* prefix,
    const char* fmt,
    ...
);

// format as log message
#define LOGFMT(_buf, _n, _fmt, ...)                    \
    log_fmt(                                           \
        (_buf), (_n), __FILE__, __LINE__, __FUNCTION__,\
        "LOG", (_fmt), ##__VA_ARGS__)

// debug log message
#define LOG(_fmt, ...)                                 \
    log_write(                                         \
        __FILE__, __LINE__, __FUNCTION__,              \
        stdout, "LOG", (_fmt), ##__VA_ARGS__)

// debug warn
#define WARN(_fmt, ...)                                \
    log_write(                                         \
        __FILE__, __LINE__, __FUNCTION__,              \
        stderr, "WRN", (_fmt), ##__VA_ARGS__)

// debug error
#define ERROR(_fmt, ...)                               \
    log_write(                                         \
        __FILE__, __LINE__, __FUNCTION__,              \
        stderr, "ERR", (_fmt), ##__VA_ARGS__)