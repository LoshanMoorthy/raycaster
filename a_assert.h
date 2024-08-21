#pragma once

#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>

#include "l_log.h"

#ifndef ERROR
#define ERROR(fmt, ...) fprintf(stderr, fmt "\n", __VA_ARGS__)
#endif

static inline void _assert_impl(int _e, const char* _fmt, ...) {
    if (!_e) {
        printf("Assertion failed in file %s at line %d with message format: %s\n", __FILE__, __LINE__, _fmt);
        va_list args;
        va_start(args, _fmt);
        vfprintf(stderr, _fmt, args);
        va_end(args);
        exit(1);
    }
}

#define ASSERT1(_e) _assert_impl((_e), "%s", "assertion failed")
#define ASSERT2(_e, _fmt, ...) _assert_impl((_e), (_fmt), ##__VA_ARGS__)

#define GET_ASSERT_MACRO(_1, _2, NAME, ...) NAME

#define ASSERT(...) GET_ASSERT_MACRO(__VA_ARGS__, ASSERT2, ASSERT1)(__VA_ARGS__)

static inline int _vassert_impl(int _v, const char* _fmt, ...) {
    if (!_v) {
        va_list args;
        va_start(args, _fmt);
        vfprintf(stderr, _fmt, args);
        va_end(args);
        exit(1); 
    }
    return _v;
}

#define VASSERT1(_v) _vassert_impl((_v), "%s", "assertion failed")
#define VASSERT2(_v, _fmt, ...) _vassert_impl((_v), (_fmt), ##__VA_ARGS__)

#define GET_VASSERT_MACRO(_1, _2, NAME, ...) NAME

#define VASSERT(...) GET_VASSERT_MACRO(__VA_ARGS__, VASSERT2, VASSERT1)(__VA_ARGS__)
