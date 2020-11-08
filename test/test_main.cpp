#pragma clang diagnostic push
#pragma ide diagnostic ignored "OCUnusedMacroInspection"
#define CATCH_CONFIG_MAIN
#pragma clang diagnostic pop

#undef DLIB_USE_LAPACK

#include <catch.hpp>

#include "../src/stdafx.h"
#include "../src/utils.hpp"
#include "../src/solnp.hpp"

#include "test_utils.cpp"
#include "test_solnp.cpp"