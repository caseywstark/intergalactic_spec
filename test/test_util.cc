/*
Test util functions.
*/

#include "util.h"

#include "catch.hpp"


TEST_CASE("periodic indexing", "[util]" ) {
    const int n = 128;
    int i;

    SECTION("positive, < size") {
        i = index_wrap(1, n);
        REQUIRE(i == 1);
    }
    SECTION("positive, >= size") {
        i = index_wrap(128, n);
        REQUIRE(i == 0);
        i = index_wrap(150, n);
        REQUIRE(i == 22);
    }
    SECTION("negative, < size") {
        i = index_wrap(-1, n);
        REQUIRE(i == 127);
        i = index_wrap(-11, n);
        REQUIRE(i == 117);
    }
    SECTION("negative, >= size") {
        i = index_wrap(-128, n);
        REQUIRE(i == 0);
        i = index_wrap(-129, n);
        REQUIRE(i == 127);
    }
}
