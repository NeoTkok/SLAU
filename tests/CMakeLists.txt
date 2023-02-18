include(FetchContent)
FetchContent_Declare(
        googletest
        URL https://github.com/google/googletest/archive/refs/tags/release-1.11.0.zip
)
FetchContent_MakeAvailable(googletest)
add_library(GTest::GTest INTERFACE IMPORTED)
target_link_libraries(GTest::GTest INTERFACE gtest_main)


add_executable(
        testSolve_Mat3Diag
        testSolve_Mat3Dig.cpp)


target_link_libraries(
        testSolve_Mat3Diag PRIVATE
        GTest::GTest
)
