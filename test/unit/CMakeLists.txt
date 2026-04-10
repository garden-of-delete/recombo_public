add_executable(unitTest
        clk_test.cpp
        newsud_test.cpp
        pseudorandom_test.cpp
        recombo_criteria_test.cpp)
target_link_libraries(unitTest recombo_core gtest gtest_main)
add_test(NAME unitTest COMMAND $<TARGET_FILE:unitTest>)
