cmake_minimum_required(VERSION 3.12)

project(SofaOrigami_test)

find_package(SofaOrigami REQUIRED)

set(SOURCE_FILES
    MyBehaviorModel_test.cpp
    )

add_executable(${PROJECT_NAME} ${SOURCE_FILES})
target_link_libraries(${PROJECT_NAME} PUBLIC Sofa.Testing SofaOrigami)

add_test(NAME ${PROJECT_NAME} COMMAND ${PROJECT_NAME})
